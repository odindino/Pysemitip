#!/usr/bin/env python3
"""
調試run_multint中的Pot0計算問題
"""
import logging
import sys
import os

# 設置debug級別日志
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(name)s - %(message)s')

# 添加項目路徑
sys.path.insert(0, os.path.abspath('.'))

from src.simulation.multint import MultIntSimulation
from src.core.filereader import YamlConfigReader

def debug_run_multint():
    """調試run_multint的Pot0計算"""
    print("調試run_multint中的Pot0計算")
    print("="*80)
    
    config_file = "data/input/examples/test/quick_test.yaml"
    
    # 讀取配置
    filereader = YamlConfigReader()
    try:
        config = filereader.read_and_validate_yaml(config_file)
        print(f"✅ 配置讀取成功")
    except Exception as e:
        print(f"❌ 配置讀取失敗: {e}")
        return
    
    # 創建模擬
    try:
        simulation = MultIntSimulation(config)
        print(f"✅ 模擬初始化成功")
    except Exception as e:
        print(f"❌ 模擬初始化失敗: {e}")
        return
    
    # 執行一次SCF迭代來檢查Pot0
    print(f"\n開始SCF循環...")
    try:
        # 手動執行一次來查看詳細信息
        V_tip = -2.07
        V_sample = 0.0
        system_fermi = 1.4187
        
        # 先做Laplace求解
        print(f"\n1. Laplace求解階段:")
        potential, iters, error = simulation.poisson_solver.solve_laplace(
            V_tip, V_sample, max_iterations=200, tolerance=1e-4
        )
        
        # 檢查Laplace結果
        pot0_laplace_raw = simulation.poisson_solver._calculate_pot0_fortran_style(potential, apply_scaling_correction=False)
        pot0_laplace_scaled = simulation.poisson_solver._calculate_pot0_fortran_style(potential, apply_scaling_correction=True)
        
        print(f"  Laplace Pot0 (原始): {pot0_laplace_raw:.6f} V")
        print(f"  Laplace Pot0 (縮放): {pot0_laplace_scaled:.6f} V")
        
        # 檢查VSINT完整求解
        print(f"\n2. VSINT完整求解階段:")
        potential_vsint, iters_vs, converged = simulation.poisson_solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=simulation.charge_density_calculator,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=200,
            tolerance_Volts=1e-4
        )
        
        # 手動計算所有變體
        vsint_array = simulation.poisson_solver._initialize_vsint_array()
        vsint_array = simulation.poisson_solver._update_vsint_with_surface_charge(
            vsint_array, potential_vsint, simulation.charge_density_calculator,
            system_fermi, V_tip)
        
        pot0_regular_raw = simulation.poisson_solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=False, apply_scaling_correction=False)
        pot0_regular_scaled = simulation.poisson_solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=False, apply_scaling_correction=True)
        pot0_vsint_raw = simulation.poisson_solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=False)
        pot0_vsint_scaled = simulation.poisson_solver._calculate_pot0_fortran_style(
            potential_vsint, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
        
        print(f"  VSINT結果分析:")
        print(f"    Regular (原始):     {pot0_regular_raw:.6f} V")
        print(f"    Regular (縮放):     {pot0_regular_scaled:.6f} V")
        print(f"    VSINT (原始):      {pot0_vsint_raw:.6f} V")
        print(f"    🎯 VSINT (縮放):   {pot0_vsint_scaled:.6f} V")
        
        print(f"  VSINT陣列值:")
        print(f"    VSINT[0,0] = {vsint_array[0,0]:.6f} V")
        print(f"    VSINT[1,0] = {vsint_array[1,0]:.6f} V")
        
        print(f"\n3. 與Fortran目標(-0.08V)比較:")
        fortran_target = -0.08
        print(f"    Laplace (縮放) 差異: {abs(pot0_laplace_scaled - fortran_target):.6f} V")
        print(f"    VSINT (縮放) 差異:   {abs(pot0_vsint_scaled - fortran_target):.6f} V")
        
        if abs(pot0_vsint_scaled - fortran_target) < 0.01:
            print(f"    ✅ VSINT達到高精度!")
        elif abs(pot0_vsint_scaled - fortran_target) < abs(pot0_laplace_scaled - fortran_target):
            print(f"    📈 VSINT比Laplace更好")
        else:
            print(f"    ⚠️  VSINT沒有改善")
        
        # 檢查為什麼run_multint看到的是-0.16V
        print(f"\n4. run_multint顯示差異分析:")
        print(f"    run_multint大約顯示: -0.16V")
        print(f"    這接近: Laplace (縮放) = {pot0_laplace_scaled:.6f} V")
        print(f"    說明run_multint主要顯示的是Laplace階段的結果")
        print(f"    而不是完整VSINT求解的結果")
        
    except Exception as e:
        print(f"❌ SCF執行失敗: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_run_multint()