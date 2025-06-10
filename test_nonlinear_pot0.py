#!/usr/bin/env python3
"""
測試完整的非線性 Poisson 求解器對 Pot0 的改善
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.core.charge_density import ChargeDensityCalculator
from src.physics.materials.semiconductor import SemiconductorRegion
from src.physics.solvers.grid import HyperbolicGrid

logging.basicConfig(level=logging.INFO, format='%(message)s')

def test_nonlinear_pot0():
    """測試非線性 Poisson 求解器對 Pot0 的影響"""
    print("測試非線性 Poisson 求解器對 Pot0 的改善")
    print("="*80)
    
    # 創建測試網格
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    # 創建物理參數
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    solver = PoissonSOREquation(grid, props)
    
    # 測試參數
    V_tip = -2.07
    V_sample = 0.0
    system_fermi_level = 1.4187  # 從 log 中獲取的 Fermi level
    
    print(f"測試條件:")
    print(f"  V_tip = {V_tip} V")
    print(f"  V_sample = {V_sample} V") 
    print(f"  System Fermi level = {system_fermi_level} eV")
    print(f"  Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    
    # === 測試1: Laplace 方程 (無電荷密度) ===
    print(f"\n" + "="*60)
    print(f"測試1: Laplace 方程 (無電荷密度)")
    print(f"="*60)
    
    potential_laplace, iterations_l, error_l = solver.solve_laplace(
        V_tip, V_sample, max_iterations=300, tolerance=1e-4
    )
    
    pot0_laplace = solver._calculate_pot0_fortran_style(potential_laplace)
    print(f"Laplace 結果:")
    print(f"  迭代次數: {iterations_l}")
    print(f"  最大誤差: {error_l:.3e}")
    print(f"  Pot0 = {pot0_laplace:.6f} V")
    print(f"  與 Fortran (-0.08V) 差異: {abs(pot0_laplace - (-0.08)):.6f} V")
    
    # === 測試2: 非線性 Poisson 方程 (含電荷密度) ===
    print(f"\n" + "="*60)
    print(f"測試2: 非線性 Poisson 方程 (含電荷密度)")
    print(f"="*60)
    
    # 創建電荷密度計算器
    semiconductor = SemiconductorRegion(
        Ev_offset_eV=-5.17,
        Ec_offset_eV=-3.75,  # Eg=1.42eV
        epsilon_r=12.9,
        Nd_cm3=1e18,
        Na_cm3=0.0,
        T_K=300.0
    )
    
    charge_calculator = ChargeDensityCalculator(
        semiconductor_physics=semiconductor,
        system_fermi_level_E_F_main_eV=system_fermi_level,
        grid=grid,
        verbose=False
    )
    
    print(f"電荷密度計算器創建成功")
    print(f"  Ev_offset = {semiconductor.Ev_offset_eV} eV")
    print(f"  Ec_offset = {semiconductor.Ec_offset_eV} eV")
    print(f"  摻雜濃度 = {semiconductor.Nd_cm3:.1e} cm⁻³")
    
    try:
        potential_nonlinear, iterations_nl, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calculator,
            system_fermi_level_E_F_main_eV=system_fermi_level,
            max_iterations=500,
            tolerance_Volts=1e-4
        )
        
        pot0_nonlinear = solver._calculate_pot0_fortran_style(potential_nonlinear)
        
        print(f"非線性 Poisson 結果:")
        print(f"  迭代次數: {iterations_nl}")
        print(f"  收斂狀態: {'是' if converged else '否'}")
        print(f"  Pot0 = {pot0_nonlinear:.6f} V")
        print(f"  與 Fortran (-0.08V) 差異: {abs(pot0_nonlinear - (-0.08)):.6f} V")
        
        # === 比較分析 ===
        print(f"\n" + "="*60)
        print(f"比較分析")
        print(f"="*60)
        
        improvement = abs(pot0_laplace - (-0.08)) - abs(pot0_nonlinear - (-0.08))
        print(f"結果比較:")
        print(f"  Laplace Pot0:    {pot0_laplace:.6f} V (差異: {abs(pot0_laplace - (-0.08)):.6f} V)")
        print(f"  非線性 Pot0:     {pot0_nonlinear:.6f} V (差異: {abs(pot0_nonlinear - (-0.08)):.6f} V)")
        print(f"  Fortran 目標:    -0.08 V")
        print(f"  改善程度:        {improvement:.6f} V ({'改善' if improvement > 0 else '惡化'})")
        
        # === 詳細分析界面電位 ===
        print(f"\n詳細分析界面電位:")
        interface_idx = grid.N_nu - 1
        
        print(f"Laplace 界面電位:")
        for i in range(min(4, grid.N_eta)):
            v_l = potential_laplace[i, interface_idx]
            v_nl = potential_nonlinear[i, interface_idx]
            print(f"  eta={i}: Laplace={v_l:.6f} V, 非線性={v_nl:.6f} V")
        
        # 手動驗證 PCENT 計算
        print(f"\n手動驗證 PCENT 計算:")
        v1_l = potential_laplace[0, interface_idx]
        v2_l = potential_laplace[1, interface_idx]
        pcent_manual_l = (9.0 * v1_l - v2_l) / 8.0
        
        v1_nl = potential_nonlinear[0, interface_idx]
        v2_nl = potential_nonlinear[1, interface_idx]
        pcent_manual_nl = (9.0 * v1_nl - v2_nl) / 8.0
        
        print(f"  Laplace 手動:    (9*{v1_l:.4f} - {v2_l:.4f})/8 = {pcent_manual_l:.6f} V")
        print(f"  非線性手動:      (9*{v1_nl:.4f} - {v2_nl:.4f})/8 = {pcent_manual_nl:.6f} V")
        
    except Exception as e:
        print(f"❌ 非線性求解失敗: {e}")
        potential_nonlinear = None
        pot0_nonlinear = None
    
    return {
        'pot0_laplace': pot0_laplace,
        'pot0_nonlinear': pot0_nonlinear,
        'potential_laplace': potential_laplace,
        'potential_nonlinear': potential_nonlinear
    }

if __name__ == "__main__":
    results = test_nonlinear_pot0()
    
    print(f"\n" + "="*80)
    print(f"總結:")
    print(f"="*80)
    
    if results['pot0_nonlinear'] is not None:
        diff_l = abs(results['pot0_laplace'] - (-0.08))
        diff_nl = abs(results['pot0_nonlinear'] - (-0.08))
        
        print(f"Pot0 結果:")
        print(f"  Laplace:     {results['pot0_laplace']:.6f} V (差異: {diff_l:.6f} V)")
        print(f"  非線性:      {results['pot0_nonlinear']:.6f} V (差異: {diff_nl:.6f} V)")
        print(f"  Fortran:     -0.08 V")
        
        if diff_nl < diff_l:
            print(f"✅ 非線性求解改善了 Pot0 準確性")
        else:
            print(f"⚠️  非線性求解沒有改善 Pot0")
        
        if diff_nl < 0.1:
            print(f"🎉 Pot0 已接近 Fortran 結果 (差異 < 0.1V)")
        elif diff_nl < 0.5:
            print(f"👍 Pot0 合理接近 Fortran 結果 (差異 < 0.5V)")
        else:
            print(f"⚠️  Pot0 仍需進一步改善")
    else:
        print(f"❌ 非線性求解失敗，僅有 Laplace 結果")
        print(f"Laplace Pot0: {results['pot0_laplace']:.6f} V")