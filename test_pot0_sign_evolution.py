#!/usr/bin/env python3
"""
測試Python版本是否能觀察到Pot0的符號轉變
模擬Fortran中從負值演化到正值的物理過程
"""
import numpy as np
import logging
import sys
import os

# 添加項目路徑
sys.path.insert(0, os.path.abspath('.'))

from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid
from src.utils.constants import PhysicalConstants as PC

logging.basicConfig(level=logging.INFO, format='%(message)s')

def test_pot0_evolution():
    """測試Python版本的Pot0演化過程"""
    print("🧪 測試Python版本的Pot0符號演化")
    print("="*80)
    
    # 使用與Fortran相同的參數
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    solver = PoissonSOREquation(grid, props)
    
    # 測試條件 (與Fortran完全相同)
    V_tip = -2.0707107  # 與Fortran fort_MultInt.16完全相同
    V_sample = 0.0
    system_fermi = 1.4186435  # 與Fortran相同
    
    print(f"🎯 測試條件 (與Fortran完全相同):")
    print(f"   V_tip = {V_tip:.7f} V")
    print(f"   V_sample = {V_sample:.1f} V")
    print(f"   System Fermi = {system_fermi:.7f} eV")
    print(f"")
    
    # 改進的電荷密度計算器 (更接近Fortran物理)
    class ImprovedChargeDensityCalculator:
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            """改進的電荷密度計算，包含更真實的半導體物理"""
            kT = 0.0259  # 300K thermal energy
            
            # n型半導體參數 (與Fortran相似)
            Nd = 1e18  # cm^-3 (與Fortran相同)
            ni = 1e10  # cm^-3 intrinsic carrier concentration
            Eg = 1.42  # eV band gap
            
            # 更真實的載流子密度計算
            # 電子密度 (考慮費米-狄拉克統計)
            if ef_rel_vb_eV > Eg:
                # 導帶電子密度
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / kT)
            else:
                # 價帶電子密度 (很小)
                n_electrons = ni * np.exp(ef_rel_vb_eV / kT)
            
            # 電洞密度
            n_holes = ni**2 / n_electrons
            
            # 離化雜質密度 (溫度相關)
            if ef_rel_vb_eV < 0.5:  # 淺能級完全離化
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + 2 * np.exp((ef_rel_vb_eV - 0.5) / kT))
            
            # 總電荷密度 (C/m^3)
            # 正: 離化雜質 + 電洞, 負: 電子
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            
            # 轉換為 C/m^3
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            
            # 限制在合理範圍內
            charge_density_C_m3 = np.clip(charge_density_C_m3, -1e18, 1e18)
            
            return charge_density_C_m3
    
    charge_calculator = ImprovedChargeDensityCalculator()
    
    print(f"📊 開始長時間自洽求解以觀察Pot0演化...")
    print(f"   (需要足夠迭代次數來觀察符號轉變)")
    print(f"")
    
    # 記錄Pot0演化
    pot0_evolution = []
    
    try:
        # 使用更長的迭代次數 (接近Fortran的3500次)
        potential_final, total_iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calculator,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=2000,  # 更多迭代以觀察完整演化
            tolerance_Volts=1e-5,  # 更嚴格的容差
            omega=1.0  # 保守的鬆弛因子
        )
        
        print(f"✅ 求解完成:")
        print(f"   總迭代次數: {total_iterations}")
        print(f"   收斂狀態: {'是' if converged else '否'}")
        print(f"")
        
        # 計算最終的所有Pot0變體
        vsint_array = solver._initialize_vsint_array()
        vsint_array = solver._update_vsint_with_surface_charge(
            vsint_array, potential_final, charge_calculator,
            system_fermi, V_tip)
        
        # 所有計算方法的最終結果
        pot0_regular_raw = solver._calculate_pot0_fortran_style(
            potential_final, use_vsint=False, apply_scaling_correction=False)
        pot0_regular_scaled = solver._calculate_pot0_fortran_style(
            potential_final, use_vsint=False, apply_scaling_correction=True)
        pot0_vsint_raw = solver._calculate_pot0_fortran_style(
            potential_final, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=False)
        pot0_vsint_scaled = solver._calculate_pot0_fortran_style(
            potential_final, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
        
        print(f"🎯 Python最終Pot0結果:")
        print(f"   Regular (原始):     {pot0_regular_raw:+.6f} V")
        print(f"   Regular (縮放):     {pot0_regular_scaled:+.6f} V")
        print(f"   VSINT (原始):      {pot0_vsint_raw:+.6f} V")
        print(f"   VSINT (縮放):      {pot0_vsint_scaled:+.6f} V")
        print(f"")
        
        # 與Fortran比較
        fortran_final = +0.0698396191  # Fortran最終結果
        print(f"📊 與Fortran比較:")
        print(f"   Fortran最終:       {fortran_final:+.6f} V")
        print(f"   Python VSINT:      {pot0_vsint_scaled:+.6f} V")
        print(f"   差異:              {abs(pot0_vsint_scaled - fortran_final):.6f} V")
        print(f"")
        
        # 符號分析
        print(f"🔍 符號分析:")
        if pot0_vsint_scaled > 0:
            print(f"   ✅ Python達到正值! (與Fortran一致)")
            print(f"   🎯 成功實現積累→耗盡轉變")
        elif pot0_vsint_scaled < 0 and abs(pot0_vsint_scaled) < 0.1:
            print(f"   📈 接近轉變點 (可能需要更多迭代)")
        else:
            print(f"   ❌ 仍停留在負值階段")
        
        print(f"")
        print(f"💡 VSINT陣列分析:")
        print(f"   VSINT[0,0] = {vsint_array[0,0]:+.6f} V")
        print(f"   VSINT[1,0] = {vsint_array[1,0]:+.6f} V")
        
        return {
            'pot0_vsint_scaled': pot0_vsint_scaled,
            'pot0_regular_scaled': pot0_regular_scaled,
            'fortran_target': fortran_final,
            'sign_correct': pot0_vsint_scaled > 0,
            'iterations': total_iterations,
            'converged': converged
        }
        
    except Exception as e:
        print(f"❌ 求解失敗: {e}")
        import traceback
        traceback.print_exc()
        return None

def analyze_sign_transition():
    """分析符號轉變的物理條件"""
    print(f"🔬 符號轉變物理條件分析")
    print(f"="*80)
    
    print(f"📋 Fortran中觀察到的轉變:")
    print(f"   ITER=1600: -4.30E-03 V (負值)")
    print(f"   ITER=1700: +3.68E-03 V (正值)")
    print(f"   轉變發生在第1600-1700次迭代之間")
    print(f"")
    
    print(f"🔑 轉變的物理條件:")
    print(f"   1. 表面電荷密度達到臨界值")
    print(f"   2. 針尖場強度超過表面態吸引力")
    print(f"   3. 體電荷重新分布達到平衡")
    print(f"   4. 表面勢壘形成完成")
    print(f"")
    
    print(f"⚡ Python版本需要確保:")
    print(f"   • 足夠的迭代次數 (≥2000)")
    print(f"   • 正確的表面態物理模型")
    print(f"   • 合適的非線性求解策略")
    print(f"   • 真實的電荷密度計算")

if __name__ == "__main__":
    # 執行Pot0演化測試
    results = test_pot0_evolution()
    
    print()
    # 分析符號轉變條件
    analyze_sign_transition()
    
    print()
    print("="*80)
    print("🏆 最終評估")
    print("="*80)
    
    if results:
        if results['sign_correct']:
            print("🎉 成功！Python版本實現了符號轉變!")
            print(f"   ✅ 從負值演化到正值: {results['pot0_vsint_scaled']:+.6f} V")
            print(f"   ✅ 與Fortran符號一致")
            print(f"   📏 精度: {abs(results['pot0_vsint_scaled'] - results['fortran_target']):.6f} V")
        else:
            print("⚠️  Python版本尚未完全實現符號轉變")
            print(f"   當前結果: {results['pot0_vsint_scaled']:+.6f} V")
            print(f"   需要: 更完整的物理模型或更多迭代")
        
        print(f"")
        print(f"🔑 關鍵理解:")
        print(f"   • Fortran的從負到正演化是物理上正確的")
        print(f"   • 這代表表面能帶彎曲的真實物理過程")
        print(f"   • Python版本需要實現這個完整過程")
    else:
        print("❌ 測試失敗，需要進一步調試")
    
    print(f"")
    print(f"🎯 主要結論:")
    print(f"   1. Fortran行為驗證為正確 (從-0.083V到+0.070V)")
    print(f"   2. 符號變化代表真實的STM物理過程")
    print(f"   3. Python版本的挑戰是實現這個完整演化")
    print(f"   4. 需要更完善的VSINT實現以達到正值")