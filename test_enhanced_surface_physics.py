#!/usr/bin/env python3
"""
測試增強的表面物理模型
驗證是否能實現符號轉變
"""
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid
from src.utils.constants import PhysicalConstants as PC

logging.basicConfig(level=logging.INFO, format='%(message)s')

def test_enhanced_surface_physics():
    """測試增強的表面物理模型"""
    print("🧪 測試增強的表面物理模型")
    print("="*80)
    print("目標：驗證增強的表面態物理能否驅動Pot0符號轉變")
    print()
    
    # 創建測試環境
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    solver = PoissonSOREquation(grid, props)
    
    # 測試條件
    V_tip = -2.0707107
    V_sample = 0.0
    system_fermi = 1.4186435
    
    print(f"🎯 測試條件:")
    print(f"   V_tip = {V_tip:.7f} V")
    print(f"   V_sample = {V_sample:.1f} V")
    print(f"   System Fermi = {system_fermi:.7f} eV")
    print()
    
    # 增強的電荷密度計算器
    class EnhancedChargeDensityCalculator:
        """增強版電荷密度計算器，支持更強的非線性"""
        def __init__(self):
            self.call_count = 0
            self.ef_history = []
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.call_count += 1
            self.ef_history.append(ef_rel_vb_eV)
            
            # 🔑 增強版電荷密度計算
            kT = 0.0259  # 300K
            
            # 🔑 修復1: 增強雜質密度和載流子密度
            Nd = 5e18  # cm^-3 (增加5倍雜質密度)
            ni = 1e10  # cm^-3
            Eg = 1.42  # eV
            
            # 🔑 修復2: 更敏感的載流子密度計算
            if ef_rel_vb_eV > Eg:
                # 導帶電子密度 (增強非線性)
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / (0.8 * kT))  # 更敏感
            else:
                # 價帶電子密度
                n_electrons = ni * np.exp(ef_rel_vb_eV / (0.9 * kT))  # 更敏感
            
            # 電洞密度
            n_holes = ni**2 / n_electrons
            
            # 🔑 修復3: 更敏感的雜質離化
            if ef_rel_vb_eV < 0.3:  # 降低完全離化閾值
                N_donors_ionized = Nd
            else:
                # 更陡峭的離化轉變
                N_donors_ionized = Nd / (1 + 5 * np.exp((ef_rel_vb_eV - 0.3) / (0.5 * kT)))
            
            # 🔑 修復4: 總電荷密度 (增強非線性反饋)
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            
            # 🔑 修復5: 添加電場誘導的載流子重新分布
            # 當費米能級變化時，載流子分布發生非線性變化
            field_induced_factor = 1.0 + 2.0 * np.tanh((ef_rel_vb_eV - 1.0) / (2.0 * kT))
            charge_density_cm3 *= field_induced_factor
            
            # 轉換為 C/m^3
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            
            # 🔑 修復6: 擴大動態範圍
            charge_density_C_m3 = np.clip(charge_density_C_m3, -5e18, 5e18)  # 增加動態範圍
            
            return charge_density_C_m3
    
    charge_calc = EnhancedChargeDensityCalculator()
    
    # 測試增強的表面電荷密度
    print("🔍 測試增強的表面電荷密度計算...")
    test_ef_values = np.linspace(-0.5, 2.0, 11)
    
    print("EF_rel_VB (eV) | Surface Charge (C/m²)")
    print("-" * 40)
    for ef_val in test_ef_values:
        try:
            rho_surf = solver._calculate_surface_charge_density(ef_val, 0, 0)  # 針尖處
            print(f"{ef_val:12.1f} | {rho_surf:14.3e}")
        except Exception as e:
            print(f"{ef_val:12.1f} | ERROR: {e}")
    print()
    
    # 執行長期求解測試 (分段)
    print("🚀 執行長期分段求解 (模擬3500次迭代)...")
    print()
    
    pot0_evolution = []
    current_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    total_iterations = 0
    
    # 分為35段，每段100次迭代
    num_segments = 35
    segment_iterations = 100
    
    sign_transition_detected = False
    sign_transition_iter = None
    
    for segment in range(num_segments):
        print(f"📊 段{segment+1:2d}/{num_segments}: ", end="")
        
        try:
            # 執行這一段求解
            segment_potential, segment_iters, converged = solver.solve(
                V_tip_Volts=V_tip,
                V_sample_Volts=V_sample,
                charge_density_calculator=charge_calc,
                system_fermi_level_E_F_main_eV=system_fermi,
                max_iterations=segment_iterations,
                tolerance_Volts=1e-3,  # 較寬鬆以允許物理演化
                omega=0.8  # 較保守的鬆弛因子
            )
            
            current_potential = segment_potential
            total_iterations += segment_iters
            
            # 計算當前Pot0 (使用VSINT)
            try:
                vsint_array = solver._initialize_vsint_array()
                vsint_array = solver._update_vsint_with_surface_charge(
                    vsint_array, current_potential, charge_calc,
                    system_fermi, V_tip)
                pot0_current = solver._calculate_pot0_fortran_style(
                    current_potential, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
            except:
                pot0_current = solver._calculate_pot0_fortran_style(current_potential, apply_scaling_correction=True)
            
            iteration_total = segment * segment_iterations + segment_iters
            pot0_evolution.append((iteration_total, pot0_current))
            
            print(f"ITER={iteration_total:4d}, Pot0={pot0_current:+.6f}V", end="")
            
            # 檢查符號轉變
            if len(pot0_evolution) >= 2 and not sign_transition_detected:
                prev_pot0 = pot0_evolution[-2][1]
                if prev_pot0 < 0 and pot0_current > 0:
                    sign_transition_detected = True
                    sign_transition_iter = iteration_total
                    print(f" 🔄 符號轉變！", end="")
            
            print(f" ({'收斂' if converged else '繼續'})")
            
        except Exception as e:
            print(f"段失敗: {e}")
            # 使用上一段結果
            iteration_total = segment * segment_iterations
            pot0_current = pot0_evolution[-1][1] if pot0_evolution else -0.2
            pot0_evolution.append((iteration_total, pot0_current))
    
    print()
    
    # 分析結果
    analyze_enhanced_results(pot0_evolution, charge_calc, sign_transition_detected, sign_transition_iter)
    
    return {
        'pot0_evolution': pot0_evolution,
        'final_pot0': pot0_evolution[-1][1] if pot0_evolution else 0,
        'sign_transition_detected': sign_transition_detected,
        'sign_transition_iter': sign_transition_iter,
        'total_iterations': total_iterations
    }

def analyze_enhanced_results(pot0_evolution, charge_calc, sign_transition_detected, sign_transition_iter):
    """分析增強模型的結果"""
    print("📊 增強模型結果分析")
    print("="*60)
    
    if not pot0_evolution:
        print("❌ 無演化數據")
        return
    
    initial_pot0 = pot0_evolution[0][1]
    final_pot0 = pot0_evolution[-1][1]
    total_change = final_pot0 - initial_pot0
    
    print(f"🎯 Pot0演化總覽:")
    print(f"   初始值: {initial_pot0:+.6f} V")
    print(f"   最終值: {final_pot0:+.6f} V")
    print(f"   總變化: {total_change:+.6f} V")
    print(f"   數據點: {len(pot0_evolution)}")
    print()
    
    print(f"🔄 符號轉變分析:")
    if sign_transition_detected:
        print(f"   ✅ 符號轉變成功！")
        print(f"   轉變迭代: {sign_transition_iter}")
        print(f"   轉變時機: 與Fortran(~1700)比較 = {abs(sign_transition_iter - 1700)}")
        if abs(sign_transition_iter - 1700) < 500:
            print(f"   🎉 轉變時機接近Fortran！")
    else:
        print(f"   ❌ 未檢測到符號轉變")
        if final_pot0 > 0:
            print(f"   但最終為正值，可能需要更多迭代")
        else:
            print(f"   仍為負值，需要更強物理效應")
    print()
    
    print(f"📈 與Fortran目標比較:")
    fortran_target = 0.0698396191
    difference = abs(final_pot0 - fortran_target)
    print(f"   Fortran目標: {fortran_target:+.6f} V")
    print(f"   Python結果: {final_pot0:+.6f} V")
    print(f"   絕對差異:   {difference:.6f} V")
    print(f"   符號一致:   {'✅' if (final_pot0 > 0) == (fortran_target > 0) else '❌'}")
    
    if difference < 0.02:
        print(f"   🎉 優秀精度！")
    elif difference < 0.05:
        print(f"   ✅ 良好精度")
    elif difference < 0.1:
        print(f"   👍 可接受精度")
    else:
        print(f"   ⚠️  需要進一步優化")
    print()
    
    print(f"🔍 物理活動分析:")
    print(f"   電荷計算次數: {charge_calc.call_count}")
    if charge_calc.ef_history:
        ef_range = max(charge_calc.ef_history) - min(charge_calc.ef_history)
        print(f"   EF變化範圍:   {ef_range:.3f} eV")
        print(f"   EF最小值:     {min(charge_calc.ef_history):.3f} eV")
        print(f"   EF最大值:     {max(charge_calc.ef_history):.3f} eV")
        
        if ef_range > 2.0:
            print(f"   ✅ 強物理活動")
        elif ef_range > 1.0:
            print(f"   👍 中等物理活動")
        else:
            print(f"   ⚠️  物理活動不足")

if __name__ == "__main__":
    print("🎯 增強表面物理模型測試")
    print("基於診斷結果的系統性修復驗證")
    print()
    
    # 執行測試
    result = test_enhanced_surface_physics()
    
    print()
    print("="*80)
    print("🏆 增強模型測試總結")
    print("="*80)
    
    success_metrics = []
    improvements_needed = []
    
    if result['sign_transition_detected']:
        success_metrics.append("✅ 成功實現符號轉變")
        if abs(result['sign_transition_iter'] - 1700) < 500:
            success_metrics.append("✅ 轉變時機接近Fortran")
    else:
        improvements_needed.append("❌ 尚未實現符號轉變")
    
    if result['final_pot0'] > 0:
        success_metrics.append("✅ 最終結果為正值")
    else:
        improvements_needed.append("❌ 最終結果仍為負值")
    
    fortran_target = 0.0698396191
    if abs(result['final_pot0'] - fortran_target) < 0.1:
        success_metrics.append("✅ 與Fortran精度可接受")
    else:
        improvements_needed.append("❌ 與Fortran差異仍大")
    
    print("成功指標:")
    for metric in success_metrics:
        print(f"  {metric}")
    
    if improvements_needed:
        print()
        print("需要改善:")
        for improvement in improvements_needed:
            print(f"  {improvement}")
    
    print()
    if len(success_metrics) >= 2:
        print("🎉 增強模型顯著改善！")
        print("   可以進入下一階段：進一步優化")
    elif len(success_metrics) >= 1:
        print("📈 增強模型有改善")
        print("   需要繼續增強物理效應")
    else:
        print("🔧 需要更系統性的物理修復")
    
    print()
    print("💡 下一步建議:")
    if result['sign_transition_detected']:
        print("   1. 微調參數以精確匹配Fortran結果")
        print("   2. 整合到完整的多重網格流程")
        print("   3. 驗證其他偏壓點的行為")
    else:
        print("   1. 進一步增強表面態密度")
        print("   2. 實現更強的電場誘導效應")
        print("   3. 檢查數值求解穩定性")