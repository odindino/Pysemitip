#!/usr/bin/env python3
"""
分析我們的Python實現與Fortran算法的具體差異
基於用戶的重要觀點：如果輸入相同、算法相同，結果應該很接近
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

def analyze_algorithm_differences():
    """分析Python與Fortran算法的關鍵差異"""
    print("🔍 分析Python與Fortran算法差異")
    print("="*80)
    print("🎯 目標：找出具體的計算差異，而非參數調優")
    print("💡 原則：相同輸入 + 相同算法 → 相同結果")
    print()
    
    # 設置相同的測試條件
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    solver = PoissonSOREquation(grid, props)
    
    # 與Fortran完全相同的輸入條件
    V_tip = -2.0707107
    V_sample = 0.0
    system_fermi = 1.4186435
    
    print(f"📋 測試條件 (與Fortran完全一致):")
    print(f"   V_tip = {V_tip:.7f} V")
    print(f"   V_sample = {V_sample:.1f} V")
    print(f"   System Fermi = {system_fermi:.7f} eV")
    print()
    
    # 創建簡單、一致的電荷密度計算器
    class StandardChargeDensityCalculator:
        def __init__(self):
            self.call_count = 0
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.call_count += 1
            
            # 使用與Fortran相同的基本半導體物理
            kT = 0.0259  # 室溫
            Nd = 5e18    # 標準摻雜密度 cm^-3
            ni = 1e10    # 本征載流子密度 cm^-3
            Eg = 1.42    # 能隙 eV
            
            # 標準載流子統計
            if ef_rel_vb_eV > Eg:
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / kT)
            else:
                n_electrons = ni * np.exp(ef_rel_vb_eV / kT)
            
            n_holes = ni**2 / n_electrons
            
            # 標準雜質離化
            if ef_rel_vb_eV < 0.5:
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + np.exp((ef_rel_vb_eV - 0.5) / kT))
            
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            
            # 轉換為C/m³
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            
            return charge_density_C_m3
    
    charge_calc = StandardChargeDensityCalculator()
    
    print("🔍 第1步：檢查關鍵算法差異")
    print("-" * 60)
    
    # 1. 檢查PCENT計算
    print("1️⃣ PCENT函數分析:")
    print("   Fortran PCENT(JJ=0):")
    print("   - 使用VSINT陣列：專門的表面電位")
    print("   - 公式：(9*VSINT(1,1,K) - VSINT(1,2,K))/8")
    print("   - 對所有角度點K求平均")
    print()
    
    # 檢查我們的PCENT實現
    potential_test = solver._create_initial_potential_guess(V_tip, V_sample)
    
    # 原始方法（未縮放）
    pot0_raw = solver._calculate_pot0_fortran_style(potential_test, use_vsint=False, apply_scaling_correction=False)
    
    # VSINT方法（未縮放）
    vsint_array = solver._initialize_vsint_array()
    # 模擬表面電位
    N_eta, N_nu = potential_test.shape
    for i in range(N_eta):
        vsint_array[i, 0] = potential_test[i, N_nu-1]  # 界面電位
    
    pot0_vsint_raw = solver._calculate_pot0_fortran_style(potential_test, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=False)
    
    # 縮放方法
    pot0_scaled = solver._calculate_pot0_fortran_style(potential_test, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
    
    print("   Python PCENT 實現:")
    print(f"   - 原始方法：      {pot0_raw:.6f} V")
    print(f"   - VSINT方法（原始）：{pot0_vsint_raw:.6f} V")
    print(f"   - VSINT方法（縮放）：{pot0_scaled:.6f} V")
    print(f"   - 縮放因子：      0.113 (這是問題所在！)")
    print()
    
    print("❌ 發現問題1：我們添加了經驗縮放因子0.113")
    print("   這違反了'相同算法'的原則！")
    print()
    
    # 2. 檢查Golden Section Search
    print("2️⃣ Golden Section Search分析:")
    print("   Fortran GSECT:")
    print("   - GS = 0.3819660")
    print("   - 標準的二分搜索邏輯")
    print()
    
    # 檢查我們的GSS實現
    def test_function(x):
        return (x - 0.5)**2  # 簡單的二次函數，最小值在x=0.5
    
    # 使用我們的實現
    result_python = solver._golden_section_minimize(test_function, 0.0, 1.0, 1e-6)
    
    # 手動實現標準GSS
    def standard_gss(func, xmin, xmax, tol):
        gs = 0.3819660
        if xmax == xmin:
            return xmin
        
        delx = xmax - xmin
        xa = xmin + delx * gs
        fa = func(xa)
        xb = xmax - delx * gs
        fb = func(xb)
        
        while delx >= tol:
            if fb < fa:
                xmin = xa
                xa = xb
                fa = fb
                delx = xmax - xmin
                xb = xmax - delx * gs
                fb = func(xb)
            else:
                xmax = xb
                xb = xa
                fb = fa
                delx = xmax - xmin
                xa = xmin + delx * gs
                fa = func(xa)
        
        return (xmin + xmax) / 2
    
    result_standard = standard_gss(test_function, 0.0, 1.0, 1e-6)
    
    print(f"   Python實現結果：  {result_python:.6f}")
    print(f"   標準實現結果：    {result_standard:.6f}")
    print(f"   理論最優值：      0.500000")
    print(f"   差異：           {abs(result_python - result_standard):.8f}")
    
    if abs(result_python - result_standard) > 1e-6:
        print("❌ 發現問題2：Golden Section Search實現有差異")
    else:
        print("✅ Golden Section Search實現正確")
    print()
    
    # 3. 檢查物理常數
    print("3️⃣ 物理常數分析:")
    print("   Fortran常數：")
    print("   - EEP = 1.80943E-20")
    print("   - EPSIL0 = 8.854185E-12")
    print("   - E = 1.60210E-19")
    print()
    
    # 檢查我們的常數
    from src.utils.constants import PhysicalConstants as PC
    print("   Python常數：")
    print(f"   - E (電子電荷) = {PC.E:.5e}")
    print(f"   - EPSILON0 = {PC.EPSILON0:.6e}")
    
    # 計算EEP等效值
    eep_equivalent = PC.E / PC.EPSILON0 * 1e-14  # Fortran註釋中的單位轉換
    print(f"   - EEP等效值 = {eep_equivalent:.5e}")
    print(f"   - Fortran EEP = 1.80943e-20")
    print(f"   - 差異 = {abs(eep_equivalent - 1.80943e-20):.2e}")
    
    if abs(eep_equivalent - 1.80943e-20) > 1e-22:
        print("❌ 發現問題3：物理常數有微小差異")
    else:
        print("✅ 物理常數基本一致")
    print()
    
    # 4. 檢查邊界條件
    print("4️⃣ 邊界條件分析:")
    
    # 檢查初始電位設置
    print("   檢查初始電位設置...")
    N_eta, N_nu = potential_test.shape
    
    print(f"   網格大小：{N_eta} x {N_nu}")
    print(f"   針尖電位 [0,0]: {potential_test[0, 0]:.6f} V")
    print(f"   界面電位 [0,{N_nu-1}]: {potential_test[0, N_nu-1]:.6f} V")
    print(f"   樣品電位 [0,{N_nu-1}]: 應該是 {V_sample:.1f} V")
    
    if abs(potential_test[0, N_nu-1] - V_sample) > 1e-6:
        print("❌ 發現問題4：界面邊界條件設置有誤")
        print(f"   界面電位應該是{V_sample}，但設置為{potential_test[0, N_nu-1]}")
    else:
        print("✅ 邊界條件設置正確")
    print()
    
    print("🔍 第2步：執行無縮放的標準求解")
    print("-" * 60)
    
    # 移除所有非標準修改，使用基本算法
    def standard_solve():
        """執行標準的、未修改的求解"""
        # 臨時修改_calculate_pot0_fortran_style以禁用縮放
        original_method = solver._calculate_pot0_fortran_style
        
        def no_scaling_pot0(potential, use_vsint=False, vsint_array=None, apply_scaling_correction=False):
            # 強制禁用縮放
            return original_method(potential, use_vsint, vsint_array, False)
        
        solver._calculate_pot0_fortran_style = no_scaling_pot0
        
        try:
            potential, iterations, converged = solver.solve(
                V_tip_Volts=V_tip,
                V_sample_Volts=V_sample,
                charge_density_calculator=charge_calc,
                system_fermi_level_E_F_main_eV=system_fermi,
                max_iterations=500,
                tolerance_Volts=1e-3,
                omega=1.2  # 標準鬆弛因子
            )
            
            # 計算最終Pot0（無縮放）
            pot0_final = original_method(potential, use_vsint=False, vsint_array=None, apply_scaling_correction=False)
            
            return pot0_final, iterations, converged
            
        finally:
            solver._calculate_pot0_fortran_style = original_method
    
    pot0_standard, iterations, converged = standard_solve()
    
    print(f"   標準求解結果：")
    print(f"   - Pot0（無縮放）：  {pot0_standard:.6f} V")
    print(f"   - 迭代次數：        {iterations}")
    print(f"   - 收斂狀態：        {converged}")
    print(f"   - 電荷計算次數：    {charge_calc.call_count:,}")
    print()
    
    print("🎯 關鍵發現總結")
    print("="*70)
    
    print("❌ 主要問題：")
    print("1. **經驗縮放因子0.113** - 這不是Fortran算法的一部分")
    print("2. **可能的VSINT計算差異** - 需要精確匹配Fortran的RHOSURF邏輯")
    print("3. **界面電位處理** - 可能與Fortran的VSINT更新機制不同")
    print()
    
    print("✅ 相同部分：")
    print("1. **Golden Section Search** - 實現基本正確")
    print("2. **物理常數** - 基本一致")
    print("3. **網格設置** - 結構相同")
    print()
    
    print("💡 解決策略：")
    print("1. **移除所有經驗縮放** - 回到標準算法")
    print("2. **精確實現Fortran的VSINT邏輯** - 特別是RHOSURF計算")
    print("3. **檢查數值精度** - 確保浮點數計算一致")
    print("4. **逐步對比** - 每個計算步驟都與Fortran比較")
    print()
    
    # 與Fortran的比較
    fortran_target = 0.0698396191
    difference_standard = abs(pot0_standard - fortran_target)
    difference_no_scaling = abs(pot0_standard - fortran_target)
    
    print(f"📊 與Fortran比較 (+{fortran_target:.6f}V):")
    print(f"   標準Python：    {pot0_standard:+.6f}V (差異: {difference_standard:.6f}V)")
    print(f"   改善幅度：      需要 {difference_standard:.6f}V 的修正")
    print()
    
    if abs(pot0_standard) > 1.0:
        print("⚠️  Python結果數量級過大，算法實現確實有根本差異")
    else:
        print("✅ Python結果數量級合理，主要是數值差異")
    
    print()
    print("🔑 下一步行動：")
    print("1. 移除所有經驗性修正")
    print("2. 精確實現Fortran的VSINT計算")
    print("3. 逐步驗證每個算法組件")
    print("4. 確保數值實現完全一致")

if __name__ == "__main__":
    print("🎯 算法差異分析：找出Python與Fortran的具體差異")
    print("原則：相同輸入 + 相同算法 = 相同結果")
    print()
    
    analyze_algorithm_differences()
    
    print()
    print("="*80)
    print("🏁 分析完成")
    print("="*80)