#!/usr/bin/env python3
"""
深入調試PCENT計算的差異
重點分析為什麼求解過程中顯示-0.061V，但最終計算得到-0.542V
"""
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

logging.basicConfig(level=logging.INFO, format='%(message)s')

def debug_pcent_calculation():
    """調試PCENT計算的具體差異"""
    print("🔍 深入調試PCENT計算")
    print("="*80)
    print("🎯 問題：求解過程顯示-0.061V，但最終計算-0.542V")
    print("💡 分析：PCENT函數實現與Fortran的細微差異")
    print()
    
    # 設置測試環境
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
    
    print(f"📋 測試條件:")
    print(f"   V_tip = {V_tip:.7f} V")
    print(f"   V_sample = {V_sample:.1f} V")
    print(f"   System Fermi = {system_fermi:.7f} eV")
    print()
    
    # 創建標準電荷計算器
    class StandardChargeDensityCalculator:
        def __init__(self):
            self.call_count = 0
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.call_count += 1
            
            kT = 0.0259
            Nd = 5e18
            ni = 1e10
            Eg = 1.42
            
            if ef_rel_vb_eV > Eg:
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / kT)
            else:
                n_electrons = ni * np.exp(ef_rel_vb_eV / kT)
            
            n_holes = ni**2 / n_electrons
            
            if ef_rel_vb_eV < 0.5:
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + np.exp((ef_rel_vb_eV - 0.5) / kT))
            
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            charge_density_C_m3 = charge_density_cm3 * 1e6 * 1.60210e-19  # 使用Fortran常數
            
            return charge_density_C_m3
    
    charge_calc = StandardChargeDensityCalculator()
    
    # 正確的邊界條件
    def correct_initial_guess(V_tip, V_sample):
        N_eta, N_nu = grid.N_eta, grid.N_nu
        potential = np.zeros((N_eta, N_nu))
        
        for i in range(N_eta):
            for j in range(N_nu):
                nu_fraction = j / max(N_nu - 1, 1)
                potential[i, j] = V_tip * (1 - nu_fraction) + V_sample * nu_fraction
        
        # 確保界面為V_sample
        for i in range(N_eta):
            potential[i, N_nu - 1] = V_sample
        
        return potential
    
    # 手動計算PCENT的不同方法
    def analyze_pcent_methods(potential, verbose=True):
        """分析不同PCENT計算方法的差異"""
        N_eta, N_nu = potential.shape
        interface_nu_idx = N_nu - 1
        
        if verbose:
            print(f"🔍 PCENT計算分析:")
            print(f"   網格大小: {N_eta} x {N_nu}")
            print(f"   界面索引: nu = {interface_nu_idx}")
            print()
        
        # 方法1：Fortran PCENT邏輯 (JJ=0, 使用界面電位)
        I = 0  # Fortran I=1 對應 Python I=0
        
        if I + 1 < N_eta:
            v1 = potential[I, interface_nu_idx]      # 界面電位 [0, 7]
            v2 = potential[I + 1, interface_nu_idx]  # 下一個徑向點 [1, 7]
            pcent_fortran_formula = (9.0 * v1 - v2) / 8.0
        else:
            pcent_fortran_formula = potential[0, interface_nu_idx]
        
        # 方法2：簡單界面電位
        pcent_simple_interface = potential[0, interface_nu_idx]
        
        # 方法3：中心軸電位 (如果有的話)
        pcent_center_axis = potential[0, 0]  # 針尖中心
        
        # 方法4：平均界面電位
        pcent_average_interface = np.mean(potential[:, interface_nu_idx])
        
        if verbose:
            print(f"   界面電位值:")
            for i in range(min(4, N_eta)):
                print(f"     [i={i}, nu={interface_nu_idx}]: {potential[i, interface_nu_idx]:.6f} V")
            print()
            
            print(f"   PCENT計算方法比較:")
            print(f"     Fortran公式 (9*v1-v2)/8: {pcent_fortran_formula:.6f} V")
            print(f"     簡單界面電位:           {pcent_simple_interface:.6f} V")
            print(f"     中心軸電位:             {pcent_center_axis:.6f} V")
            print(f"     平均界面電位:           {pcent_average_interface:.6f} V")
            print()
        
        return {
            'fortran_formula': pcent_fortran_formula,
            'simple_interface': pcent_simple_interface,
            'center_axis': pcent_center_axis,
            'average_interface': pcent_average_interface
        }
    
    print("🔍 第1步：分析初始電位的PCENT計算")
    print("-" * 60)
    
    # 分析初始電位
    potential_init = correct_initial_guess(V_tip, V_sample)
    pcent_init = analyze_pcent_methods(potential_init, verbose=True)
    
    print("🔍 第2步：執行求解並追蹤PCENT變化")
    print("-" * 60)
    
    # 替換solver的初始猜測
    original_initial_guess = solver._create_initial_potential_guess
    solver._create_initial_potential_guess = correct_initial_guess
    
    # 修改solver的PCENT計算以進行調試
    original_pot0_calc = solver._calculate_pot0_fortran_style
    
    def debug_pot0_calculation(potential, use_vsint=False, vsint_array=None, apply_scaling_correction=False):
        """調試版本的Pot0計算"""
        pcent_results = analyze_pcent_methods(potential, verbose=False)
        
        # 使用Fortran公式（無縮放）
        result = pcent_results['fortran_formula']
        
        # 記錄調試信息
        if hasattr(debug_pot0_calculation, 'call_count'):
            debug_pot0_calculation.call_count += 1
        else:
            debug_pot0_calculation.call_count = 1
        
        if debug_pot0_calculation.call_count % 100 == 0:  # 每100次迭代記錄一次
            print(f"     調試 (迭代 {debug_pot0_calculation.call_count}): Fortran公式={result:.6f}V, 簡單界面={pcent_results['simple_interface']:.6f}V")
        
        return result
    
    solver._calculate_pot0_fortran_style = debug_pot0_calculation
    
    try:
        print("🚀 開始求解 (將顯示調試信息)...")
        
        potential, iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calc,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=500,
            tolerance_Volts=1e-3,
            omega=1.2
        )
        
        print(f"✅ 求解完成:")
        print(f"   迭代次數: {iterations}")
        print(f"   收斂狀態: {converged}")
        print()
        
        print("🔍 第3步：詳細分析最終結果")
        print("-" * 60)
        
        # 分析最終電位
        pcent_final = analyze_pcent_methods(potential, verbose=True)
        
        # 檢查電位分布
        N_eta, N_nu = potential.shape
        print(f"   最終電位分布檢查:")
        print(f"     針尖區域 [0,0]: {potential[0, 0]:.6f} V")
        print(f"     中間區域 [0,{N_nu//2}]: {potential[0, N_nu//2]:.6f} V")
        print(f"     界面區域 [0,{N_nu-1}]: {potential[0, N_nu-1]:.6f} V")
        print()
        
        print(f"     徑向分布 (界面處):")
        for i in range(min(4, N_eta)):
            print(f"       [i={i}, nu={N_nu-1}]: {potential[i, N_nu-1]:.6f} V")
        print()
        
        # 與求解過程中顯示的值比較
        print(f"🔍 關鍵發現:")
        print(f"   求解過程顯示的值: 約 -0.061V")
        print(f"   Fortran公式計算:    {pcent_final['fortran_formula']:+.6f} V")
        print(f"   簡單界面電位:       {pcent_final['simple_interface']:+.6f} V")
        print()
        
        # 檢查是否有數值問題
        v1 = potential[0, N_nu-1]
        v2 = potential[1, N_nu-1] if N_eta > 1 else v1
        
        print(f"   Fortran公式詳細計算:")
        print(f"     v1 (界面電位) = {v1:.6f} V")
        print(f"     v2 (下一點)   = {v2:.6f} V")
        print(f"     9*v1 = {9*v1:.6f} V")
        print(f"     9*v1 - v2 = {9*v1 - v2:.6f} V")
        print(f"     (9*v1 - v2)/8 = {(9*v1 - v2)/8:.6f} V")
        print()
        
        # 檢查是否求解過程和最終計算使用了不同的數據
        print(f"🤔 可能的問題:")
        if abs(pcent_final['fortran_formula'] - (-0.061)) > 0.1:
            print(f"   1. 求解過程中的Pot0計算與最終PCENT計算使用了不同邏輯")
            print(f"   2. 可能存在數據更新時序問題")
            print(f"   3. VSINT陣列與電位矩陣不同步")
        else:
            print(f"   ✅ 求解過程與最終計算基本一致")
        
        # 與Fortran目標比較
        fortran_target = 0.0698396191
        best_match = min(pcent_final.values(), key=lambda x: abs(x - fortran_target))
        best_method = [k for k, v in pcent_final.items() if v == best_match][0]
        
        print(f"   最接近Fortran的方法: {best_method}")
        print(f"   差異: {abs(best_match - fortran_target):.6f} V")
        
    except Exception as e:
        print(f"❌ 求解失敗: {e}")
        
    finally:
        # 恢復原始方法
        solver._create_initial_potential_guess = original_initial_guess
        solver._calculate_pot0_fortran_style = original_pot0_calc
    
    print()
    print("🎯 PCENT調試總結")
    print("="*70)
    print("🔍 關鍵發現:")
    print("1. 需要確認求解過程中的Pot0計算邏輯")
    print("2. 檢查VSINT陣列的更新機制")
    print("3. 可能需要完全重新實現PCENT函數")
    print("4. 確保電位矩陣的數據一致性")

if __name__ == "__main__":
    print("🎯 深入調試PCENT計算差異")
    print("分析為什麼求解過程和最終計算結果不同")
    print()
    
    debug_pcent_calculation()
    
    print()
    print("="*80)
    print("🏁 PCENT調試完成")
    print("="*80)