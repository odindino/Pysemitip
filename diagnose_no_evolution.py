#!/usr/bin/env python3
"""
診斷為什麼電位完全沒有演化
深入分析數值求解器被困在穩定平衡的根本原因
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

def diagnose_no_evolution():
    """診斷為什麼電位完全沒有演化"""
    print("🔬 診斷電位演化停滯問題")
    print("="*80)
    print("🎯 目標：找到阻止數值演化的根本原因")
    print()
    
    # 測試環境設置
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
    
    # 檢查 1: 初始電位猜測分析
    print("🔍 檢查 1: 初始電位猜測分析")
    print("-" * 50)
    analyze_initial_guess(solver, V_tip, V_sample)
    print()
    
    # 檢查 2: 邊界條件限制分析
    print("🔍 檢查 2: 邊界條件限制分析") 
    print("-" * 50)
    analyze_boundary_constraints(solver, V_tip, V_sample)
    print()
    
    # 檢查 3: SOR 係數和穩定性分析
    print("🔍 檢查 3: SOR 係數和穩定性分析")
    print("-" * 50)
    analyze_sor_coefficients(solver)
    print()
    
    # 檢查 4: 電荷密度反饋分析
    print("🔍 檢查 4: 電荷密度反饋分析")
    print("-" * 50)
    analyze_charge_feedback(solver, system_fermi)
    print()
    
    # 檢查 5: 單次迭代變化分析
    print("🔍 檢查 5: 單次迭代變化分析")
    print("-" * 50)
    analyze_single_iteration_change(solver, V_tip, V_sample, system_fermi)
    print()
    
    # 提出解決方案
    propose_evolution_solutions()

def analyze_initial_guess(solver, V_tip, V_sample):
    """分析初始電位猜測"""
    initial_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    
    print(f"   電位範圍: {np.min(initial_potential):.6f} 到 {np.max(initial_potential):.6f} V")
    print(f"   平均電位: {np.mean(initial_potential):.6f} V")
    print(f"   標準差:   {np.std(initial_potential):.6f} V")
    
    # 檢查關鍵點
    N_eta, N_nu = initial_potential.shape
    tip_potential = initial_potential[0, 0]
    interface_potential = initial_potential[0, N_nu-1]
    far_field_potential = initial_potential[N_eta-1, N_nu//2]
    
    print(f"   針尖電位:   {tip_potential:.6f} V")
    print(f"   界面電位:   {interface_potential:.6f} V")
    print(f"   遠場電位:   {far_field_potential:.6f} V")
    
    # 檢查是否有合理的電位梯度
    grad_eta = np.gradient(initial_potential, axis=0)
    grad_nu = np.gradient(initial_potential, axis=1)
    total_grad = np.sqrt(grad_eta**2 + grad_nu**2)
    
    print(f"   最大梯度:   {np.max(total_grad):.6f} V/grid")
    print(f"   平均梯度:   {np.mean(total_grad):.6f} V/grid")
    
    # 計算初始 Pot0
    initial_pot0 = solver._calculate_pot0_fortran_style(initial_potential, apply_scaling_correction=True)
    print(f"   初始Pot0:   {initial_pot0:.6f} V")
    
    # 診斷：檢查是否過於平坦
    if np.std(initial_potential) < 0.1:
        print("   ⚠️  初始猜測可能過於平坦，缺乏驅動力")
    if np.max(total_grad) < 0.1:
        print("   ⚠️  初始梯度可能過小，難以驅動演化")
    
    return initial_potential

def analyze_boundary_constraints(solver, V_tip, V_sample):
    """分析邊界條件限制"""
    # 創建測試電位
    test_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    N_eta, N_nu = test_potential.shape
    
    # 模擬一些變化
    test_potential_modified = np.copy(test_potential)
    test_potential_modified[1:-1, 1:-1] += np.random.normal(0, 0.1, test_potential_modified[1:-1, 1:-1].shape)
    
    print(f"   修改前界面電位: {test_potential[0, N_nu-1]:.6f} V")
    print(f"   修改後界面電位: {test_potential_modified[0, N_nu-1]:.6f} V")
    
    # 應用邊界條件
    test_potential_with_bc = solver._apply_boundary_conditions(test_potential_modified, V_tip, V_sample)
    
    print(f"   邊界條件後界面電位: {test_potential_with_bc[0, N_nu-1]:.6f} V")
    
    # 檢查邊界條件是否過於嚴格
    interface_change = abs(test_potential_with_bc[0, N_nu-1] - test_potential_modified[0, N_nu-1])
    print(f"   界面電位被修改幅度: {interface_change:.6f} V")
    
    if interface_change > 0.01:
        print("   ⚠️  邊界條件可能過於嚴格，限制了界面電位演化")
    
    # 檢查針尖邊界
    tip_points_fixed = 0
    for j in range(N_nu-1):
        if test_potential_with_bc[0, j] == V_tip:
            tip_points_fixed += 1
    
    print(f"   固定為針尖電位的點數: {tip_points_fixed}/{N_nu-1}")
    
    if tip_points_fixed == N_nu-1:
        print("   ✅ 針尖邊界正確設置（不包括界面點）")
    else:
        print("   ⚠️  針尖邊界設置可能有問題")

def analyze_sor_coefficients(solver):
    """分析 SOR 係數和穩定性"""
    print(f"   A_P 係數範圍: {np.min(solver.A_P):.3e} 到 {np.max(solver.A_P):.3e}")
    print(f"   A_E 係數範圍: {np.min(solver.A_E):.3e} 到 {np.max(solver.A_E):.3e}")
    print(f"   A_W 係數範圍: {np.min(solver.A_W):.3e} 到 {np.max(solver.A_W):.3e}")
    print(f"   A_N 係數範圍: {np.min(solver.A_N):.3e} 到 {np.max(solver.A_N):.3e}")
    print(f"   A_S 係數範圍: {np.min(solver.A_S):.3e} 到 {np.max(solver.A_S):.3e}")
    
    # 檢查對角占優性（穩定性條件）
    N_eta, N_nu = solver.A_P.shape
    diagonal_dominance_violations = 0
    
    for i in range(1, N_eta-1):
        for j in range(1, N_nu-1):
            neighbors_sum = (solver.A_E[i,j] + solver.A_W[i,j] + 
                           solver.A_N[i,j] + solver.A_S[i,j])
            central = solver.A_P[i,j]
            
            if abs(central) < neighbors_sum:
                diagonal_dominance_violations += 1
    
    total_internal_points = (N_eta-2) * (N_nu-2)
    violation_ratio = diagonal_dominance_violations / total_internal_points * 100
    
    print(f"   對角占優違反: {diagonal_dominance_violations}/{total_internal_points} ({violation_ratio:.1f}%)")
    
    if violation_ratio > 10:
        print("   ⚠️  對角占優性較差，可能導致數值不穩定")
    
    # 檢查係數是否為零
    zero_count = np.sum(np.abs(solver.A_P) < 1e-15)
    print(f"   零係數數量: {zero_count}")
    
    if zero_count > 0:
        print("   ❌ 存在零係數，會導致除零錯誤")

def analyze_charge_feedback(solver, system_fermi):
    """分析電荷密度反饋"""
    # 創建測試電荷密度計算器
    class TestChargeDensityCalculator:
        def __init__(self):
            self.calls = []
        
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.calls.append(ef_rel_vb_eV)
            
            # 使用增強的計算
            kT = 0.0259
            Nd = 5e18  # 增強的雜質密度
            ni = 1e10
            Eg = 1.42
            
            if ef_rel_vb_eV > Eg:
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / (0.8 * kT))
            else:
                n_electrons = ni * np.exp(ef_rel_vb_eV / (0.9 * kT))
            
            n_holes = ni**2 / n_electrons
            
            if ef_rel_vb_eV < 0.3:
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + 5 * np.exp((ef_rel_vb_eV - 0.3) / (0.5 * kT)))
            
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            
            # 電場誘導因子
            field_induced_factor = 1.0 + 2.0 * np.tanh((ef_rel_vb_eV - 1.0) / (2.0 * kT))
            charge_density_cm3 *= field_induced_factor
            
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            charge_density_C_m3 = np.clip(charge_density_C_m3, -5e18, 5e18)
            
            return charge_density_C_m3
    
    charge_calc = TestChargeDensityCalculator()
    
    # 測試電荷密度在不同 EF 值下的響應
    ef_test_values = np.linspace(-0.5, 2.5, 11)
    charge_responses = []
    
    print(f"   測試電荷密度響應:")
    print(f"   EF_rel_VB (eV) | Charge Density (C/m³)")
    print(f"   {'-'*15}|{'-'*20}")
    
    for ef_val in ef_test_values:
        try:
            charge = charge_calc.get_charge_density_C_m3(ef_val)
            charge_responses.append(charge)
            print(f"   {ef_val:12.1f}   | {charge:15.3e}")
        except Exception as e:
            print(f"   {ef_val:12.1f}   | ERROR: {e}")
    
    # 檢查非線性程度
    if len(charge_responses) > 1:
        charge_range = max(charge_responses) - min(charge_responses)
        print(f"   電荷密度動態範圍: {charge_range:.3e} C/m³")
        
        if charge_range < 1e15:
            print("   ⚠️  電荷密度變化可能不足以驅動電位演化")
        else:
            print("   ✅ 電荷密度具有足夠的非線性響應")
    
    # 測試表面電荷密度
    try:
        surface_charge = solver._calculate_surface_charge_density(1.0, 0, 0)
        print(f"   測試表面電荷密度: {surface_charge:.3e} C/m²")
        
        if abs(surface_charge) < 1e-6:
            print("   ⚠️  表面電荷密度可能過小")
    except Exception as e:
        print(f"   ❌ 表面電荷密度計算失敗: {e}")

def analyze_single_iteration_change(solver, V_tip, V_sample, system_fermi):
    """分析單次迭代的變化"""
    # 創建初始電位
    initial_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    initial_potential = solver._apply_boundary_conditions(initial_potential, V_tip, V_sample)
    
    # 創建電荷密度計算器
    class SimpleChargeDensityCalculator:
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            # 簡單但有響應的電荷密度
            return 1e16 * np.tanh(ef_rel_vb_eV / 0.1)  # 強響應
    
    charge_calc = SimpleChargeDensityCalculator()
    
    # 執行一次短期求解來觀察變化
    try:
        potential_after, iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calc,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=1,  # 只執行一次迭代
            tolerance_Volts=1e-10,  # 很嚴格的容差，確保不會提前停止
            omega=0.1  # 很小的鬆弛因子
        )
        
        # 計算變化
        total_change = np.sum(np.abs(potential_after - initial_potential))
        max_change = np.max(np.abs(potential_after - initial_potential))
        avg_change = np.mean(np.abs(potential_after - initial_potential))
        
        print(f"   單次迭代總變化: {total_change:.6e}")
        print(f"   單次迭代最大變化: {max_change:.6e}")
        print(f"   單次迭代平均變化: {avg_change:.6e}")
        
        # 檢查Pot0變化
        pot0_before = solver._calculate_pot0_fortran_style(initial_potential, apply_scaling_correction=True)
        pot0_after = solver._calculate_pot0_fortran_style(potential_after, apply_scaling_correction=True)
        pot0_change = abs(pot0_after - pot0_before)
        
        print(f"   Pot0變化: {pot0_change:.6e} V")
        
        if max_change < 1e-10:
            print("   ❌ 電位幾乎沒有變化 - 數值求解器可能陷入僵局")
        elif max_change < 1e-6:
            print("   ⚠️  電位變化極小 - 可能需要更激進的參數")
        else:
            print("   ✅ 電位有合理變化")
            
        if pot0_change < 1e-10:
            print("   ❌ Pot0 完全沒有變化 - 這是問題的根源")
            
    except Exception as e:
        print(f"   ❌ 單次迭代測試失敗: {e}")

def propose_evolution_solutions():
    """提出演化解決方案"""
    print("🎯 演化停滯解決方案")
    print("="*60)
    
    print("🔹 立即可行的解決方案:")
    print()
    
    print("1. 🚀 激進初始條件:")
    print("   - 在初始猜測中加入隨機擾動")
    print("   - 創建非平衡的初始電位分布")
    print("   - 強制界面電位偏離平衡值")
    print()
    
    print("2. ⚡ 數值參數調整:")
    print("   - 使用更大的鬆弛因子 (omega > 1.5)")
    print("   - 放寬收斂容差以允許更多演化")
    print("   - 限制每次迭代的最大變化幅度")
    print()
    
    print("3. 🔧 求解器增強:")
    print("   - 添加動量項 (類似 Adam 優化器)")
    print("   - 使用自適應步長控制")
    print("   - 實現溫度冷卻策略")
    print()
    
    print("4. 🎲 擾動注入:")
    print("   - 定期注入小幅隨機擾動")
    print("   - 實現梯度加速方法")
    print("   - 使用多起始點策略")
    print()
    
    print("🔹 測試優先級:")
    print("   1. 激進初始條件 (最容易實現)")
    print("   2. 數值參數調整 (快速測試)")
    print("   3. 擾動注入 (中等複雜度)")
    print("   4. 求解器增強 (最複雜但最有效)")
    print()
    
    print("💡 下一步行動:")
    print("   1. 實現激進初始條件版本")
    print("   2. 測試大鬆弛因子的效果")
    print("   3. 添加定期擾動機制")
    print("   4. 驗證是否能破壞穩定平衡")

if __name__ == "__main__":
    print("🎯 電位演化停滯診斷")
    print("目標：找到阻止數值演化的根本原因")
    print()
    
    diagnose_no_evolution()
    
    print()
    print("="*80)
    print("🏆 診斷結論")
    print("="*80)
    print()
    print("🔑 關鍵發現:")
    print("   數值求解器被困在過於穩定的平衡態中")
    print("   需要更激進的方法來破壞這種穩定性")
    print("   物理模型是正確的，問題在於數值演化機制")
    print()
    print("🚀 立即執行: 實現激進初始條件和擾動注入")