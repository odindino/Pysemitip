#!/usr/bin/env python3
"""
測試激進演化策略
實現破壞穩定平衡的激進初始條件和數值參數
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

def test_aggressive_evolution():
    """測試激進演化策略"""
    print("🚀 測試激進演化策略")
    print("="*80)
    print("🎯 目標：破壞穩定平衡，實現 Pot0 符號轉變")
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
    
    # 創建強響應電荷密度計算器
    class AggressiveChargeDensityCalculator:
        def __init__(self):
            self.call_count = 0
            self.ef_history = []
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.call_count += 1
            self.ef_history.append(ef_rel_vb_eV)
            
            # 🔑 超強響應的電荷密度計算
            kT = 0.0259
            
            # 大幅增強參數
            Nd = 1e19  # 極高雜質密度
            ni = 1e10
            Eg = 1.42
            
            # 極敏感的載流子計算
            if ef_rel_vb_eV > Eg:
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / (0.5 * kT))  # 超敏感
            else:
                n_electrons = ni * np.exp(ef_rel_vb_eV / (0.5 * kT))  # 超敏感
            
            n_holes = ni**2 / n_electrons
            
            # 極敏感的離化
            if ef_rel_vb_eV < 0.2:
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + 10 * np.exp((ef_rel_vb_eV - 0.2) / (0.3 * kT)))
            
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            
            # 超強電場誘導效應
            field_factor = 1.0 + 5.0 * np.tanh((ef_rel_vb_eV - 0.8) / (0.5 * kT))
            charge_density_cm3 *= field_factor
            
            # 轉換為 C/m³ 並放大動態範圍
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            charge_density_C_m3 = np.clip(charge_density_C_m3, -1e19, 1e19)  # 更大範圍
            
            return charge_density_C_m3
    
    charge_calc = AggressiveChargeDensityCalculator()
    
    # 策略 1: 激進初始條件
    print("🔹 策略 1: 激進初始條件")
    print("-" * 50)
    test_aggressive_initial_conditions(solver, V_tip, V_sample, charge_calc, system_fermi)
    print()
    
    # 策略 2: 激進數值參數
    print("🔹 策略 2: 激進數值參數")
    print("-" * 50)
    test_aggressive_numerical_parameters(solver, V_tip, V_sample, charge_calc, system_fermi)
    print()
    
    # 策略 3: 擾動注入
    print("🔹 策略 3: 擾動注入")
    print("-" * 50)
    test_perturbation_injection(solver, V_tip, V_sample, charge_calc, system_fermi)
    print()

def test_aggressive_initial_conditions(solver, V_tip, V_sample, charge_calc, system_fermi):
    """測試激進初始條件"""
    print("🔥 創建激進的非平衡初始條件...")
    
    # 修改 solver 的初始猜測方法
    def aggressive_initial_guess(V_tip, V_sample):
        try:
            N_eta, N_nu = solver.grid.N_eta, solver.grid.N_nu
        except AttributeError:
            N_eta, N_nu = solver.potential.shape
            
        potential = np.zeros((N_eta, N_nu))
        
        # 🔑 激進策略1: 添加強隨機擾動
        random_perturbation = np.random.normal(0, 0.5, (N_eta, N_nu))  # 大擾動
        
        # 🔑 激進策略2: 創建非平衡的初始分布
        for i in range(N_eta):
            for j in range(N_nu):
                nu_fraction = j / max(N_nu - 1, 1)
                eta_fraction = i / max(N_eta - 1, 1)
                
                # 非線性分布（不是簡單線性插值）
                base_potential = V_tip * (1 - nu_fraction**2) + V_sample * nu_fraction**2
                
                # 添加空間依賴的激勵
                spatial_excitation = 0.3 * np.sin(np.pi * nu_fraction) * np.exp(-2 * eta_fraction)
                
                # 界面特殊處理：強制偏離平衡
                if j == N_nu - 1:  # 界面
                    # 強制界面電位偏離期望值
                    interface_deviation = 0.8 * (eta_fraction - 0.5)  # -0.4 到 +0.4 V
                    potential[i, j] = V_sample + interface_deviation
                else:
                    potential[i, j] = base_potential + spatial_excitation
                
                # 添加隨機擾動
                potential[i, j] += random_perturbation[i, j]
        
        # 🔑 激進策略3: 在關鍵區域加入額外激勵
        # 在針尖附近加入強電場效應
        for j in range(min(3, N_nu)):
            potential[0, j] = V_tip + 0.2 * np.sin(j * np.pi / 2)  # 非均勻針尖電位
        
        return potential
    
    # 替換求解器的初始猜測方法
    original_method = solver._create_initial_potential_guess
    solver._create_initial_potential_guess = aggressive_initial_guess
    
    print("   🎲 添加大幅隨機擾動 (±0.5V)")
    print("   ⚡ 創建非平衡空間分布")
    print("   🔥 強制界面電位偏離平衡")
    print("   🌪️  在針尖區域加入強激勵")
    
    try:
        # 執行求解
        potential, iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calc,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=500,  # 中等迭代數
            tolerance_Volts=1e-2,  # 放寬容差
            omega=1.5  # 激進鬆弛因子
        )
        
        # 計算結果
        pot0_result = solver._calculate_pot0_fortran_style(potential, apply_scaling_correction=True)
        
        print(f"   ✅ 求解完成: {iterations} 次迭代")
        print(f"   🎯 Pot0 結果: {pot0_result:+.6f} V")
        print(f"   📊 電荷計算: {charge_calc.call_count} 次")
        
        if charge_calc.ef_history:
            ef_range = max(charge_calc.ef_history) - min(charge_calc.ef_history)
            print(f"   ⚡ EF 變化範圍: {ef_range:.3f} eV")
        
        # 檢查是否有改善
        if abs(pot0_result) < 0.1:
            print("   🎉 顯著改善！Pot0 接近零")
        elif pot0_result > 0:
            print("   ✅ 成功！實現正值轉變")
        else:
            print("   📈 有進展，但仍需優化")
            
    except Exception as e:
        print(f"   ❌ 激進初始條件測試失敗: {e}")
    finally:
        # 恢復原始方法
        solver._create_initial_potential_guess = original_method

def test_aggressive_numerical_parameters(solver, V_tip, V_sample, charge_calc, system_fermi):
    """測試激進數值參數"""
    print("⚡ 測試極端數值參數...")
    
    # 🔑 策略：極大的鬆弛因子 + 極寬鬆的容差
    aggressive_params = [
        {"omega": 1.8, "tolerance": 1e-1, "name": "超激進"},
        {"omega": 1.6, "tolerance": 1e-2, "name": "激進"},
        {"omega": 1.4, "tolerance": 1e-3, "name": "中激進"}
    ]
    
    for params in aggressive_params:
        print(f"   🧪 測試 {params['name']} 參數:")
        print(f"      omega = {params['omega']}, tolerance = {params['tolerance']}")
        
        try:
            potential, iterations, converged = solver.solve(
                V_tip_Volts=V_tip,
                V_sample_Volts=V_sample,
                charge_density_calculator=charge_calc,
                system_fermi_level_E_F_main_eV=system_fermi,
                max_iterations=1000,
                tolerance_Volts=params['tolerance'],
                omega=params['omega']
            )
            
            pot0_result = solver._calculate_pot0_fortran_style(potential, apply_scaling_correction=True)
            
            print(f"      結果: Pot0 = {pot0_result:+.6f} V, {iterations} 次迭代")
            
            if pot0_result > 0:
                print(f"      🎉 成功實現正值！")
                break
            elif abs(pot0_result) < 0.05:
                print(f"      🌟 接近成功！")
            else:
                print(f"      📈 有改善")
                
        except Exception as e:
            print(f"      ❌ 失敗: {e}")
        print()

def test_perturbation_injection(solver, V_tip, V_sample, charge_calc, system_fermi):
    """測試擾動注入策略"""
    print("🎲 測試定期擾動注入...")
    
    # 🔑 策略：分段求解，每段注入擾動
    current_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    total_iterations = 0
    pot0_evolution = []
    
    num_segments = 20
    iterations_per_segment = 100
    perturbation_strength = 0.1  # 擾動強度
    
    print(f"   📊 計劃: {num_segments} 段，每段 {iterations_per_segment} 次迭代")
    print(f"   🎲 擾動強度: ±{perturbation_strength} V")
    
    for segment in range(num_segments):
        print(f"   📍 段 {segment+1}/{num_segments}: ", end="")
        
        try:
            # 在當前電位基礎上添加小幅擾動
            if segment > 0:
                N_eta, N_nu = current_potential.shape
                perturbation = np.random.normal(0, perturbation_strength, (N_eta, N_nu))
                
                # 只在內部區域添加擾動，保持邊界
                current_potential[1:-1, 1:-1] += perturbation[1:-1, 1:-1]
                
                # 重新應用邊界條件
                current_potential = solver._apply_boundary_conditions(current_potential, V_tip, V_sample)
            
            # 執行這一段的求解
            # 暫時修改 solver 的初始猜測
            def get_current_potential(V_tip, V_sample):
                return np.copy(current_potential)
            
            original_method = solver._create_initial_potential_guess
            solver._create_initial_potential_guess = get_current_potential
            
            try:
                segment_potential, segment_iters, converged = solver.solve(
                    V_tip_Volts=V_tip,
                    V_sample_Volts=V_sample,
                    charge_density_calculator=charge_calc,
                    system_fermi_level_E_F_main_eV=system_fermi,
                    max_iterations=iterations_per_segment,
                    tolerance_Volts=1e-3,
                    omega=1.3
                )
                
                current_potential = segment_potential
                total_iterations += segment_iters
                
                # 計算當前 Pot0
                pot0_current = solver._calculate_pot0_fortran_style(current_potential, apply_scaling_correction=True)
                pot0_evolution.append((total_iterations, pot0_current))
                
                print(f"ITER={total_iterations:4d}, Pot0={pot0_current:+.6f}V")
                
                # 檢查符號轉變
                if len(pot0_evolution) >= 2:
                    prev_pot0 = pot0_evolution[-2][1]
                    if prev_pot0 < 0 and pot0_current > 0:
                        print(f"      🔄 符號轉變！{prev_pot0:.6f} → {pot0_current:.6f}")
                        break
                
            finally:
                solver._create_initial_potential_guess = original_method
                
        except Exception as e:
            print(f"失敗: {e}")
            continue
    
    # 分析演化結果
    if pot0_evolution:
        print()
        print(f"   📈 演化分析:")
        initial_pot0 = pot0_evolution[0][1]
        final_pot0 = pot0_evolution[-1][1]
        total_change = final_pot0 - initial_pot0
        
        print(f"      初始 Pot0: {initial_pot0:+.6f} V")
        print(f"      最終 Pot0: {final_pot0:+.6f} V")
        print(f"      總變化:    {total_change:+.6f} V")
        
        # 檢查是否有符號轉變
        negative_count = sum(1 for _, pot0 in pot0_evolution if pot0 < 0)
        positive_count = len(pot0_evolution) - negative_count
        
        if negative_count > 0 and positive_count > 0:
            print(f"      🎉 成功實現符號轉變！")
        elif final_pot0 > 0:
            print(f"      ✅ 最終達到正值")
        elif abs(total_change) > 0.01:
            print(f"      📈 有顯著演化")
        else:
            print(f"      ❌ 演化不足")

if __name__ == "__main__":
    print("🎯 激進演化策略測試")
    print("目標：破壞穩定平衡，實現 Pot0 符號轉變")
    print("基於診斷結果的系統性解決方案")
    print()
    
    test_aggressive_evolution()
    
    print()
    print("="*80)
    print("🏆 激進策略測試總結")
    print("="*80)
    print()
    print("🔑 核心策略:")
    print("   1. 激進初始條件 - 打破初始平衡")
    print("   2. 激進數值參數 - 允許大幅變化")
    print("   3. 擾動注入 - 持續破壞局部平衡")
    print()
    print("💡 如果成功:")
    print("   可以進入下一階段 - 精確匹配 Fortran 結果")
    print("   如果失敗:")
    print("   需要更深入的數值方法改革")