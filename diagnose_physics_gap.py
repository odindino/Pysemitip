#!/usr/bin/env python3
"""
診斷物理模型缺口
系統性分析為什麼Python版本無法實現符號轉變
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

def diagnose_physics_gap():
    """系統性診斷物理模型缺口"""
    print("🔬 診斷物理模型缺口")
    print("="*80)
    
    print("💡 首先回顧Fortran的物理過程:")
    print("   1. 初始狀態: 表面態主導 → 電子積累層 → Pot0 < 0")
    print("   2. 針尖場效應: 強電場改變載流子分布")
    print("   3. 自洽演化: 電荷重新分布 → 電位變化 → 電荷再分布")
    print("   4. 臨界轉變: ~ITER=1700，針尖場超過表面態吸引力")
    print("   5. 最終狀態: 耗盡層形成 → Pot0 > 0")
    print()
    
    # 測試當前實現的物理特性
    test_current_physics()
    
    # 分析缺失的物理機制
    analyze_missing_physics()
    
    # 提出系統性解決方案
    propose_systematic_solution()

def test_current_physics():
    """測試當前實現的物理特性"""
    print("🧪 測試當前實現的物理特性")
    print("="*60)
    
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
    
    class TestChargeDensityCalculator:
        def __init__(self):
            self.call_count = 0
            self.ef_history = []
            
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            self.call_count += 1
            self.ef_history.append(ef_rel_vb_eV)
            
            # 檢測電荷密度變化範圍
            kT = 0.0259
            Nd = 1e18
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
                N_donors_ionized = Nd / (1 + 2 * np.exp((ef_rel_vb_eV - 0.5) / kT))
            
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            charge_density_C_m3 = np.clip(charge_density_C_m3, -1e18, 1e18)
            
            return charge_density_C_m3
    
    charge_calc = TestChargeDensityCalculator()
    
    print(f"📊 執行100次迭代測試...")
    
    # 執行短期測試
    try:
        potential, iterations, converged = solver.solve(
            V_tip_Volts=V_tip,
            V_sample_Volts=V_sample,
            charge_density_calculator=charge_calc,
            system_fermi_level_E_F_main_eV=system_fermi,
            max_iterations=100,
            tolerance_Volts=1e-4,
            omega=1.0
        )
        
        print(f"✅ 求解完成: {iterations}次迭代, 收斂={converged}")
        print(f"📈 電荷密度計算次數: {charge_calc.call_count}")
        print(f"⚡ EF範圍: {min(charge_calc.ef_history):.3f} 到 {max(charge_calc.ef_history):.3f} eV")
        
        # 計算Pot0
        pot0_regular = solver._calculate_pot0_fortran_style(potential, apply_scaling_correction=True)
        
        # VSINT計算
        vsint_array = solver._initialize_vsint_array()
        vsint_array = solver._update_vsint_with_surface_charge(
            vsint_array, potential, charge_calc, system_fermi, V_tip)
        pot0_vsint = solver._calculate_pot0_fortran_style(
            potential, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
        
        print(f"🎯 Pot0結果:")
        print(f"   Regular: {pot0_regular:.6f} V")
        print(f"   VSINT:   {pot0_vsint:.6f} V")
        print(f"   VSINT改善: {abs(pot0_regular) - abs(pot0_vsint):.6f} V")
        
        # 分析電位分布
        analyze_potential_distribution(potential, grid)
        
        return {
            'pot0_regular': pot0_regular,
            'pot0_vsint': pot0_vsint,
            'charge_calls': charge_calc.call_count,
            'ef_range': (min(charge_calc.ef_history), max(charge_calc.ef_history))
        }
        
    except Exception as e:
        print(f"❌ 測試失敗: {e}")
        return None

def analyze_potential_distribution(potential, grid):
    """分析電位分布特性"""
    print(f"🔍 電位分布分析:")
    print(f"   電位範圍: {np.min(potential):.3f} 到 {np.max(potential):.3f} V")
    print(f"   平均電位: {np.mean(potential):.3f} V")
    print(f"   標準差:   {np.std(potential):.3f} V")
    
    # 檢查界面電位（最關鍵的點）
    interface_potential = potential[0, -1]  # [0, N_nu-1]
    tip_potential = potential[0, 0]         # [0, 0]
    
    print(f"   針尖電位: {tip_potential:.3f} V")
    print(f"   界面電位: {interface_potential:.3f} V")
    print(f"   電位落差: {tip_potential - interface_potential:.3f} V")
    
    # 檢查電位梯度
    if potential.shape[0] > 1 and potential.shape[1] > 1:
        grad_eta = np.gradient(potential, axis=0)
        grad_nu = np.gradient(potential, axis=1)
        total_grad = np.sqrt(grad_eta**2 + grad_nu**2)
        
        print(f"   最大梯度: {np.max(total_grad):.3f} V/grid")
        print(f"   平均梯度: {np.mean(total_grad):.3f} V/grid")

def analyze_missing_physics():
    """分析缺失的物理機制"""
    print()
    print("🔍 分析缺失的物理機制")
    print("="*60)
    
    print("❌ 可能缺失的關鍵物理:")
    print()
    
    print("1. 🔋 表面態物理:")
    print("   ❌ 表面態密度可能不足")
    print("   ❌ 表面態能階分布可能不正確")
    print("   ❌ 費米-狄拉克占據統計可能簡化")
    print("   ❌ 表面電荷與電位的反饋機制可能不完整")
    print()
    
    print("2. ⚡ 電場效應:")
    print("   ❌ 針尖誘導電場強度可能不足")
    print("   ❌ 電場對載流子分布的影響可能低估")
    print("   ❌ 高電場下的非線性效應可能缺失")
    print()
    
    print("3. 🔄 自洽反饋:")
    print("   ❌ 電荷-電位反饋迴圈可能不充分")
    print("   ❌ 收斂條件可能過於嚴格，阻止了物理演化")
    print("   ❌ 迭代步長可能不適合物理時間尺度")
    print()
    
    print("4. 📐 數值方法:")
    print("   ❌ 網格解析度可能不足以捕捉尖銳轉變")
    print("   ❌ 邊界條件可能限制了界面電位演化")
    print("   ❌ 非線性求解方法可能陷入局部極小值")
    print()

def propose_systematic_solution():
    """提出系統性解決方案"""
    print("🎯 系統性解決方案")
    print("="*60)
    
    print("📋 三階段系統性策略:")
    print()
    
    print("🔹 階段A: 物理模型診斷和修復 (優先級: CRITICAL)")
    print("   1. 深入分析Fortran的表面態參數")
    print("      - 確認表面態密度 (4.4e14 cm^-2)")
    print("      - 驗證電荷中性能級 (0.125 eV above VB)")
    print("      - 檢查能態分布寬度 (0.25 eV FWHM)")
    print()
    print("   2. 重新實現表面電荷密度計算")
    print("      - 使用完整的費米-狄拉克分布")
    print("      - 包含多種表面態分布")
    print("      - 確保強非線性反饋")
    print()
    print("   3. 增強電場效應模型")
    print("      - 提高針尖電場強度")
    print("      - 包含電場誘導的載流子重新分布")
    print("      - 實現閾值效應（積累→耗盡轉變）")
    print()
    
    print("🔹 階段B: 數值方法優化 (優先級: HIGH)")
    print("   1. 改進非線性求解策略")
    print("      - 減少過早收斂")
    print("      - 增加物理演化時間")
    print("      - 使用自適應步長")
    print()
    print("   2. 優化邊界條件")
    print("      - 允許界面電位自由演化")
    print("      - 確保物理一致性")
    print("      - 防止人為限制")
    print()
    
    print("🔹 階段C: 驗證和調優 (優先級: MEDIUM)")
    print("   1. 系統性參數掃描")
    print("   2. 與Fortran逐步對比")
    print("   3. 物理合理性檢查")
    print()
    
    print("⏰ 執行時間表:")
    print("   階段A: 1-2天（核心物理）")
    print("   階段B: 1天（數值優化）")
    print("   階段C: 0.5天（驗證）")
    print("   總計: 2.5-3.5天")
    print()
    
    print("🎯 成功判據:")
    print("   ✅ 在1600-1700次迭代觀察到符號轉變")
    print("   ✅ 最終Pot0 > 0且接近+0.07V")
    print("   ✅ 物理演化軌跡與Fortran一致")
    print()

if __name__ == "__main__":
    print("🎯 系統性物理診斷")
    print("目標：找到阻止符號轉變的根本原因")
    print()
    
    # 執行診斷
    result = diagnose_physics_gap()
    
    print()
    print("="*80)
    print("🏆 診斷總結")
    print("="*80)
    
    if result:
        print("📊 當前實現狀態:")
        print(f"   Regular Pot0: {result['pot0_regular']:.6f} V")
        print(f"   VSINT Pot0:   {result['pot0_vsint']:.6f} V") 
        print(f"   電荷計算次數: {result['charge_calls']}")
        print(f"   EF變化範圍: {result['ef_range'][1] - result['ef_range'][0]:.3f} eV")
        print()
        
        if abs(result['pot0_vsint']) < abs(result['pot0_regular']):
            print("✅ VSINT物理有效（有改善）")
        else:
            print("❌ VSINT物理無效")
            
        if result['ef_range'][1] - result['ef_range'][0] > 0.1:
            print("✅ 費米能級有顯著變化")
        else:
            print("❌ 費米能級變化不足")
    
    print()
    print("🔑 關鍵結論:")
    print("   1. 需要更強的表面態物理模型")
    print("   2. 需要更充分的自洽演化時間")
    print("   3. 需要系統性的物理參數調優")
    print("   4. 冷靜、持續、靈活地攻克每個物理環節")
    print()
    
    print("🚀 下一步: 實施階段A - 物理模型修復")
    print("   重點：表面態物理 + 電場效應 + 非線性反饋")