#!/usr/bin/env python3
"""
測試多重網格實現
驗證能否觀察到Pot0符號轉變
"""
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

import numpy as np
import logging
from src.physics.solvers.multigrid import MultiGridPoissonSolver
from src.utils.constants import PhysicalConstants as PC

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(name)s - %(message)s')

def test_multigrid_solver():
    """測試多重網格求解器"""
    print("🧪 測試多重網格Poisson求解器")
    print("="*80)
    
    # 物理參數 (與Fortran相同)
    class Props:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = Props()
    
    # 針尖參數
    R = 1.0  # nm
    Z_TS = 1.0  # nm
    shank_slope = 1.0
    
    # 創建多重網格求解器
    multigrid_solver = MultiGridPoissonSolver(R, Z_TS, shank_slope, props)
    
    # 測試條件 (與Fortran完全相同)
    V_tip = -2.0707107
    V_sample = 0.0
    system_fermi = 1.4186435
    
    print(f"🎯 測試條件:")
    print(f"   V_tip = {V_tip:.7f} V")
    print(f"   V_sample = {V_sample:.1f} V")
    print(f"   System Fermi = {system_fermi:.7f} eV")
    print(f"   網格階段: {len(multigrid_solver.stages)}")
    print()
    
    for i, stage in enumerate(multigrid_solver.stages):
        print(f"   階段{i+1} ({stage.name}): {stage.grid_size[0]}×{stage.grid_size[1]}, "
              f"{stage.max_iterations}次迭代")
    print()
    
    # 改進的電荷密度計算器
    class ImprovedChargeDensityCalculator:
        def get_charge_density_C_m3(self, ef_rel_vb_eV):
            """更真實的電荷密度計算"""
            kT = 0.0259  # 300K
            
            # n型半導體參數
            Nd = 1e18  # cm^-3
            ni = 1e10  # cm^-3
            Eg = 1.42  # eV
            
            # 載流子密度計算
            if ef_rel_vb_eV > Eg:
                n_electrons = ni * np.exp((ef_rel_vb_eV - Eg) / kT)
            else:
                n_electrons = ni * np.exp(ef_rel_vb_eV / kT)
            
            n_holes = ni**2 / n_electrons
            
            # 離化雜質
            if ef_rel_vb_eV < 0.5:
                N_donors_ionized = Nd
            else:
                N_donors_ionized = Nd / (1 + 2 * np.exp((ef_rel_vb_eV - 0.5) / kT))
            
            # 總電荷密度
            charge_density_cm3 = N_donors_ionized + n_holes - n_electrons
            charge_density_C_m3 = charge_density_cm3 * 1e6 * PC.E
            
            # 限制範圍
            charge_density_C_m3 = np.clip(charge_density_C_m3, -1e18, 1e18)
            
            return charge_density_C_m3
    
    charge_calculator = ImprovedChargeDensityCalculator()
    
    print("🚀 開始多重網格求解...")
    print()
    
    try:
        # 執行多重網格求解
        result = multigrid_solver.solve_with_multigrid(
            V_tip=V_tip,
            V_sample=V_sample,
            charge_density_calculator=charge_calculator,
            system_fermi_level_E_F_main_eV=system_fermi
        )
        
        print("✅ 多重網格求解完成!")
        print()
        
        # 分析結果
        analyze_multigrid_results(result)
        
        return result
        
    except Exception as e:
        print(f"❌ 多重網格求解失敗: {e}")
        import traceback
        traceback.print_exc()
        return None

def analyze_multigrid_results(result):
    """分析多重網格結果"""
    print("📊 多重網格結果分析")
    print("="*60)
    
    if not result['success']:
        print("❌ 求解未成功")
        return
    
    print(f"🎯 總體結果:")
    print(f"   總迭代次數: {result['total_iterations']}")
    print(f"   最終Pot0: {result['final_pot0']:+.6f} V")
    print(f"   符號轉變: {'✅ 是' if result['sign_transition_achieved'] else '❌ 否'}")
    print()
    
    # 階段詳細分析
    print(f"📋 各階段詳細結果:")
    print()
    
    for i, stage_result in enumerate(result['stage_results']):
        stage_name = stage_result['stage_name']
        grid_size = stage_result['grid_size']
        iterations = stage_result['iterations']
        final_pot0 = stage_result['final_pot0']
        
        print(f"🔹 階段{i+1} ({stage_name}):")
        print(f"   網格大小: {grid_size[0]}×{grid_size[1]}")
        print(f"   迭代次數: {iterations}")
        print(f"   最終Pot0: {final_pot0:+.6f} V")
        
        # 第一階段的特別分析
        if i == 0 and 'pot0_evolution' in stage_result:
            analyze_coarse_stage_evolution(stage_result)
        
        print()
    
    # 與Fortran比較
    if 'fortran_comparison' in result:
        analyze_fortran_comparison(result['fortran_comparison'])

def analyze_coarse_stage_evolution(stage_result):
    """分析粗網格階段的Pot0演化"""
    pot0_evolution = stage_result.get('pot0_evolution', [])
    
    if not pot0_evolution:
        print("   ⚠️  無Pot0演化數據")
        return
    
    print(f"   📈 Pot0演化 ({len(pot0_evolution)}個數據點):")
    
    # 顯示關鍵演化點
    key_points = []
    if len(pot0_evolution) >= 1:
        key_points.append(("初始", pot0_evolution[0]))
    if len(pot0_evolution) >= 5:
        key_points.append(("中期", pot0_evolution[len(pot0_evolution)//2]))
    if len(pot0_evolution) >= 1:
        key_points.append(("最終", pot0_evolution[-1]))
    
    for label, (iter_num, pot0_val) in key_points:
        print(f"     {label}: ITER={iter_num:4d}, Pot0={pot0_val:+.6f} V")
    
    # 檢查符號轉變
    negative_points = [(iter_num, pot0) for iter_num, pot0 in pot0_evolution if pot0 < 0]
    positive_points = [(iter_num, pot0) for iter_num, pot0 in pot0_evolution if pot0 > 0]
    
    if negative_points and positive_points:
        last_negative = negative_points[-1]
        first_positive = positive_points[0]
        print(f"   🔄 符號轉變檢測:")
        print(f"     最後負值: ITER={last_negative[0]:4d}, Pot0={last_negative[1]:+.6f} V")
        print(f"     首個正值: ITER={first_positive[0]:4d}, Pot0={first_positive[1]:+.6f} V")
        
        # 估算轉變迭代
        estimated_transition = (last_negative[0] + first_positive[0]) // 2
        print(f"     估算轉變點: ~ITER={estimated_transition}")
    elif positive_points:
        print(f"   ✅ 全程正值 (可能從上階段繼承)")
    elif negative_points:
        print(f"   ❌ 全程負值 (未發生轉變)")

def analyze_fortran_comparison(comparison):
    """分析與Fortran的比較"""
    print(f"🔍 與Fortran比較:")
    print(f"   Fortran結果: {comparison['fortran_final']:+.6f} V")
    print(f"   Python結果:  {comparison['python_final']:+.6f} V") 
    print(f"   絕對差異:    {comparison['absolute_difference']:.6f} V")
    print(f"   相對誤差:    {comparison['relative_error']:.1f}%")
    print(f"   符號正確:    {'✅ 是' if comparison['sign_correct'] else '❌ 否'}")
    
    # 精度評估
    if comparison['absolute_difference'] < 0.01:
        print(f"   🎉 優秀精度! (<0.01V)")
    elif comparison['absolute_difference'] < 0.05:
        print(f"   ✅ 良好精度! (<0.05V)")
    elif comparison['absolute_difference'] < 0.1:
        print(f"   👍 可接受精度 (<0.1V)")
    else:
        print(f"   ⚠️  需要改善 (>0.1V)")

if __name__ == "__main__":
    print("🎯 多重網格Poisson求解器測試")
    print("目標：驗證能否觀察到Pot0從負值到正值的轉變")
    print()
    
    # 執行測試
    result = test_multigrid_solver()
    
    if result:
        print()
        print("="*80)
        print("🏆 測試總結")
        print("="*80)
        
        success_indicators = []
        issues = []
        
        # 檢查成功指標
        if result['sign_transition_achieved']:
            success_indicators.append("✅ 成功實現符號轉變")
        else:
            issues.append("❌ 未實現符號轉變")
        
        if result['final_pot0'] > 0:
            success_indicators.append("✅ 最終結果為正值")
        else:
            issues.append("❌ 最終結果仍為負值")
        
        if 'fortran_comparison' in result:
            comp = result['fortran_comparison']
            if comp['sign_correct']:
                success_indicators.append("✅ 與Fortran符號一致")
            else:
                issues.append("❌ 與Fortran符號不一致")
                
            if comp['absolute_difference'] < 0.1:
                success_indicators.append("✅ 與Fortran精度可接受")
            else:
                issues.append("❌ 與Fortran差異過大")
        
        print("成功指標:")
        for indicator in success_indicators:
            print(f"  {indicator}")
        
        if issues:
            print()
            print("需要改善:")
            for issue in issues:
                print(f"  {issue}")
        
        print()
        if len(success_indicators) >= 3:
            print("🎉 多重網格實現基本成功!")
            print("   可以進入下一階段：完善和優化")
        elif len(success_indicators) >= 2:
            print("📈 多重網格實現部分成功")
            print("   需要調整參數或物理模型")
        else:
            print("🔧 多重網格實現需要重大改進")
            print("   需要檢查核心算法實現")
    else:
        print("❌ 測試失敗，需要修復基礎實現")