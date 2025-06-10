#!/usr/bin/env python3
"""
手動分析 Pot0 計算過程，對比 Fortran PCENT 函數實現
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

logging.basicConfig(level=logging.INFO, format='%(message)s')

def manual_pot0_analysis():
    """手動追蹤 Pot0 計算過程，對比 Fortran PCENT 實現"""
    print("手動分析 Pot0 計算過程")
    print("="*80)
    
    # 創建測試網格
    grid = HyperbolicGrid(N_eta=16, N_nu=8, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    # 創建 Poisson 求解器
    class MinimalProps:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = MinimalProps()
    solver = PoissonSOREquation(grid, props)
    
    # 測試參數
    V_tip = -2.07  # Fortran 調整後的偏壓
    V_sample = 0.0
    
    print(f"測試條件:")
    print(f"  V_tip = {V_tip} V")
    print(f"  V_sample = {V_sample} V")
    print(f"  Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    print(f"  Grid params: f={grid.f:.4f} nm, eta_tip={grid.eta_tip:.4f}")
    
    # 運行 Laplace 求解
    potential, iterations, error = solver.solve_laplace(
        V_tip, V_sample, max_iterations=500, tolerance=1e-4
    )
    
    print(f"\nLaplace 求解完成: {iterations} 次迭代, 誤差={error:.3e}")
    
    # === 分析 Fortran PCENT 函數邏輯 ===
    print(f"\n" + "="*80)
    print(f"Fortran PCENT 函數分析:")
    print(f"="*80)
    
    print(f"Fortran PCENT(0) 邏輯:")
    print(f"  JJ=0 時使用 VSINT 陣列")
    print(f"  I=1 (固定)")
    print(f"  公式: SUM = Σ(K=1 to NP) [(9*VSINT(1,1,K) - VSINT(1,2,K))/8]")
    print(f"  結果: PCENT = SUM/NP")
    
    # === 模擬 Fortran 的座標系統 ===
    print(f"\n分析座標映射:")
    
    # Fortran 使用的是：
    # I=1,2,... 對應徑向方向 (我們的 eta 方向)
    # K=1,2,... 對應角度方向 (我們的部分 nu 方向？)
    # VSINT(1,I,K) 是半導體表面電位
    
    # 在我們的實現中，半導體表面對應 nu=N_nu-1
    interface_nu_idx = grid.N_nu - 1
    
    print(f"  半導體表面位置: nu_idx = {interface_nu_idx}")
    print(f"  Fortran I=1 對應我們的 eta_idx = 0")
    print(f"  Fortran I=2 對應我們的 eta_idx = 1")
    
    # 提取"VSINT"等價值
    print(f"\n提取等價 VSINT 值:")
    vsint_1_1 = potential[0, interface_nu_idx]  # VSINT(1,1,K) 等價
    vsint_1_2 = potential[1, interface_nu_idx]  # VSINT(1,2,K) 等價
    
    print(f"  VSINT(1,1,K) 等價: potential[0,{interface_nu_idx}] = {vsint_1_1:.6f} V")
    print(f"  VSINT(1,2,K) 等價: potential[1,{interface_nu_idx}] = {vsint_1_2:.6f} V")
    
    # 應用 Fortran PCENT 公式
    print(f"\n應用 Fortran PCENT 公式:")
    
    # 由於我們是 2D 軸對稱，NP=1 (只有一個角度分量)
    NP = 1  # 軸對稱情況
    weighted_value = (9.0 * vsint_1_1 - vsint_1_2) / 8.0
    pot0_fortran_style = weighted_value / NP
    
    print(f"  加權值: (9*{vsint_1_1:.6f} - {vsint_1_2:.6f})/8 = {weighted_value:.6f}")
    print(f"  Fortran PCENT: {weighted_value:.6f}/1 = {pot0_fortran_style:.6f} V")
    
    # 對比我們當前的計算
    pot0_current = solver._calculate_pot0_fortran_style(potential)
    print(f"  我們的計算: {pot0_current:.6f} V")
    
    print(f"\n" + "="*80)
    print(f"關鍵差異分析:")
    print(f"="*80)
    
    print(f"1. 數值對比:")
    print(f"   Fortran 參考值: -0.08 V")
    print(f"   手動 PCENT 計算: {pot0_fortran_style:.6f} V")
    print(f"   我們當前計算: {pot0_current:.6f} V")
    print(f"   與 Fortran 差異: {abs(pot0_fortran_style - (-0.08)):.6f} V")
    
    print(f"\n2. 正負號分析:")
    if pot0_fortran_style * (-0.08) > 0:
        print(f"   ✅ 正負號一致")
    else:
        print(f"   ❌ 正負號相反！這是主要問題")
    
    print(f"\n3. 可能的問題:")
    print(f"   a) 座標系統定義不同")
    print(f"   b) VSINT 在 Fortran 中是專門計算的表面電位")
    print(f"   c) 我們使用的是整體電位矩陣，不是專門的表面電位")
    print(f"   d) 邊界條件或初始條件不同")
    print(f"   e) 電荷密度效應的影響")
    
    # === 測試不同的 VSINT 解釋 ===
    print(f"\n" + "="*80)
    print(f"測試不同的 VSINT 解釋:")
    print(f"="*80)
    
    # 可能性1: VSINT 是相對於某個參考值的
    print(f"\n可能性1: VSINT 是相對值")
    ref_potential = V_sample  # 或其他參考值
    vsint_rel_1 = vsint_1_1 - ref_potential
    vsint_rel_2 = vsint_1_2 - ref_potential
    weighted_rel = (9.0 * vsint_rel_1 - vsint_rel_2) / 8.0
    print(f"   相對 VSINT(1,1): {vsint_rel_1:.6f} V")
    print(f"   相對 VSINT(1,2): {vsint_rel_2:.6f} V")
    print(f"   相對 PCENT: {weighted_rel:.6f} V")
    
    # 可能性2: 座標軸定義不同
    print(f"\n可能性2: 座標軸或層序不同")
    # 試試其他位置
    for eta_test in range(min(3, grid.N_eta)):
        for nu_test in [grid.N_nu-1, grid.N_nu-2]:
            v_test = potential[eta_test, nu_test]
            print(f"   potential[{eta_test},{nu_test}] = {v_test:.6f} V")
    
    # 可能性3: 符號定義不同
    print(f"\n可能性3: 符號定義不同")
    pot0_negative = -pot0_fortran_style
    print(f"   負號版本: {pot0_negative:.6f} V")
    print(f"   與 Fortran 差異: {abs(pot0_negative - (-0.08)):.6f} V")
    
    return {
        'pot0_manual': pot0_fortran_style,
        'pot0_current': pot0_current,
        'vsint_1_1': vsint_1_1,
        'vsint_1_2': vsint_1_2,
        'potential': potential
    }

if __name__ == "__main__":
    results = manual_pot0_analysis()
    
    print(f"\n" + "="*80)
    print(f"總結:")
    print(f"="*80)
    print(f"手動 PCENT 計算: {results['pot0_manual']:.6f} V")
    print(f"與 Fortran 目標 (-0.08V) 差異: {abs(results['pot0_manual'] - (-0.08)):.6f} V")
    
    if abs(results['pot0_manual'] - (-0.08)) > abs(results['pot0_current'] - (-0.08)):
        print(f"❌ 手動計算更差，說明問題不在公式本身")
    else:
        print(f"✅ 手動計算更好，問題在實現細節")