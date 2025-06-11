#!/usr/bin/env python3
"""
系統性分析剩餘0.63V Pot0差異的根源
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid
from src.utils.constants import PhysicalConstants as PC

logging.basicConfig(level=logging.INFO, format='%(message)s')

def analyze_remaining_difference():
    """系統性分析Pot0計算差異"""
    print("系統性分析剩餘Pot0差異")
    print("="*80)
    
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
    
    # 測試參數
    V_tip = -2.07
    V_sample = 0.0
    system_fermi_level = 1.4187
    fortran_target = -0.08
    
    print(f"測試條件:")
    print(f"  V_tip = {V_tip} V")
    print(f"  V_sample = {V_sample} V")
    print(f"  System Fermi level = {system_fermi_level} eV")
    print(f"  Fortran 目標 = {fortran_target} V")
    
    # 基準：Laplace求解
    potential_laplace, _, _ = solver.solve_laplace(V_tip, V_sample, max_iterations=200, tolerance=1e-4)
    pot0_laplace = solver._calculate_pot0_fortran_style(potential_laplace)
    
    print(f"\n" + "="*60)
    print(f"基準結果 (Laplace):")
    print(f"  Pot0 = {pot0_laplace:.6f} V")
    print(f"  差異 = {abs(pot0_laplace - fortran_target):.6f} V")
    
    # 分析1: 測試不同的PCENT公式變體
    print(f"\n" + "="*60)
    print(f"分析1: 測試不同PCENT公式變體")
    print(f"="*60)
    
    interface_idx = grid.N_nu - 1
    v1 = potential_laplace[0, interface_idx]
    v2 = potential_laplace[1, interface_idx]
    
    # Fortran原始公式: (9*V1 - V2)/8
    pcent_original = (9.0 * v1 - v2) / 8.0
    
    # 變體1: 不同權重
    pcent_var1 = (8.0 * v1 - v2) / 7.0
    pcent_var2 = (10.0 * v1 - v2) / 9.0
    
    # 變體2: 添加常數偏移
    pcent_with_offset1 = pcent_original + 0.63  # 直接補償差異
    pcent_with_offset2 = pcent_original * 0.113  # 縮放到Fortran範圍 (-0.08/-0.707)
    
    # 變體3: 使用不同網格點
    if grid.N_eta > 2:
        v3 = potential_laplace[2, interface_idx]
        pcent_3points = (16.0 * v1 - 9.0 * v2 + v3) / 8.0  # 三點公式
    else:
        pcent_3points = pcent_original
    
    print(f"  界面電位值:")
    print(f"    V[0,{interface_idx}] = {v1:.6f} V")
    print(f"    V[1,{interface_idx}] = {v2:.6f} V")
    if grid.N_eta > 2:
        print(f"    V[2,{interface_idx}] = {potential_laplace[2, interface_idx]:.6f} V")
    
    print(f"  PCENT公式變體:")
    print(f"    原始 (9*V1-V2)/8:     {pcent_original:.6f} V (差異: {abs(pcent_original-fortran_target):.6f})")
    print(f"    變體1 (8*V1-V2)/7:     {pcent_var1:.6f} V (差異: {abs(pcent_var1-fortran_target):.6f})")
    print(f"    變體2 (10*V1-V2)/9:    {pcent_var2:.6f} V (差異: {abs(pcent_var2-fortran_target):.6f})")
    print(f"    偏移1 +0.63:          {pcent_with_offset1:.6f} V (差異: {abs(pcent_with_offset1-fortran_target):.6f})")
    print(f"    縮放 *0.113:          {pcent_with_offset2:.6f} V (差異: {abs(pcent_with_offset2-fortran_target):.6f})")
    print(f"    三點公式:            {pcent_3points:.6f} V (差異: {abs(pcent_3points-fortran_target):.6f})")
    
    # 分析2: 檢查物理參數和單位轉換
    print(f"\n" + "="*60)
    print(f"分析2: 物理參數和單位一致性")
    print(f"="*60)
    
    # 檢查網格參數
    print(f"  網格參數:")
    print(f"    N_eta = {grid.N_eta}, N_nu = {grid.N_nu}")
    print(f"    f = {grid.f:.6f} nm")
    print(f"    eta_tip = {grid.eta_tip:.6f}")
    print(f"    a_nm = {getattr(grid, 'a_nm', 'N/A')}")
    
    # 檢查電介常數
    print(f"  物理參數:")
    print(f"    epsilon_r = {props.semiconductor_props.epsilon_r}")
    print(f"    Ev_offset = {props.semiconductor_props.Ev_offset_eV} eV")
    
    # 分析3: 測試座標系統映射
    print(f"\n" + "="*60)
    print(f"分析3: 座標系統映射")
    print(f"="*60)
    
    # 檢查雙曲座標
    eta_vals = np.linspace(grid.eta_tip, 2.0, grid.N_eta)
    nu_vals = np.linspace(0, np.pi/2, grid.N_nu)
    
    print(f"  雙曲座標範圍:")
    print(f"    eta: {eta_vals[0]:.4f} to {eta_vals[-1]:.4f}")
    print(f"    nu: {nu_vals[0]:.4f} to {nu_vals[-1]:.4f}")
    
    # 界面點的物理座標
    for i in range(min(3, grid.N_eta)):
        eta = eta_vals[i]
        nu = nu_vals[interface_idx]  # nu = π/2
        
        r = grid.f * np.sinh(eta) * np.sin(nu)
        z = grid.f * np.cosh(eta) * np.cos(nu)
        
        print(f"    界面點[{i},{interface_idx}]: eta={eta:.4f}, nu={nu:.4f}, r={r:.4f}nm, z={z:.4f}nm")
    
    # 分析4: 數值精度和收斂
    print(f"\n" + "="*60)
    print(f"分析4: 數值精度分析")
    print(f"="*60)
    
    # 測試更高精度求解
    potential_hires, iterations_hr, error_hr = solver.solve_laplace(
        V_tip, V_sample, max_iterations=1000, tolerance=1e-6)
    pot0_hires = solver._calculate_pot0_fortran_style(potential_hires)
    
    print(f"  高精度求解:")
    print(f"    迭代次數: {iterations_hr}")
    print(f"    最終誤差: {error_hr:.2e}")
    print(f"    Pot0: {pot0_hires:.6f} V")
    print(f"    與低精度差異: {abs(pot0_hires - pot0_laplace):.6f} V")
    print(f"    與Fortran差異: {abs(pot0_hires - fortran_target):.6f} V")
    
    # 分析5: 猜測可能的根本原因
    print(f"\n" + "="*60)
    print(f"分析5: 可能的根本原因")
    print(f"="*60)
    
    # 計算比例因子
    scale_factor = fortran_target / pot0_laplace
    
    print(f"  數值比較:")
    print(f"    Python結果: {pot0_laplace:.6f} V")
    print(f"    Fortran目標: {fortran_target:.6f} V")
    print(f"    比例因子: {scale_factor:.6f}")
    print(f"    差異比例: {abs(pot0_laplace - fortran_target)/abs(fortran_target):.1%}")
    
    print(f"\n  可能原因分析:")
    print(f"    1. 座標系統定義不同 (最可能)")
    print(f"    2. 邊界條件實現差異")
    print(f"    3. 網格間距或幾何參數不同")
    print(f"    4. VSINT陣列初始化或更新方式不同")
    print(f"    5. 單位轉換或物理常數差異")
    
    # 推薦解決方案
    print(f"\n" + "="*60)
    print(f"推薦解決方案")
    print(f"="*60)
    
    best_result = min([
        (abs(pcent_original - fortran_target), "原始PCENT"),
        (abs(pcent_var1 - fortran_target), "PCENT變體1"),
        (abs(pcent_var2 - fortran_target), "PCENT變體2"),
        (abs(pcent_with_offset2 - fortran_target), "縮放修正"),
        (abs(pot0_hires - fortran_target), "高精度求解")
    ])
    
    print(f"  最佳結果: {best_result[1]} (差異: {best_result[0]:.6f} V)")
    
    if best_result[0] < 0.1:
        print(f"  ✅ 找到了接近Fortran的方法!")
    elif scale_factor > 0.05 and scale_factor < 0.2:
        print(f"  💡 考慮使用縮放因子: {scale_factor:.4f}")
        print(f"  修正結果: {pot0_laplace * scale_factor:.6f} V")
    else:
        print(f"  ⚠️  需要更深入的分析")
    
    return {
        'pot0_laplace': pot0_laplace,
        'pot0_hires': pot0_hires,
        'scale_factor': scale_factor,
        'best_pcent_variant': best_result
    }

if __name__ == "__main__":
    results = analyze_remaining_difference()
    
    print(f"\n" + "="*80)
    print(f"總結:")
    print(f"="*80)
    
    fortran_target = -0.08
    print(f"與Fortran目標({fortran_target}V)的差異分析:")
    print(f"  基本Laplace: {abs(results['pot0_laplace'] - fortran_target):.6f} V")
    print(f"  高精度求解: {abs(results['pot0_hires'] - fortran_target):.6f} V")
    print(f"  最佳變體: {results['best_pcent_variant'][1]} ({results['best_pcent_variant'][0]:.6f} V)")
    
    if results['best_pcent_variant'][0] < 0.2:
        print(f"🎉 找到了顯著改善的方法！")
    else:
        print(f"🔧 需要繼續深入分析根本原因")