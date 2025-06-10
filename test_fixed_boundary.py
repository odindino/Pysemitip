#!/usr/bin/env python3
"""
測試修正的邊界條件
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

logging.basicConfig(level=logging.INFO, format='%(message)s')

def test_fixed_boundary():
    """測試修正後的邊界條件是否允許界面電位演化"""
    print("測試修正的邊界條件...")
    print("="*60)
    
    # 創建網格
    grid = HyperbolicGrid(N_eta=8, N_nu=6, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    # 創建物理參數
    class MinimalProps:
        def __init__(self):
            class SemiconductorProps:
                epsilon_r = 12.9
                Ev_offset_eV = -5.17
            self.semiconductor_props = SemiconductorProps()
    
    props = MinimalProps()
    solver = PoissonSOREquation(grid, props)
    
    # 測試參數
    V_tip = -2.07
    V_sample = 0.0
    
    # 創建初始電位
    potential = solver._create_initial_potential_guess(V_tip, V_sample)
    
    print(f"初始電位 (應該是線性插值):")
    print(f"  [0,0] = {potential[0,0]:.4f} V (針尖頂點)")
    print(f"  [0,5] = {potential[0,5]:.4f} V (界面中心，應該不等於 V_tip)")
    
    # 應用邊界條件
    potential_with_bc = solver._apply_boundary_conditions(potential.copy(), V_tip, V_sample)
    
    print(f"\n應用邊界條件後:")
    print(f"  [0,0] = {potential_with_bc[0,0]:.4f} V (應該 = {V_tip:.4f})")
    print(f"  [0,5] = {potential_with_bc[0,5]:.4f} V (應該保持不變)")
    
    # 檢查是否正確
    if abs(potential_with_bc[0,0] - V_tip) < 1e-10:
        print("✅ 針尖頂點正確設置為 V_tip")
    else:
        print("❌ 針尖頂點設置錯誤")
        
    if abs(potential_with_bc[0,5] - potential[0,5]) < 1e-10:
        print("✅ 界面中心保持原值，可以自由演化")
    else:
        print("❌ 界面中心被錯誤覆蓋")
    
    # 運行短時間的 Laplace 求解
    print(f"\n運行 Laplace 求解測試...")
    potential_solved, iterations, error = solver.solve_laplace(
        V_tip, V_sample, max_iterations=200, tolerance=1e-3
    )
    
    print(f"\nLaplace 求解後 (迭代 {iterations} 次):")
    print(f"  [0,0] = {potential_solved[0,0]:.4f} V")
    print(f"  [0,5] = {potential_solved[0,5]:.4f} V (應該已經演化)")
    
    # 計算 Pot0
    pot0 = solver._calculate_pot0_fortran_style(potential_solved)
    print(f"\nPot0 = {pot0:.4f} V")
    print(f"與 Fortran 差異: {abs(pot0 - (-0.08)):.4f} V")
    
    # 檢查電位分布合理性
    print(f"\n電位分布檢查:")
    for i in range(min(4, grid.N_eta)):
        v_interface = potential_solved[i, grid.N_nu-1]
        print(f"  eta={i}, nu=5: V = {v_interface:.4f} V")
    
    return pot0

if __name__ == "__main__":
    pot0 = test_fixed_boundary()
    print(f"\n{'='*60}")
    if abs(pot0 - (-0.08)) < 0.5:  # 在 0.5V 內算接近
        print("✅ Pot0 值更接近 Fortran 結果")
    else:
        print("❌ Pot0 值仍然與 Fortran 差異較大")
    print(f"{'='*60}")