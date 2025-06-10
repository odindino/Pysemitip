#!/usr/bin/env python3
"""
診斷 Pot0 計算問題，分析不同位置的電位值
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

# 設置日誌
logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')
logger = logging.getLogger(__name__)

def diagnose_pot0_calculation():
    """診斷 Pot0 計算的各個環節"""
    print("Diagnosing Pot0 calculation...")
    print("="*60)
    
    # 創建測試網格
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
    V_tip = -2.07  # 使用 Fortran 的調整後偏壓
    V_sample = 0.0
    
    print(f"Test conditions: V_tip={V_tip}V, V_sample={V_sample}V")
    print(f"Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    
    # 運行 Laplace 方程獲得初始電位分布
    potential, iterations, max_error = solver.solve_laplace(
        V_tip, V_sample, max_iterations=500, tolerance=1e-4
    )
    
    print(f"\nLaplace solution: {iterations} iterations, error={max_error:.3e}")
    
    # 分析不同位置的電位
    print("\n電位分布分析:")
    print("-"*60)
    
    # 1. 針尖處 (eta=0)
    print(f"1. 針尖電位 (eta=0):")
    for j in range(min(4, grid.N_nu)):
        print(f"   nu={j}: V = {potential[0, j]:.6f} V")
    
    # 2. 樣品表面 (nu=N_nu-1)
    print(f"\n2. 樣品表面電位 (nu={grid.N_nu-1}):")
    for i in range(min(4, grid.N_eta)):
        print(f"   eta={i}: V = {potential[i, grid.N_nu-1]:.6f} V")
    
    # 3. 計算不同版本的 Pot0
    print("\n3. 不同方法計算的 Pot0:")
    
    # 方法1: 使用修正後的 PCENT 方法
    pot0_pcent = solver._calculate_pot0_fortran_style(potential)
    print(f"   PCENT method: Pot0 = {pot0_pcent:.6f} V")
    
    # 方法2: 直接使用界面中心點
    interface_center = potential[0, grid.N_nu-1]
    print(f"   Interface center: V = {interface_center:.6f} V")
    
    # 方法3: 界面平均值
    interface_avg = np.mean(potential[:, grid.N_nu-1])
    print(f"   Interface average: V = {interface_avg:.6f} V")
    
    # 方法4: 使用 Fortran 公式的詳細計算
    if grid.N_eta >= 2:
        v1 = potential[0, grid.N_nu-1]
        v2 = potential[1, grid.N_nu-1]
        weighted = (9.0 * v1 - v2) / 8.0
        print(f"   Weighted formula: (9*{v1:.4f} - {v2:.4f})/8 = {weighted:.6f} V")
    
    # 4. 檢查電位梯度
    print("\n4. 電位梯度分析:")
    # 徑向梯度
    if grid.N_eta >= 2:
        grad_eta = potential[1, grid.N_nu-1] - potential[0, grid.N_nu-1]
        print(f"   徑向梯度 (at interface): {grad_eta:.6f} V/step")
    
    # 角度梯度
    if grid.N_nu >= 2:
        grad_nu = potential[0, grid.N_nu-1] - potential[0, grid.N_nu-2]
        print(f"   角度梯度 (at eta=0): {grad_nu:.6f} V/step")
    
    # 5. 邊界條件檢查
    print("\n5. 邊界條件檢查:")
    print(f"   Tip (eta=0, all nu): {potential[0, 0]:.6f} V (應該 = {V_tip:.6f} V)")
    print(f"   Sample nominal: {V_sample:.6f} V")
    print(f"   Far field (eta={grid.N_eta-1}): {potential[grid.N_eta-1, grid.N_nu-1]:.6f} V")
    
    # 6. 物理參數影響
    print("\n6. 關鍵差異分析:")
    print(f"   Fortran Pot0: ~-0.08 V")
    print(f"   Python Pot0: {pot0_pcent:.6f} V")
    print(f"   差異: {abs(pot0_pcent - (-0.08)):.6f} V")
    
    # 可能的原因
    print("\n可能的差異原因:")
    print("   1. VSINT 在 Fortran 中是獨立計算和存儲的")
    print("   2. 界面電位可能需要特殊的邊界條件處理")
    print("   3. 電荷密度對界面電位的影響可能不同")
    print("   4. 初始猜測或收斂條件的差異")
    
    return potential, pot0_pcent

if __name__ == "__main__":
    potential, pot0 = diagnose_pot0_calculation()
    print(f"\n{'='*60}")
    print(f"診斷完成 - Pot0 = {pot0:.6f} V")
    print(f"{'='*60}")