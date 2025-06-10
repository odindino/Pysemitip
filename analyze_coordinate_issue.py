#!/usr/bin/env python3
"""
分析座標系統問題和 Pot0 計算
"""
import numpy as np
from src.physics.solvers.grid import HyperbolicGrid

def analyze_coordinate_mapping():
    """分析雙曲座標系統映射"""
    print("分析雙曲座標系統...")
    print("="*60)
    
    # 創建網格
    grid = HyperbolicGrid(N_eta=8, N_nu=6, R=1.0, Z_TS=1.0, shank_slope=1.0)
    
    print(f"網格參數:")
    print(f"  N_eta = {grid.N_eta}, N_nu = {grid.N_nu}")
    print(f"  R = {grid.R} nm, Z_TS = {grid.Z_TS} nm") 
    print(f"  f = {grid.f:.4f} nm, eta_tip = {grid.eta_tip:.4f}")
    
    # 檢查關鍵點的物理位置
    print(f"\n關鍵點的物理座標:")
    
    # 創建 eta 和 nu 陣列
    eta_vals = np.linspace(grid.eta_tip, 2.0, grid.N_eta)  # 從 eta_tip 到遠場
    nu_vals = np.linspace(0, np.pi/2, grid.N_nu)  # 從 0 到 π/2
    
    # 1. 針尖頂點
    eta_idx, nu_idx = 0, 0
    eta_val = eta_vals[eta_idx]
    nu_val = nu_vals[nu_idx]
    r_tip = grid.f * np.sinh(eta_val) * np.sin(nu_val)
    z_tip = grid.f * np.cosh(eta_val) * np.cos(nu_val)
    print(f"1. 針尖頂點 (eta_idx=0, nu_idx=0): r={r_tip:.4f}, z={z_tip:.4f} nm")
    
    # 2. 樣品表面中心
    eta_idx, nu_idx = 0, grid.N_nu-1
    eta_val = eta_vals[eta_idx]
    nu_val = nu_vals[nu_idx]
    r_center = grid.f * np.sinh(eta_val) * np.sin(nu_val)
    z_center = grid.f * np.cosh(eta_val) * np.cos(nu_val)
    print(f"2. 樣品表面中心 (eta_idx=0, nu_idx={grid.N_nu-1}): r={r_center:.4f}, z={z_center:.4f} nm")
    
    # 3. 檢查 nu 方向
    print(f"\n沿 eta_idx=0 的 nu 方向 (從針尖到樣品):")
    for j in range(grid.N_nu):
        eta_val = eta_vals[0]
        nu_val = nu_vals[j]
        r = grid.f * np.sinh(eta_val) * np.sin(nu_val)
        z = grid.f * np.cosh(eta_val) * np.cos(nu_val)
        print(f"  nu_idx={j}: r={r:.4f}, z={z:.4f} nm")
    
    # 4. 檢查界面位置
    print(f"\n界面位置分析 (nu_idx={grid.N_nu-1}):")
    for i in range(min(4, grid.N_eta)):
        eta_val = eta_vals[i]
        nu_val = nu_vals[grid.N_nu-1]
        r = grid.f * np.sinh(eta_val) * np.sin(nu_val)
        z = grid.f * np.cosh(eta_val) * np.cos(nu_val)
        print(f"  eta_idx={i}: r={r:.4f}, z={z:.4f} nm")
    
    # 5. 關鍵發現
    print(f"\n關鍵發現:")
    print(f"1. eta={eta_vals[0]:.4f} 對應的是一條線，不是一個點")
    print(f"2. 在 eta_idx=0, nu_idx={grid.N_nu-1} 處，z = {z_center:.4f} nm")
    print(f"3. 這個點在樣品表面 (z=0) {'之上' if z_center > 0 else '之下'}")
    
    # 6. 真正的樣品表面在哪？
    print(f"\n尋找真正的樣品表面 (z ≈ 0):")
    # 對於每個 eta，找到 z 最接近 0 的 nu
    for i in range(min(4, grid.N_eta)):
        min_z = float('inf')
        best_j = -1
        for j in range(grid.N_nu):
            eta_val = eta_vals[i]
            nu_val = nu_vals[j]
            z = grid.f * np.cosh(eta_val) * np.cos(nu_val)
            if abs(z) < abs(min_z):
                min_z = z
                best_j = j
        print(f"  eta_idx={i}: 最接近 z=0 的是 nu_idx={best_j}, z={min_z:.4f} nm")
    
    # 7. Fortran 的 VSINT 可能對應什麼？
    print(f"\nFortran VSINT 分析:")
    print(f"VSINT 可能對應的是真正的半導體/真空界面")
    print(f"而不是簡單的 nu = N_nu-1")
    print(f"這解釋了為什麼 Pot0 差異如此大")
    
    return grid

if __name__ == "__main__":
    grid = analyze_coordinate_mapping()
    print(f"\n{'='*60}")
    print("結論：座標映射可能是 Pot0 差異的主要原因")
    print("需要確認 Fortran 中 VSINT 的確切物理位置")
    print(f"{'='*60}")