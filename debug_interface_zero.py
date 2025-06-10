#!/usr/bin/env python3
"""
調試界面電位為 0 的問題
"""
import numpy as np
import logging
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

logging.basicConfig(level=logging.INFO, format='%(message)s')

def debug_interface_zero():
    """調試為什麼界面電位是 0"""
    print("調試界面電位為 0 的問題")
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
    
    V_tip = -2.07
    V_sample = 0.0
    
    print(f"Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    print(f"V_tip={V_tip}, V_sample={V_sample}")
    
    # 步驟1: 檢查初始猜測
    print(f"\n" + "="*50)
    print(f"步驟1: 檢查初始電位猜測")
    print(f"="*50)
    
    initial_potential = solver._create_initial_potential_guess(V_tip, V_sample)
    interface_idx = grid.N_nu - 1
    
    print(f"初始猜測 - 界面電位 (nu={interface_idx}):")
    for i in range(min(4, grid.N_eta)):
        v = initial_potential[i, interface_idx]
        print(f"  eta={i}: V = {v:.6f} V")
    
    print(f"初始猜測 - 其他重要點:")
    print(f"  針尖頂點 [0,0]: {initial_potential[0,0]:.6f} V")
    print(f"  中心軸 [0,0]: {initial_potential[0,0]:.6f} V")
    
    # 步驟2: 檢查邊界條件應用
    print(f"\n" + "="*50)
    print(f"步驟2: 檢查邊界條件應用")
    print(f"="*50)
    
    potential_bc = initial_potential.copy()
    potential_bc = solver._apply_boundary_conditions(potential_bc, V_tip, V_sample)
    
    print(f"應用邊界條件後 - 界面電位:")
    for i in range(min(4, grid.N_eta)):
        v_before = initial_potential[i, interface_idx]
        v_after = potential_bc[i, interface_idx]
        print(f"  eta={i}: {v_before:.6f} → {v_after:.6f} V (變化: {v_after-v_before:.6f})")
    
    print(f"應用邊界條件後 - 檢查針尖邊界:")
    for j in range(min(4, grid.N_nu)):
        v = potential_bc[0, j]
        should_be_tip = j < grid.N_nu - 1  # 除了界面點
        print(f"  [0,{j}]: {v:.6f} V ({'應該=V_tip' if should_be_tip else '應該自由'})")
    
    # 步驟3: 檢查界面更新
    print(f"\n" + "="*50)
    print(f"步驟3: 檢查界面電位更新")
    print(f"="*50)
    
    potential_updated = potential_bc.copy()
    potential_updated = solver._update_interface_potential(potential_updated)
    
    print(f"界面更新後:")
    for i in range(min(4, grid.N_eta)):
        v_before = potential_bc[i, interface_idx]
        v_after = potential_updated[i, interface_idx]
        print(f"  eta={i}: {v_before:.6f} → {v_after:.6f} V (變化: {v_after-v_before:.6f})")
    
    # 步驟4: 檢查完整 Laplace 求解
    print(f"\n" + "="*50)
    print(f"步驟4: 檢查完整 Laplace 求解")
    print(f"="*50)
    
    # 運行短時間求解
    potential_solved, iterations, error = solver.solve_laplace(
        V_tip, V_sample, max_iterations=100, tolerance=1e-2
    )
    
    print(f"Laplace 求解 ({iterations} 次迭代):")
    for i in range(min(4, grid.N_eta)):
        v = potential_solved[i, interface_idx]
        print(f"  eta={i}: V = {v:.6f} V")
    
    # 步驟5: 檢查是否有特殊的邊界處理
    print(f"\n" + "="*50)
    print(f"步驟5: 檢查中心軸邊界條件")
    print(f"="*50)
    
    print(f"中心軸 (nu=0) 應該使用 Neumann 邊界:")
    for i in range(min(4, grid.N_eta)):
        v_center = potential_solved[i, 0]
        v_next = potential_solved[i, 1]
        print(f"  eta={i}: V[nu=0]={v_center:.6f}, V[nu=1]={v_next:.6f}")
    
    # 步驟6: 分析問題
    print(f"\n" + "="*50)
    print(f"步驟6: 問題分析")
    print(f"="*50)
    
    interface_center = potential_solved[0, interface_idx]
    
    if abs(interface_center) < 1e-6:
        print(f"❌ 界面中心電位確實是 0")
        print(f"可能原因:")
        print(f"  1. 初始猜測就設為 0")
        print(f"  2. 邊界條件強制設為 0")
        print(f"  3. 對稱性導致的數值結果")
        print(f"  4. 座標系統原點定義問題")
    else:
        print(f"✅ 界面中心電位不是 0: {interface_center:.6f} V")
    
    # 步驟7: 檢查座標系統
    print(f"\n" + "="*50)
    print(f"步驟7: 檢查座標系統定義")
    print(f"="*50)
    
    print(f"檢查物理座標:")
    # 創建座標陣列
    eta_vals = np.linspace(grid.eta_tip, 2.0, grid.N_eta)
    nu_vals = np.linspace(0, np.pi/2, grid.N_nu)
    
    # 檢查界面點的物理位置
    for i in range(min(3, grid.N_eta)):
        eta = eta_vals[i]
        nu = nu_vals[interface_idx]  # nu = π/2
        
        # 計算物理座標
        r = grid.f * np.sinh(eta) * np.sin(nu)
        z = grid.f * np.cosh(eta) * np.cos(nu)
        
        print(f"  eta_idx={i} (eta={eta:.4f}): r={r:.4f}, z={z:.4f} nm")
    
    print(f"\n在 nu=π/2 時，cos(π/2)=0，所以 z=0")
    print(f"這確實是半導體表面！")
    
    return {
        'initial': initial_potential,
        'after_bc': potential_bc, 
        'after_update': potential_updated,
        'final': potential_solved,
        'interface_idx': interface_idx
    }

if __name__ == "__main__":
    results = debug_interface_zero()
    
    print(f"\n" + "="*80)
    print(f"總結:")
    print(f"="*80)
    
    interface_idx = results['interface_idx']
    final_interface = results['final'][0, interface_idx]
    
    print(f"最終界面電位: {final_interface:.6f} V")
    
    if abs(final_interface) < 1e-6:
        print(f"❌ 確實存在界面電位為 0 的問題")
        print(f"主要原因可能是座標系統或邊界條件設定")
    else:
        print(f"✅ 界面電位正常")