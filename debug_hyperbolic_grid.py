import numpy as np
import matplotlib.pyplot as plt
from src.physics.solvers.grid import HyperbolicGrid

def main():
    """
    主函式，用於生成並視覺化一個雙曲面網格。
    """
    # --- 設定 ---
    # 【重要修正】調整參數以滿足新模型的物理約束 (Z_TS > R)
    # 原始 Fortran 程式碼通常在 Z_TS 較大的情況下運作。
    R = 5.0          # 針尖曲率半徑 (nm)
    Z_TS = 10.0      # 針尖-樣品距離 (nm)
    
    # 網格維度 (點數)
    N_eta = 100
    N_nu = 150

    # --- 網格生成 ---
    try:
        # 實例化網格
        print(f"正在生成網格，參數: R={R} nm, Z_TS={Z_TS} nm")
        hyperbolic_grid = HyperbolicGrid(N_eta, N_nu, R, Z_TS)
        
        # 獲取座標陣列
        # 【次要修正】從 .eta 改為 .eta_grid 以匹配新類別的屬性名
        eta_grid, nu = hyperbolic_grid.eta_grid, hyperbolic_grid.nu
        r, z = hyperbolic_grid.r, hyperbolic_grid.z

        print("成功生成雙曲面網格。")
        print(f"網格維度 (eta_grid x nu): {eta_grid.shape}")
        print(f"笛卡爾座標維度 (r x z): {r.shape}")

    except ValueError as ve:
        print(f"網格生成失敗，發生數值或參數錯誤: {ve}")
        return
    except Exception as e:
        print(f"網格生成過程中發生未知錯誤: {e}")
        return

    # --- 視覺化 ---
    # 這部分的程式碼無需修改，它的邏輯是正確的
    try:
        fig, ax = plt.subplots(figsize=(10, 8))

        # 繪製恆定 nu 和恆定 eta 的網格線
        ax.plot(r, z, 'b-', lw=0.5)
        ax.plot(r.T, z.T, 'b-', lw=0.5)

        # 高亮針尖表面 (eta_grid = 0)
        ax.plot(r[0, :], z[0, :], 'r-', lw=2, label=f'針尖表面 (η_grid=0, R={R}nm)')
        
        # 高亮樣品表面 (nu = pi/2)
        ax.plot(r[:, -1], z[:, -1], 'g-', lw=2, label=f'樣品表面 (ν=π/2, z=0)')
        
        ax.set_title('雙曲面座標網格視覺化')
        ax.set_xlabel('r (nm)')
        ax.set_ylabel('z (nm)')
        ax.set_aspect('equal', adjustable='box')
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        
        # 設定繪圖範圍，更好地聚焦於針尖-樣品區域
        ax.set_xlim(0, 3 * R)
        ax.set_ylim(-0.5 * Z_TS, 1.5 * Z_TS)

        print("正在顯示圖形...")
        plt.show()

    except Exception as e:
        print(f"繪圖過程中發生錯誤: {e}")

if __name__ == '__main__':
    main()