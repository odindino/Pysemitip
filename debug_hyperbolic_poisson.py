import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace # 用於建立一個簡單的物件來模擬 props

from src.physics.solvers.grid import HyperbolicGrid
from src.physics.core.poisson import PoissonSOREquation

def run_and_plot_case(grid, props, V_tip, V_sample, omega, tolerance, max_iterations):
    """
    為一組給定的邊界條件，執行 SOR 求解並將結果視覺化。

    Args:
        grid (HyperbolicGrid): 網格物件。
        props (SimpleNamespace): 模擬的物理屬性物件。
        V_tip (float): 針尖電位。
        V_sample (float): 樣品電位。
        omega (float): SOR 鬆弛因子。
        tolerance (float): 收斂容忍度。
        max_iterations (int): 最大迭代次數。
    """
    print("-" * 50)
    print(f"開始新案例: 針尖電位 = {V_tip}V, 樣品電位 = {V_sample}V")
    print("-" * 50)

    try:
        # 1. 初始化解法器
        sor_solver = PoissonSOREquation(grid, props)
        
        # 2. 求解拉普拉斯方程式
        print("正在求解拉普拉斯方程式...")
        potential, iterations, error = sor_solver.solve_laplace(
            V_tip=V_tip, 
            V_sample=V_sample,
            omega=omega,
            tolerance=tolerance,
            max_iterations=max_iterations
        )

        print(f"SOR 解法器在 {iterations} 次迭代後完成。")
        print(f"最終誤差: {error:.2e}")

        # 3. 視覺化結果
        fig, ax = plt.subplots(figsize=(10, 8))
        
        contour = ax.contourf(grid.r, grid.z, potential, levels=50, cmap='viridis')
        fig.colorbar(contour, ax=ax, label='電位 (V)')
        
        ax.plot(grid.r, grid.z, 'k-', lw=0.3, alpha=0.5)
        ax.plot(grid.r.T, grid.z.T, 'k-', lw=0.3, alpha=0.5)
        
        ax.plot(grid.r[0, :], grid.z[0, :], 'r-', lw=2, label=f'針尖 (V={V_tip}V)')
        ax.plot(grid.r[:, -1], grid.z[:, -1], 'g-', lw=2, label=f'樣品 (V={V_sample}V)')
        
        title = f'拉普拉斯方程式解 (V_tip={V_tip}V, V_sample={V_sample}V)'
        ax.set_title(title)
        ax.set_xlabel('r (nm)')
        ax.set_ylabel('z (nm)')
        ax.set_aspect('equal', adjustable='box')
        ax.legend()
        ax.set_xlim(0, 3 * grid.R)
        ax.set_ylim(-0.5 * grid.Z_TS, 1.5 * grid.Z_TS)

        print(f"正在顯示圖形: {title}")
        plt.show()

    except Exception as e:
        print(f"案例執行過程中發生錯誤 (V_tip={V_tip}, V_sample={V_sample}): {e}")


def main():
    """
    主函式，設定共享資源並執行多個測試案例。
    """
    # --- 通用設定 ---
    # 1. 網格參數
    R = 5.0
    Z_TS = 10.0
    N_eta = 100
    N_nu = 150

    # 2. SOR 解法器參數
    omega = 1.8
    tolerance = 1e-6
    max_iterations = 20000

    # --- 建立共享資源 (只需一次) ---
    # 模擬一個最小化的 props 物件
    mock_props = SimpleNamespace(
        semiconductor=SimpleNamespace(
            epsilon_r=12.9  # GaAs
        )
    )
    # 建立網格
    try:
        grid = HyperbolicGrid(N_eta, N_nu, R, Z_TS)
        print("雙曲面網格建立成功。")
    except Exception as e:
        print(f"網格建立失敗: {e}")
        return

    # --- 執行測試案例 ---
    # 案例一：針尖 -1V，樣品 0V
    run_and_plot_case(grid, mock_props, V_tip=-1.0, V_sample=0.0,
                      omega=omega, tolerance=tolerance, max_iterations=max_iterations)

    # 【新增測試】案例二：針尖 0V，樣品 2V
    run_and_plot_case(grid, mock_props, V_tip=0.0, V_sample=2.0,
                      omega=omega, tolerance=tolerance, max_iterations=max_iterations)
    
    # 【新增測試】案例三：針尖 2V，樣品 2V，兩者等電位
    run_and_plot_case(grid, mock_props, V_tip=2.0, V_sample=2.0,
                      omega=omega, tolerance=tolerance, max_iterations=max_iterations)


if __name__ == '__main__':
    main()