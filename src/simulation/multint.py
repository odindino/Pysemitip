import numpy as np
import logging
import time

# 【修改】導入新的核心模組
from src.physics.solvers.grid import HyperbolicGrid
from src.physics.core.poisson import PoissonSOREquation
# 舊的導入將不再使用： from src.physics.solvers.grid import CylindricalGrid
# 舊的導入將不再使用： from src.physics.core.poisson import PoissonSolver

from src.physics.core.charge_density import ChargeDensityCalculator

# ... (如果 MultInt 中有其他導入，保持不變)

logger = logging.getLogger(__name__)

class MultInt:
    """
    執行自洽計算以模擬 STM 與半導體樣品間的交互作用。
    【此版本已更新，以使用雙曲面網格和新的 SOR 解法器】
    """
    def __init__(self, config):
        self.config = config
        # self.props = SimulationProperties(config)
        
        # 【修改】建立新的雙曲面網格
        self.grid = self._create_grid()
        
        # 電荷密度計算模組保持不變，但現在會在雙曲面網格上運作
        self.charge_density = ChargeDensityCalculator(self.grid, self.props)

        # 【修改】初始化新的 SOR 解法器
        self.poisson_solver = PoissonSOREquation(self.grid, self.props)
        
        self.results = {}

    def _create_grid(self):
        """
        【修改】根據物理參數建立雙曲面網格。
        """
        grid_params = self.config.grid
        tip_params = self.config.tip
        
        logger.info("正在建立雙曲面網格...")
        try:
            grid = HyperbolicGrid(
                N_eta=grid_params.nx,           # 使用 nx 作為 eta 方向的點數
                N_nu=grid_params.nz,            # 使用 nz 作為 nu 方向的點數
                R=tip_params.radius,            # 傳入針尖曲率半徑
                Z_TS=tip_params.separation      # 傳入針尖-樣品距離
            )
            logger.info("雙曲面網格建立成功。")
            return grid
        except ValueError as e:
            logger.error(f"建立網格失敗: {e}")
            raise

    def mix_potential(self, old_potential, new_potential):
        """混合新舊電位以穩定收斂過程。"""
        alpha = self.config.convergence.mixing_alpha
        return (1 - alpha) * old_potential + alpha * new_potential

    def run_self_consistent_loop(self):
        """
        執行核心的自洽迴圈。
        """
        logger.info("開始執行自洽迴圈...")
        start_time = time.time()

        # 從設定檔獲取參數
        max_iterations = self.config.convergence.max_iterations
        tolerance = self.config.convergence.tolerance_potential
        V_tip = self.config.voltage_scan.initial_voltage
        V_sample = 0.0 # 假設樣品接地

        # --- 1. 初始解：求解拉普拉斯方程式 ---
        logger.info(f"正在求解初始電位 (拉普拉斯方程式) V_tip={V_tip}V...")
        try:
            # 【修改】呼叫新的 solve_laplace 方法
            potential, iters, err = self.poisson_solver.solve_laplace(
                V_tip=V_tip, 
                V_sample=V_sample,
                omega=self.config.sor.omega,
                tolerance=self.config.sor.tolerance
            )
            logger.info(f"拉普拉斯方程式求解完成，耗時 {iters} 次迭代。")
        except Exception as e:
            logger.error(f"拉普拉斯方程式求解失敗: {e}")
            raise

        # --- 2. 自洽迴圈 ---
        for i in range(max_iterations):
            logger.info(f"--- 自洽迭代: {i + 1}/{max_iterations} ---")
            
            # a. 根據當前電位，計算半導體中的電荷密度
            # 注意：charge_density 模組可能也需要適配雙曲面網格，但目前我們先假設其介面不變
            rho = self.charge_density.calculate(potential)
            
            # b. 求解帕松方程式，得到新的電位分佈
            # 【修改】呼叫新的 solve 方法，傳入電荷密度和當前電位作為初始猜測
            new_potential, sor_iters, sor_err = self.poisson_solver.solve(
                charge_density_rho=rho,
                potential_guess=potential, # 使用上一次的結果作為初始猜測
                omega=self.config.sor.omega,
                tolerance=self.config.sor.tolerance
            )
            logger.info(f"帕松方程式求解完成，耗時 {sor_iters} 次迭代，誤差 {sor_err:.2e}")

            # c. 檢查收斂性
            potential_diff = np.max(np.abs(new_potential - potential))
            logger.info(f"電位最大變化量: {potential_diff:.4e} V")
            if potential_diff < tolerance:
                logger.info(f"自洽迴圈在第 {i + 1} 次迭代後收斂！")
                self.results['potential'] = new_potential
                self.results['converged'] = True
                break

            # d. 混合電位，準備下一次迭代
            potential = self.mix_potential(potential, new_potential)
        else:
            # 如果 for 迴圈正常結束（未被 break）
            logger.warning("自洽迴圈達到最大迭代次數，未收斂。")
            self.results['potential'] = potential
            self.results['converged'] = False
        
        end_time = time.time()
        logger.info(f"自洽計算結束，總耗時: {end_time - start_time:.2f} 秒。")