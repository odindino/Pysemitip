import numpy as np
from src.physics.solvers.grid import HyperbolicGrid
# 為了獲取物理常數和半導體材料屬性，我們需要 SimulationProperties。
# 由於該檔案的確切路徑和名稱可能不同，這裡先註解掉，但在 multint.py 中傳入時是有效的。
# from src.simulation.properties import SimulationProperties 
from src.utils.constants import E, EPSILON0

class PoissonSOREquation:
    """
    使用連續過鬆弛法 (SOR) 在雙曲面網格上求解帕松方程式。
    
    該實作旨在忠實復現原始 Fortran Semitip 的物理模型與演算法，
    並設計為可整合進 pysemitip 的自洽迴圈中。
    """

    def __init__(self, grid: HyperbolicGrid, props):
        """
        初始化解法器。

        Args:
            grid (HyperbolicGrid): 已經初始化的雙曲面網格物件。
            props (SimulationProperties): 包含所有物理參數的物件 (例如介電常數)。
        """
        if not isinstance(grid, HyperbolicGrid):
            raise TypeError("grid 參數必須是 HyperbolicGrid 的一個實例。")
            
        self.grid = grid
        self.props = props
        
        self.potential = np.zeros((grid.N_eta, grid.N_nu))
        
        # 預先計算幾何係數，以加速迭代過程
        self._precompute_coefficients()
        
    def _precompute_coefficients(self):
        """
        根據網格幾何，預先計算 SOR 公式中所有與位置相關的係數。
        這些係數源於雙曲面座標系下拉普拉斯算子的有限差分法。
        """
        eta = self.grid.eta_grid + self.grid.eta_tip
        nu = self.grid.nu
        N_eta, N_nu = self.grid.N_eta, self.grid.N_nu

        # 檢查網格是否至少為 3x3，以便計算間距
        if N_eta < 3 or N_nu < 3:
            raise ValueError("網格維度必須至少為 3x3 才能進行 SOR 計算。")

        d_eta = eta[1, 0] - eta[0, 0]
        d_nu = nu[0, 1] - nu[0, 0]
        d_eta_sq = d_eta**2
        d_nu_sq = d_nu**2

        # 建立係數矩陣 (內部點)
        # 係數 A^E, A^W, A^N, A^S 分別對應東, 西, 北, 南四個相鄰點
        self.A_E = np.zeros_like(self.potential)
        self.A_W = np.zeros_like(self.potential)
        self.A_N = np.zeros_like(self.potential)
        self.A_S = np.zeros_like(self.potential)

        # 根據 "Scaling and Convergence of the Potential.pdf" 中的公式
        # 計算內部點 (1..N_eta-2, 1..N_nu-2) 的係數
        eta_int = eta[1:-1, 1:-1]
        nu_int = nu[1:-1, 1:-1]
        
        # sinh 和 sin 在半步長節點上的值
        sinh_eta_p_half = np.sinh(eta_int + d_eta / 2)
        sinh_eta_m_half = np.sinh(eta_int - d_eta / 2)
        sin_nu_p_half = np.sin(nu_int + d_nu / 2)
        sin_nu_m_half = np.sin(nu_int - d_nu / 2)

        self.A_E[1:-1, 1:-1] = (sinh_eta_p_half / np.sinh(eta_int)) / d_eta_sq
        self.A_W[1:-1, 1:-1] = (sinh_eta_m_half / np.sinh(eta_int)) / d_eta_sq
        self.A_N[1:-1, 1:-1] = (sin_nu_p_half / np.sin(nu_int)) / d_nu_sq
        self.A_S[1:-1, 1:-1] = (sin_nu_m_half / np.sin(nu_int)) / d_nu_sq

        # 中心點的係數 C 是四個方向係數之和
        self.C = self.A_E + self.A_W + self.A_N + self.A_S
        
        # 計算電荷項的係數
        # RHS = (rho / epsilon) * f^2 * (sinh^2(eta) + sin^2(nu))
        epsilon = self.props.semiconductor.epsilon_r * EPSILON0
        self.charge_factor = (self.grid.f**2 / epsilon) * \
                             (np.sinh(eta)**2 + np.sin(nu)**2)


    def _apply_boundary_conditions(self, potential, V_tip, V_sample):
        """在給定的電位矩陣上強制設定邊界條件。"""
        # 針尖邊界 (Dirichlet)
        potential[0, :] = V_tip
        
        # 樣品邊界 (Dirichlet)
        potential[:, -1] = V_sample
        
        # 遠端邊界 (eta_max, Neumann, 電場法向分量為0)
        potential[-1, :] = potential[-2, :]
        
        # 中心軸邊界 (nu=0, Neumann, 電場法向分量為0)
        potential[:, 0] = potential[:, 1]
        
        # 修正左上角和右下角的點，避免衝突
        potential[0, 0] = V_tip # 針尖尖端
        potential[-1, -1] = V_sample # 樣品遠處

    def solve(self, charge_density_rho: np.ndarray, potential_guess: np.ndarray,
              omega: float = 1.8, tolerance: float = 1e-7, max_iterations: int = 20000) -> tuple:
        """
        求解帕松方程式 V'' = -rho/epsilon。

        Args:
            charge_density_rho (np.ndarray): 與網格相同形狀的電荷密度陣列 (C/m^3)。
            potential_guess (np.ndarray): 用於開始迭代的初始電位猜測。
            omega (float): SOR 鬆弛因子。
            tolerance (float): 收斂容忍度。
            max_iterations (int): 最大迭代次數。

        Returns:
            tuple[np.ndarray, int, float]: (計算出的電位, 實際迭代次數, 最終誤差)
        """
        potential = np.copy(potential_guess)
        N_eta, N_nu = self.grid.N_eta, self.grid.N_nu
        
        # 計算右側項 (RHS), 即電荷項
        rhs = charge_density_rho * self.charge_factor

        for iteration in range(max_iterations):
            max_error = 0.0
            
            # 遍歷所有內部網格點
            for i in range(1, N_eta - 1):
                for j in range(1, N_nu - 1):
                    # 儲存舊值以計算誤差
                    old_phi = potential[i, j]
                    
                    # 根據五點差分公式計算新值
                    new_phi_unrelaxed = (self.A_E[i,j] * potential[i+1, j] +
                                         self.A_W[i,j] * potential[i-1, j] +
                                         self.A_N[i,j] * potential[i, j+1] +
                                         self.A_S[i,j] * potential[i, j-1] -
                                         rhs[i,j]) / self.C[i,j]
                    
                    # 應用 SOR
                    new_phi = (1 - omega) * old_phi + omega * new_phi_unrelaxed
                    
                    potential[i, j] = new_phi
                    
                    # 更新最大誤差
                    error = abs(new_phi - old_phi)
                    if error > max_error:
                        max_error = error
            
            # 檢查是否收斂
            if max_error < tolerance:
                # print(f"SOR 收斂於第 {iteration + 1} 次迭代。")
                return potential, iteration + 1, max_error
                
        # print(f"警告：SOR達到最大迭代次數 {max_iterations} 未收斂。")
        return potential, max_iterations, max_error

    def solve_laplace(self, V_tip: float, V_sample: float,
                      omega: float = 1.8, tolerance: float = 1e-7, max_iterations: int = 20000) -> tuple:
        """
        求解拉普拉斯方程式 V'' = 0，用於獲取初始電位。

        Args:
            V_tip (float): 針尖電位。
            V_sample (float): 樣品電位。
        
        Returns:
            tuple[np.ndarray, int, float]: (計算出的電位, 實際迭代次數, 最終誤差)
        """
        # 建立一個零電荷密度陣列
        zero_rho = np.zeros_like(self.potential)
        
        # 建立一個初始電位猜測（例如，從針尖到樣品的線性梯度）
        initial_potential = np.linspace(V_tip, V_sample, self.grid.N_nu).reshape(1, -1)
        initial_potential = np.tile(initial_potential, (self.grid.N_eta, 1))

        # 在初始猜測上應用嚴格的邊界條件
        self._apply_boundary_conditions(initial_potential, V_tip, V_sample)
        
        # 呼叫主求解器
        return self.solve(zero_rho, initial_potential, omega, tolerance, max_iterations)