import numpy as np
from scipy.optimize import brentq # Added for non-linear solver
from src.physics.solvers.grid import HyperbolicGrid
# Ensure ChargeDensityCalculator can be imported for type hinting and use
# This path might need adjustment based on actual project structure if run as a module
try:
    from .charge_density import ChargeDensityCalculator
except ImportError:
    # Fallback for cases where the script might be run in a different context or for testing
    from src.physics.core.charge_density import ChargeDensityCalculator

from src.utils.constants import E, EPSILON0

class PoissonSOREquation:
    """
    使用連續過鬆弛法 (SOR) 在雙曲面網格上求解帕松方程式。
    
    該實作旨在忠實復現原始 Fortran Semitip 的物理模型與演算法，
    並設計為可整合進 pysemitip 的自洽迴圈中。
    """

    def __init__(self, grid: HyperbolicGrid, props): # props is SemitipConfig
        """
        初始化解法器。

        Args:
            grid (HyperbolicGrid): 已經初始化的雙曲面網格物件。
            props (SimulationProperties): 包含所有物理參數的物件 (例如介電常數)。
                                         Actually, props is SemitipConfig
        """
        if not isinstance(grid, HyperbolicGrid):
            raise TypeError("grid 參數必須是 HyperbolicGrid 的一個實例。")
            
        self.grid = grid
        self.props = props # SemitipConfig
        
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
        
        # 計算幾何因子，不包含 epsilon
        # RHS term will be (rho / epsilon_eff) * geom_factor
        self.geom_factor = (self.grid.f**2) * (np.sinh(eta)**2 + np.sin(nu)**2)


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

    def solve(self, charge_density_rho_initial: np.ndarray, 
              potential_guess: np.ndarray,
              charge_density_calculator: 'ChargeDensityCalculator', 
              region_id_map: np.ndarray, # Map of region IDs for each grid point
              system_fermi_level_E_F_main: float, # System Fermi level in eV
              omega: float = 1.8, tolerance: float = 1e-7, 
              max_iterations: int = 20000) -> tuple:
        """
        求解帕松方程式 V'' = -rho/epsilon。
        For semiconductor regions, rho is a function of V, requiring a non-linear solve.
        For other regions, rho is taken from charge_density_rho_initial.

        Args:
            charge_density_rho_initial (np.ndarray): Initial charge density (C/m^3) based on potential_guess.
                                                     Used for non-semiconductor regions.
            potential_guess (np.ndarray): Initial potential guess (Volts).
            charge_density_calculator (ChargeDensityCalculator): Instance for rho(V) calculations.
            region_id_map (np.ndarray): Integer map of region ID for each grid point.
            system_fermi_level_E_F_main (float): System Fermi level (E_F_main) in eV.
            omega (float): SOR鬆弛因子。
            tolerance (float): 收斂容忍度。
            max_iterations (int): 最大迭代次數。

        Returns:
            tuple[np.ndarray, int, float]: (計算出的電位, 實際迭代次數, 最終誤差)
        """
        potential = np.copy(potential_guess)
        N_eta, N_nu = self.grid.N_eta, self.grid.N_nu

        # Determine local epsilon_r and absolute epsilon for each grid point
        epsilon_r_map = np.ones_like(region_id_map, dtype=float) # Default to vacuum (epsilon_r = 1)
        if charge_density_calculator.semiconductor_regions_params:
            for k, params in enumerate(charge_density_calculator.semiconductor_regions_params):
                mask = (region_id_map == k)
                epsilon_r_map[mask] = params['dielectric_constant']
        
        # For regions not in semiconductor_regions_params (e.g. id < 0 or > len-1), epsilon_r remains 1.
        # This assumes other material types (if any) also have epsilon_r=1 or are handled by region_id_map giving appropriate semiconductor_region_params indices.
        # If specific dielectric constants are needed for non-semiconductor regions, region_id_map and props.insulator/vacuum handling would need to provide them.
        # For now, any point not explicitly a semiconductor region gets epsilon_r=1.

        local_epsilon_abs_map = epsilon_r_map * EPSILON0

        for iteration in range(max_iterations):
            max_error = 0.0
            
            for i in range(1, N_eta - 1):
                for j in range(1, N_nu - 1):
                    old_phi = potential[i, j]
                    
                    neighbor_sum = (self.A_E[i,j] * potential[i+1, j] +
                                    self.A_W[i,j] * potential[i-1, j] +
                                    self.A_N[i,j] * potential[i, j+1] +
                                    self.A_S[i,j] * potential[i, j-1])
                    
                    current_region_id = int(region_id_map[i, j])
                    current_local_epsilon_abs = local_epsilon_abs_map[i,j]

                    is_semiconductor = (charge_density_calculator.semiconductor_regions_params and
                                       0 <= current_region_id < len(charge_density_calculator.semiconductor_regions_params))

                    if is_semiconductor:
                        region_params = charge_density_calculator.semiconductor_regions_params[current_region_id]

                        def nonlinear_eq_to_solve(phi_trial_volt):
                            # rho_C_per_cm3 from calculator (takes potential in eV, which is numerically phi_trial_volt)
                            rho_C_per_cm3 = charge_density_calculator.calculate_total_charge_density_at_point_direct(
                                current_region_id, phi_trial_volt, system_fermi_level_E_F_main
                            )
                            rho_C_per_m3 = rho_C_per_cm3 * 1e6
                            # Equation: phi_trial * C - neighbor_sum + rho(phi_trial) * geom_factor / epsilon_abs = 0
                            return (phi_trial_volt * self.C[i,j] - neighbor_sum +
                                    rho_C_per_m3 * self.geom_factor[i,j] / current_local_epsilon_abs)

                        # Define search bounds for brentq for phi_trial_volt
                        # Based on system Fermi level, region's DelVB, and table energy range
                        # V_local_eV = E_F_main - DelVB - E_f_prime_local_rel_VB
                        # Min V_local_eV (max E_f_prime)
                        v_brentq_min = system_fermi_level_E_F_main - region_params['del_vb'] - charge_density_calculator.table_energy_end_eV
                        # Max V_local_eV (min E_f_prime)
                        v_brentq_max = system_fermi_level_E_F_main - region_params['del_vb'] - charge_density_calculator.table_energy_start_eV
                        
                        # Ensure bounds are not identical and are reasonably ordered
                        if v_brentq_min >= v_brentq_max:
                             # Fallback if bounds are problematic, e.g. use old_phi +/- range
                            v_brentq_min = old_phi - 2.0 # Heuristic range
                            v_brentq_max = old_phi + 2.0
                        
                        # Ensure old_phi is not one of the bounds to avoid issues if root is at old_phi
                        if abs(old_phi - v_brentq_min) < 1e-9: v_brentq_min -= 1e-3
                        if abs(old_phi - v_brentq_max) < 1e-9: v_brentq_max += 1e-3
                        # Ensure min < max
                        if v_brentq_min > v_brentq_max: v_brentq_min, v_brentq_max = v_brentq_max, v_brentq_min


                        try:
                            val_at_vmin = nonlinear_eq_to_solve(v_brentq_min)
                            val_at_vmax = nonlinear_eq_to_solve(v_brentq_max)

                            if np.sign(val_at_vmin) == np.sign(val_at_vmax):
                                # Bracketing failed. Try a small interval around old_phi.
                                # This indicates that the root is not within the wide E_f_prime-derived bounds,
                                # or the function is too flat, or has multiple roots.
                                # A more robust strategy might be needed here, e.g. Newton step or smaller search.
                                # print(f"Warning: Brentq bracketing failed at ({i},{j}). Signs: {val_at_vmin:.2e}, {val_at_vmax:.2e}. Bounds: {v_brentq_min:.2f}, {v_brentq_max:.2f}. Using linear step.")
                                raise ValueError("Brentq bracketing failed") # Force fallback

                            solved_phi = brentq(nonlinear_eq_to_solve, v_brentq_min, v_brentq_max, 
                                                xtol=tolerance*0.01, rtol=tolerance*0.01) # Tighter tolerance for inner solve
                            new_phi_unrelaxed = solved_phi
                        except ValueError: # Handles brentq failure (e.g. bracketing, non-convergence)
                            # Fallback to linearized update using rho from old_phi
                            rho_C_per_cm3_old = charge_density_calculator.calculate_total_charge_density_at_point_direct(
                                current_region_id, old_phi, system_fermi_level_E_F_main
                            )
                            rho_C_per_m3_old = rho_C_per_cm3_old * 1e6
                            rhs_val_linearized = rho_C_per_m3_old * self.geom_factor[i,j] / current_local_epsilon_abs
                            new_phi_unrelaxed = (neighbor_sum - rhs_val_linearized) / self.C[i,j]
                            # if iteration % 100 == 0 and i == N_eta // 2 and j == N_nu // 2: # Debug print for a central point
                                # print(f"Iter {iteration}, Point ({i},{j}): Brentq fallback. OldPhi={old_phi:.3f}, NewPhiUnrelaxed={new_phi_unrelaxed:.3f}")


                    else: # Not a semiconductor region (e.g., vacuum, insulator)
                        # Use charge_density_rho_initial (C/m^3)
                        rhs_val = charge_density_rho_initial[i,j] * self.geom_factor[i,j] / current_local_epsilon_abs
                        new_phi_unrelaxed = (neighbor_sum - rhs_val) / self.C[i,j]
                    
                    new_phi = (1 - omega) * old_phi + omega * new_phi_unrelaxed
                    potential[i, j] = new_phi
                    
                    error = abs(new_phi - old_phi)
                    if error > max_error:
                        max_error = error
            
            if max_error < tolerance:
                return potential, iteration + 1, max_error
                
        return potential, max_iterations, max_error

    def solve_laplace(self, V_tip: float, V_sample: float,
                      charge_density_calculator: 'ChargeDensityCalculator', 
                      region_id_map: np.ndarray, 
                      system_fermi_level_E_F_main: float,
                      omega: float = 1.8, tolerance: float = 1e-7, 
                      max_iterations: int = 20000) -> tuple:
        """
        求解拉普拉斯方程式 V'' = 0，用於獲取初始電位。
        This now calls the main 'solve' method, passing necessary parameters.
        """
        zero_rho = np.zeros_like(self.potential) # Charge density is zero for Laplace
        
        initial_potential = np.linspace(V_tip, V_sample, self.grid.N_nu).reshape(1, -1)
        initial_potential = np.tile(initial_potential, (self.grid.N_eta, 1))
        self._apply_boundary_conditions(initial_potential, V_tip, V_sample)
        
        return self.solve(zero_rho, initial_potential, 
                          charge_density_calculator, region_id_map, system_fermi_level_E_F_main,
                          omega, tolerance, max_iterations)

# ... (rest of the file, if any) ...