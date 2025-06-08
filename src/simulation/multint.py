import numpy as np
import logging
import time

# 【修改】導入新的核心模組
from src.physics.solvers.grid import HyperbolicGrid
from src.physics.core.poisson import PoissonSOREquation
# 舊的導入將不再使用： from src.physics.solvers.grid import CylindricalGrid
# 舊的導入將不再使用： from src.physics.core.poisson import PoissonSolver

from src.physics.core.charge_density import ChargeDensityCalculator
from src.physics.materials.semiconductor import SemiconductorRegion
from src.physics.materials.surface_states import SurfaceRegion
from dataclasses import dataclass
from typing import List

# ... (如果 MultInt 中有其他導入，保持不變)

logger = logging.getLogger(__name__)


@dataclass
class SimulationProperties:
    """Simple container for simulation properties needed by physics modules."""
    semiconductor: 'SemiconductorConfig'
    surface: 'SurfaceConfig'
    tip: 'TipConfig'
    
    @property
    def fermi_level(self):
        """Get Fermi level from first semiconductor region."""
        if self.semiconductor.regions:
            return self.semiconductor.regions[0].fermi_level()
        return 0.0

class MultInt:
    """
    執行自洽計算以模擬 STM 與半導體樣品間的交互作用。
    【此版本已更新，以使用雙曲面網格和新的 SOR 解法器】
    """
    def __init__(self, config):
        self.config = config
        
        # Create props object with necessary attributes
        self.props = self._create_simulation_properties(config)
        
        # 【修改】建立新的雙曲面網格
        self.grid = self._create_grid()
        
        # 電荷密度計算模組保持不變，但現在會在雙曲面網格上運作
        self.charge_density = ChargeDensityCalculator(self.grid, self.props)

        # 【修改】初始化新的 SOR 解法器
        self.poisson_solver = PoissonSOREquation(self.grid, self.props)
        
        self.results = {}

    def _create_simulation_properties(self, config):
        """
        Create a SimulationProperties object from config.
        
        This converts the config data classes into the physics material classes
        needed by the physics modules.
        """
        # Convert semiconductor regions from config to physics objects
        semiconductor_regions = []
        for region_config in config.semiconductor.regions:
            semiconductor_region = SemiconductorRegion(
                region_id=region_config.id,
                donor_concentration=region_config.donor_concentration,
                acceptor_concentration=region_config.acceptor_concentration,
                band_gap=region_config.band_gap,
                valence_band_offset=region_config.valence_band_offset,
                electron_affinity=region_config.affinity,
                donor_binding_energy=region_config.donor_binding_energy,
                acceptor_binding_energy=region_config.acceptor_binding_energy,
                cb_effective_mass=region_config.effective_mass.conduction_band,
                vb_effective_mass_heavy=region_config.effective_mass.valence_band_heavy,
                vb_effective_mass_light=region_config.effective_mass.valence_band_light,
                vb_effective_mass_so=region_config.effective_mass.split_off,
                spin_orbit_splitting=region_config.spin_orbit_splitting,
                permittivity=region_config.permittivity,
                allow_degeneracy=region_config.allow_degeneracy,
                allow_inversion=region_config.allow_inversion,
                temperature=config.environment.temperature
            )
            semiconductor_regions.append(semiconductor_region)
        
        # Convert surface regions from config to physics objects
        surface_regions = []
        for region_config in config.surface.regions:
            # Get position from config if available, otherwise default to (0, 0)
            position = (getattr(region_config, 'x_position', 0.0), 
                       getattr(region_config, 'y_position', 0.0))
            
            # Convert distributions from config to physics objects
            from src.physics.materials.surface_states import SurfaceStateDistribution
            
            dist1 = None
            if region_config.first_distribution:
                d = region_config.first_distribution
                dist1 = SurfaceStateDistribution(
                    density=d.density,
                    neutrality_level=d.neutrality_level,
                    fwhm=d.fwhm,
                    center_energy=d.center_energy
                )
            
            dist2 = None  
            if region_config.second_distribution:
                d = region_config.second_distribution
                dist2 = SurfaceStateDistribution(
                    density=d.density,
                    neutrality_level=d.neutrality_level,
                    fwhm=d.fwhm,
                    center_energy=d.center_energy
                )
            
            surface_region = SurfaceRegion(
                region_id=region_config.id,
                position=position,
                distribution1=dist1,
                distribution2=dist2,
                temperature=config.environment.temperature
            )
            surface_regions.append(surface_region)
        
        # Create the properties object with proper structure
        # We'll create a simple object that mimics the config structure
        # but contains the physics objects
        class Props:
            def __init__(self):
                self.semiconductor = type('obj', (object,), {
                    'regions': semiconductor_regions,
                    'epsilon_r': config.environment.dielectric_constant,
                    'fermi_level': semiconductor_regions[0].fermi_level() if semiconductor_regions else 0.0
                })()
                self.surface = type('obj', (object,), {
                    'regions': surface_regions
                })()
                self.tip = config.tip
        
        return Props()
    
    def _create_grid(self):
        """
        【修改】根據物理參數建立雙曲面網格。
        """
        grid_params = self.config.grid
        tip_params = self.config.tip
        
        logger.info("正在建立雙曲面網格...")
        try:
            # Ensure separation > radius for hyperbolic grid model
            separation = tip_params.separation
            if separation <= tip_params.radius:
                logger.warning(f"調整針尖-樣品距離從 {separation} nm 到 {tip_params.radius * 1.1} nm 以滿足雙曲面網格要求")
                separation = tip_params.radius * 1.1
                
            grid = HyperbolicGrid(
                N_eta=grid_params.radial_points,      # 使用 radial_points 作為 eta 方向的點數
                N_nu=grid_params.angular_points,      # 使用 angular_points 作為 nu 方向的點數
                R=tip_params.radius,                  # 傳入針尖曲率半徑
                Z_TS=separation                       # 傳入針尖-樣品距離
            )
            logger.info("雙曲面網格建立成功。")
            return grid
        except ValueError as e:
            logger.error(f"建立網格失敗: {e}")
            raise

    def mix_potential(self, old_potential, new_potential):
        """混合新舊電位以穩定收斂過程。"""
        # Use default mixing parameter if not in config
        alpha = getattr(getattr(self.config, 'convergence', None), 'mixing_alpha', 0.3)
        return (1 - alpha) * old_potential + alpha * new_potential

    def run_self_consistent_loop(self):
        """
        執行核心的自洽迴圈。
        """
        logger.info("開始執行自洽迴圈...")
        start_time = time.time()

        # 從設定檔獲取參數
        # Use computation.max_iterations if available, otherwise default
        max_iterations = self.config.computation.max_iterations[0] if hasattr(self.config, 'computation') else 10000
        tolerance = self.config.computation.convergence_parameters[0] if hasattr(self.config, 'computation') else 1e-3
        V_tip = self.config.voltage_scan.start if hasattr(self.config, 'voltage_scan') else 0.0
        V_sample = 0.0 # 假設樣品接地

        # --- 1. 初始解：求解拉普拉斯方程式 ---
        logger.info(f"正在求解初始電位 (拉普拉斯方程式) V_tip={V_tip}V...")
        try:
            # 【修改】呼叫新的 solve_laplace 方法
            # Use default SOR parameters if not in config
            omega = 1.5  # Default SOR relaxation parameter
            sor_tolerance = 1e-6  # Default SOR tolerance
            
            potential, iters, err = self.poisson_solver.solve_laplace(
                V_tip=V_tip, 
                V_sample=V_sample,
                omega=omega,
                tolerance=sor_tolerance
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
                omega=omega,
                tolerance=sor_tolerance
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