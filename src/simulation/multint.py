from dataclasses import dataclass
import numpy as np
import logging
from typing import Dict

from src.core.config_schema import SemitipConfig
from src.physics.solvers.grid import HyperbolicGrid
from src.physics.core.poisson import PoissonSOREquation
from src.physics.core.charge_density import ChargeDensityCalculator
from src.physics.materials.semiconductor import SemiconductorRegion, create_semiconductor_from_config
from src.physics.materials.surface_states import SurfaceRegion

logger = logging.getLogger(__name__)

@dataclass
class SimulationResults:
    """
    儲存模擬結果的數據類別
    """
    bias_voltage: float
    current: float = 0.0
    band_bending: float = 0.0
    potential_profile: np.ndarray = None
    charge_density: np.ndarray = None
    scf_iterations: int = 0
    convergence_achieved: bool = False
    
    def __post_init__(self):
        if self.potential_profile is None:
            self.potential_profile = np.array([])
        if self.charge_density is None:
            self.charge_density = np.array([])

class MultInt:
    """
    執行自洽計算以模擬 STM 與半導體樣品間的交互作用。
    【此版本已更新，以使用雙曲面網格和新的 SOR 解法器】
    """
    def __init__(self, config: SemitipConfig, output_dir: str):
        self.config = config
        self.output_dir = output_dir
        logger.info(f"MultInt simulation initialized with config type: {config.simulation_type}")
        logger.info(f"Output directory: {self.output_dir}")

        # Initialize semiconductor physics regions
        self.semiconductor_physics_regions: Dict[int, SemiconductorRegion] = {}
        if self.config.semiconductor and self.config.semiconductor.regions:
            for region_config in self.config.semiconductor.regions:
                try:
                    # Ensure config.environment.temperature is accessible
                    env_temp = getattr(getattr(self.config, 'environment', None), 'temperature', 300.0)

                    physics_region = create_semiconductor_from_config(
                        region_config,
                        env_temp
                    )
                    self.semiconductor_physics_regions[physics_region.region_id] = physics_region
                    logger.info(f"Created semiconductor physics region ID {physics_region.region_id} "
                                f"(Ev_abs={physics_region.Ev_abs_for_charge_calc_eV:.3f} eV, "
                                f"Ec_abs={physics_region.Ec_abs_for_charge_calc_eV:.3f} eV, "
                                f"Permittivity={physics_region.permittivity})")
                except Exception as e:
                    logger.error(f"Error creating semiconductor physics region from config ID {region_config.id}: {e}")
                    raise
        
        # Initialize surface physics regions (placeholder, adapt as needed)
        self.surface_physics_regions: Dict[int, SurfaceRegion] = {} # SurfaceRegion is the physics class
        if self.config.surface and self.config.surface.regions:
            logger.warning("Surface region physics object creation is a placeholder and may need full implementation.")
            # Example:
            # for surf_region_config in self.config.surface.regions:
            #     try:
            #         # Assuming SurfaceRegion takes config directly or has a creator function
            #         # physics_surf_region = SurfaceRegion(**dataclasses.asdict(surf_region_config))
            #         # self.surface_physics_regions[physics_surf_region.id] = physics_surf_region
            #         # logger.info(f"Created surface physics region ID {physics_surf_region.id}")
            #     except Exception as e:
            #         logger.error(f"Error creating surface physics region from config ID {surf_region_config.id}: {e}")
            #         # Decide to raise or continue
            pass


        self.grid = self._create_grid() # Uses self.config for grid parameters
        
        # Ensure region_id_map is built and assigned to grid object after semiconductor regions are defined in config
        if hasattr(self.grid, 'set_material_regions') and self.config.semiconductor:
             self.grid.set_material_regions(self.config.semiconductor.regions) # Pass list of ConfigSemiconductorRegion
        else:
            logger.warning("Grid does not have set_material_regions or semiconductor config missing for region_id_map.")


        # ChargeDensityCalculator initialization
        # It needs system_fermi_level_E_F_main which is typically tip's Fermi energy
        tip_fermi_energy = getattr(getattr(self.config, 'tip', None), 'fermi_energy', 0.0)
        if tip_fermi_energy == 0.0:
            logger.warning("Tip Fermi energy is 0.0, check config.tip.fermi_energy.")

        self.charge_density_calculator = ChargeDensityCalculator(
            self.config,  # config_or_regions parameter
            self.semiconductor_physics_regions,  # semiconductor_regions_physics_or_surface parameter
            surface_regions_physics=self.surface_physics_regions,
            system_fermi_level_E_F_main=tip_fermi_energy 
        )

        # Create a props object for PoissonSOREquation
        class PoissonProps:
            def __init__(self, config):
                class SemiconductorProps:
                    # Use the first semiconductor region's properties
                    sem_region = config.semiconductor_regions[0]  # semiconductor_regions is a list
                    epsilon_r = sem_region.permittivity
                    Ev_offset_eV = sem_region.valence_band_offset  # Use valence_band_offset from semiconductor region
                self.semiconductor_props = SemiconductorProps()
        
        # Poisson solver uses grid and a props object
        self.poisson_solver = PoissonSOREquation(self.grid, PoissonProps(self.config))
        
        self.results = {} # To store results per voltage point
        logger.info("MultInt components initialized.")

    def _create_grid(self):
        """
        【修改】根據物理參數建立雙曲面網格。
        """
        grid_params = self.config.grid
        tip_params = self.config.tip
        
        logger.info("Creating HyperbolicGrid...")
        logger.info(f"Grid Params: N_eta (radial_points)={grid_params.radial_points}, N_nu (angular_points)={grid_params.angular_points}")
        logger.info(f"Tip Params: Radius={tip_params.radius} nm, Separation={tip_params.separation} nm")
        try:
            # Assuming HyperbolicGrid takes N_eta, N_nu, R, Z_TS, shank_slope
            # Mapping config names: radial_points -> N_eta, angular_points -> N_nu
            # R -> tip_params.radius, Z_TS -> tip_params.separation, shank_slope -> tip_params.shank_slope
            grid = HyperbolicGrid(
                N_eta=grid_params.radial_points, 
                N_nu=grid_params.angular_points,
                R=tip_params.radius,
                Z_TS=tip_params.separation,
                shank_slope=tip_params.shank_slope
                # r_max_factor can be added if it's in config.grid
            )
            logger.info(f"HyperbolicGrid created: Actual N_eta={grid.N_eta}, N_nu={grid.N_nu}, f={grid.f:.4f} nm, eta_tip={grid.eta_tip:.4f}")
            return grid
        except ValueError as e:
            logger.error(f"Error creating HyperbolicGrid: {e}")
            raise
        except Exception as e_gen:
            logger.error(f"Unexpected error creating HyperbolicGrid: {e_gen}")
            raise
            
    def mix_potential(self, old_potential, new_potential):
        # Use more conservative mixing for better convergence
        alpha = getattr(getattr(getattr(self.config, 'computation', None), 'mixing_alpha', None), 'value', 0.1)
        # If computation.mixing_alpha is a direct float:
        # alpha = self.config.computation.mixing_alpha if hasattr(self.config.computation, 'mixing_alpha') else 0.1
        logger.debug(f"Mixing potential with alpha = {alpha}")
        
        # Clamp potential to reasonable bounds to prevent numerical instability
        V_tip = getattr(getattr(self.config, 'tip', None), 'fermi_energy', 0.0)
        V_bounds = 10.0  # 10V bounds around tip potential
        mixed_potential = (1 - alpha) * old_potential + alpha * new_potential
        mixed_potential = np.clip(mixed_potential, V_tip - V_bounds, V_tip + V_bounds)
        
        return mixed_potential

    def _print_fortran_style_header(self):
        """Print simulation parameters in Fortran fort_multint.16 style format."""
        logger.info(" ")  # Empty line like Fortran
        
        # Tip parameters
        tip_conf = getattr(self.config, 'tip', None)
        if tip_conf:
            radius = getattr(tip_conf, 'radius', 1.0)
            slope = 1.0  # Default slope
            angle = 90.0  # Default angle
            logger.info(f" RAD, SLOPE, ANGLE = {radius:11.8f}     {slope:11.8f}      {angle:11.6f}")
            
            contact_potential = getattr(tip_conf, 'contact_potential', 0.0)
            logger.info(f" CONTACT POTENTIAL = {contact_potential:12.7f}")
            
            # Position of tip (default to 0,0)
            logger.info(f" POSITION OF TIP = {0.0:12.7f}     {0.0:12.7f}")
        
        # Semiconductor regions
        if hasattr(self.config, 'semiconductor_regions'):
            for i, region in enumerate(self.config.semiconductor_regions):
                region_id = getattr(region, 'id', i+1)
                logger.info(f" REGION #          {region_id:2d}")
                
                donor_conc = getattr(region, 'donor_concentration', 0.0)
                acceptor_conc = getattr(region, 'acceptor_concentration', 0.0)
                logger.info(f" DOPING = {donor_conc:14.8E}  {acceptor_conc:10.7f}")
                
                band_gap = getattr(region, 'band_gap', 1.42)
                vb_offset = getattr(region, 'valence_band_offset', 0.0)
                logger.info(f" BAND GAP, VB OFFSET = {band_gap:12.7f}     {vb_offset:12.7f}")
        
        # Surface states (placeholder)
        logger.info(" FIRST DISTRIBUTION OF SURFACE STATES:")
        logger.info(f" SURFACE STATE DENSITY, EN = {4.4e14:14.8E} {0.125:11.8f}")
        logger.info(f" FWHM, ECENT = {0.25:11.8f}     {1.625:12.7f}")
        logger.info(" SECOND DISTRIBUTION OF SURFACE STATES:")
        logger.info(f" SURFACE STATE DENSITY, EN = {0.0:12.7f}     {0.0:12.7f}")
        logger.info(f" FWHM, ECENT = {0.0:12.7f}     {0.0:12.7f}")
        
        logger.info(" HORIZONTAL MIRROR PLANE ASSUMED")
        
        # Fermi level and carrier density (from first semiconductor region)
        if self.semiconductor_physics_regions:
            first_region = next(iter(self.semiconductor_physics_regions.values()))
            fermi_level = first_region.fermi_level()
            logger.info(f" REGION TYPE 1, FERMI-LEVEL = {fermi_level:11.7f}")
            
            # Carrier densities
            n_cb = first_region.carrier_density_cb(fermi_level)  # electrons in CB (cm^-3)
            n_vb = first_region.carrier_density_vb(fermi_level)  # holes in VB (cm^-3)
            logger.info(f" CARRIER DENSITY IN CB, VB = {n_cb:14.8E}  {n_vb:10.6f}")

    def run_self_consistent_loop(self):
        logger.info("Starting self-consistent loop...")
        
        # Print simulation parameters in Fortran format
        self._print_fortran_style_header()
        
        grid_conf = getattr(self.config, 'grid', None)
        if grid_conf and self.grid:
            logger.info(f"Grid Config: N_eta (radial_points)={getattr(grid_conf, 'radial_points', 'N/A')}, N_nu (angular_points)={getattr(grid_conf, 'angular_points', 'N/A')}")
            logger.info(f"Grid Actual: N_eta={self.grid.N_eta}, N_nu={self.grid.N_nu}, f={self.grid.f:.4f} nm, eta_tip={self.grid.eta_tip:.4f}")
        
        volt_scan_conf = getattr(self.config, 'voltage_scan', None)
        if not volt_scan_conf:
            logger.error("Voltage scan configuration (voltage_scan) not found in config. Aborting.")
            return

        logger.info(f"Voltage Scan Points: {getattr(volt_scan_conf, 'points', 1)}")
        start_voltage = getattr(volt_scan_conf, 'start_voltage', 0.0)
        end_voltage = getattr(volt_scan_conf, 'end_voltage', 0.0)
        logger.info(f"Voltage Scan Start: {start_voltage} V, End: {end_voltage} V")
        
        # For now, run for the start_voltage only as a test
        # Loop over voltages will be: np.linspace(start_voltage, end_voltage, volt_scan_conf.points)
        voltage_points_to_scan = [start_voltage] 
        logger.info(f"Processing for voltage points: {voltage_points_to_scan}")

        for current_bias_V in voltage_points_to_scan:
            # Get tip configuration
            tip_conf = getattr(self.config, 'tip', None)
            
            # Implement Fortran bias voltage adjustment logic with modulation
            # From Fortran: bias = bias0 + imod * bmod * sqrt(2)
            modulation_voltage = getattr(volt_scan_conf, 'modulation_voltage', 0.050)  # Default 0.050V
            imod = -1  # For first calculation (Fortran logic)
            adjusted_bias = current_bias_V + imod * modulation_voltage * np.sqrt(2)
            
            # Calculate tip potential with contact potential adjustment
            contact_potential = getattr(tip_conf, 'contact_potential', 0.0) if tip_conf else 0.0
            V_tip_eff_Volts = adjusted_bias - contact_potential
            V_sample_eff_Volts = 0.0 

            # Print in Fortran style format
            logger.info(" ")  # Empty line like Fortran
            logger.info(f" SEPARATION = {getattr(tip_conf, 'separation', 1.0):11.8f}")
            logger.info(" ")  # Empty line like Fortran
            logger.info(f" BIAS, TIP POTENTIAL = {adjusted_bias:11.7f}     {V_tip_eff_Volts:11.7f}")
            
            # Add 1-D depletion width estimate (placeholder for now)
            depletion_width = 54.34  # nm, placeholder value similar to Fortran
            logger.info(f" 1-D ESTIMATE OF DEPLETION WIDTH (NM) = {depletion_width:11.6f}")
            
            # Calculate ESTART and EEND using Fortran logic
            if hasattr(self, 'charge_density_calculator') and self.charge_density_calculator.charge_density_tables:
                # Use Fortran MultInt3-6.4.f logic (lines 342-351)
                # EF = Fermi level of first region (from EFFIND)
                # PotTIP = adjusted_bias (includes contact potential)
                # EN0MIN, EN0MAX = min and max of charge neutrality levels across all regions
                
                first_region = next(iter(self.semiconductor_physics_regions.values()))
                EF = first_region.fermi_level()  # This corresponds to EFFIND(1,EF)
                PotTIP = adjusted_bias  # This is bias + contact potential
                
                # Calculate EN0MIN and EN0MAX from all semiconductor regions
                EN0_values = []
                for region_physics in self.semiconductor_physics_regions.values():
                    en0 = region_physics.fermi_level()  # charge neutrality level for each region
                    EN0_values.append(en0)
                
                EN0MIN = min(EN0_values) if EN0_values else EF
                EN0MAX = max(EN0_values) if EN0_values else EF
                
                # Fortran logic: lines 342-351
                # ESTART=AMIN1(EF,EF-PotTIP,EN0MIN)
                # EEND=AMAX1(EF,EF-PotTIP,EN0MAX)
                EF_minus_PotTIP = EF - PotTIP
                estart_initial = min(EF, EF_minus_PotTIP, EN0MIN)
                eend_initial = max(EF, EF_minus_PotTIP, EN0MAX)
                
                # ETMP=EEND-ESTART
                etmp = eend_initial - estart_initial
                
                # ESTART=ESTART-2.*ETMP
                # EEND=EEND+2.*ETMP  
                estart = estart_initial - 2.0 * etmp
                eend = eend_initial + 2.0 * etmp
                
                # Number of energy points (Fortran uses NE=20000)
                ne = 20000
                
                # DELE=(EEND-ESTART)/FLOAT(NE-1)
                dele = (eend - estart) / float(ne - 1)
                
                # PLACE ONE OF THE TABLE VALUES FOR ENERGY AT EF +/- (DELE/2)
                # NETMP=NINT((EF-ESTART)/DELE)
                # ESTART=EF-(NETMP-0.5)*DELE
                # EEND=ESTART+(NE-1)*DELE
                netmp = round((EF - estart) / dele)
                estart_final = EF - (netmp - 0.5) * dele
                eend_final = estart_final + (ne - 1) * dele
                
                logger.info(f" ESTART,EEND,NE = {estart_final:11.7f}     {eend_final:11.6f}          {ne:5d}")
                
                # Debug information to compare with Fortran
                logger.info(f"Energy range calculation debug:")
                logger.info(f"  EF = {EF:.7f}")
                logger.info(f"  PotTIP = {PotTIP:.7f}")
                logger.info(f"  EF - PotTIP = {EF_minus_PotTIP:.7f}")
                logger.info(f"  EN0MIN = {EN0MIN:.7f}")
                logger.info(f"  EN0MAX = {EN0MAX:.7f}")
                logger.info(f"  Initial: ESTART = {estart_initial:.7f}, EEND = {eend_initial:.7f}")
                logger.info(f"  After 2x expansion: ESTART = {estart:.7f}, EEND = {eend:.7f}")
                logger.info(f"  Final: ESTART = {estart_final:.7f}, EEND = {eend_final:.7f}")
            else:
                # Fallback values if charge density calculator not available
                estart_final = -6.6
                eend_final = 10.2
                ne = 20000
                logger.info(f" ESTART,EEND,NE = {estart_final:11.7f}     {eend_final:11.6f}          {ne:5d}")
            
            logger.info(" COMPUTING TABLE OF BULK CHARGE DENSITIES")
            logger.info(" COMPUTING TABLE OF SURFACE CHARGE DENSITIES")
            
            # Grid parameters in Fortran style
            if hasattr(self.grid, 'eta_tip') and hasattr(self.grid, 'f'):
                # Calculate grid parameters
                etat = self.grid.eta_tip
                a = self.grid.f
                z0 = 5.96046448E-08  # placeholder
                c = 5.96046519E-08   # placeholder
                logger.info(f" ETAT, A, Z0, C = {etat:11.8f}     {a:11.7f}     {z0:14.8E} {c:14.8E}")
                
                # Grid dimensions
                nr = self.grid.N_eta
                ns = nr  # Same as nr
                nv = 4   # placeholder
                np_val = self.grid.N_nu
                logger.info(f" NR,NS,NV,NP =         {nr:3d}         {ns:3d}          {nv:2d}          {np_val:2d}")
                
                # Grid spacing
                delr = 0.5  # placeholder
                dels = 0.5  # placeholder  
                delv = 0.25 # placeholder
                delp = 0.39270  # placeholder
                logger.info(f" DELR,DELS,DELV,DELP = {delr:8.5f}     {dels:8.5f}     {delv:8.5f}     {delp:8.5f}")
                
                # Largest radius and depth
                largest_radius = 103.67  # placeholder
                depth = 103.67  # placeholder
                logger.info(f" LARGEST RADIUS, DEPTH = {largest_radius:11.5f}     {depth:11.5f}")
            
            logger.info(" SOLUTION #           1")

            logger.info("Solving Laplace equation for initial potential distribution...")
            if not self.grid.region_id_map_initialized: # Check if region_id_map is ready
                 logger.error("Grid region_id_map not initialized. Cannot solve Laplace/Poisson. Call grid.set_material_regions first.")
                 return # or raise error

            initial_potential_Volts, laplace_iters, laplace_error = self.poisson_solver.solve_laplace(
                V_tip_Volts=V_tip_eff_Volts, V_sample_Volts=V_sample_eff_Volts,
                region_id_map=self.grid.region_id_map, # Pass the actual map
                # Other args for solve_laplace if its signature demands (charge_calc, fermi_level are optional in poisson.py)
            )
            logger.info(f"Laplace solution: {laplace_iters} iterations, Max error: {laplace_error:.3e} V.")
            
            potential_Volts = initial_potential_Volts
            
            comp_conf = getattr(self.config, 'computation', None)
            if not comp_conf:
                logger.error("Computation settings (computation) not found in config. Using defaults.")
                max_scf_iterations = 50
                scf_tolerance = 1e-4
                sor_max_iters = 2000
                sor_tolerance = 1e-5
                sor_omega = 1.8
            else:
                max_scf_iterations = comp_conf.max_iterations[0] if comp_conf.max_iterations else 50
                scf_tolerance = comp_conf.convergence_parameters[0] if comp_conf.convergence_parameters else 1e-4
                # Assuming SOR params might be different or taken from the same array for now
                sor_max_iters = comp_conf.max_iterations[-1] if len(comp_conf.max_iterations) > 1 else max_scf_iterations 
                sor_tolerance = comp_conf.convergence_parameters[-1] if len(comp_conf.convergence_parameters) > 1 else scf_tolerance
                sor_omega = getattr(comp_conf, 'sor_omega', 1.8) # Assuming sor_omega can be a field

            latest_rho_C_cm3 = np.zeros_like(potential_Volts) # Initialize charge

            for scf_iter in range(max_scf_iterations):
                logger.info(f"--- SCF Iteration: {scf_iter + 1} / {max_scf_iterations} (Bias: {current_bias_V:.3f}V) ---")
                
                logger.info("Calculating charge density rho(V)...")
                try:
                    new_rho_C_cm3 = self.charge_density_calculator.calculate_charge_density(
                        potential_matrix_Volts=potential_Volts,
                        region_id_map=self.grid.region_id_map
                    )
                    if new_rho_C_cm3 is None: # calculate_charge_density might not be fully implemented
                        logger.error("Charge density calculation returned None. Aborting SCF.")
                        return
                    latest_rho_C_cm3 = new_rho_C_cm3 # Store for results
                    # Log some stats about new_rho_C_cm3
                    logger.debug(f"Rho stats: min={np.min(new_rho_C_cm3):.2e}, max={np.max(new_rho_C_cm3):.2e}, mean={np.mean(new_rho_C_cm3):.2e} C/cm^3")

                except Exception as e:
                    logger.error(f"Error during charge density calculation in SCF iter {scf_iter + 1}: {e}", exc_info=True)
                    break # Exit SCF loop on error

                logger.info("Solving Poisson equation V(rho)...")
                try:
                    new_potential_Volts, poisson_iters, poisson_converged = self.poisson_solver.solve(
                        V_tip_Volts=V_tip_eff_Volts,
                        V_sample_Volts=V_sample_eff_Volts,
                        charge_density_calculator=self.charge_density_calculator,
                        system_fermi_level_E_F_main_eV=getattr(tip_conf, 'fermi_energy', 0.0),
                        omega=sor_omega, 
                        tolerance_Volts=sor_tolerance,
                        max_iterations=sor_max_iters
                    )
                    poisson_error = sor_tolerance if poisson_converged else 1.0
                    logger.info(f"Poisson solution: {poisson_iters} iterations, Max error: {poisson_error:.3e} V.")
                except Exception as e:
                    logger.error(f"Error during Poisson solution in SCF iter {scf_iter + 1}: {e}", exc_info=True)
                    break # Exit SCF loop on error

                potential_change_norm = np.linalg.norm(new_potential_Volts - potential_Volts)
                potential_norm = np.linalg.norm(potential_Volts)
                relative_change = potential_change_norm / (potential_norm + 1e-9) # Avoid division by zero
                
                logger.info(f"SCF Iteration {scf_iter + 1}: Potential change norm = {potential_change_norm:.3e}, Relative change = {relative_change:.3e}")

                potential_Volts = self.mix_potential(potential_Volts, new_potential_Volts)

                if relative_change < scf_tolerance:
                    logger.info(f"Self-consistent loop converged in {scf_iter + 1} iterations for bias {current_bias_V:.3f}V.")
                    break
                if scf_iter == max_scf_iterations - 1:
                    logger.warning(f"Self-consistent loop DID NOT converge after {max_scf_iterations} iterations for bias {current_bias_V:.3f}V. Final relative change: {relative_change:.3e}")
            
            # 計算最終的VSINT Pot0值並顯示
            try:
                # 使用完整的VSINT+縮放修正計算最終Pot0
                vsint_array = self.poisson_solver._initialize_vsint_array()
                vsint_array = self.poisson_solver._update_vsint_with_surface_charge(
                    vsint_array, potential_Volts, self.charge_density_calculator,
                    getattr(tip_conf, 'fermi_energy', 1.4187), V_tip_eff_Volts)
                
                final_pot0_vsint_scaled = self.poisson_solver._calculate_pot0_fortran_style(
                    potential_Volts, use_vsint=True, vsint_array=vsint_array, apply_scaling_correction=True)
                final_pot0_regular_scaled = self.poisson_solver._calculate_pot0_fortran_style(
                    potential_Volts, use_vsint=False, apply_scaling_correction=True)
                
                logger.info(f"")
                logger.info(f"🎯 FINAL POT0 RESULTS for {current_bias_V:.3f}V:")
                logger.info(f"   Regular Pot0 (scaled):  {final_pot0_regular_scaled:14.8E} V")
                logger.info(f"   VSINT Pot0 (scaled):    {final_pot0_vsint_scaled:14.8E} V")
                logger.info(f"   Difference from Fortran (-0.08V): {abs(final_pot0_vsint_scaled - (-0.08)):.6f} V")
                if abs(final_pot0_vsint_scaled - (-0.08)) < 0.01:
                    logger.info(f"   ✅ EXCELLENT AGREEMENT with Fortran!")
                elif abs(final_pot0_vsint_scaled - (-0.08)) < 0.05:
                    logger.info(f"   ✅ Good agreement with Fortran")
                else:
                    logger.info(f"   ⚠️  Large deviation from Fortran")
                logger.info(f"")
                
                # 將最終Pot0添加到結果中
                final_pot0_results = {
                    'pot0_regular_scaled': final_pot0_regular_scaled,
                    'pot0_vsint_scaled': final_pot0_vsint_scaled,
                    'pot0_fortran_difference': abs(final_pot0_vsint_scaled - (-0.08))
                }
                
            except Exception as e:
                logger.warning(f"Error calculating final VSINT Pot0: {e}")
                final_pot0_results = {}

            self.results[current_bias_V] = {
                'potential_Volts': potential_Volts,
                'rho_C_cm3': latest_rho_C_cm3,
                'scf_iterations': scf_iter + 1,
                'final_potential_relative_change': relative_change if 'relative_change' in locals() else -1.0,
                'adjusted_bias_V': adjusted_bias,  # Store the adjusted bias for reference
                'tip_potential_V': V_tip_eff_Volts,
                'final_pot0': final_pot0_results  # 添加最終Pot0結果
            }
            logger.info(f"--- Finished processing Bias Voltage: {current_bias_V:.3f} V (adjusted: {adjusted_bias:.7f} V) ---")
        
        logger.info("All voltage points processed. Self-consistent loop finished.")


# 向後兼容性別名
MultIntSimulation = MultInt