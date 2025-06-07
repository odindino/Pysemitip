"""
MultInt main simulation program.

This module implements the main MultInt simulation flow, corresponding to
the MultInt3-6.3.f Fortran program. It coordinates all the physics modules
to perform STM simulations with multiple semiconductor regions.
"""

import numpy as np
import time
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import logging

from ..core.config_schema import SemitipConfig
from ..physics.materials import (SemiconductorRegion, SurfaceRegion, TipModel,
                               create_semiconductor_from_config,
                               create_surface_region_from_config,
                               create_tip_from_config)
from ..physics.solvers import Grid3D, create_grid_from_config
from ..physics.core import (ChargeDensityCalculator, ChargeDensityTables,
                          PoissonSolver, PoissonSolverParameters,
                          PotentialProcessor, estimate_energy_range,
                          SchrodingerSolver, create_schrodinger_solver)
from ..utils.constants import PhysicalConstants as PC


# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class SimulationResults:
    """Container for simulation results at a single bias voltage."""
    bias_voltage: float
    tip_potential: float
    band_bending: float
    depletion_width: float
    potential_3d: np.ndarray
    convergence_info: dict
    current: float = 0.0  # Will be calculated by current module
    
    # Optional detailed results
    potential_profile: Optional[np.ndarray] = None
    charge_density_profile: Optional[np.ndarray] = None


class MultIntSimulation:
    """
    Main MultInt simulation class.
    
    Coordinates the solution of Poisson equation and calculation of
    tunneling current for multiple bias voltages.
    """
    
    def __init__(self, config: SemitipConfig):
        """
        Initialize simulation from configuration.
        
        Args:
            config: Validated SEMITIP configuration
        """
        self.config = config
        
        # Set up logging
        self._setup_logging()
        
        # Initialize components
        self._initialize_materials()
        self._initialize_grid()
        self._initialize_solvers()
        
        # Results storage
        self.results: List[SimulationResults] = []
    
    def _setup_logging(self):
        """Set up logging with file output matching Fortran fort.16."""
        # Create output file handler
        self.output_file = open('output_MultInt.log', 'w')
        
        # Log initial parameters
        self._log_parameters()
    
    def _log_parameters(self):
        """Log simulation parameters matching Fortran output format."""
        self.output_file.write("\n")
        self.output_file.write(f"RAD, SLOPE, ANGLE = {self.config.tip.radius} "
                             f"{getattr(self.config.tip, 'slope', 1.0)} "
                             f"{90.0}\n")  # Angle calculation
        self.output_file.write(f"CONTACT POTENTIAL = {self.config.contact_potential}\n")
        self.output_file.write(f"POSITION OF TIP = {self.config.tip.position.x} "
                             f"{self.config.tip.position.y}\n")
        
        # Log semiconductor regions
        for i, region in enumerate(self.config.semiconductor_regions):
            self.output_file.write(f"REGION # {i+1}\n")
            self.output_file.write(f"DOPING = {region.donor_concentration:.6e} "
                                 f"{region.acceptor_concentration:.6e}\n")
            self.output_file.write(f"BAND GAP, VB OFFSET = {region.band_gap} "
                                 f"{region.valence_band_offset}\n")
        
        # Log surface regions
        for region in self.config.surface_regions:
            self.output_file.write("FIRST DISTRIBUTION OF SURFACE STATES:\n")
            self.output_file.write(f"SURFACE STATE DENSITY, EN = "
                                 f"{region.first_distribution.density:.6e} "
                                 f"{region.first_distribution.neutrality_level}\n")
            self.output_file.write(f"FWHM, ECENT = {region.first_distribution.fwhm} "
                                 f"{region.first_distribution.center_energy}\n")
            
            if region.second_distribution:
                self.output_file.write("SECOND DISTRIBUTION OF SURFACE STATES:\n")
                self.output_file.write(f"SURFACE STATE DENSITY, EN = "
                                     f"{region.second_distribution.density:.6e} "
                                     f"{region.second_distribution.neutrality_level}\n")
        
        if self.config.mirror_symmetry:
            self.output_file.write("HORIZONTAL MIRROR PLANE ASSUMED\n")
        
        self.output_file.flush()
    
    def _initialize_materials(self):
        """Initialize material models from configuration."""
        # Create semiconductor regions
        self.semiconductor_regions = []
        for config_region in self.config.semiconductor_regions:
            region = create_semiconductor_from_config(config_region)
            region.temperature = self.config.temperature
            self.semiconductor_regions.append(region)
        
        # Create surface regions
        self.surface_regions = []
        for i, config_region in enumerate(self.config.surface_regions):
            region = create_surface_region_from_config(config_region, i+1)
            region.temperature = self.config.temperature
            self.surface_regions.append(region)
        
        # Create tip model
        self.tip = create_tip_from_config(self.config)
        
        # Calculate Fermi level for primary region
        self.fermi_level = self._calculate_fermi_level()
        
        logger.info(f"Fermi level: {self.fermi_level:.6f} eV")
        self.output_file.write(f"REGION TYPE 1, FERMI-LEVEL = {self.fermi_level:.6f}\n")
        
        # Calculate carrier densities
        region1 = self.semiconductor_regions[0]
        n = region1.carrier_density_cb(self.fermi_level)
        p = region1.carrier_density_vb(self.fermi_level)
        self.output_file.write(f"CARRIER DENSITY IN CB, VB = {n:.6e} {p:.6e}\n")
        self.output_file.flush()
    
    def _calculate_fermi_level(self) -> float:
        """
        Calculate Fermi level for the primary semiconductor region.
        
        This implements the EFFIND functionality.
        """
        # For now, use simplified calculation
        # Full implementation would solve charge neutrality equation
        region = self.semiconductor_regions[0]
        return region.fermi_level()
    
    def _initialize_grid(self):
        """Initialize computational grid."""
        self.grid = create_grid_from_config(self.config, self.tip)
        
        logger.info(f"Grid initialized: nr={self.grid.params.nr}, "
                   f"nv={self.grid.params.nv}, ns={self.grid.params.ns}, "
                   f"np={self.grid.params.np}")
    
    def _initialize_solvers(self):
        """Initialize numerical solvers."""
        # Charge density calculator
        self.charge_calculator = ChargeDensityCalculator(
            self.semiconductor_regions,
            self.surface_regions,
            self.fermi_level
        )
        
        # Poisson solver
        solver_params = PoissonSolverParameters(
            tolerance=1e-6,
            max_iterations=2000,
            omega=0.8,
            adaptive_omega=True,
            verbose=True
        )
        self.poisson_solver = PoissonSolver(self.grid, self.tip, solver_params)
        
        # Potential processor
        self.potential_processor = PotentialProcessor(self.grid)
        
        # Schrödinger solver for current calculation
        self.schrodinger_solver = create_schrodinger_solver()
        
        logger.info("Solvers initialized")
    
    def run(self, bias_voltages: Optional[List[float]] = None,
            save_interval: int = 1, max_points: Optional[int] = None) -> List[SimulationResults]:
        """
        Run simulation for specified bias voltages.
        
        Args:
            bias_voltages: List of bias voltages (V). If None, use config values.
            save_interval: Save results every N bias points
            max_points: Maximum number of voltage points to process (for testing)
            
        Returns:
            List of SimulationResults
        """
        if bias_voltages is None:
            bias_voltages = self.config.voltage_scan.voltages
        
        # Limit number of points for testing
        if max_points is not None:
            bias_voltages = bias_voltages[:max_points]
        
        logger.info(f"Starting simulation for {len(bias_voltages)} bias voltages")
        self.output_file.write("\n")
        
        # Main bias voltage loop
        for i, bias in enumerate(bias_voltages):
            logger.info(f"\nBias voltage {i+1}/{len(bias_voltages)}: {bias:.4f} V")
            
            # Run single bias point
            result = self._solve_single_bias(bias)
            self.results.append(result)
            
            # Save intermediate results
            if (i + 1) % save_interval == 0:
                self._save_results()
        
        # Final save
        self._save_results()
        
        logger.info("Simulation completed")
        return self.results
    
    def _solve_single_bias(self, bias_voltage: float) -> SimulationResults:
        """
        Solve for a single bias voltage with modulation (following Fortran logic).
        
        Args:
            bias_voltage: Base bias voltage (V)
            
        Returns:
            SimulationResults object
        """
        start_time = time.time()
        
        # Get modulation voltage from config
        modulation_voltage = self.config.voltage_scan.modulation_voltage
        
        # Implement Fortran modulation logic: imod loop from -1 to 1 step 2
        # For now, use imod = -1 (first iteration) to match Fortran output
        if modulation_voltage > 0:
            # bias = bias0 + imod * bmod * sqrt(2.)
            import math
            imod = -1  # First iteration (like Fortran)
            bias_actual = bias_voltage + imod * modulation_voltage * math.sqrt(2.0)
        else:
            bias_actual = bias_voltage
        
        # Update tip potential with actual bias
        self.tip.bias_voltage = bias_actual
        tip_potential = self.tip.tip_potential
        
        self.output_file.write(f"\nSEPARATION = {self.tip.separation}\n")
        self.output_file.write(f"\nBIAS, TIP POTENTIAL = {bias_actual} {tip_potential}\n")
        
        # Estimate depletion width (1D approximation)
        depl_width_1d = self._estimate_depletion_width_1d(bias_actual)
        self.output_file.write(f"1-D ESTIMATE OF DEPLETION WIDTH (NM) = {depl_width_1d:.6f}\n")
        
        # Create charge density tables
        energy_range = self._calculate_energy_range(bias_actual)
        self.output_file.write(f"ESTART,EEND,NE = {energy_range[0]:.6f} "
                             f"{energy_range[1]:.6f} {self.config.charge_table_points}\n")
        
        logger.info("Computing charge density tables...")
        self.output_file.write("COMPUTING TABLE OF BULK CHARGE DENSITIES\n")
        self.output_file.write("COMPUTING TABLE OF SURFACE CHARGE DENSITIES\n")
        
        # Add hyperboloid parameters (ETAT, A, Z0, C) like Fortran
        eta, a, z0, c = self.tip.hyperboloid_parameters()
        self.output_file.write(f"ETAT, A, Z0, C = {eta:.8f} {a:.8f} {z0:.8e} {c:.8e}\n")
        
        self.output_file.flush()
        
        charge_tables = self.charge_calculator.create_charge_tables(
            energy_range, 
            self.config.charge_table_points
        )
        
        # Store charge tables for use in charge density functions
        self._current_charge_tables = charge_tables
        
        # Multi-grid solution sequence
        # solution_3d_coarse stores the potential from the previous (coarser) grid stage
        solution_3d_coarse = None
        previous_stage_grid_params = None # Stores GridParameters of the coarser grid

        grid_sequence = self._get_grid_sequence()
        
        # Base delr, dels for the tangent formula are from the initial grid configuration
        # These are used as DELR0, DELS0 in the Fortran tangent formula for r and s coordinates.
        # The Grid3D constructor takes these base values.
        # For multi-grid, these base values themselves are scaled.
        base_delr_config = self.grid.params.delr
        base_dels_config = self.grid.params.dels
        base_delv_config = self.grid.params.delv # This is for linear vacuum grid; hyperboloidal is handled by Grid3D
        base_delp_config = self.grid.params.delp
        # Base extents from the initial full grid configuration
        base_rmax_config = self.grid.params.rmax
        base_vmax_config = self.grid.params.vmax
        base_smax_config = self.grid.params.smax

        for grid_level, (nr_stage, nv_stage, ns_stage, np_stage) in enumerate(grid_sequence):
            self.output_file.write(f"NR,NS,NV,NP = {nr_stage} {ns_stage} {nv_stage} {np_stage}\n")

            # Calculate grid spacing and extents for this level
            # spacing_factor = 2.0 ** grid_level # This was used to scale from finest to current
            # For multi-grid, delr, dels, etc. usually refer to the actual spacing of the current grid.
            # However, our Grid3D uses delr/dels as input to tangent formula (like DELR0, DELS0)
            # The Fortran code effectively scales these DELR0, DELS0 by 1/(2^level)

            # Effective DELR0, DELS0 for the tangent formula at this grid level
            # If grid_level 0 is coarsest, then factor is 2.0**(num_levels - 1 - grid_level)
            # If grid_level 0 is finest (as current self.grid is), then factor is 1.0 / (2.0**grid_level_from_finest)
            # The current grid_sequence starts coarse and goes fine.
            # So, for grid_level 0 (coarsest in sequence), spacing is largest.
            # Let's assume grid_sequence[0] is the coarsest.
            # The scaling factor should be relative to the *initial* finest grid parameters.
            # The `grid_sequence` gives absolute number of points for current stage.
            # `self.grid.params` refers to the *initial* (usually finest) grid.

            # Let's assume the delr, dels, delv, delp in self.grid.params are for the *finest* grid.
            # The scaling factor for DELR0, DELS0, etc. needs to be determined based on
            # how many "doubling" steps coarser the current stage is from the finest.
            # Example: finest is 64x64. Stage 0 is 16x16 (2 doublings away). Factor = 2^2 = 4.
            # Stage 1 is 32x32 (1 doubling away). Factor = 2^1 = 2.
            # Stage 2 is 64x64 (0 doublings away). Factor = 2^0 = 1.
            num_stages = len(grid_sequence)
            # grid_level goes from 0 (coarsest) to num_stages-1 (finest)
            # Let's say finest_nr = self.grid.params.nr
            # current_nr = nr_stage
            # scaling_factor_nr = finest_nr / current_nr (assuming powers of 2)
            # This scaling_factor is for the number of points.
            # The DELR0, DELS0 for tangent formula scale inversely with number of points for same extent
            # OR, if extent also scales, it's more complex.

            # The existing delr, dels, delv, delp calculations in the loop header are correct for *actual spacing*.
            # We need to pass these to GridParameters.
            spacing_factor_for_actual_spacing = 2.0 ** (num_stages - 1 - grid_level) # if grid_level 0 is coarsest

            # Recalculate delr, dels, delv, delp based on the *initial grid's parameters*
            # These are the DELR0, DELS0, etc. for the current stage.
            # The current `grid_sequence` starts with the coarsest grid.
            # `self.grid` is the initial grid, usually the finest.
            # Let's find the scaling factor from the current stage to the reference (finest) grid.
            # If nr_stage is nr_finest / K, then delr_for_tangent_formula_stage = delr_finest / K

            # The parameters `delr`, `dels`, `delv`, `delp` calculated in the loop were the *actual* cell sizes.
            # GridParameters expects the DELR0, DELS0 type values for tangent formula.
            # The existing code was: `delr = self.grid.params.delr / spacing_factor` where spacing_factor
            # was relative to the current `grid_level` assuming `self.grid.params.delr` was for level 0.
            # This interpretation needs to be fixed. `self.grid` is the full grid.

            # Let's use the number of points to derive the scaling for DELR0, DELS0.
            # If nr_stage = self.grid.params.nr / K, then current_delr0 = self.grid.params.delr / K
            # This assumes rmax, smax, vmax are kept constant across stages for the tangent formula scaling.
            # The rmax, smax, vmax in the loop are *effective extents* which grow for finer actual spacing.

            # Simplification: Assume the `delr, dels, delv, delp` calculated in the loop are correct *inputs*
            # for the GridParameters of the current stage.
            # This means Grid3D's tangent formula will use these as the "base" spacings for that stage.

            # Calculate actual cell spacings for the current stage
            # This logic appears to assume self.grid.params are for the *coarsest* grid in the sequence.
            # Let's correct this. self.grid.params should be for the *reference* (usually finest) grid.
            # The grid_sequence defines N_stage. We need DEL_stage.
            # DEL_stage = RMAX_stage / N_stage (for linear)
            # For tangent, DEL0_stage is such that R(N_stage, DEL0_stage) = RMAX_stage.
            # R(I) = (2*N*DEL0/PI)*TAN(PI*(I-0.5)/(2.*N))
            # RMAX = R(N) = (2*N*DEL0/PI)*TAN(PI*(N-0.5)/(2.*N))
            # So DEL0 = (RMAX * PI) / (2*N*TAN(PI*(N-0.5)/(2.*N)))

            # The existing code for rmax, smax in the loop:
            # rmax_stage = self.grid.r[-1] * spacing_factor_from_coarsest
            # smax_stage = self.grid.zs_positive[-1] * spacing_factor_from_coarsest
            # This implies self.grid.r[-1] is max radius of coarsest grid.
            # This is confusing. Let's use the config's rmax, smax, vmax as fixed physical extents.

            current_rmax = self.config.grid.radial_extent
            current_smax = self.config.grid.semiconductor_extent
            current_vmax = self.config.grid.vacuum_extent

            # Calculate DELR0, DELS0 for the current stage such that the extent is matched
            # Using the inverse of the tangent formula for the last point
            # For R(nr_stage) = current_rmax:
            # current_delr0 = (current_rmax * np.pi) / (2 * nr_stage * np.tan(np.pi * (nr_stage - 0.5) / (2 * nr_stage)))
            # current_dels0 = (current_smax * np.pi) / (2 * ns_stage * np.tan(np.pi * (ns_stage - 0.5) / (2 * ns_stage)))
            # This is complex. Let's use the simpler Fortran approach: DELR0, DELS0 are scaled.
            # Assume self.grid.params.delr and .dels are the DELR0, DELS0 for the *finest* grid.
            # If current stage is K times coarser in points (e.g. nr_finest / nr_stage = K)
            # Then DELR0_stage = DELR0_finest * K (to maintain same extent with fewer points in tangent formula)
            # OR DELR0_stage = DELR0_finest / (2**level_diff_from_finest)

            # Let's use the scaling factor method like Fortran:
            # The delr, dels, delv, delp in self.grid.params are the *actual* cell sizes of the *initial* grid.
            # The multi-grid stages have N points, and we need to set their DELR0, DELS0, DELV0, DELP0.
            # The Fortran code scales DELR, DELS by 1/2 at each stage. These are the actual dx, dy.
            # It seems GridParameters.delr is more like DELR0 for tangent formula.

            # Re-evaluate delr, dels, delv, delp for current stage GridParameters:
            # These are the 'DELR0', 'DELS0' for the tangent formula for *this specific grid stage*.
            # The physical extent (rmax, smax, vmax) is assumed constant.
            # DELR0_stage = (RMAX * PI) / (2*nr_stage*TAN(PI*(nr_stage-0.5)/(2*nr_stage)))
            # This is how GridParameters should be populated if rmax, smax are fixed.
            if nr_stage > 1:
                delr_param = (current_rmax * np.pi) / (2 * nr_stage * np.tan(np.pi * (nr_stage - 0.5) / (2.0 * nr_stage)))
            else: # nr_stage = 1
                delr_param = current_rmax # Single point, spacing is extent

            if ns_stage > 1:
                dels_param = (current_smax * np.pi) / (2 * ns_stage * np.tan(np.pi * (ns_stage - 0.5) / (2.0 * ns_stage)))
            else: # ns_stage = 1
                dels_param = current_smax

            if nv_stage > 1:
                delv_param = current_vmax / (nv_stage -1) # Linear spacing for vacuum extent
            else: # nv_stage = 1
                delv_param = current_vmax

            if np_stage > 0:
                delp_param = (2 * np.pi if not self.config.mirror_symmetry else np.pi) / np_stage
            else:
                delp_param = np.pi # Should not happen

            self.output_file.write(f"PARAM DELR,DELS,DELV,DELP = {delr_param:.5f} {dels_param:.5f} {delv_param:.5f} {delp_param:.5f}\n")
            self.output_file.write(f"PARAM RMAX,SMAX,VMAX = {current_rmax:.5f} {current_smax:.5f} {current_vmax:.5f}\n")

            current_stage_grid_params = self.grid.params.__class__( # Use same type as self.grid.params
                nr=nr_stage, nv=nv_stage, ns=ns_stage, np=np_stage,
                delr=delr_param, dels=dels_param, delv=delv_param, delp=delp_param,
                rmax=current_rmax, vmax=current_vmax, smax=current_smax,
                mirror_symmetry=self.config.mirror_symmetry
            )
            current_stage_grid = Grid3D(current_stage_grid_params, self.tip)
            
            # Update Poisson solver parameters for this grid level
            if hasattr(self.config.computation, 'convergence_parameters'):
                tolerance = self.config.computation.convergence_parameters[grid_level]
            else:
                default_tol = [1e-3, 1e-3, 1e-4] # Default tolerances matching Fortran
                tolerance = default_tol[min(grid_level, len(default_tol)-1)]
            
            stage_max_iter = self.config.computation.max_iterations[min(grid_level, len(self.config.computation.max_iterations)-1)]
            if grid_level == 0: # Coarsest grid
                stage_max_iter = max(stage_max_iter, 4000)
            
            solver_params_for_stage = PoissonSolverParameters(
                tolerance=tolerance,
                max_iterations=stage_max_iter,
                omega=0.8, # TODO: Consider if omega should vary per stage
                adaptive_omega=True, # TODO: Review adaptive omega strategy for multigrid
                verbose=True
            )
            
            # Create new Poisson solver for the current grid stage
            current_stage_poisson_solver = PoissonSolver(current_stage_grid, self.tip, solver_params_for_stage)
            
            self.output_file.write(f"SOLUTION # {grid_level + 1}\n")
            self.output_file.flush()
            
            # Pass previous stage solution (solution_3d_coarse) and its grid_params
            initial_guess_data = (solution_3d_coarse, previous_stage_grid_params)

            solution_3d_fine, conv_info = current_stage_poisson_solver.solve(
                self._bulk_charge_func,
                self._surface_charge_func,
                initial_guess_data=initial_guess_data
            )
            
            # Update for next iteration
            solution_3d_coarse = solution_3d_fine
            previous_stage_grid_params = current_stage_grid_params

            # Log convergence (using the solver for the current stage)
            bb = current_stage_poisson_solver.get_band_bending()
            self.output_file.write(f"NUMBER OF ITERATIONS = {conv_info['iterations']}\n")
            self.output_file.write(f"BAND BENDING AT MIDPOINT = {bb:.8f}\n")
            self.output_file.flush()
        
        # After the loop, solution_3d_coarse holds the solution from the finest grid
        solution_3d = solution_3d_coarse

        # The final poisson_solver instance available is the one from the finest grid.
        # Depletion width and potential profile should be calculated using this finest grid solver/processor.
        # Ensure potential_processor is updated if its grid changes, or make it stateless.
        # For now, assume self.potential_processor can work with solution_3d from finest grid.
        # If PotentialProcessor is grid-dependent, it needs to be updated/recreated for the finest grid.
        # Let's assume the last `current_stage_poisson_solver` is what we need.
        final_poisson_solver = current_stage_poisson_solver
        final_potential_processor = PotentialProcessor(final_poisson_solver.grid) # Use grid from final solver

        depl_width = final_poisson_solver.get_depletion_width()
        profile = final_potential_processor.extract_profile(solution_3d) # solution_3d is from finest grid
        
        # Calculate tunneling current
        self.output_file.write("\nCOMPUTATION OF CURRENT:\n")
        
        # Prepare band parameters
        region = self.semiconductor_regions[0]  # Primary region
        band_params = {
            'band_gap': region.band_gap,
            'cb_effective_mass': region.cb_effective_mass,
            'vb_effective_mass_heavy': region.vb_effective_mass_heavy,
            'vb_effective_mass_light': region.vb_effective_mass_light,
            'fermi_level': self.fermi_level,
            'tip_fermi_level': self.tip.fermi_level
        }
        
        # Calculate current
        current_result = self.schrodinger_solver.solve_tunneling_current(
            profile, band_params, bias_voltage, self.config.temperature
        )
        
        # Log current results
        self.output_file.write(f"number of VB light-hole localized states = "
                             f"{len([s for s in current_result.localized_states if s.band == 'vb'])}\n")
        self.output_file.write(f"number of CB localized states = "
                             f"{len([s for s in current_result.localized_states if s.band == 'cb'])}\n")
        self.output_file.write(f"valence band current ext,loc = "
                             f"{current_result.vb_current - current_result.vb_localized:.6e} "
                             f"{current_result.vb_localized:.6e}\n")
        self.output_file.write(f"conduction band current ext,loc = "
                             f"{current_result.cb_current - current_result.cb_localized:.6e} "
                             f"{current_result.cb_localized:.6e}\n")
        self.output_file.flush()
        
        logger.info(f"Total current: {current_result.total_current:.3e} A")
        
        # Create results
        result = SimulationResults(
            bias_voltage=bias_actual,  # Use actual bias with modulation
            tip_potential=tip_potential,
            band_bending=bb,
            depletion_width=depl_width,
            potential_3d=solution_3d.copy(),
            convergence_info=conv_info,
            potential_profile=profile,
            current=current_result.total_current
        )
        
        # Log computation time
        comp_time = time.time() - start_time
        logger.info(f"Bias point completed in {comp_time:.2f} s")
        
        return result
    
    def _estimate_depletion_width_1d(self, bias_voltage: float) -> float:
        """
        Estimate depletion width using 1D approximation.
        
        Args:
            bias_voltage: Applied bias (V)
            
        Returns:
            Estimated depletion width (nm)
        """
        # Simple depletion approximation for primary region
        region = self.semiconductor_regions[0]
        
        # For depletion width estimation, use only bias voltage (matches Fortran)
        # Built-in potential contributes to band bending but not depletion width calc
        total_potential = abs(bias_voltage)
        
        # Depletion width formula
        eps = region.permittivity * PC.EPSILON0
        doping = abs(region.net_doping) * 1e6  # Convert to m^-3
        
        if doping > 0:
            w = np.sqrt(2 * eps * total_potential / (PC.E * doping)) * 1e9  # Convert to nm
        else:
            w = 100.0  # Default for intrinsic
        
        return w
    
    def _calculate_energy_range(self, bias_voltage: float) -> Tuple[float, float]:
        """
        Calculate appropriate energy range for charge tables following Fortran logic.
        
        Based on Fortran MultInt3-6.4.f lines 342-351:
        ESTART=AMIN1(EF,EF-PotTIP,EN0MIN)
        EEND=AMAX1(EF,EF-PotTIP,EN0MAX)
        ETMP=EEND-ESTART
        ESTART=ESTART-2.*ETMP
        EEND=EEND+2.*ETMP
        
        Key insight: EN0 is the charge neutrality level, not the energy spread.
        For single surface distributions, EN0 = EN (neutrality level parameter).
        """
        # Get tip potential (BIAS + contact potential)
        tip_potential = self.tip.tip_potential
        
        # Calculate charge neutrality levels (EN0) following Fortran logic
        # In Fortran, EN0 is found by ENFIND routine for each surface area
        en0_values = []
        
        for region in self.surface_regions:
            # For single distribution case, EN0 = EN (neutrality level)
            if hasattr(region, 'distribution1') and region.distribution1:
                d1 = region.distribution1
                if d1.density > 0:
                    # If only first distribution is non-zero
                    if (not hasattr(region, 'distribution2') or 
                        region.distribution2 is None or 
                        region.distribution2.density == 0):
                        en0_values.append(d1.neutrality_level)
                    else:
                        # Both distributions non-zero: would need ENFIND calculation
                        # For now, use first distribution neutrality level as approximation
                        en0_values.append(d1.neutrality_level)
            
            # Second distribution only
            elif (hasattr(region, 'distribution2') and region.distribution2 and 
                  region.distribution2.density > 0):
                en0_values.append(region.distribution2.neutrality_level)
        
        # Handle case where no surface states are defined
        if not en0_values:
            # If en0_values is empty, initialize en0_min to float('inf')
            # and en0_max to float('-inf')
            en0_min = float('inf')
            en0_max = float('-inf')
        else:
            # Find min and max of charge neutrality levels across all areas
            en0_min = min(en0_values)
            en0_max = max(en0_values)
        
        # Fortran logic: ESTART=AMIN1(EF,EF-PotTIP,EN0MIN)
        ef_minus_tip = self.fermi_level - tip_potential
        estart = min(self.fermi_level, ef_minus_tip, en0_min)
        
        # Fortran logic: EEND=AMAX1(EF,EF-PotTIP,EN0MAX)
        eend = max(self.fermi_level, ef_minus_tip, en0_max)
        
        # Fortran logic: expand range by factor of 2
        etmp = eend - estart
        estart = estart - 2.0 * etmp
        eend = eend + 2.0 * etmp

        # Final adjustment of ESTART and EEND to align grid with EF
        NE = self.config.charge_table_points

        if NE <= 1:
            # If NE=1, the table has one point, no step.
            # The initial estart, eend are fine.
            # dele is undefined or could be considered eend-estart.
            # The Fortran code doesn't explicitly handle NE=1 for this adjustment part.
            # Let's assume NE > 1 based on typical usage.
            pass # Keep estart and eend as they are
        else: # NE > 1
            dele = (eend - estart) / (NE - 1)

            # Ensure dele is not zero to prevent division errors
            if dele == 0:
                # If dele is zero (e.g. estart == eend),
                # this adjustment doesn't make sense.
                # Keep previous estart, eend. This case should be rare.
                pass # Keep estart and eend as they are
            else:
                netmp = np.round((self.fermi_level - estart) / dele) # NINT equivalent
                estart_new = self.fermi_level - (netmp - 0.5) * dele
                eend_new = estart_new + (NE - 1) * dele
                estart = estart_new
                eend = eend_new

        return (estart, eend)
    
    def _get_grid_sequence(self) -> List[Tuple[int, int, int, int]]:
        """
        Get multi-grid sequence following Fortran SEMITIP3 logic.
        
        The Fortran code implements progressive grid doubling:
        1. Start with initial grid from config
        2. Double all dimensions (NR, NS, NV, NP) each stage
        3. If memory limit reached, try doubling only NS
        4. Continue until max stages (IPMAX=3) reached
        
        Returns:
            List of (nr, nv, ns, np) tuples
        """
        # Initial grid from configuration (matches Fortran NRIN, NSIN, NVIN, NPIN)
        nr_start = self.grid.params.nr
        nv_start = self.grid.params.nv  # This should be 4 from config
        ns_start = self.grid.params.ns
        np_start = self.grid.params.np
        
        # Maximum dimensions (Fortran limits from PARAMETER statement)
        # NRDIM=512, NVDIM=64, NSDIM=512, NPDIM=64
        max_nr = 512
        max_nv = 64
        max_ns = 512
        max_np = 64
        
        sequence = []
        nr, nv, ns, np = nr_start, nv_start, ns_start, np_start
        
        # Stage 1: Initial grid
        sequence.append((nr, nv, ns, np))
        
        # Stage 2: First doubling (if possible)
        if (nr*2 <= max_nr and nv*2 <= max_nv and 
            ns*2 <= max_ns and np*2 <= max_np):
            nr, nv, ns, np = nr*2, nv*2, ns*2, np*2
            sequence.append((nr, nv, ns, np))
            
            # Stage 3: Second doubling (if possible)  
            if (nr*2 <= max_nr and nv*2 <= max_nv and 
                ns*2 <= max_ns and np*2 <= max_np):
                nr, nv, ns, np = nr*2, nv*2, ns*2, np*2
                sequence.append((nr, nv, ns, np))
            else:
                # Try doubling only semiconductor grid (NS)
                if ns*2 <= max_ns:
                    sequence.append((nr, nv, ns*2, np))
        else:
            # If can't double everything, try doubling only semiconductor grid
            if ns*2 <= max_ns:
                sequence.append((nr, nv, ns*2, np))
        
        return sequence
    
    def _get_region_id(self, r: float, z: float, phi: float) -> int:
        """
        Determine semiconductor region ID at given position.
        
        For now, simple implementation - full version would use
        region boundaries from configuration.
        """
        # Example: region 1 for n-type, region 2 for p-type
        # This needs to be implemented based on actual geometry
        if z > -50:  # Near surface
            return 1
        else:
            return 2
    
    def _get_area_id(self, r: float, phi: float) -> int:
        """
        Determine surface area ID at given position.
        
        For now, returns area 1. Full implementation would check
        area boundaries.
        """
        return 1
    
    def _bulk_charge_func(self, r, z, phi, pot):
        """
        Bulk charge density function for Poisson solver.
        
        Following Fortran RHOBULK logic:
        ENER = EF - Pot  (energy represents local Fermi level)
        IENER = NINT((ENER-ESTART)/DELE) + 1
        RHO = RHOBTAB(IREG,IENER)
        
        Args:
            r, z, phi: Spatial coordinates
            pot: Electric potential (V)
            
        Returns:
            Charge density (C/cm^3)
        """
        try:
            # Determine which region
            region_id = self._get_region_id(r, z, phi)
            
            # Convert pot to float - handle all numeric types safely
            pot_val = float(pot)
            
            # Calculate local Fermi level: ENER = EF - Pot
            # This is the key insight from Fortran RHOBULK
            local_fermi_level = self.fermi_level - pot_val
            
            # Safety checks - check if values are finite
            if not (np.isfinite(pot_val) and np.isfinite(local_fermi_level)):
                return 0.0
            
            # Get charge tables
            charge_tables = self._current_charge_tables
            
            # Check if local Fermi level is within table range
            energy_min = charge_tables.energy_start
            energy_max = charge_tables.energy_start + (charge_tables.num_points - 1) * charge_tables.energy_step
            
            if energy_min <= local_fermi_level <= energy_max:
                # Use table interpolation - this is the primary path
                density = charge_tables.interpolate_bulk_density(region_id, local_fermi_level)
                
                # Check for valid density
                if np.isfinite(density):
                    return density
            
            # Fallback to direct calculation for out-of-range energies
            # Use local_fermi_level as the Fermi level, zero potential in direct calculation
            return self._calculate_bulk_charge_direct(region_id, local_fermi_level)
            
        except Exception as e:
            # Return zero for any error
            return 0.0
    
    def _surface_charge_func(self, r, phi, pot):
        """
        Surface charge density function for Poisson solver.
        
        Args:
            r, phi: Surface coordinates  
            pot: Surface potential
            
        Returns:
            Surface charge density (C/cm^2)
        """
        try:
            # Determine which area
            area_id = self._get_area_id(r, phi)
            
            # Convert pot to float
            pot_val = float(pot)
            
            # Calculate energy
            energy = self.fermi_level - pot_val
            
            # Safety checks
            if not (np.isfinite(pot_val) and np.isfinite(energy)):
                return 0.0
            
            # Get charge tables
            charge_tables = self._current_charge_tables
            
            # Limit energy to table range
            energy_min = charge_tables.energy_start
            energy_max = charge_tables.energy_start + (charge_tables.num_points - 1) * charge_tables.energy_step
            
            if energy_min <= energy <= energy_max:
                density = charge_tables.interpolate_surface_density(area_id, energy)
                
                if np.isfinite(density):
                    return density
            
            # Fallback calculation
            result = self._calculate_surface_charge_direct(area_id, energy)
            
            # Debug surface charge at origin
            if abs(r) < 0.1 and abs(phi) < 0.1:
                print(f"    Surface charge debug: r={r:.3f}, phi={phi:.3f}, pot={pot_val:.6f}, energy={energy:.3f}, density={result:.2e}")
            
            return result
            
        except Exception as e:
            return 0.0
    
    def _calculate_bulk_charge_direct(self, region_id: int, local_fermi_level: float) -> float:
        """
        Direct calculation of bulk charge density for out-of-range energies.
        
        Args:
            region_id: Semiconductor region ID
            local_fermi_level: Local Fermi level (EF - potential)
            
        Returns:
            Charge density (C/cm^3)
        """
        try:
            # Use charge calculator for direct computation with zero potential
            # since the local_fermi_level already accounts for the potential shift
            return self.charge_calculator.calculate_bulk_density(region_id, local_fermi_level, 0.0)
        except:
            return 0.0
    
    def _calculate_surface_charge_direct(self, area_id: int, energy: float) -> float:
        """
        Direct calculation of surface charge density for out-of-range energies.
        """
        try:
            # Calculate potential from energy
            potential = self.fermi_level - energy
            return self.charge_calculator.calculate_surface_density(area_id, self.fermi_level, potential)
        except:
            return 0.0
    
    def _save_results(self):
        """Save simulation results to file."""
        # Save in format compatible with analysis tools
        import pickle
        
        with open('multint_results.pkl', 'wb') as f:
            pickle.dump({
                'config': self.config,
                'results': self.results,
                'fermi_level': self.fermi_level
            }, f)
        
        logger.info("Results saved to multint_results.pkl")
    
    def __del__(self):
        """Clean up resources."""
        if hasattr(self, 'output_file'):
            self.output_file.close()


def run_multint_simulation(config_file: str, max_points: Optional[int] = None) -> List[SimulationResults]:
    """
    Run MultInt simulation from configuration file.
    
    Args:
        config_file: Path to YAML configuration file
        max_points: Maximum number of voltage points to process (for testing)
        
    Returns:
        List of simulation results
    """
    from ..core.filereader import load_yaml_config
    
    # Load configuration
    config = load_yaml_config(config_file)
    
    # Create and run simulation
    simulation = MultIntSimulation(config)
    results = simulation.run(max_points=max_points)
    
    return results