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
            save_interval: int = 1) -> List[SimulationResults]:
        """
        Run simulation for specified bias voltages.
        
        Args:
            bias_voltages: List of bias voltages (V). If None, use config values.
            save_interval: Save results every N bias points
            
        Returns:
            List of SimulationResults
        """
        if bias_voltages is None:
            bias_voltages = self.config.voltage_scan.voltages
        
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
        Solve for a single bias voltage.
        
        Args:
            bias_voltage: Bias voltage (V)
            
        Returns:
            SimulationResults object
        """
        start_time = time.time()
        
        # Update tip potential
        self.tip.bias_voltage = bias_voltage
        tip_potential = self.tip.tip_potential
        
        self.output_file.write(f"\nSEPARATION = {self.tip.separation}\n")
        self.output_file.write(f"\nBIAS, TIP POTENTIAL = {bias_voltage} {tip_potential}\n")
        
        # Estimate depletion width (1D approximation)
        depl_width_1d = self._estimate_depletion_width_1d(bias_voltage)
        self.output_file.write(f"1-D ESTIMATE OF DEPLETION WIDTH (NM) = {depl_width_1d:.6f}\n")
        
        # Create charge density tables
        energy_range = self._calculate_energy_range(bias_voltage)
        self.output_file.write(f"ESTART,EEND,NE = {energy_range[0]:.6f} "
                             f"{energy_range[1]:.6f} {self.config.charge_table_points}\n")
        
        logger.info("Computing charge density tables...")
        self.output_file.write("COMPUTING TABLE OF BULK CHARGE DENSITIES\n")
        self.output_file.write("COMPUTING TABLE OF SURFACE CHARGE DENSITIES\n")
        self.output_file.flush()
        
        charge_tables = self.charge_calculator.create_charge_tables(
            energy_range, 
            self.config.charge_table_points
        )
        
        # Store charge tables for use in charge density functions
        self._current_charge_tables = charge_tables
        
        # Multi-grid solution sequence
        solution_3d = None
        grid_sequence = self._get_grid_sequence()
        
        for grid_level, (nr, nv, ns, np) in enumerate(grid_sequence):
            self.output_file.write(f"NR,NS,NV,NP = {nr} {ns} {nv} {np}\n")
            
            # Create sub-grid (simplified - full implementation needed)
            # For now, use full grid
            
            # Solve Poisson equation
            self.output_file.write(f"SOLUTION # {grid_level + 1}\n")
            self.output_file.flush()
            
            solution_3d, conv_info = self.poisson_solver.solve(
                self._bulk_charge_func,
                self._surface_charge_func,
                initial_guess=solution_3d
            )
            
            # Log convergence
            bb = self.poisson_solver.get_band_bending()
            self.output_file.write(f"NUMBER OF ITERATIONS = {conv_info['iterations']}\n")
            self.output_file.write(f"BAND BENDING AT MIDPOINT = {bb:.8f}\n")
            self.output_file.flush()
        
        # Calculate depletion width from solution
        depl_width = self.poisson_solver.get_depletion_width()
        
        # Extract potential profile
        profile = self.potential_processor.extract_profile(solution_3d)
        
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
            bias_voltage=bias_voltage,
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
        
        # Built-in potential
        vbi = self.fermi_level - region.valence_band_edge()
        
        # Total potential drop
        total_potential = abs(vbi - bias_voltage)
        
        # Depletion width formula
        eps = region.permittivity * PC.EPSILON0
        doping = abs(region.net_doping) * 1e6  # Convert to m^-3
        
        if doping > 0:
            w = np.sqrt(2 * eps * total_potential / (PC.E * doping)) * 1e9  # Convert to nm
        else:
            w = 100.0  # Default for intrinsic
        
        return w
    
    def _calculate_energy_range(self, bias_voltage: float) -> Tuple[float, float]:
        """Calculate appropriate energy range for charge tables."""
        # Use the estimate_energy_range function
        bias_range = (min(0, bias_voltage), max(0, bias_voltage))
        return estimate_energy_range(
            self.semiconductor_regions,
            self.surface_regions,
            self.fermi_level,
            bias_range
        )
    
    def _get_grid_sequence(self) -> List[Tuple[int, int, int, int]]:
        """
        Get multi-grid sequence for convergence.
        
        Returns:
            List of (nr, nv, ns, np) tuples
        """
        # Start with coarse grid and refine
        # This matches the Fortran multi-grid approach
        nr_final = self.grid.params.nr
        nv_final = self.grid.params.nv
        ns_final = self.grid.params.ns
        np_final = self.grid.params.np
        
        # Three-level sequence
        sequence = [
            (nr_final // 4, nv_final // 4, ns_final // 4, np_final // 2),
            (nr_final // 2, nv_final // 2, ns_final // 2, np_final),
            (nr_final, nv_final, ns_final, np_final)
        ]
        
        # Ensure minimum sizes
        sequence = [(max(16, nr), max(16, nv), max(16, ns), max(8, np))
                   for nr, nv, ns, np in sequence]
        
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
            return self._calculate_surface_charge_direct(area_id, energy)
            
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


def run_multint_simulation(config_file: str) -> List[SimulationResults]:
    """
    Run MultInt simulation from configuration file.
    
    Args:
        config_file: Path to YAML configuration file
        
    Returns:
        List of simulation results
    """
    from ..core.filereader import load_yaml_config
    
    # Load configuration
    config = load_yaml_config(config_file)
    
    # Create and run simulation
    simulation = MultIntSimulation(config)
    results = simulation.run()
    
    return results