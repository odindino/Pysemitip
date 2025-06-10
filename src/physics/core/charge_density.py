"""
Charge density calculations for SEMITIP simulations.

This module implements the charge density calculations from SEMIRHOMULT
and SURFRHOMULT, including tabulation of charge densities for efficient
evaluation during Poisson equation solving.
"""
import numpy as np # Ensure numpy is imported
from scipy import special, optimize
 # Added for type hinting if potential_matrix_Volts is np.ndarray

from typing import Dict, List, Tuple, Optional, TYPE_CHECKING, Any # Ensure Dict, Optional, Any are imported
from dataclasses import dataclass, field

from ...utils.constants import PhysicalConstants as PC
# from ..materials.semiconductor import SemiconductorRegion # Physics material class
# from ..materials.surface_states import SurfaceRegion # Physics material class
from ...core.config_schema import SemitipConfig, TipConfig, VoltageScanConfig

if TYPE_CHECKING:
    from ..materials.semiconductor import SemiconductorRegion as SemiconductorRegionPhysics
    from ..materials.surface_states import SurfaceRegion as SurfaceRegionPhysics # Assuming this is the physics class for surface states
    # from ...core.config_schema import SemiconductorRegion as ConfigSemiconductorRegion


@dataclass
class ChargeDensityTable:
    """
    Stores the tabulated charge density for a single semiconductor region.
    The potential axis represents the local Fermi level relative to the local
    valence band maximum (ef_local_rel_local_vb).
    """
    region_id: int
    # Potential axis relative to the region's valence band maximum (ef_val_rel_vb in Fortran)
    # This is E_F - E_V_local where E_V_local = Ev_abs - q*Potential_local
    # So, ef_val_rel_vb = E_F_system - (Ev_abs_region - q*Potential_local)
    potential_axis_rel_vb: np.ndarray = field(default_factory=lambda: np.array([])) # Units: eV
    # Corresponding charge densities (rho in C/cm^3)
    charge_density_axis: np.ndarray = field(default_factory=lambda: np.array([])) # Units: C/cm^3

    # Store band edges and Fermi level used for this table for context/debugging
    Ev_abs_eV: Optional[float] = None # Absolute VB max energy (at V=0) for this region
    Ec_abs_eV: Optional[float] = None # Absolute CB min energy (at V=0) for this region
    system_fermi_level_eV: Optional[float] = None # System Fermi level E_F_main


class ChargeDensityCalculator:
    """
    Manages charge density calculations for all semiconductor and surface regions.
    """

    def __init__(
        self,
        config_or_regions,  # Can be SemitipConfig or List of regions (backward compatibility)
        semiconductor_regions_physics_or_surface=None, # Can be Dict[int, SemiconductorRegionPhysics] or List (backward compatibility)
        surface_regions_physics: Optional[Dict[int, 'SurfaceRegionPhysics']] = None, # Use SurfaceRegionPhysics type hint
        system_fermi_level_E_F_main: float = 0.0,  # System Fermi level in eV
        fermi_level: Optional[float] = None  # Backward compatibility parameter
    ):
        """
        Initializes the ChargeDensityCalculator.
        
        Supports two calling patterns:
        1. New: ChargeDensityCalculator(config, semiconductor_regions_physics, surface_regions_physics, system_fermi_level)
        2. Old: ChargeDensityCalculator([regions], [], fermi_level=fermi_level)

        Args:
            config_or_regions: Either SemitipConfig object or List of semiconductor regions (backward compatibility)
            semiconductor_regions_physics_or_surface: Either Dict[int, SemiconductorRegionPhysics] or empty list (backward compatibility)
            surface_regions_physics: (Optional) A dictionary mapping region_id to SurfaceRegion (physics) objects.
            system_fermi_level_E_F_main: The system's main Fermi level (EF_main) in eV.
            fermi_level: Backward compatibility parameter, used when called with old signature
        """
        # Detect calling pattern
        if isinstance(config_or_regions, list) and fermi_level is not None:
            # Old calling pattern: ChargeDensityCalculator([regions], [], fermi_level=fermi_level)
            regions_list = config_or_regions
            self.config = None  # No config in old pattern
            self.semiconductor_regions_physics = {region.region_id: region for region in regions_list}
            self.surface_regions_physics = {}
            self.system_fermi_level_E_F_main_eV = fermi_level
        else:
            # New calling pattern
            self.config = config_or_regions
            self.semiconductor_regions_physics = semiconductor_regions_physics_or_surface
            self.surface_regions_physics = surface_regions_physics or {}
            self.system_fermi_level_E_F_main_eV = system_fermi_level_E_F_main

        # Ev_abs_for_charge_calc_eV and Ec_abs_for_charge_calc_eV are now attributes
        # of the SemiconductorRegionPhysics objects themselves, calculated in their __post_init__.
        # No need to calculate and setattr here. They can be accessed via:
        # region_physics = self.semiconductor_regions_physics[region_id]
        # ev_abs = region_physics.Ev_abs_for_charge_calc_eV
        # ec_abs = region_physics.Ec_abs_for_charge_calc_eV
        
        # Log to confirm access if needed during debugging:
        # for r_id, r_phys in self.semiconductor_regions_physics.items():
        #     logger.debug(f"Region {r_id}: Ev_abs={r_phys.Ev_abs_for_charge_calc_eV}, Ec_abs={r_phys.Ec_abs_for_charge_calc_eV}")

        self.charge_density_tables: Dict[int, ChargeDensityTable] = {}
        
        # Only build tables if we have config (new calling pattern)
        if self.config is not None:
            self._build_charge_density_tables() # This method will populate the tables

    def _estimate_energy_range_for_region(
        self,
        # region_config: SemiconductorRegionConfig, # This was config class, now we use physics class
        region_physics: 'SemiconductorRegionPhysics',
        tip_config: TipConfig, # from main config
        voltage_scan_config: Optional[VoltageScanConfig] # from main config
    ) -> Tuple[float, float]:
        """
        Estimates the energy range (min and max ef_local_rel_local_vb) for a given
        semiconductor region's charge density table. This estimation is crucial for
        ensuring the table covers all relevant local Fermi level variations relative
        to the local valence band maximum during the simulation.

        The method aims to replicate the core logic of Fortran's ESTART/EEND estimation
        but adapted for per-region tables. The table axis, potential_axis_rel_vb,
        represents ef_local_rel_local_vb = E_F_local - Ev_local.
        Based on the problem description, E_F_local is typically E_F_main (system Fermi level),
        and Ev_local = Ev_abs_zero_potential - V_local (where V_local is local electrostatic potential).
        Thus, ef_local_rel_local_vb = E_F_main - (Ev_abs_zero_potential - V_local)
                                   = (E_F_main - Ev_abs_zero_potential) + V_local
                                   = region_physics.ef_vb_equilibrium + V_local.
        Alternatively, as stated in the original docstring for this function:
        ef_local_rel_local_vb = region_physics.ef_vb_equilibrium - V_local_alt,
        where V_local_alt is the local electrostatic potential relative to E_F_main, defined such
        that positive V_local_alt lowers electron energies. If V_local_alt = -V_local, the forms are consistent.
        The implementation uses ef_vb_equilibrium - V_local_potential_wrt_EF_main.

        Args:
            region_config: Configuration for the semiconductor region.
            region_physics: Physics properties for the semiconductor region.
            tip_config: Tip configuration.
            voltage_scan_config: Voltage scan configuration.

        Returns:
            A tuple (min_energy_ev, max_energy_ev) for the potential_axis_rel_vb in eV.
        """
        # en0_region is the equilibrium Fermi level for this specific region,
        # relative to its own valence band maximum (E_F_equil_region - Ev_bulk_region).
        # This is a key reference point for this region's table.
        en0_region = region_physics.fermi_level(potential=0.0)  # Calculate equilibrium Fermi level

        # 1. Determine the range of applied bias voltages.
        min_bias = 0.0
        max_bias = 0.0
        if voltage_scan_config and voltage_scan_config.num_steps > 0:
            # Ensure num_steps is at least 1 for linspace if start and end are different,
            # or handle single point scan (num_steps=1 means start_voltage only).
            if voltage_scan_config.num_steps == 1:
                voltages = np.array([voltage_scan_config.start_voltage_V])
            else:
                voltages = np.linspace(voltage_scan_config.start_voltage_V,
                                       voltage_scan_config.end_voltage_V,
                                       voltage_scan_config.num_steps)
            min_bias = np.min(voltages)
            max_bias = np.max(voltages)
        elif self.config.computation: # Single bias from computation config
            min_bias = self.config.computation.bias_voltage_V
            max_bias = self.config.computation.bias_voltage_V
        
        # 2. Calculate Contact Potential (CPot) for this region.
        # CPot = tip_work_function - semiconductor_bulk_work_function
        # semiconductor_bulk_work_function = Chi_s + (Ec_bulk - EF_bulk)
        #                                = electron_affinity + Eg - ef_vb_equilibrium
        semiconductor_wf = region_physics.electron_affinity + region_physics.Eg - en0_region
        contact_potential = tip_config.work_function_eV - semiconductor_wf

        # 3. Calculate PotTIP: Tip potential relative to the system Fermi level (E_F_main).
        # PotTIP = applied_bias + contact_potential. This is the effective potential
        # difference driving band bending near the tip.
        min_pot_tip = min_bias + contact_potential
        max_pot_tip = max_bias + contact_potential

        # 4. Estimate the range of ef_local_rel_local_vb.
        # ef_local_rel_local_vb = en0_region - V_local, where V_local is the local
        # electrostatic potential at a point in the semiconductor, relative to E_F_main.
        # V_local can be approximated by:
        #   - 0 (deep in the bulk, far from the tip influence)
        #   - min_pot_tip (potential at semiconductor surface under tip influence, one extreme)
        #   - max_pot_tip (potential at semiconductor surface under tip influence, other extreme)
        # These approximations help define the initial span of ef_local_rel_local_vb.
        
        val_at_bulk_potential = en0_region  # V_local = 0
        # If V_local = max_pot_tip, ef_local_rel_local_vb = en0_region - max_pot_tip
        val_at_min_tip_potential_effect = en0_region - max_pot_tip
        # If V_local = min_pot_tip, ef_local_rel_local_vb = en0_region - min_pot_tip
        val_at_max_tip_potential_effect = en0_region - min_pot_tip
        
        # Initial estimates for the range of ef_local_rel_local_vb, similar to Fortran's
        # use of AMIN1/AMAX1 with relevant energy levels for the specific region.
        estart_initial = min(val_at_bulk_potential, val_at_min_tip_potential_effect, val_at_max_tip_potential_effect)
        eend_initial = max(val_at_bulk_potential, val_at_min_tip_potential_effect, val_at_max_tip_potential_effect)

        # 5. Add initial fixed padding (Fortran often uses +/- 0.5 eV before expansion).
        # This value can be configured in the simulation settings.
        initial_padding_ev = self.config.computation.charge_density_table_initial_padding_eV \
            if self.config.computation and self.config.computation.charge_density_table_initial_padding_eV is not None \
            else 0.5  # Default based on common Fortran practice

        estart_padded = estart_initial - initial_padding_ev
        eend_padded = eend_initial + initial_padding_ev
        
        # 6. Fortran-style range expansion.
        # ETMP = EEND_padded - ESTART_padded
        # ESTART_new = ESTART_padded - expansion_factor * ETMP
        # EEND_new = EEND_padded + expansion_factor * ETMP
        # A common expansion_factor (N in Fortran) is 2, leading to a ~5x wider range.
        # This expansion provides a safety margin for interpolation accuracy and
        # to cover unexpected potential variations during the Poisson solver's iterations.
        # This factor can be configured in the simulation settings.
        fortran_expansion_factor = self.config.computation.charge_density_table_expansion_factor \
            if self.config.computation and self.config.computation.charge_density_table_expansion_factor is not None \
            else 2.0  # Default to N=2 for 5x expansion, a common Fortran heuristic

        etmp = eend_padded - estart_padded
        
        # Handle cases where etmp might be zero or very small (e.g., intrinsic material at zero bias).
        # This ensures a minimum width for ETMP before expansion.
        min_etmp_width = 1e-3 # eV
        if etmp < min_etmp_width:
            # Fallback: ensure a minimal sensible range.
            # If bandgap is significant, use it as a basis for ETMP.
            # Otherwise, use a default span (e.g., 1 eV).
            if region_physics.Eg > 0.1: # eV
                etmp_fallback = region_physics.Eg 
            else: # For metals or very small gap materials
                etmp_fallback = 1.0 # Default 1 eV span
            
            # Recenter the padded range around its midpoint using the new etmp_fallback
            mid_point = (estart_padded + eend_padded) / 2.0
            estart_padded = mid_point - etmp_fallback / 2.0
            eend_padded = mid_point + etmp_fallback / 2.0
            etmp = eend_padded - estart_padded # Recalculate etmp, should now be etmp_fallback

        final_min_energy = estart_padded - fortran_expansion_factor * etmp
        final_max_energy = eend_padded + fortran_expansion_factor * etmp
        
        # 7. Sanity checks and constraints.
        # Ensure the range is physically sensible and adequately covers band edges.
        # ef_local_rel_local_vb = 0 corresponds to Ev_local.
        # ef_local_rel_local_vb = Eg corresponds to Ec_local.
        # The table should extend somewhat into the valence and conduction bands.
        # These bounds ensure the table covers at least a certain range around Ev and Ec.
        # For example, from well within the valence band to well within the conduction band.
        # `max(0.5 * Eg, 1.0)` ensures a margin of at least 1eV or 0.5*Eg.
        margin_factor = 0.5 # Multiplier for Eg
        min_abs_margin = 1.0 # eV, absolute minimum margin
        
        # Sensible lower bound: Should be below Ev (0) by a margin.
        # (e.g., Ev - (0.5*Eg + padding) or Ev - (1eV + padding))
        min_sensible_bound = -max(margin_factor * region_physics.Eg, min_abs_margin) - initial_padding_ev
        
        # Sensible upper bound: Should be above Ec (Eg) by a margin.
        # (e.g., Ec + (0.5*Eg + padding) or Ec + (1eV + padding))
        max_sensible_bound = region_physics.Eg + max(margin_factor * region_physics.Eg, min_abs_margin) + initial_padding_ev
        
        # Ensure the calculated final range is at least as wide as these sensible bounds.
        # If expanded range is narrower, widen it to sensible bounds.
        # If expanded range is wider, keep the wider (more conservative) range.
        final_min_energy = min(final_min_energy, min_sensible_bound)
        final_max_energy = max(final_max_energy, max_sensible_bound)

        # 8. Final fallback: Ensure min is strictly less than max.
        # This handles pathological cases (e.g., Eg is zero/negative, or extreme parameters).
        if final_max_energy <= final_min_energy:
            # Fallback to an arbitrary sensible range around the region's equilibrium Fermi level.
            final_min_energy = en0_region - 3.0 # eV
            final_max_energy = en0_region + region_physics.Eg + 3.0 # eV
            if final_max_energy <= final_min_energy: # Ultimate fallback for extreme cases
                 final_min_energy = -4.0 # eV
                 final_max_energy = 4.0  # eV
        
        # Debug printing (can be enabled if needed by uncommenting print statements in the original code)
        # print(f"Region {region_config.id} ({region_physics.name}):")
        # print(f"  ef_vb_equilibrium (en0_region): {en0_region:.3f} eV")
        # print(f"  Bias: [{min_bias:.3f}, {max_bias:.3f}] V, CPot: {contact_potential:.3f} eV")
        # print(f"  PotTIP range: [{min_pot_tip:.3f}, {max_pot_tip:.3f}] eV (relative to E_F_main)")
        # print(f"  Initial ef_local_rel_local_vb estimates (val_bulk, val_minPotTIP, val_maxPotTIP): {val_at_bulk_potential:.3f}, {val_at_min_tip_potential:.3f}, {val_at_max_tip_potential:.3f}")
        # print(f"  estart_initial: {estart_initial:.3f}, eend_initial: {eend_initial:.3f}")
        # print(f"  Padded estart/eend: {estart_padded:.3f}, {eend_padded:.3f} (padding: {initial_padding_ev:.2f} eV)")
        # print(f"  ETMP: {etmp:.3f} eV")
        # print(f"  Fortran expansion factor: {fortran_expansion_factor}")
        # print(f"  Final energy range for table (ef_local_rel_local_vb): [{final_min_energy:.3f}, {final_max_energy:.3f}] eV")
        # print(f"  Sensible bounds check: min_sensible={min_sensible_bound:.3f}, max_sensible={max_sensible_bound:.3f}")

        return final_min_energy, final_max_energy

    def _build_charge_density_tables(self):
        """
        Builds the charge density lookup tables for each semiconductor region.
        The table maps ef_local_rel_local_vb to charge density.
        """
        if not self.config.computation:
            # Should not happen if config is validated
            return

        num_energy_points = self.config.computation.charge_density_table_size
        if num_energy_points <= 1:
            num_energy_points = 200 # Default fallback

        for region_id, region_physics in self.semiconductor_regions_physics.items():
            min_energy, max_energy = self._estimate_energy_range_for_region(
                region_physics,
                self.config.tip,
                self.config.voltage_scan
            )
            
            # potential_axis_rel_vb is ef_local relative to local Ev
            potential_axis_rel_vb = np.linspace(min_energy, max_energy, num_energy_points)
            
            rho_values = np.array([
                self._calculate_rho_for_table_entry(ef_val_rel_vb, region_id)
                for ef_val_rel_vb in potential_axis_rel_vb
            ])

            self.charge_density_tables[region_id] = ChargeDensityTable(
                region_id=region_id,
                potential_axis_rel_vb=potential_axis_rel_vb,
                charge_density_axis=rho_values
            )
            # print(f"Built charge density table for region {region_id}:")
            # print(f"  Potential axis (ef_local_rel_local_vb): {potential_axis_rel_vb.min():.2f} to {potential_axis_rel_vb.max():.2f} eV, {len(potential_axis_rel_vb)} points")
            # print(f"  Rho values: {rho_values.min():.2e} to {rho_values.max():.2e} C/m^3")


    def _calculate_rho_for_table_entry(
        self,
        ef_val_rel_vb: float, # This is the local Fermi level relative to the local valence band max (E_F_local - Ev_local)
        region_id: int
    ) -> float:
        """
        Calculates the total charge density for a given local Fermi level relative
        to the local valence band maximum (ef_val_rel_vb) for a specific region.
        This is used to populate the charge density tables.

        Args:
            ef_val_rel_vb: The local Fermi level relative to the local valence band maximum (E_F_local - Ev_local), in eV.
            region_id: The ID of the semiconductor region.

        Returns:
            Total charge density (rho_total) in C/m^3.
        """
        region_physics = self.semiconductor_regions_physics[region_id]
        
        # Calculate electron and hole concentrations
        # n = Nc * F_1/2((E_F_local - Ec_local) / kT)
        # p = Nv * F_1/2((Ev_local - E_F_local) / kT)
        #
        # We have ef_val_rel_vb = E_F_local - Ev_local
        # So, E_F_local - Ec_local = E_F_local - (Ev_local + Eg) = (E_F_local - Ev_local) - Eg = ef_val_rel_vb - Eg
        # And Ev_local - E_F_local = -(E_F_local - Ev_local) = -ef_val_rel_vb

        kT_eV = PC.KB_EV * region_physics.temperature

        eta_c = (ef_val_rel_vb - region_physics.band_gap) / kT_eV
        eta_v = -ef_val_rel_vb / kT_eV # Note: F_1/2 argument is (Ev-EF)/kT

        n = region_physics.Nc * region_physics._simple_fermi_integral(eta_c)  # cm^-3
        p = region_physics.Nv * region_physics._simple_fermi_integral(eta_v)  # cm^-3
        
        # Ionized dopants
        # For donors: N_D_plus = N_D / (1 + g_D * exp((E_F_local - E_D_local) / kT))
        # E_D_local is the donor energy level. Ed_binding is Ec_local - E_D_local.
        # So E_D_local = Ec_local - Ed_binding = (Ev_local + Eg) - Ed_binding
        # E_F_local - E_D_local = E_F_local - (Ev_local + Eg - Ed_binding)
        #                       = (E_F_local - Ev_local) - Eg + Ed_binding
        #                       = ef_val_rel_vb - Eg + Ed_binding
        N_D_plus = 0.0
        if region_physics.donor_concentration > 1e6 and region_physics.donor_binding_energy is not None: # Check if donors are present
            # donor_binding_energy is the ionization energy (usually positive, Ec - Ed)
            # The energy level Ed is band_gap - donor_binding_energy relative to Ev
            # Ef - Ed = ef_val_rel_vb - (band_gap - donor_binding_energy)
            donor_exp_arg = (ef_val_rel_vb - region_physics.band_gap + region_physics.donor_binding_energy) / kT_eV
            N_D_plus = region_physics.donor_concentration / (1.0 + 2.0 * np.exp(donor_exp_arg))

        # For acceptors: N_A_minus = N_A / (1 + g_A * exp((E_A_local - E_F_local) / kT))
        # E_A_local is the acceptor energy level. Ea_binding is E_A_local - Ev_local.
        # So E_A_local = Ev_local + Ea_binding
        # E_A_local - E_F_local = (Ev_local + Ea_binding) - E_F_local
        #                       = Ea_binding - (E_F_local - Ev_local)
        #                       = Ea_binding - ef_val_rel_vb
        N_A_minus = 0.0
        if region_physics.acceptor_concentration > 1e6 and region_physics.acceptor_binding_energy is not None: # Check if acceptors are present
            # acceptor_binding_energy is the ionization energy (usually positive, Ea - Ev)
            # The energy level Ea is acceptor_binding_energy relative to Ev
            # Ea - Ef = acceptor_binding_energy - ef_val_rel_vb
            acceptor_exp_arg = (region_physics.acceptor_binding_energy - ef_val_rel_vb) / kT_eV
            N_A_minus = region_physics.acceptor_concentration / (1.0 + 4.0 * np.exp(acceptor_exp_arg))
            
        rho_total_SI = PC.E * (p - n + N_D_plus - N_A_minus) # C/m^3
        return rho_total_SI


    def calculate_total_charge_density_at_point_direct(
        self,
        potential_V: float, # Electrostatic potential at the point (Volts)
        region_id: int
    ) -> float:
        """
        Calculates the total charge density at a specific point given the local
        electrostatic potential, by directly evaluating the physical formulas.
        This is used by the Poisson solver during its non-linear iteration.

        Args:
            potential_V: Local electrostatic potential at the point in Volts.
                         (V=0 corresponds to the system's main Fermi level E_F_main).
            region_id: The ID of the semiconductor region.

        Returns:
            Total charge density (rho_total) in C/m^3.
        """
        # ef_val_for_rho is the local Fermi level relative to the local valence band maximum (E_F_local - Ev_local)
        # E_F_local = E_F_main - q*potential_V  (if q is positive elementary charge)
        # Ev_local = Ev_abs_for_charge_calc[region_id] - e*potential_V
        # So, E_F_local - Ev_local = (E_F_main - e*potential_V) - (Ev_abs_for_charge_calc[region_id] - e*potential_V)
        #                          = E_F_main - Ev_abs_for_charge_calc[region_id]
        # This is region_physics.ef_vb_equilibrium. This is only if potential_V is defined differently.

        # Let's use the definition from `calculate` method's `ef_query_values`:
        # ef_query_values = self.E_F_main - potential_matrix_Volts - self.Ev_abs_for_charge_calc[region_id]
        # This ef_query_value is ef_local_rel_local_vb.
        # Here, potential_V is equivalent to potential_matrix_Volts at a point.
        ef_local_rel_local_vb = self.E_F_main - potential_V - self.Ev_abs_for_charge_calc[region_id]
        
        return self._calculate_rho_for_table_entry(ef_local_rel_local_vb, region_id)


    def calculate_charge_density(
        self,
        potential_matrix_Volts: np.ndarray, # 2D or 3D array of potentials
        region_id_map: np.ndarray,          # Same shape as potential_matrix_Volts
                                            # Identifies the region for each grid point
        charge_output_data: Optional[Dict] = None # For debugging output
    ) -> np.ndarray:
        """
        Calculates the total charge density for all points in the grid using
        interpolation from the pre-built tables.

        Args:
            potential_matrix_Volts: Grid of local electrostatic potential values (Volts).
                                    V=0 corresponds to the system's main Fermi level E_F_main.
            region_id_map: Grid identifying the semiconductor region ID for each point.
                           Points with ID <= 0 are typically not semiconductor regions.
            charge_output_data: Optional dictionary to store intermediate data for debugging.


        Returns:
            Grid of total charge densities (rho_total) in C/m^3.
        """
        rho_matrix_SI = np.zeros_like(potential_matrix_Volts, dtype=float)

        for region_id, table in self.charge_density_tables.items():
            mask = (region_id_map == region_id)
            if not np.any(mask):
                continue

            # ef_query_values: local Fermi level relative to local valence band max (E_F_local - Ev_local)
            # E_F_local = E_F_main - potential_matrix_Volts  (if potential_V is energy in eV, or e*potential_V if V is volts)
            # Let's be careful with units. potential_matrix_Volts is in Volts.
            # E_F_local_eV = self.E_F_main (eV) - potential_matrix_Volts (V)
            # Ev_local_eV  = self.Ev_abs_for_charge_calc[region_id] (eV) - potential_matrix_Volts (V)
            # So, ef_query_values = E_F_local_eV - Ev_local_eV
            #                     = (self.E_F_main - potential_matrix_Volts[mask]) - (self.Ev_abs_for_charge_calc[region_id] - potential_matrix_Volts[mask])
            #                     = self.E_F_main - self.Ev_abs_for_charge_calc[region_id]
            # This is ef_vb_equilibrium for the region, which is constant. This cannot be right.
            # The potential_matrix_Volts *does* shift the local bands.

            # Correct calculation for ef_query_values (ef_local relative to local Ev):
            # E_F_local = E_F_main - potential_matrix_Volts[mask] (energy of Fermi level at point, relative to E_F_main=0 if potential is 0)
            # Ev_local  = self.Ev_abs_for_charge_calc[region_id] - potential_matrix_Volts[mask] (energy of VB at point)
            # ef_query_values = E_F_local - Ev_local
            #                 = (self.E_F_main - potential_matrix_Volts[mask]) - (self.Ev_abs_for_charge_calc[region_id] - potential_matrix_Volts[mask])
            #                 = self.E_F_main - self.Ev_abs_for_charge_calc[region_id]

            # Let's re-derive E_F_local and Ev_local carefully.
            # E_F_main is the absolute energy of the system's Fermi level.
            # potential_matrix_Volts (V) is the local electrostatic potential. V=0 means energy is E_F_main.
            # If V > 0, local energy levels are shifted down by V (electron energy = -e * V).
            # If V < 0, local energy levels are shifted up by |V|.
            # So, E_F_local_abs = E_F_main - potential_matrix_Volts[mask] (if we treat potential as energy shift for e=1)
            # This is what was used in `calculate_total_charge_density_at_point_direct`.
            # And Ev_abs_local = self.Ev_abs_for_charge_calc[region_id] - potential_matrix_Volts[mask]
            # So, ef_local_rel_local_vb = E_F_local_abs - Ev_abs_local
            #                            = (self.E_F_main - potential_matrix_Volts[mask]) - (self.Ev_abs_for_charge_calc[region_id] - potential_matrix_Volts[mask])
            #                            = self.E_F_main - self.Ev_abs_for_charge_calc[region_id]
            # This is ef_vb_equilibrium for that region. This is constant. This is the error from before.

            # The `potential_axis_rel_vb` of the table IS `ef_local_rel_local_vb`.
            # The value we need to query this axis with is the actual `ef_local_rel_local_vb` at each point.
            # `ef_local_rel_local_vb` = E_F_local - Ev_local
            # E_F_local is the local Fermi energy. In self-consistency, this is E_F_main everywhere.
            # Ev_local is the local valence band energy: Ev_local = Ev_abs_zero_potential - potential_V
            # So, ef_local_rel_local_vb = E_F_main - (Ev_abs_zero_potential - potential_V)
            #                           = (E_F_main - Ev_abs_zero_potential) + potential_V
            #                           = ef_vb_equilibrium_region + potential_V
            # This was the correction made in `calculate_total_charge_density_at_point_direct` (implicitly).
            # Let's verify `calculate_total_charge_density_at_point_direct` again:
            # `ef_local_rel_local_vb = self.E_F_main - potential_V - self.Ev_abs_for_charge_calc[region_id]`
            # This is (E_F_main - Ev_abs_for_charge_calc[region_id]) - potential_V
            # = ef_vb_equilibrium_region - potential_V. This seems correct.

            ef_query_values = (self.E_F_main - self.Ev_abs_for_charge_calc[region_id]) - potential_matrix_Volts[mask]

            # Interpolate
            # np.interp requires xp to be increasing. Table axis should be.
            rho_matrix_SI[mask] = np.interp(
                ef_query_values,
                table.potential_axis_rel_vb,
                table.charge_density_axis,
                left=table.charge_density_axis[0],  # Extrapolate with edge values
                right=table.charge_density_axis[-1]
            )
            
            if charge_output_data is not None and region_id not in charge_output_data:
                charge_output_data[region_id] = {}
            if charge_output_data is not None:
                charge_output_data[region_id]['ef_query_values'] = ef_query_values
                charge_output_data[region_id]['potentials_at_mask'] = potential_matrix_Volts[mask]
                charge_output_data[region_id]['rho_at_mask'] = rho_matrix_SI[mask]
                charge_output_data[region_id]['table_axis'] = table.potential_axis_rel_vb
                charge_output_data[region_id]['table_rho'] = table.charge_density_axis


        # Handle surface charges if any (simplified for now)
        # This part would need actual SurfaceRegion physics and tables
        for region_id, surf_region_physics in self.surface_regions_physics.items():
            mask = (region_id_map == region_id)
            if not np.any(mask) or not hasattr(surf_region_physics, 'calculate_surface_charge_density'):
                continue
            
            # Example: surface charge might depend on local Fermi level at surface
            # This needs a proper implementation based on surfrhomult
            # For now, assuming it's a fixed charge or simple model
            # ef_surface = self.E_F_main - potential_matrix_Volts[mask] # Simplified local E_F relative to E_F_main
            # rho_matrix_SI[mask] += surf_region_physics.calculate_surface_charge_density(ef_surface)
            pass # Placeholder

        return rho_matrix_SI

    def _calculate_equilibrium_fermi_level_for_region(self, region_id: int) -> float:
        """
        Calculates the equilibrium Fermi level (relative to its own valence band maximum)
        for a given semiconductor region by enforcing charge neutrality.

        Args:
            region_id: The ID of the semiconductor region.

        Returns:
            Equilibrium Fermi level (ef_vb_equilibrium) in eV, relative to Ev_region.
        """
        region_physics = self.semiconductor_regions_physics[region_id]

        # Objective function: total charge density, should be zero at equilibrium
        def charge_neutrality_condition(ef_rel_vb_trial):
            # ef_rel_vb_trial is E_F - Ev for this region
            return self._calculate_rho_for_table_entry(ef_rel_vb_trial, region_id) / PC.e_C # in carriers/m^3

        # Estimate search range for ef_rel_vb_trial
        # Fermi level should be within the band gap, or slightly outside for degenerate cases
        # Search from slightly below Ev (e.g., Ev - 0.5 eV) to slightly above Ec (e.g., Ec + 0.5 eV)
        # Relative to Ev, this is -0.5 eV to Eg + 0.5 eV
        search_min_ev = -0.5 
        search_max_ev = region_physics.Eg + 0.5
        
        try:
            # Check if the function crosses zero in the interval
            rho_min = charge_neutrality_condition(search_min_ev)
            rho_max = charge_neutrality_condition(search_max_ev)
            if np.sign(rho_min) == np.sign(rho_max) and rho_min !=0 :
                # Try expanding the search range if no sign change
                # This can happen for highly intrinsic or very heavily doped materials
                # print(f"Warning: Charge neutrality for region {region_id} might not bracket zero in [{search_min_ev:.2f}, {search_max_ev:.2f}] eV. Rho: [{rho_min:.2e}, {rho_max:.2e}]. Expanding search.")
                search_min_ev -= region_physics.Eg # Expand down
                search_max_ev += region_physics.Eg # Expand up
                rho_min_expanded = charge_neutrality_condition(search_min_ev)
                rho_max_expanded = charge_neutrality_condition(search_max_ev)
                if np.sign(rho_min_expanded) == np.sign(rho_max_expanded) and rho_min_expanded != 0:
                    # print(f"Error: Still no zero crossing for region {region_id} in expanded range [{search_min_ev:.2f}, {search_max_ev:.2f}] ev. Rho: [{rho_min_expanded:.2e}, {rho_max_expanded:.2e}]. Using mid-gap as fallback.")
                    # Fallback: intrinsic Fermi level or mid-gap (approximation)
                    # For simplicity, if heavily N-doped, EF is near Ec. If P-doped, near Ev.
                    if region_physics.Nd > region_physics.Na:
                        return region_physics.Eg - 0.05 # Close to Ec
                    elif region_physics.Na > region_physics.Nd:
                        return 0.05 # Close to Ev
                    else: # Intrinsic or balanced
                        return region_physics.Eg / 2.0 
                else: # Use expanded range
                    rho_min, rho_max = rho_min_expanded, rho_max_expanded


            ef_vb_equilibrium, r = optimize.brentq(
                charge_neutrality_condition,
                search_min_ev,
                search_max_ev,
                xtol=1e-7, rtol=1e-7, # Tighter tolerances
                full_output=True
            )
            if not r.converged:
                print(f"Warning: Brentq for equilibrium Fermi level (region {region_id}) did not converge. Result: {ef_vb_equilibrium:.4f}")
            return ef_vb_equilibrium
        except ValueError as e:
            # This typically means f(a) and f(b) must have different signs for brentq
            print(f"Error finding equilibrium Fermi level for region {region_id} with brentq: {e}")
            print(f"  Search range [{search_min_ev:.2f}, {search_max_ev:.2f}] e")
            print(f"  Rho at bounds: rho({search_min_ev:.2f}) = {charge_neutrality_condition(search_min_ev):.2e}, rho({search_max_ev:.2f}) = {charge_neutrality_condition(search_max_ev):.2e}")
            # Fallback to a simpler estimate, e.g., intrinsic level or mid-gap
            # This should be handled more robustly, perhaps by adjusting search range or method
            return region_physics.Eg / 2.0 # Simplistic fallback

    def _electron_density_direct(self, region_physics: 'SemiconductorRegionPhysics', fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate electron density directly using Fermi-Dirac statistics.
        
        Args:
            region_physics: The semiconductor region physics object
            fermi_level: Fermi level position (eV) relative to valence band maximum  
            potential: Electrostatic potential (V)
            
        Returns:
            Electron density (cm^-3)
        """
        # Convert fermi_level and potential to ef_val_rel_vb
        ef_val_rel_vb = fermi_level - potential  # Convert V to eV assuming e=1
        
        kT_eV = PC.KB_EV * region_physics.temperature
        eta_c = (ef_val_rel_vb - region_physics.band_gap) / kT_eV
        
        # Calculate electron density in cm^-3 using the custom Fermi-Dirac integral
        fermi_integral = region_physics._simple_fermi_integral(eta_c)
        n_cm3 = region_physics.Nc * fermi_integral
        
        return n_cm3

    def _hole_density_direct(self, region_physics: 'SemiconductorRegionPhysics', fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate hole density directly using Fermi-Dirac statistics.
        
        Args:
            region_physics: The semiconductor region physics object
            fermi_level: Fermi level position (eV) relative to valence band maximum
            potential: Electrostatic potential (V)
            
        Returns:
            Hole density (cm^-3)
        """
        # Convert fermi_level and potential to ef_val_rel_vb
        ef_val_rel_vb = fermi_level - potential  # Convert V to eV assuming e=1
        
        kT_eV = PC.KB_EV * region_physics.temperature
        eta_v = -ef_val_rel_vb / kT_eV  # Note: F_1/2 argument is (Ev-EF)/kT
        
        # Calculate hole density in cm^-3 using the custom Fermi-Dirac integral
        fermi_integral = region_physics._simple_fermi_integral(eta_v)
        p_cm3 = region_physics.Nv * fermi_integral
        
        return p_cm3

    def calculate_bulk_density(self, region_id: int, energy: float, potential: float = 0.0) -> float:
        """
        Calculate bulk charge density for a region at a given energy and potential.
        This method is used by tests and provides a simple interface to calculate
        charge density at a specific point.
        
        Args:
            region_id: The ID of the semiconductor region
            energy: The Fermi energy level (eV) relative to valence band maximum
            potential: The electrostatic potential (V)
            
        Returns:
            Bulk charge density (C/m^3)
        """
        # Use the table-based calculation if available
        if region_id in self.charge_density_tables:
            # Convert energy and potential to ef_val_rel_vb
            ef_val_rel_vb = energy - potential  # Assuming 1V = 1eV for simplicity
            
            table = self.charge_density_tables[region_id]
            # Interpolate from the table
            rho_SI = np.interp(
                ef_val_rel_vb,
                table.potential_axis_rel_vb,
                table.charge_density_axis,
                left=table.charge_density_axis[0],
                right=table.charge_density_axis[-1]
            )
            return rho_SI
        else:
            # Fallback to direct calculation
            return self._calculate_rho_for_table_entry(energy - potential, region_id)

# Example usage (for testing purposes, if run directly)
if __name__ == '__main__':
    # This section would require mock objects for Config, SemiconductorRegionPhysics etc.
    # For example:
    # mock_region_physics = SemiconductorRegionPhysics(...)
    # mock_config = SemitipConfig(...)
    # calculator = ChargeDensityCalculator(mock_config, {1: mock_region_physics})
    # print(f"Table for region 1: {calculator.charge_density_tables[1].potential_axis_rel_vb}")
    print("ChargeDensityCalculator module loaded. Run through simulation scripts for full testing.")