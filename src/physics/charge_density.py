"""
Charge Density Calculation Module

This module implements comprehensive charge density calculations for semiconductor
systems, including bulk charge density from free carriers and ionized dopants,
and surface charge density from surface states.

Based on the physics models in semirhomult-6.0.f and surfrhomult-6.2.f

Author: odindino
"""

import numpy as np
from typing import Dict, Optional, Tuple, Union, List, Callable
from dataclasses import dataclass
import warnings

from .materials import (
    SemiconductorMaterial, SurfaceStateParameters, MaterialParameters,
    PhysicalConstants
)

# Import numerical functions with proper path handling
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.numerical import fermi_dirac_integral, fermi_dirac_occupation


@dataclass
class ChargeDensityConfig:
    """Configuration for charge density calculations"""
    
    # Numerical integration parameters
    energy_resolution: float = 0.001  # eV - energy step for integration
    integration_range: float = 10.0  # kT - integration range around Fermi level
    lookup_table_size: int = 50000  # Size of charge density lookup table (NEDIM)
    
    # Physical model switches
    include_ionized_dopants: bool = True
    include_band_gap_narrowing: bool = False
    include_quantum_corrections: bool = False
    
    # Convergence parameters
    charge_neutrality_tolerance: float = 1e12  # cm^-3
    max_iterations: int = 1000


class ChargeDensityCalculator:
    """
    Comprehensive charge density calculator for semiconductor systems.
    
    This class implements the charge density models from the original Fortran
    code, including:
    - Free electron density in conduction band (RHOCB function)
    - Free hole density in valence band (RHOVB function) 
    - Ionized dopant densities
    - Surface state charge densities (RHOS function)
    
    Maintains numerical consistency with the original implementation.
    """
    
    def __init__(self, config: Optional[ChargeDensityConfig] = None):
        self.config = config or ChargeDensityConfig()
        self.constants = PhysicalConstants()
        
        # Charge density lookup tables (like RHOBTAB, RHOSTAB in Fortran)
        self._bulk_charge_table: Dict[int, np.ndarray] = {}
        self._surface_charge_table: Dict[int, np.ndarray] = {}
        self._energy_grid: Optional[np.ndarray] = None
        
    def calculate_electron_density(self, 
                                 material: SemiconductorMaterial,
                                 fermi_level: float,
                                 potential: Union[float, np.ndarray],
                                 region_id: int = 1) -> Union[float, np.ndarray]:
        """
        Calculate free electron density in conduction band.
        
        Implements RHOCB function from semirhomult-6.0.f
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
            potential: Electrostatic potential [eV]
            region_id: Material region identifier
            
        Returns:
            Electron density [cm^-3]
        """
        # Conduction band edge relative to valence band
        conduction_band_edge = material.bandgap + material.valence_band_offset
        
        # Effective density of states constant (from Fortran: C = 6.815E21)
        # C = (2/√π) * 2 * (m*kT/(2πℏ²))^(3/2) in eV^-1.5 cm^-3
        dos_constant = self.constants.CHARGE_DENSITY_CONSTANT
        
        # Energy difference: E_F - E_c - qV
        energy_diff = fermi_level - conduction_band_edge - potential
        
        if material.temperature == 0.0:
            # Zero temperature case - step function
            n = np.zeros_like(energy_diff) if hasattr(energy_diff, '__len__') else 0.0
            positive_mask = energy_diff > 0
            
            if hasattr(energy_diff, '__len__'):
                n[positive_mask] = (2.0 * dos_constant / 3.0) * \
                    np.sqrt((material.conduction_band_mass * energy_diff[positive_mask])**3)
            else:
                if positive_mask:
                    n = (2.0 * dos_constant / 3.0) * \
                        np.sqrt((material.conduction_band_mass * energy_diff)**3)
            
            return n
        else:
            # Finite temperature case - Fermi-Dirac integral
            thermal_energy = material.thermal_energy
            effective_mass_term = (material.conduction_band_mass * thermal_energy)**1.5
            
            # Calculate Fermi-Dirac integral F_{1/2}(η)
            eta = energy_diff / thermal_energy
            fd_integral = fermi_dirac_integral(0.5, eta)  # j=1 corresponds to F_{1/2}
            
            n = dos_constant * np.sqrt(effective_mass_term**3) * fd_integral
            
            return np.maximum(n, 0.0)  # Ensure non-negative
    
    def calculate_hole_density(self,
                             material: SemiconductorMaterial, 
                             fermi_level: float,
                             potential: Union[float, np.ndarray],
                             region_id: int = 1) -> Union[float, np.ndarray]:
        """
        Calculate free hole density in valence band.
        
        Implements RHOVB function from semirhomult-6.0.f
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
            potential: Electrostatic potential [eV] 
            region_id: Material region identifier
            
        Returns:
            Hole density [cm^-3]
        """
        # Valence band edge (reference point)
        valence_band_edge = material.valence_band_offset
        
        # Effective density of states constant
        dos_constant = self.constants.CHARGE_DENSITY_CONSTANT
        
        # Energy difference: E_v - E_F + qV (note sign convention for holes)
        energy_diff = valence_band_edge - fermi_level + potential
        
        if material.temperature == 0.0:
            # Zero temperature case
            p = np.zeros_like(energy_diff) if hasattr(energy_diff, '__len__') else 0.0
            positive_mask = energy_diff > 0
            
            if hasattr(energy_diff, '__len__'):
                p[positive_mask] = (2.0 * dos_constant / 3.0) * \
                    np.sqrt((material.valence_band_mass_eff * energy_diff[positive_mask])**3)
            else:
                if positive_mask:
                    p = (2.0 * dos_constant / 3.0) * \
                        np.sqrt((material.valence_band_mass_eff * energy_diff)**3)
            
            return p
        else:
            # Finite temperature case
            thermal_energy = material.thermal_energy
            effective_mass_term = (material.valence_band_mass_eff * thermal_energy)**1.5
            
            # Calculate Fermi-Dirac integral F_{1/2}(η)
            eta = energy_diff / thermal_energy
            fd_integral = fermi_dirac_integral(0.5, eta)
            
            p = dos_constant * np.sqrt(effective_mass_term**3) * fd_integral
            
            return np.maximum(p, 0.0)
    
    def calculate_ionized_donor_density(self,
                                      material: SemiconductorMaterial,
                                      fermi_level: float,
                                      potential: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate ionized donor density N_D^+.
        
        Uses Fermi-Dirac statistics for donor ionization.
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
            potential: Electrostatic potential [eV]
            
        Returns:
            Ionized donor density [cm^-3]
        """
        if material.donor_concentration == 0.0:
            return np.zeros_like(potential) if hasattr(potential, '__len__') else 0.0
        
        # Donor energy level
        donor_level = material.bandgap + material.valence_band_offset - material.donor_binding_energy
        
        # Ionization probability: f_D = 1 / (1 + 2*exp((E_D - E_F + qV)/(kT)))
        # Factor of 2 accounts for spin degeneracy
        thermal_energy = material.thermal_energy
        
        if thermal_energy > 0:
            energy_diff = (donor_level - fermi_level + potential) / thermal_energy
            # Avoid overflow in exponential
            energy_diff = np.clip(energy_diff, -50, 50)
            ionization_prob = 1.0 / (1.0 + 2.0 * np.exp(energy_diff))
        else:
            # Zero temperature: step function
            ionization_prob = (fermi_level - potential) > donor_level
        
        return material.donor_concentration * ionization_prob
    
    def calculate_ionized_acceptor_density(self,
                                         material: SemiconductorMaterial,
                                         fermi_level: float, 
                                         potential: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate ionized acceptor density N_A^-.
        
        Uses Fermi-Dirac statistics for acceptor ionization.
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
            potential: Electrostatic potential [eV]
            
        Returns:
            Ionized acceptor density [cm^-3]
        """
        if material.acceptor_concentration == 0.0:
            return np.zeros_like(potential) if hasattr(potential, '__len__') else 0.0
        
        # Acceptor energy level
        acceptor_level = material.valence_band_offset + material.acceptor_binding_energy
        
        # Ionization probability: f_A = 1 / (1 + 4*exp((E_F - E_A - qV)/(kT)))
        # Factor of 4 accounts for spin degeneracy  
        thermal_energy = material.thermal_energy
        
        if thermal_energy > 0:
            energy_diff = (fermi_level - acceptor_level - potential) / thermal_energy
            # Avoid overflow
            energy_diff = np.clip(energy_diff, -50, 50)
            ionization_prob = 1.0 / (1.0 + 4.0 * np.exp(energy_diff))
        else:
            # Zero temperature: step function
            ionization_prob = (fermi_level - potential) < acceptor_level
        
        return material.acceptor_concentration * ionization_prob
    
    def calculate_bulk_charge_density(self,
                                    material: SemiconductorMaterial,
                                    fermi_level: float,
                                    potential: Union[float, np.ndarray],
                                    region_id: int = 1) -> Union[float, np.ndarray]:
        """
        Calculate total bulk charge density.
        
        ρ_bulk = -n + p + N_D^+ - N_A^-
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
            potential: Electrostatic potential [eV]
            region_id: Material region identifier
            
        Returns:
            Bulk charge density [e·cm^-3] (positive for positive charge)
        """
        # Free carriers
        n = self.calculate_electron_density(material, fermi_level, potential, region_id)
        p = self.calculate_hole_density(material, fermi_level, potential, region_id)
        
        # Ionized dopants
        if self.config.include_ionized_dopants:
            n_d_plus = self.calculate_ionized_donor_density(material, fermi_level, potential)
            n_a_minus = self.calculate_ionized_acceptor_density(material, fermi_level, potential)
        else:
            n_d_plus = np.zeros_like(potential) if hasattr(potential, '__len__') else 0.0
            n_a_minus = np.zeros_like(potential) if hasattr(potential, '__len__') else 0.0
        
        # Total bulk charge density (in units of elementary charge)
        rho_bulk = -n + p + n_d_plus - n_a_minus
        
        return rho_bulk
    
    def calculate_surface_charge_density(self,
                                       surface_params: SurfaceStateParameters,
                                       fermi_level: float,
                                       surface_potential: float,
                                       region_id: int = 1) -> float:
        """
        Calculate surface charge density from surface states.
        
        Implements RHOS function from surfrhomult-6.2.f
        
        Args:
            surface_params: Surface state parameters
            fermi_level: Fermi level relative to valence band [eV]
            surface_potential: Surface potential [eV]
            region_id: Surface region identifier
            
        Returns:
            Surface charge density [e·cm^-2]
        """
        if surface_params is None:
            return 0.0
        
        # Energy reference: Fermi level at surface
        ef_surface = fermi_level - surface_potential
        
        # Calculate donor-like and acceptor-like surface state charges
        rho_donor = self._calculate_surface_state_charge(
            surface_params.donor_density,
            surface_params.donor_energy,
            surface_params.donor_center_energy,
            surface_params.donor_width,
            ef_surface,
            is_donor=True,
            temperature_dependent=surface_params.temperature_dependent
        )
        
        rho_acceptor = self._calculate_surface_state_charge(
            surface_params.acceptor_density,
            surface_params.acceptor_energy,
            surface_params.acceptor_center_energy,
            surface_params.acceptor_width,
            ef_surface,
            is_donor=False,
            temperature_dependent=surface_params.temperature_dependent
        )
        
        return rho_donor + rho_acceptor
    
    def _calculate_surface_state_charge(self,
                                      density: float,
                                      characteristic_energy: float,
                                      center_energy: float,
                                      width: float,
                                      fermi_level: float,
                                      is_donor: bool,
                                      temperature_dependent: bool,
                                      thermal_energy: float = 0.026) -> float:
        """
        Calculate charge from a single type of surface state.
        
        Implements the integration in RHOS1 and RHOS2 functions.
        """
        if density == 0.0:
            return 0.0
        
        # Energy integration parameters
        energy_step = self.config.energy_resolution
        
        if not temperature_dependent:
            # Zero temperature case - simple integration up to Fermi level
            if is_donor:
                # Donor states: positively charged when empty (E > E_F)
                if fermi_level < characteristic_energy:
                    # All states above Fermi level are empty (positively charged)
                    return density
                else:
                    # States below Fermi level are filled (neutral)
                    return 0.0
            else:
                # Acceptor states: negatively charged when filled (E < E_F)
                if fermi_level > characteristic_energy:
                    # All states below Fermi level are filled (negatively charged)
                    return -density
                else:
                    # States above Fermi level are empty (neutral)  
                    return 0.0
        else:
            # Finite temperature case - integrate with Fermi function
            # Use Gaussian distribution of surface states
            
            # Integration limits
            e_min = center_energy - 5 * width
            e_max = center_energy + 5 * width
            n_points = int((e_max - e_min) / energy_step) + 1
            energies = np.linspace(e_min, e_max, n_points)
            
            # Gaussian density of states
            dos = density * np.exp(-(energies - center_energy)**2 / (2 * width**2))
            dos /= (np.sqrt(2 * np.pi) * width)  # Normalize
            
            # Fermi-Dirac occupation
            occupation = fermi_dirac_occupation(energies, fermi_level, thermal_energy)
            
            if is_donor:
                # Donor charge: +e when empty, 0 when filled
                charge_integrand = dos * (1 - occupation)
                return np.trapz(charge_integrand, energies)
            else:
                # Acceptor charge: -e when filled, 0 when empty
                charge_integrand = -dos * occupation
                return np.trapz(charge_integrand, energies)
    
    def build_charge_density_lookup_table(self,
                                        material: SemiconductorMaterial,
                                        potential_range: Tuple[float, float],
                                        fermi_level: float,
                                        region_id: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build lookup table for charge density vs potential.
        
        Equivalent to building RHOBTAB in the Fortran code.
        
        Args:
            material: Semiconductor material parameters
            potential_range: (min_potential, max_potential) in eV
            fermi_level: Fermi level relative to valence band [eV]
            region_id: Material region identifier
            
        Returns:
            Tuple of (potential_grid, charge_density_grid)
        """
        v_min, v_max = potential_range
        n_points = self.config.lookup_table_size
        
        potential_grid = np.linspace(v_min, v_max, n_points)
        charge_density_grid = np.zeros(n_points)
        
        for i, potential in enumerate(potential_grid):
            charge_density_grid[i] = self.calculate_bulk_charge_density(
                material, fermi_level, potential, region_id
            )
        
        # Store in lookup table
        self._bulk_charge_table[region_id] = charge_density_grid
        if self._energy_grid is None:
            self._energy_grid = potential_grid
        
        return potential_grid, charge_density_grid
    
    def interpolate_charge_density(self,
                                 potential: Union[float, np.ndarray],
                                 region_id: int = 1) -> Union[float, np.ndarray]:
        """
        Interpolate charge density from lookup table.
        
        Args:
            potential: Electrostatic potential(s) [eV]
            region_id: Material region identifier
            
        Returns:
            Interpolated charge density [e·cm^-3]
        """
        if region_id not in self._bulk_charge_table or self._energy_grid is None:
            raise ValueError(f"No lookup table found for region {region_id}. "
                           "Call build_charge_density_lookup_table first.")
        
        charge_grid = self._bulk_charge_table[region_id]
        
        # Linear interpolation
        return np.interp(potential, self._energy_grid, charge_grid)
    
    def calculate_charge_neutrality_condition(self,
                                            material: SemiconductorMaterial,
                                            potential: float = 0.0) -> Callable[[float], float]:
        """
        Create charge neutrality condition function for finding equilibrium Fermi level.
        
        The charge neutrality condition is:
        n - p - N_D^+ + N_A^- = 0
        
        Args:
            material: Semiconductor material parameters
            potential: Fixed potential [eV] (default: flat band)
            
        Returns:
            Function that takes Fermi level and returns charge imbalance
        """
        def charge_balance(fermi_level: float) -> float:
            """Charge balance function for root finding"""
            n = self.calculate_electron_density(material, fermi_level, potential)
            p = self.calculate_hole_density(material, fermi_level, potential)
            n_d_plus = self.calculate_ionized_donor_density(material, fermi_level, potential)
            n_a_minus = self.calculate_ionized_acceptor_density(material, fermi_level, potential)
            
            return n - p - n_d_plus + n_a_minus
        
        return charge_balance
    
    def find_equilibrium_fermi_level(self,
                                   material: SemiconductorMaterial,
                                   potential: float = 0.0) -> float:
        """
        Find equilibrium Fermi level using charge neutrality condition.
        
        Args:
            material: Semiconductor material parameters
            potential: Fixed potential [eV]
            
        Returns:
            Equilibrium Fermi level [eV]
        """
        from scipy.optimize import brentq
        
        balance_func = self.calculate_charge_neutrality_condition(material, potential)
        
        # Search range: typically within bandgap ± some margin
        ef_min = -0.5  # Below valence band
        ef_max = material.bandgap + 0.5  # Above conduction band
        
        try:
            fermi_level = brentq(balance_func, ef_min, ef_max, 
                               xtol=1e-12, maxiter=self.config.max_iterations)
        except ValueError as e:
            warnings.warn(f"Failed to find charge neutrality: {e}. "
                         f"Using intrinsic Fermi level.")
            # Fallback: intrinsic Fermi level
            fermi_level = material.bandgap / 2.0
        
        return fermi_level


# Convenience functions for common use cases

def calculate_intrinsic_fermi_level(material: SemiconductorMaterial) -> float:
    """
    Calculate intrinsic Fermi level.
    
    For parabolic bands: E_Fi = E_v + E_g/2 + (3/4)kT ln(m_v*/m_c*)
    """
    thermal_energy = material.thermal_energy
    mass_ratio = material.valence_band_mass_eff / material.conduction_band_mass
    
    if thermal_energy > 0:
        correction = 0.75 * thermal_energy * np.log(mass_ratio)
    else:
        correction = 0.0
    
    return material.bandgap / 2.0 + correction


def calculate_debye_length(material: SemiconductorMaterial, 
                         carrier_concentration: float) -> float:
    """
    Calculate Debye screening length.
    
    L_D = sqrt(ε*kT / (e²*n))
    
    Args:
        material: Semiconductor material parameters
        carrier_concentration: Total carrier concentration [cm^-3]
        
    Returns:
        Debye length [nm]
    """
    if carrier_concentration <= 0:
        return np.inf
    
    constants = PhysicalConstants()
    
    # Permittivity in F/cm
    permittivity = constants.VACUUM_PERMITTIVITY * material.relative_permittivity * 1e-2
    
    # Thermal energy in J
    thermal_energy_j = constants.BOLTZMANN_CONSTANT * material.temperature
    
    # Debye length in cm
    debye_length_cm = np.sqrt(permittivity * thermal_energy_j / 
                             (constants.ELEMENTARY_CHARGE**2 * carrier_concentration))
    
    # Convert to nm
    return debye_length_cm * 1e7


if __name__ == "__main__":
    # Demo usage
    from .materials import default_materials
    
    # Get silicon n-type material
    si_n = default_materials.get_material("Si_n")
    
    # Create charge density calculator
    calculator = ChargeDensityCalculator()
    
    # Find equilibrium Fermi level
    ef_eq = calculator.find_equilibrium_fermi_level(si_n)
    print(f"Equilibrium Fermi level: {ef_eq:.3f} eV")
    
    # Calculate charge densities at equilibrium
    potential = 0.0  # Flat band
    n = calculator.calculate_electron_density(si_n, ef_eq, potential)
    p = calculator.calculate_hole_density(si_n, ef_eq, potential)
    
    print(f"Electron density: {n:.2e} cm⁻³")
    print(f"Hole density: {p:.2e} cm⁻³")
    print(f"n·p product: {n*p:.2e} cm⁻⁶")
    print(f"Intrinsic n_i²: {si_n.intrinsic_concentration**2:.2e} cm⁻⁶")
    
    # Calculate charge density vs potential
    potentials = np.linspace(-0.5, 0.5, 100)
    charge_densities = [calculator.calculate_bulk_charge_density(si_n, ef_eq, v) 
                       for v in potentials]
    
    print(f"\nCharge density range: [{min(charge_densities):.2e}, {max(charge_densities):.2e}] cm⁻³")
    
    # Calculate Debye length
    debye_length = calculate_debye_length(si_n, si_n.donor_concentration)
    print(f"Debye length: {debye_length:.1f} nm")
