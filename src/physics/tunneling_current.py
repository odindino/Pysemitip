"""
Tunneling Current Calculation Module

This module implements the quantum mechanical tunneling current calculation
for scanning tunneling microscopy based on Bardeen's transfer Hamiltonian
approach. It includes both extended and localized state contributions.

Based on the tunneling current models in intcurr-6.2.f

Author: odindino
"""

import numpy as np
from typing import Dict, Optional, Tuple, Union, List, Callable
from dataclasses import dataclass
from enum import Enum
import warnings

from .materials import SemiconductorMaterial, PhysicalConstants
from .charge_density import ChargeDensityCalculator

# Import numerical functions with proper path handling
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.numerical import fermi_dirac_occupation


class StateType(Enum):
    """Types of electronic states"""
    EXTENDED = "extended"      # Bulk band states
    LOCALIZED = "localized"    # Surface/defect states
    BOTH = "both"              # Combined calculation


@dataclass
class TunnelingConfig:
    """Configuration for tunneling current calculations"""
    
    # Energy integration parameters
    energy_points: int = 50     # NEE in Fortran
    energy_range: float = 2.0   # eV around Fermi level
    
    # k-space integration parameters  
    k_points: int = 20          # NWK in Fortran - parallel momentum points
    k_max_factor: float = 3.0   # Maximum k as factor of Fermi k
    
    # Wavefunction calculation
    wavefunction_points: int = 2048  # Points for wavefunction calculation
    expansion_factor: float = 20.0   # EXPANI - wavefunction expansion factor
    barrier_points: int = 1024       # Points in tunnel barrier
    
    # Physical model options
    include_image_potential: bool = True   # IMPOT flag
    include_valence_bands: bool = True     # Include all valence subbands
    include_conduction_band: bool = True   # Include conduction band
    state_types: StateType = StateType.BOTH
    
    # Numerical parameters
    integration_tolerance: float = 1e-6
    wavefunction_tolerance: float = 1e-8
    max_localized_states: int = 100
    
    # Temperature parameters
    tip_temperature: float = 4.2      # K - tip temperature
    sample_temperature: float = 4.2   # K - sample temperature


class TunnelingCurrentCalculator:
    """
    Quantum mechanical tunneling current calculator for STM.
    
    Implements Bardeen's transfer Hamiltonian method to calculate tunneling
    current between tip and sample, including contributions from both
    extended and localized electronic states.
    
    Current formula: I = (e/ħ) ∫ |M|² ρ_tip(E) ρ_sample(E) [f_tip(E) - f_sample(E)] dE
    """
    
    def __init__(self, config: Optional[TunnelingConfig] = None):
        self.config = config or TunnelingConfig()
        self.constants = PhysicalConstants()
        
        # Calculated quantities
        self.barrier_profile: Optional[np.ndarray] = None
        self.transmission_function: Optional[Callable] = None
        self.tip_density_of_states: float = 1.0  # Normalized tip DOS
        
        # Current components
        self.conduction_current: float = 0.0
        self.valence_current_light: float = 0.0
        self.valence_current_heavy: float = 0.0
        self.valence_current_split_off: float = 0.0
        self.localized_current: float = 0.0
        
    def calculate_barrier_profile(self,
                                potential_1d: np.ndarray,
                                z_grid: np.ndarray,
                                tip_position: float,
                                sample_position: float,
                                tip_work_function: float = 4.5,
                                sample_work_function: float = 4.0) -> np.ndarray:
        """
        Calculate tunnel barrier profile along z-direction.
        
        Args:
            potential_1d: 1D potential profile [eV]
            z_grid: Z-coordinate grid [nm]
            tip_position: Tip position [nm]
            sample_position: Sample surface position [nm]
            tip_work_function: Tip work function [eV]
            sample_work_function: Sample work function [eV]
            
        Returns:
            Barrier potential profile [eV]
        """
        # Base barrier is vacuum level plus electrostatic potential
        barrier = np.zeros_like(z_grid)
        
        # Find indices for tip and sample positions
        tip_idx = np.argmin(np.abs(z_grid - tip_position))
        sample_idx = np.argmin(np.abs(z_grid - sample_position))
        
        # Set barrier between tip and sample
        for i, z in enumerate(z_grid):
            if tip_position <= z <= sample_position:
                # Linear interpolation of work functions
                fraction = (z - tip_position) / (sample_position - tip_position)
                local_work_function = (1 - fraction) * tip_work_function + fraction * sample_work_function
                barrier[i] = local_work_function + potential_1d[i]
            elif z < tip_position:
                barrier[i] = tip_work_function + potential_1d[i]
            else:
                barrier[i] = sample_work_function + potential_1d[i]
        
        # Add image potential correction if enabled
        if self.config.include_image_potential:
            barrier = self._add_image_potential(barrier, z_grid, tip_position, sample_position)
        
        self.barrier_profile = barrier
        return barrier
    
    def _add_image_potential(self,
                           barrier: np.ndarray,
                           z_grid: np.ndarray,
                           tip_position: float,
                           sample_position: float) -> np.ndarray:
        """
        Add image potential correction to barrier.
        
        Implements the image potential model from potexpand-6.1.f
        """
        corrected_barrier = barrier.copy()
        
        # Image potential parameter (from Fortran: lambda factor)
        image_strength = 1.15  # eV⋅nm
        
        barrier_width = sample_position - tip_position
        
        for i, z in enumerate(z_grid):
            if tip_position < z < sample_position:
                # Distance from surfaces
                d_tip = z - tip_position
                d_sample = sample_position - z
                
                # Image potential contribution
                image_correction = -image_strength * barrier_width**2 / (4 * d_tip * d_sample)
                corrected_barrier[i] += image_correction
        
        return corrected_barrier
    
    def calculate_transmission_coefficient(self,
                                        energy: float,
                                        k_parallel: float,
                                        barrier_profile: np.ndarray,
                                        z_grid: np.ndarray) -> float:
        """
        Calculate transmission coefficient through tunnel barrier.
        
        Uses WKB approximation for quantum tunneling.
        
        Args:
            energy: Electron energy [eV]
            k_parallel: Parallel momentum [nm^-1]
            barrier_profile: Barrier potential profile [eV]
            z_grid: Z-coordinate grid [nm]
            
        Returns:
            Transmission coefficient (0 ≤ T ≤ 1)
        """
        # Kinetic energy conversion factor: 2m/ħ² in convenient units
        # C = 26.254 nm^-2⋅eV^-1 from Fortran
        kinetic_factor = self.constants.C_KINETIC
        
        # Perpendicular kinetic energy
        perpendicular_energy = energy - k_parallel**2 / kinetic_factor
        
        # Calculate local wave vector in barrier
        local_k = np.zeros_like(z_grid, dtype=complex)  # Allow complex values
        
        for i, (z, potential) in enumerate(zip(z_grid, barrier_profile)):
            local_energy = perpendicular_energy - potential
            
            if local_energy > 0:
                # Propagating wave
                local_k[i] = np.sqrt(kinetic_factor * local_energy)
            else:
                # Evanescent wave (imaginary k)
                local_k[i] = 1j * np.sqrt(kinetic_factor * abs(local_energy))
        
        # WKB integral for transmission
        # T = exp(-2 * ∫ Im(k) dz)
        integrand = np.abs(np.imag(local_k))  # Ensure positive real values
        wkb_integral = np.trapz(integrand, z_grid)  # Use NumPy trapz
        
        transmission = np.exp(-2 * abs(wkb_integral))  # Ensure real result
        
        return min(float(transmission), 1.0)  # Cap at 1 and ensure float
    
    def calculate_sample_wavefunction(self,
                                    energy: float,
                                    k_parallel: float,
                                    material: SemiconductorMaterial,
                                    band_type: str = "valence") -> Tuple[float, float]:
        """
        Calculate sample wavefunction and its derivative at surface.
        
        Args:
            energy: Electron energy [eV]
            k_parallel: Parallel momentum [nm^-1]
            material: Semiconductor material
            band_type: "valence" or "conduction"
            
        Returns:
            Tuple of (wavefunction_value, wavefunction_derivative)
        """
        # Effective mass for the band
        if band_type == "conduction":
            effective_mass = material.conduction_band_mass
            band_edge = material.bandgap + material.valence_band_offset
        elif band_type == "valence_heavy":
            effective_mass = material.valence_band_mass_heavy
            band_edge = material.valence_band_offset
        elif band_type == "valence_light":
            effective_mass = material.valence_band_mass_light
            band_edge = material.valence_band_offset
        elif band_type == "valence_split_off":
            effective_mass = material.valence_band_mass_split_off
            band_edge = material.valence_band_offset - material.spin_orbit_splitting
        else:
            effective_mass = material.valence_band_mass_eff
            band_edge = material.valence_band_offset
        
        # Kinetic energy in semiconductor
        kinetic_factor = self.constants.C_KINETIC
        perpendicular_energy = energy - band_edge - k_parallel**2 / kinetic_factor
        
        if perpendicular_energy <= 0:
            return 0.0, 0.0
        
        # Wave vector in semiconductor
        k_perp = np.sqrt(kinetic_factor * effective_mass * perpendicular_energy)
        
        # Wavefunction at surface (normalized)
        wavefunction = 1.0  # Normalized plane wave amplitude
        wavefunction_derivative = k_perp  # Real part of d/dz of exp(ikz)
        
        return abs(wavefunction), abs(wavefunction_derivative)
    
    def calculate_density_of_states(self,
                                  material: SemiconductorMaterial,
                                  energy: float,
                                  band_type: str = "valence") -> float:
        """
        Calculate density of states for given band and energy.
        
        Args:
            material: Semiconductor material
            energy: Energy relative to band edge [eV]
            band_type: Band type identifier
            
        Returns:
            Density of states [eV^-1⋅cm^-3]
        """
        if band_type == "conduction":
            effective_mass = material.conduction_band_mass
            band_edge = material.bandgap + material.valence_band_offset
            degeneracy = 1  # Conduction band
        elif band_type == "valence_heavy":
            effective_mass = material.valence_band_mass_heavy
            band_edge = material.valence_band_offset
            degeneracy = 1  # Heavy hole band
        elif band_type == "valence_light":
            effective_mass = material.valence_band_mass_light
            band_edge = material.valence_band_offset
            degeneracy = 1  # Light hole band
        elif band_type == "valence_split_off":
            effective_mass = material.valence_band_mass_split_off
            band_edge = material.valence_band_offset - material.spin_orbit_splitting
            degeneracy = 1  # Split-off band
        else:
            effective_mass = material.valence_band_mass_eff
            band_edge = material.valence_band_offset
            degeneracy = 1
        
        # Energy relative to band edge
        relative_energy = energy - band_edge
        
        if relative_energy <= 0:
            return 0.0
        
        # 3D density of states: g(E) = (1/2π²) * (2m*/ħ²)^(3/2) * √E
        dos_prefactor = (1.0 / (2 * np.pi**2)) * (2 * effective_mass * 0.511e6)**1.5 / (1.97e-7)**3
        dos = dos_prefactor * np.sqrt(relative_energy) * degeneracy
        
        return dos * 1e-21  # Convert to eV^-1⋅cm^-3
    
    def calculate_extended_state_current(self,
                                       material: SemiconductorMaterial,
                                       bias_voltage: float,
                                       fermi_level: float,
                                       barrier_profile: np.ndarray,
                                       z_grid: np.ndarray) -> Dict[str, float]:
        """
        Calculate tunneling current from extended (bulk) states.
        
        Args:
            material: Semiconductor material
            bias_voltage: Applied bias [V]
            fermi_level: Sample Fermi level [eV]
            barrier_profile: Tunnel barrier profile [eV]
            z_grid: Z-coordinate grid [nm]
            
        Returns:
            Dictionary of current components [A]
        """
        current_components = {
            'conduction': 0.0,
            'valence_heavy': 0.0,
            'valence_light': 0.0,
            'valence_split_off': 0.0,
            'total': 0.0
        }
        
        # Energy integration grid
        energy_max = max(fermi_level + bias_voltage, fermi_level) + self.config.energy_range
        energy_min = min(fermi_level + bias_voltage, fermi_level) - self.config.energy_range
        energies = np.linspace(energy_min, energy_max, self.config.energy_points)
        energy_step = energies[1] - energies[0]
        
        # k-parallel integration grid
        # Estimate maximum k from Fermi energy
        thermal_energy = material.thermal_energy
        k_fermi = np.sqrt(2 * material.conduction_band_mass * 0.511e6 * abs(fermi_level)) / 1.97e-7  # nm^-1
        k_max = self.config.k_max_factor * k_fermi if k_fermi > 0 else 1.0  # nm^-1
        
        k_points = np.linspace(0, k_max, self.config.k_points)
        k_step = k_points[1] - k_points[0] if len(k_points) > 1 else 0.0
        
        # Band types to calculate
        band_types = []
        if self.config.include_conduction_band:
            band_types.append('conduction')
        if self.config.include_valence_bands:
            band_types.extend(['valence_heavy', 'valence_light', 'valence_split_off'])
        
        # Current calculation constant
        # I = (e/ħ) * (2π)² * ∫∫ |M|² ρ_tip ρ_sample [f_tip - f_sample] dE dk_parallel
        current_prefactor = self.constants.ELEMENTARY_CHARGE / self.constants.HBAR
        current_prefactor *= (2 * np.pi)**2  # Phase space factor
        
        # Convert to practical units (A)
        current_prefactor *= 1e-21  # Conversion factor
        
        for band_type in band_types:
            band_current = 0.0
            
            for energy in energies:
                for k_parallel in k_points:
                    # Occupations
                    f_tip = fermi_dirac_occupation(
                        energy - bias_voltage, fermi_level, 
                        self.constants.K_B_EV * self.config.tip_temperature
                    )
                    f_sample = fermi_dirac_occupation(
                        energy, fermi_level,
                        self.constants.K_B_EV * self.config.sample_temperature  
                    )
                    
                    occupation_diff = f_tip - f_sample
                    if abs(occupation_diff) < 1e-12:
                        continue
                    
                    # Transmission coefficient
                    transmission = self.calculate_transmission_coefficient(
                        energy, k_parallel, barrier_profile, z_grid
                    )
                    
                    if transmission < 1e-12:
                        continue
                    
                    # Sample wavefunction
                    wf_value, wf_deriv = self.calculate_sample_wavefunction(
                        energy, k_parallel, material, band_type
                    )
                    
                    if wf_value < 1e-12:
                        continue
                    
                    # Density of states
                    sample_dos = self.calculate_density_of_states(material, energy, band_type)
                    
                    if sample_dos < 1e-12:
                        continue
                    
                    # k-space degeneracy factor
                    k_degeneracy = 2 * np.pi * k_parallel if k_parallel > 0 else 1
                    
                    # Matrix element (simplified)
                    matrix_element_squared = float(wf_value**2 * transmission)
                    
                    # Current contribution
                    current_contribution = (current_prefactor * matrix_element_squared * 
                                          self.tip_density_of_states * sample_dos * 
                                          occupation_diff * k_degeneracy * 
                                          energy_step * k_step)
                    
                    band_current += float(current_contribution)
            
            current_components[band_type] = band_current
        
        # Total current
        current_components['total'] = sum(current_components[band] for band in band_types)
        
        return current_components
    
    def calculate_localized_state_current(self,
                                        surface_states: Dict,
                                        bias_voltage: float,
                                        fermi_level: float,
                                        barrier_profile: np.ndarray,
                                        z_grid: np.ndarray) -> float:
        """
        Calculate tunneling current from localized surface states.
        
        Args:
            surface_states: Dictionary of surface state parameters
            bias_voltage: Applied bias [V]
            fermi_level: Sample Fermi level [eV]
            barrier_profile: Tunnel barrier profile [eV]
            z_grid: Z-coordinate grid [nm]
            
        Returns:
            Localized state current [A]
        """
        # Simplified implementation - would need detailed surface state calculation
        # This is a placeholder for the complex localized state analysis
        return 0.0
    
    def calculate_total_current(self,
                              material: SemiconductorMaterial,
                              bias_voltage: float,
                              fermi_level: float,
                              potential_profile: np.ndarray,
                              z_grid: np.ndarray,
                              tip_position: float,
                              sample_position: float,
                              surface_states: Optional[Dict] = None) -> Dict:
        """
        Calculate total tunneling current including all contributions.
        
        Args:
            material: Semiconductor material
            bias_voltage: Applied bias voltage [V]
            fermi_level: Sample Fermi level [eV]
            potential_profile: 1D electrostatic potential [eV]
            z_grid: Z-coordinate grid [nm]
            tip_position: Tip position [nm]
            sample_position: Sample surface position [nm]
            surface_states: Optional surface state parameters
            
        Returns:
            Dictionary with current components and total current
        """
        # Calculate barrier profile
        barrier = self.calculate_barrier_profile(
            potential_profile, z_grid, tip_position, sample_position
        )
        
        # Extended state current
        extended_current = {}
        if self.config.state_types in [StateType.EXTENDED, StateType.BOTH]:
            extended_current = self.calculate_extended_state_current(
                material, bias_voltage, fermi_level, barrier, z_grid
            )
        
        # Localized state current
        localized_current = 0.0
        if (self.config.state_types in [StateType.LOCALIZED, StateType.BOTH] and 
            surface_states is not None):
            localized_current = self.calculate_localized_state_current(
                surface_states, bias_voltage, fermi_level, barrier, z_grid
            )
        
        # Store components
        self.conduction_current = extended_current.get('conduction', 0.0)
        self.valence_current_heavy = extended_current.get('valence_heavy', 0.0)
        self.valence_current_light = extended_current.get('valence_light', 0.0)
        self.valence_current_split_off = extended_current.get('valence_split_off', 0.0)
        self.localized_current = localized_current
        
        # Total current
        total_extended = extended_current.get('total', 0.0)
        total_current = total_extended + localized_current
        
        return {
            'extended_states': extended_current,
            'localized_states': localized_current,
            'total_current': total_current,
            'current_components': {
                'conduction': self.conduction_current,
                'valence_heavy': self.valence_current_heavy,
                'valence_light': self.valence_current_light,
                'valence_split_off': self.valence_current_split_off,
                'localized': self.localized_current
            },
            'barrier_profile': barrier
        }
    
    def calculate_conductance(self,
                            material: SemiconductorMaterial,
                            bias_range: Tuple[float, float],
                            fermi_level: float,
                            potential_profile: np.ndarray,
                            z_grid: np.ndarray,
                            tip_position: float,
                            sample_position: float,
                            bias_points: int = 50) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate differential conductance dI/dV vs bias.
        
        Returns:
            Tuple of (bias_array, conductance_array)
        """
        bias_min, bias_max = bias_range
        bias_array = np.linspace(bias_min, bias_max, bias_points)
        current_array = np.zeros(bias_points)
        
        # Calculate current for each bias
        for i, bias in enumerate(bias_array):
            result = self.calculate_total_current(
                material, bias, fermi_level, potential_profile, 
                z_grid, tip_position, sample_position
            )
            current_array[i] = result['total_current']
        
        # Calculate differential conductance
        conductance_array = np.gradient(current_array, bias_array)
        
        return bias_array, conductance_array
    
    def get_current_summary(self) -> Dict:
        """Get summary of current calculation"""
        return {
            'total_current': (self.conduction_current + self.valence_current_heavy + 
                            self.valence_current_light + self.valence_current_split_off + 
                            self.localized_current),
            'conduction_current': self.conduction_current,
            'valence_current_heavy': self.valence_current_heavy,
            'valence_current_light': self.valence_current_light,
            'valence_current_split_off': self.valence_current_split_off,
            'localized_current': self.localized_current,
            'valence_fraction': (self.valence_current_heavy + self.valence_current_light + 
                               self.valence_current_split_off) / 
                              max(abs(self.conduction_current + self.valence_current_heavy + 
                                    self.valence_current_light + self.valence_current_split_off), 1e-15)
        }


# Convenience functions

def calculate_simple_stm_current(material: SemiconductorMaterial,
                               bias_voltage: float,
                               separation: float = 1.0,
                               tip_work_function: float = 4.5) -> float:
    """
    Calculate STM current for simple geometry with minimal setup.
    
    Returns current in amperes.
    """
    calculator = TunnelingCurrentCalculator()
    
    # Create bias-dependent barrier model
    z_grid = np.linspace(0, separation, 100)
    
    # Find Fermi level
    charge_calc = ChargeDensityCalculator()
    fermi_level = charge_calc.find_equilibrium_fermi_level(material)
    
    # Create realistic barrier profile with bias voltage
    # Use more realistic work function values for STM
    sample_work_function = material.electron_affinity + material.bandgap/2 + fermi_level
    
    # Adjust work functions to be more realistic for STM
    tip_wf = tip_work_function
    sample_wf = 4.0  # Typical semiconductor work function
    
    # Create bias-dependent potential profile
    potential_profile = np.zeros_like(z_grid)
    
    # Simple trapezoidal barrier with bias
    vacuum_level = max(tip_wf, sample_wf + bias_voltage) + 1.0  # 1 eV above highest work function
    
    for i, z in enumerate(z_grid):
        # Linear potential drop from tip to sample due to bias
        bias_drop = bias_voltage * z / separation
        
        # Basic barrier shape (trapezoidal)
        if z < separation * 0.1:  # Near tip
            potential_profile[i] = tip_wf + 1.0  # 1 eV barrier
        elif z > separation * 0.9:  # Near sample
            potential_profile[i] = sample_wf + bias_drop + 1.0  # 1 eV barrier
        else:  # Middle region
            # Linear interpolation with smaller barrier
            tip_potential = tip_wf + 1.0
            sample_potential = sample_wf + bias_voltage + 1.0
            potential_profile[i] = tip_potential + (sample_potential - tip_potential) * z / separation
        
        # Add small image potential correction
        image_correction = -0.1 / (4 * z + 0.05) if z > 0.01 else 0  # Smaller correction
        potential_profile[i] += image_correction
    
    result = calculator.calculate_total_current(
        material, bias_voltage, fermi_level, potential_profile,
        z_grid, 0.0, separation
    )
    
    return result['total_current']


if __name__ == "__main__":
    # Demo usage
    from .materials import default_materials
    
    # Get silicon material
    si_n = default_materials.get_material("Si_n")
    
    # Create tunneling calculator
    calculator = TunnelingCurrentCalculator()
    
    print(f"Calculating STM current for {si_n.name}")
    print(f"Doping: {si_n.net_doping:.1e} cm⁻³")
    
    # Simple current calculation
    bias_voltage = 1.0  # V
    current = calculate_simple_stm_current(si_n, bias_voltage)
    
    print(f"STM current at {bias_voltage} V bias: {current:.2e} A")
    
    # Calculate I-V curve
    biases = np.linspace(-2, 2, 21)
    currents = [calculate_simple_stm_current(si_n, v) for v in biases]
    
    print(f"\nI-V curve:")
    for v, i in zip(biases[::4], currents[::4]):  # Every 4th point
        print(f"  V = {v:+.1f} V: I = {i:.2e} A")
