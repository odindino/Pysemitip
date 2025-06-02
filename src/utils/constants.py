"""
Physical constants and unit conversions for SEMITIP simulations.

This module contains fundamental physical constants used throughout the simulation,
matching the values used in the original Fortran code.
"""

import numpy as np


class PhysicalConstants:
    """Collection of physical constants in SI units."""
    
    # Fundamental constants
    E = 1.60217663e-19          # Elementary charge (C)
    KB = 1.38064852e-23         # Boltzmann constant (J/K)
    H = 6.62607015e-34          # Planck constant (J·s)
    HBAR = H / (2 * np.pi)      # Reduced Planck constant (J·s)
    M0 = 9.1093837015e-31       # Electron mass (kg)
    EPSILON0 = 8.8541878128e-12 # Vacuum permittivity (F/m)
    
    # Useful combinations
    KB_EV = KB / E              # Boltzmann constant (eV/K)
    HBAR_EV = HBAR / E          # Reduced Planck constant (eV·s)
    
    # Unit conversions
    EV_TO_J = E                 # Convert eV to Joules
    J_TO_EV = 1.0 / E           # Convert Joules to eV
    NM_TO_M = 1e-9              # Convert nanometers to meters
    M_TO_NM = 1e9               # Convert meters to nanometers
    CM3_TO_M3 = 1e-6            # Convert cm^-3 to m^-3
    CM2_TO_M2 = 1e-4            # Convert cm^-2 to m^-2
    
    # Temperature conversion factor (K to eV)
    # Used in Fortran: TK = TEM * 8.617E-5
    K_TO_EV = 8.617333262e-5    # More precise value
    
    # Mathematical constants
    PI = np.pi
    
    @classmethod
    def thermal_energy(cls, temperature_k):
        """Calculate thermal energy kT in eV."""
        return cls.KB_EV * temperature_k
    
    @classmethod
    def thermal_voltage(cls, temperature_k):
        """Calculate thermal voltage kT/e in Volts."""
        return cls.KB * temperature_k / cls.E
    
    @classmethod
    def fermi_integral_minus_half(cls, eta, temperature_k):
        """
        Calculate Fermi-Dirac integral of order -1/2.
        
        F_{-1/2}(η) = ∫[0,∞] x^{-1/2} / (1 + exp(x - η)) dx
        
        This is used for carrier density calculations.
        """
        # This will be implemented using scipy.special functions
        # For now, return a placeholder
        raise NotImplementedError("Fermi integral calculation to be implemented")


class MaterialConstants:
    """Material-specific constants and properties."""
    
    # GaAs default parameters (matching Fortran defaults)
    GAAS = {
        'band_gap': 1.42,           # eV at 300K
        'affinity': 4.07,           # eV
        'permittivity': 12.9,       # Relative permittivity
        'cb_effective_mass': 0.0635, # Conduction band effective mass (m*/m0)
        'vb_effective_mass_heavy': 0.5, # Heavy hole mass
        'vb_effective_mass_light': 0.076, # Light hole mass
        'vb_effective_mass_so': 0.145,    # Split-off band mass
        'spin_orbit_splitting': 0.341,    # eV
        'donor_binding_energy': 0.0058,   # eV (Si donor)
        'acceptor_binding_energy': 0.031, # eV (C acceptor)
    }
    
    # Silicon parameters
    SI = {
        'band_gap': 1.12,           # eV at 300K
        'affinity': 4.05,           # eV
        'permittivity': 11.7,       # Relative permittivity
        'cb_effective_mass': 0.26,  # Average CB effective mass
        'vb_effective_mass_heavy': 0.49,
        'vb_effective_mass_light': 0.16,
        'donor_binding_energy': 0.045,    # eV (P donor)
        'acceptor_binding_energy': 0.045, # eV (B acceptor)
    }


class ComputationalConstants:
    """Constants related to numerical computation."""
    
    # Default grid parameters
    DEFAULT_GRID = {
        'radial_points': 201,
        'vacuum_points': 251,
        'semiconductor_points': 401,
        'angular_points': 61,
    }
    
    # Convergence criteria
    CONVERGENCE = {
        'potential_tolerance': 1e-6,  # V
        'charge_tolerance': 1e-10,    # Relative
        'current_tolerance': 1e-12,   # A
        'max_iterations': 10000,
    }
    
    # Integration parameters
    INTEGRATION = {
        'energy_points': 20000,       # Number of energy points for DOS
        'k_parallel_points': 20,      # Number of parallel k-points
        'expansion_factor': 20,       # For wave function expansion
    }


# Convenience aliases
PC = PhysicalConstants
MC = MaterialConstants
CC = ComputationalConstants

# Commonly used values
E = PC.E
KB = PC.KB
HBAR = PC.HBAR
M0 = PC.M0
EPSILON0 = PC.EPSILON0
PI = PC.PI


def print_constants():
    """Print all physical constants for verification."""
    print("Physical Constants:")
    print(f"  Elementary charge: {PC.E:.6e} C")
    print(f"  Boltzmann constant: {PC.KB:.6e} J/K")
    print(f"  Boltzmann constant: {PC.KB_EV:.6e} eV/K")
    print(f"  Planck constant: {PC.H:.6e} J·s")
    print(f"  Reduced Planck constant: {PC.HBAR:.6e} J·s")
    print(f"  Electron mass: {PC.M0:.6e} kg")
    print(f"  Vacuum permittivity: {PC.EPSILON0:.6e} F/m")
    print(f"  Temperature conversion: {PC.K_TO_EV:.6e} eV/K")
    
    print("\nGaAs Material Constants:")
    for key, value in MC.GAAS.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    print_constants()