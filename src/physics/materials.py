"""
Materials and Parameters Management Module

This module provides a comprehensive materials database and parameter management
system for semiconductor materials used in STM simulations. It encapsulates all
the physical parameters needed for electrostatic and transport calculations.

Based on the physics analysis in Phase1.5_Physics_and_Computation_Deep_Dive.md

Author: odindino
"""

import numpy as np
from typing import Dict, Optional, List, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
import json
import warnings


class RegionType(Enum):
    """Semiconductor region types"""
    BULK = "bulk"
    SURFACE = "surface" 
    INTERFACE = "interface"


@dataclass
class SemiconductorMaterial:
    """
    Comprehensive semiconductor material properties.
    
    All parameters correspond to COMMON /SEMI/ block variables in Fortran code.
    Energy units: eV, Length units: nm, Concentration units: cm^-3
    """
    
    # Basic material properties
    name: str = "Silicon"
    relative_permittivity: float = 11.7  # EPSIL
    bandgap: float = 1.12  # EGAP [eV]
    electron_affinity: float = 4.05  # CHI [eV]
    
    # Effective masses (in units of free electron mass m0)
    conduction_band_mass: float = 0.26  # ACB [m0]
    valence_band_mass_heavy: float = 0.49  # AVBH [m0] 
    valence_band_mass_light: float = 0.16  # AVBL [m0]
    valence_band_mass_split_off: float = 0.23  # AVBSO [m0]
    valence_band_mass_eff: float = 0.81  # AVB [m0] - effective mass for density
    
    # Band structure parameters
    spin_orbit_splitting: float = 0.044  # ESO [eV]
    valence_band_offset: float = 0.0  # DELVB [eV] - for heterostructures
    
    # Doping parameters
    donor_concentration: float = 1e15  # CD [cm^-3]
    acceptor_concentration: float = 0.0  # CA [cm^-3]
    donor_binding_energy: float = 0.045  # ED [eV]
    acceptor_binding_energy: float = 0.045  # EA [eV]
    
    # Temperature
    temperature: float = 300.0  # T [K]
    
    def __post_init__(self):
        """Calculate derived properties"""
        # Thermal energy in eV
        self.thermal_energy = 8.617333e-5 * self.temperature  # kT [eV]
        
        # Net doping
        self.net_doping = self.donor_concentration - self.acceptor_concentration
        
        # Intrinsic carrier concentration (approximate for Si at 300K)
        self.intrinsic_concentration = 1.45e10 * (self.temperature / 300.0) ** 1.5 * \
                                     np.exp(-self.bandgap / (2 * self.thermal_energy))


@dataclass 
class SurfaceStateParameters:
    """
    Surface state parameters corresponding to COMMON /SURF/ block.
    
    Based on the surface charge model in surfrhomult-6.2.f
    """
    
    # Surface state densities [cm^-2 eV^-1]
    donor_density: float = 1e12  # DENS(IAR,1) - donor-like states
    acceptor_density: float = 1e12  # DENS(IAR,2) - acceptor-like states
    
    # Characteristic energies [eV] (relative to valence band)
    donor_energy: float = 0.8  # EN(IAR,1) 
    acceptor_energy: float = 0.4  # EN(IAR,2)
    charge_neutrality_level: float = 0.6  # EN0(IAR)
    
    # Gaussian distribution parameters for continuous surface states
    donor_center_energy: float = 0.8  # ECENT(IAR,1) [eV]
    acceptor_center_energy: float = 0.4  # ECENT(IAR,2) [eV]
    donor_width: float = 0.1  # FWHM(IAR,1) [eV]
    acceptor_width: float = 0.1  # FWHM(IAR,2) [eV]
    
    # Temperature dependence flag
    temperature_dependent: bool = True  # ISTK


@dataclass
class MaterialParameters:
    """
    Complete material parameter set for a semiconductor region.
    
    Combines bulk semiconductor properties with surface state parameters.
    """
    
    semiconductor: SemiconductorMaterial
    surface_states: Optional[SurfaceStateParameters] = None
    region_type: RegionType = RegionType.BULK
    region_id: int = 1
    
    def get_fermi_dirac_parameters(self) -> Dict[str, float]:
        """Get parameters needed for Fermi-Dirac calculations"""
        return {
            'temperature': self.semiconductor.temperature,
            'thermal_energy': self.semiconductor.thermal_energy,
            'bandgap': self.semiconductor.bandgap,
            'conduction_band_mass': self.semiconductor.conduction_band_mass,
            'valence_band_mass_heavy': self.semiconductor.valence_band_mass_heavy,
            'valence_band_mass_light': self.semiconductor.valence_band_mass_light,
            'valence_band_mass_split_off': self.semiconductor.valence_band_mass_split_off
        }


class MaterialDatabase:
    """
    Database of predefined semiconductor materials.
    
    Provides easy access to common semiconductor parameters used in STM simulations.
    Based on experimental data and the parameter sets used in the original Fortran code.
    """
    
    def __init__(self):
        self._materials = {}
        self._initialize_database()
    
    def _initialize_database(self):
        """Initialize database with common semiconductor materials"""
        
        # Silicon (n-type, moderate doping)
        si_n = SemiconductorMaterial(
            name="Silicon_n_type",
            relative_permittivity=11.7,
            bandgap=1.12,
            electron_affinity=4.05,
            conduction_band_mass=0.26,
            valence_band_mass_heavy=0.49,
            valence_band_mass_light=0.16,
            valence_band_mass_split_off=0.23,
            valence_band_mass_eff=0.81,
            spin_orbit_splitting=0.044,
            donor_concentration=1e15,
            acceptor_concentration=0.0,
            donor_binding_energy=0.045,
            temperature=300.0
        )
        self._materials["Si_n"] = si_n
        
        # Silicon (p-type, moderate doping)
        si_p = SemiconductorMaterial(
            name="Silicon_p_type", 
            relative_permittivity=11.7,
            bandgap=1.12,
            electron_affinity=4.05,
            conduction_band_mass=0.26,
            valence_band_mass_heavy=0.49,
            valence_band_mass_light=0.16,
            valence_band_mass_split_off=0.23,
            valence_band_mass_eff=0.81,
            spin_orbit_splitting=0.044,
            donor_concentration=0.0,
            acceptor_concentration=1e15,
            acceptor_binding_energy=0.045,
            temperature=300.0
        )
        self._materials["Si_p"] = si_p
        
        # Gallium Arsenide (n-type)
        gaas_n = SemiconductorMaterial(
            name="GaAs_n_type",
            relative_permittivity=13.1,
            bandgap=1.42,
            electron_affinity=4.07,
            conduction_band_mass=0.067,
            valence_band_mass_heavy=0.45,
            valence_band_mass_light=0.082,
            valence_band_mass_split_off=0.15,
            valence_band_mass_eff=0.53,
            spin_orbit_splitting=0.34,
            donor_concentration=1e16,
            acceptor_concentration=0.0,
            donor_binding_energy=0.0058,
            temperature=300.0
        )
        self._materials["GaAs_n"] = gaas_n
        
        # Intrinsic Silicon for reference
        si_intrinsic = SemiconductorMaterial(
            name="Silicon_intrinsic",
            relative_permittivity=11.7,
            bandgap=1.12,
            electron_affinity=4.05,
            conduction_band_mass=0.26,
            valence_band_mass_heavy=0.49,
            valence_band_mass_light=0.16,
            valence_band_mass_split_off=0.23,
            valence_band_mass_eff=0.81,
            spin_orbit_splitting=0.044,
            donor_concentration=0.0,
            acceptor_concentration=0.0,
            temperature=300.0
        )
        self._materials["Si_intrinsic"] = si_intrinsic
    
    def get_material(self, name: str) -> SemiconductorMaterial:
        """Get material by name"""
        if name not in self._materials:
            available = list(self._materials.keys())
            raise ValueError(f"Material '{name}' not found. Available: {available}")
        return self._materials[name]
    
    def list_materials(self) -> List[str]:
        """List all available materials"""
        return list(self._materials.keys())
    
    def add_material(self, name: str, material: SemiconductorMaterial):
        """Add custom material to database"""
        self._materials[name] = material
    
    def create_material_from_dict(self, material_dict: Dict) -> SemiconductorMaterial:
        """Create material from dictionary (useful for JSON loading)"""
        return SemiconductorMaterial(**material_dict)
    
    def save_to_json(self, filename: str):
        """Save database to JSON file"""
        data = {}
        for name, material in self._materials.items():
            # Convert dataclass to dict
            material_dict = {
                'name': material.name,
                'relative_permittivity': material.relative_permittivity,
                'bandgap': material.bandgap,
                'electron_affinity': material.electron_affinity,
                'conduction_band_mass': material.conduction_band_mass,
                'valence_band_mass_heavy': material.valence_band_mass_heavy,
                'valence_band_mass_light': material.valence_band_mass_light,
                'valence_band_mass_split_off': material.valence_band_mass_split_off,
                'valence_band_mass_eff': material.valence_band_mass_eff,
                'spin_orbit_splitting': material.spin_orbit_splitting,
                'valence_band_offset': material.valence_band_offset,
                'donor_concentration': material.donor_concentration,
                'acceptor_concentration': material.acceptor_concentration,
                'donor_binding_energy': material.donor_binding_energy,
                'acceptor_binding_energy': material.acceptor_binding_energy,
                'temperature': material.temperature
            }
            data[name] = material_dict
        
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
    
    def load_from_json(self, filename: str):
        """Load database from JSON file"""
        with open(filename, 'r') as f:
            data = json.load(f)
        
        for name, material_dict in data.items():
            material = self.create_material_from_dict(material_dict)
            self._materials[name] = material


class PhysicalConstants:
    """
    Physical constants used throughout the simulation.
    
    These correspond to the constants defined in the Fortran code and are
    essential for maintaining numerical consistency.
    """
    
    # Fundamental constants
    ELEMENTARY_CHARGE = 1.602176634e-19  # C (exact, 2019 SI)
    VACUUM_PERMITTIVITY = 8.8541878128e-12  # F/m (exact, 2019 SI) 
    PLANCK_CONSTANT = 6.62607015e-34  # J⋅s (exact, 2019 SI)
    HBAR = PLANCK_CONSTANT / (2 * np.pi)  # ℏ
    ELECTRON_MASS = 9.1093837015e-31  # kg (exact, 2019 SI)
    BOLTZMANN_CONSTANT = 1.380649e-23  # J/K (exact, 2019 SI)
    
    # Derived constants in convenient units
    K_B_EV = 8.617333262145e-5  # eV/K (Boltzmann constant in eV)
    
    # MultInt specific constants (from Fortran code)
    EEP = 1.80943e-20  # V⋅cm = e/ε₀ × 10⁻¹⁴ cm²/nm²
    C_KINETIC = 26.254  # nm⁻²⋅eV⁻¹ = 2m/ℏ² kinetic energy constant  
    RESISTANCE_QUANTUM = 12900.0  # Ω = h/e² resistance quantum
    
    # Numerical constants from semitip3-6.1.f
    CHARGE_DENSITY_CONSTANT = 6.815e21  # eV^-1.5 cm^-3 for Fermi-Dirac
    
    @classmethod
    def get_constants_dict(cls) -> Dict[str, float]:
        """Get all constants as dictionary"""
        return {
            'elementary_charge': cls.ELEMENTARY_CHARGE,
            'vacuum_permittivity': cls.VACUUM_PERMITTIVITY,
            'planck_constant': cls.PLANCK_CONSTANT,
            'hbar': cls.HBAR,
            'electron_mass': cls.ELECTRON_MASS,
            'boltzmann_constant': cls.BOLTZMANN_CONSTANT,
            'k_b_ev': cls.K_B_EV,
            'eep': cls.EEP,
            'c_kinetic': cls.C_KINETIC,
            'resistance_quantum': cls.RESISTANCE_QUANTUM,
            'charge_density_constant': cls.CHARGE_DENSITY_CONSTANT
        }


def create_default_surface_states() -> SurfaceStateParameters:
    """Create default surface state parameters"""
    return SurfaceStateParameters(
        donor_density=1e12,
        acceptor_density=1e12,
        donor_energy=0.8,
        acceptor_energy=0.4,
        charge_neutrality_level=0.6,
        donor_center_energy=0.8,
        acceptor_center_energy=0.4,
        donor_width=0.1,
        acceptor_width=0.1,
        temperature_dependent=True
    )


def create_multimaterial_system(materials: List[Tuple[str, SemiconductorMaterial]]) -> Dict[int, MaterialParameters]:
    """
    Create a multi-material system for heterostructures.
    
    Args:
        materials: List of (region_name, material) tuples
        
    Returns:
        Dictionary mapping region_id to MaterialParameters
    """
    system = {}
    for i, (name, material) in enumerate(materials):
        region_id = i + 1  # 1-based indexing like Fortran
        system[region_id] = MaterialParameters(
            semiconductor=material,
            region_type=RegionType.BULK,
            region_id=region_id
        )
    return system


# Pre-instantiated database for easy access
default_materials = MaterialDatabase()


if __name__ == "__main__":
    # Demo usage
    db = MaterialDatabase()
    
    print("Available materials:")
    for mat_name in db.list_materials():
        material = db.get_material(mat_name)
        print(f"  {mat_name}: {material.name}, Eg = {material.bandgap:.2f} eV")
    
    # Show Si n-type parameters
    si_n = db.get_material("Si_n")
    print(f"\nSilicon n-type parameters:")
    print(f"  Bandgap: {si_n.bandgap:.2f} eV")
    print(f"  Thermal energy (300K): {si_n.thermal_energy:.4f} eV")
    print(f"  Net doping: {si_n.net_doping:.2e} cm⁻³")
    print(f"  Intrinsic concentration: {si_n.intrinsic_concentration:.2e} cm⁻³")
    
    # Show constants
    print(f"\nPhysical constants:")
    constants = PhysicalConstants.get_constants_dict()
    for name, value in constants.items():
        print(f"  {name}: {value:.6e}")
