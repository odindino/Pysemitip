"""
Charge density calculations for SEMITIP simulations.

This module implements the charge density calculations from SEMIRHOMULT
and SURFRHOMULT, including tabulation of charge densities for efficient
evaluation during Poisson equation solving.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy import special
from scipy import integrate

from ...utils.constants import PhysicalConstants as PC
from ..materials.semiconductor import SemiconductorRegion
from ..materials.surface_states import SurfaceRegion


@dataclass
class ChargeDensityTables:
    """
    Tabulated charge densities for efficient lookup.
    
    Corresponds to the /CD/ common block in Fortran.
    """
    energy_start: float    # ESTART - starting energy (eV)
    energy_step: float     # DELE - energy step size (eV)
    num_points: int        # NE - number of energy points
    
    # Bulk charge density tables (one per semiconductor region)
    bulk_tables: Dict[int, np.ndarray]  # RHOBTAB[region_id] = density array
    
    # Surface charge density tables (one per surface area)
    surface_tables: Dict[int, np.ndarray]  # RHOSTAB[area_id] = density array
    
    def get_energy_array(self) -> np.ndarray:
        """Get the energy array for the tables."""
        return np.linspace(self.energy_start, 
                          self.energy_start + (self.num_points - 1) * self.energy_step,
                          self.num_points)
    
    def interpolate_bulk_density(self, region_id: int, energy: float) -> float:
        """
        Interpolate bulk charge density at given energy.
        
        Args:
            region_id: Semiconductor region ID
            energy: Energy value (eV)
            
        Returns:
            Charge density (C/cm^3)
        """
        if region_id not in self.bulk_tables:
            raise ValueError(f"No bulk table for region {region_id}")
        
        # Check for invalid energy values
        if not np.isfinite(energy):
            return 0.0
        
        # Find interpolation indices
        idx = (energy - self.energy_start) / self.energy_step
        
        # Check for invalid index
        if not np.isfinite(idx):
            return 0.0
            
        idx0 = int(np.floor(idx))
        idx1 = idx0 + 1
        
        # Boundary handling
        if idx0 < 0:
            return self.bulk_tables[region_id][0]
        if idx1 >= self.num_points:
            return self.bulk_tables[region_id][-1]
        
        # Linear interpolation
        w = idx - idx0
        return ((1 - w) * self.bulk_tables[region_id][idx0] + 
                w * self.bulk_tables[region_id][idx1])
    
    def interpolate_surface_density(self, area_id: int, energy: float) -> float:
        """
        Interpolate surface charge density at given energy.
        
        Args:
            area_id: Surface area ID
            energy: Energy value (eV)
            
        Returns:
            Charge density (C/cm^2)
        """
        if area_id not in self.surface_tables:
            raise ValueError(f"No surface table for area {area_id}")
        
        # Check for invalid energy values
        if not np.isfinite(energy):
            return 0.0
        
        # Similar interpolation logic
        idx = (energy - self.energy_start) / self.energy_step
        
        # Check for invalid index
        if not np.isfinite(idx):
            return 0.0
            
        idx0 = int(np.floor(idx))
        idx1 = idx0 + 1
        
        if idx0 < 0:
            return self.surface_tables[area_id][0]
        if idx1 >= self.num_points:
            return self.surface_tables[area_id][-1]
        
        w = idx - idx0
        return ((1 - w) * self.surface_tables[area_id][idx0] + 
                w * self.surface_tables[area_id][idx1])


class ChargeDensityCalculator:
    """
    Calculator for semiconductor and surface charge densities.
    
    Implements the functionality of SEMIRHO and SURFRHO subroutines.
    """
    
    def __init__(self, semiconductor_regions: List[SemiconductorRegion],
                 surface_regions: List[SurfaceRegion],
                 fermi_level: float):
        """
        Initialize the charge density calculator.
        
        Args:
            semiconductor_regions: List of semiconductor regions
            surface_regions: List of surface regions
            fermi_level: Bulk Fermi level (eV)
        """
        self.semiconductor_regions = {r.region_id: r for r in semiconductor_regions}
        self.surface_regions = {r.region_id: r for r in surface_regions}
        self.fermi_level = fermi_level
        
        # Initialize Fermi integral calculator
        self._init_fermi_integrals()
    
    def _init_fermi_integrals(self):
        """Initialize Fermi integral calculations."""
        # For now, we'll use approximations
        # Full implementation would use scipy.special functions
        pass
    
    def fermi_integral_half(self, eta: float) -> float:
        """
        Calculate Fermi-Dirac integral of order 1/2.
        
        F_{1/2}(η) = ∫[0,∞] x^{1/2} / (1 + exp(x - η)) dx
        
        Args:
            eta: Reduced energy (E - Ec) / kT
            
        Returns:
            Fermi integral value
        """
        if eta > 5:
            # Degenerate limit (Sommerfeld expansion)
            return (2 / 3) * np.sqrt(np.pi) * eta**(3/2) * (
                1 + np.pi**2 / 8 * eta**(-2) + 7 * np.pi**4 / 640 * eta**(-4)
            )
        elif eta < -5:
            # Non-degenerate limit (Boltzmann)
            return np.sqrt(np.pi) / 2 * np.exp(eta)
        else:
            # Intermediate region - use numerical integration
            from scipy import integrate
            
            def integrand(x):
                if x > 100:  # Avoid overflow
                    return 0.0
                return np.sqrt(x) / (1 + np.exp(x - eta))
            
            result, _ = integrate.quad(integrand, 0, max(50, eta + 20))
            return result
    
    def calculate_bulk_density(self, region_id: int, energy: float, 
                             potential: float = 0.0) -> float:
        """
        Calculate bulk charge density at given energy and potential.
        
        Implements SEMIRHO functionality.
        
        Args:
            region_id: Semiconductor region ID
            energy: Energy relative to vacuum (eV)
            potential: Electrostatic potential (V)
            
        Returns:
            Total charge density (C/cm^3)
        """
        if region_id not in self.semiconductor_regions:
            raise ValueError(f"Unknown semiconductor region: {region_id}")
        
        region = self.semiconductor_regions[region_id]
        
        # Adjust energy for potential
        energy_local = energy + potential
        
        # Calculate carrier densities
        n = self._electron_density(region, energy_local, self.fermi_level)
        p = self._hole_density(region, energy_local, self.fermi_level)
        
        # Ionized dopant densities (complete ionization assumed)
        n_d = region.donor_concentration
        n_a = region.acceptor_concentration
        
        # Net charge density (C/cm^3)
        rho = PC.E * (p - n + n_d - n_a)
        
        return rho
    
    def _electron_density(self, region: SemiconductorRegion, 
                         energy: float, fermi_level: float) -> float:
        """
        Calculate electron density in conduction band.
        
        Args:
            region: Semiconductor region
            energy: Local energy (including potential)
            fermi_level: Fermi level
            
        Returns:
            Electron density (cm^-3)
        """
        # Conduction band edge
        Ec = region.conduction_band_edge()
        
        # Reduced energy (note: energy here includes potential shift)
        # eta = (EF - EC) / kT
        eta = (fermi_level - (Ec - energy)) / region.kT
        
        # Effective density of states
        Nc = region.Nc
        
        if region.allow_degeneracy:
            # Fermi-Dirac statistics
            n = Nc * (2 / np.sqrt(np.pi)) * self.fermi_integral_half(eta)
        else:
            # Boltzmann statistics
            n = Nc * np.exp(eta) if eta < 0 else Nc  # Cap at Nc for positive eta
        
        return n
    
    def _hole_density(self, region: SemiconductorRegion,
                     energy: float, fermi_level: float) -> float:
        """
        Calculate hole density in valence band.
        
        Args:
            region: Semiconductor region
            energy: Local energy (including potential)
            fermi_level: Fermi level
            
        Returns:
            Hole density (cm^-3)
        """
        # Valence band edge
        Ev = region.valence_band_edge()
        
        # Reduced energy for holes
        # eta = (EV - EF) / kT
        eta = ((Ev - energy) - fermi_level) / region.kT
        
        # Effective density of states
        Nv = region.Nv
        
        if region.allow_degeneracy:
            # Fermi-Dirac statistics
            p = Nv * (2 / np.sqrt(np.pi)) * self.fermi_integral_half(eta)
        else:
            # Boltzmann statistics
            p = Nv * np.exp(eta) if eta < 0 else Nv  # Cap at Nv
        
        return p
    
    def calculate_surface_density(self, area_id: int, energy: float,
                                potential: float = 0.0) -> float:
        """
        Calculate surface charge density at given energy and potential.
        
        Implements SURFRHO functionality.
        
        Args:
            area_id: Surface area ID
            energy: Energy relative to vacuum (eV)
            potential: Surface potential (V)
            
        Returns:
            Surface charge density (C/cm^2)
        """
        if area_id not in self.surface_regions:
            raise ValueError(f"Unknown surface region: {area_id}")
        
        region = self.surface_regions[area_id]
        
        # Surface Fermi level
        ef_surface = self.fermi_level + potential
        
        # Get charge from surface states
        return region.surface_charge_density(self.fermi_level, potential)
    
    def create_charge_tables(self, energy_range: Tuple[float, float],
                           num_points: int = 20000) -> ChargeDensityTables:
        """
        Create tabulated charge densities for all regions.
        
        Args:
            energy_range: (min, max) energy range (eV)
            num_points: Number of energy points
            
        Returns:
            ChargeDensityTables object
        """
        energy_start, energy_end = energy_range
        energy_step = (energy_end - energy_start) / (num_points - 1)
        
        # Energy array
        energies = np.linspace(energy_start, energy_end, num_points)
        
        # Create bulk tables
        bulk_tables = {}
        for region_id, region in self.semiconductor_regions.items():
            print(f"Computing bulk charge density table for region {region_id}")
            densities = np.zeros(num_points)
            
            for i, energy in enumerate(energies):
                # Calculate at zero potential for table
                densities[i] = self.calculate_bulk_density(region_id, energy, 0.0)
            
            bulk_tables[region_id] = densities
        
        # Create surface tables
        surface_tables = {}
        for area_id, region in self.surface_regions.items():
            print(f"Computing surface charge density table for area {area_id}")
            densities = np.zeros(num_points)
            
            for i, energy in enumerate(energies):
                # Energy here represents surface Fermi level
                # So potential = energy - bulk_fermi_level
                potential = energy - self.fermi_level
                densities[i] = self.calculate_surface_density(area_id, 
                                                            self.fermi_level, 
                                                            potential)
            
            surface_tables[area_id] = densities
        
        return ChargeDensityTables(
            energy_start=energy_start,
            energy_step=energy_step,
            num_points=num_points,
            bulk_tables=bulk_tables,
            surface_tables=surface_tables
        )


def estimate_energy_range(semiconductor_regions: List[SemiconductorRegion],
                         surface_regions: List[SurfaceRegion],
                         fermi_level: float,
                         bias_range: Tuple[float, float]) -> Tuple[float, float]:
    """
    Estimate appropriate energy range for charge density tables.
    
    Args:
        semiconductor_regions: List of semiconductor regions
        surface_regions: List of surface regions
        fermi_level: Bulk Fermi level
        bias_range: (min, max) bias voltage range
        
    Returns:
        (min, max) energy range for tables
    """
    # Get temperature
    kT = semiconductor_regions[0].kT if semiconductor_regions else 0.026
    
    # Get band edges
    ec_min = min(r.conduction_band_edge() for r in semiconductor_regions)
    ev_max = max(r.valence_band_edge() for r in semiconductor_regions)
    
    # Add bias range and thermal broadening
    e_min = ev_max + min(bias_range) - 10 * kT
    e_max = ec_min + max(bias_range) + 10 * kT
    
    # Include surface state ranges - handle both config and physics objects
    for region in surface_regions:
        try:
            # Try physics object attributes first (distribution1, distribution2)
            if hasattr(region, 'distribution1') and region.distribution1:
                d1 = region.distribution1
                e_min = min(e_min, d1.center_energy - 3 * max(d1.fwhm, 0.1))
                e_max = max(e_max, d1.center_energy + 3 * max(d1.fwhm, 0.1))
            
            if hasattr(region, 'distribution2') and region.distribution2:
                d2 = region.distribution2
                e_min = min(e_min, d2.center_energy - 3 * max(d2.fwhm, 0.1))
                e_max = max(e_max, d2.center_energy + 3 * max(d2.fwhm, 0.1))
                
            # Try config object attributes as fallback (first_distribution, second_distribution)
            elif hasattr(region, 'first_distribution') and region.first_distribution:
                d1 = region.first_distribution
                e_min = min(e_min, d1.center_energy - 3 * max(d1.fwhm, 0.1))
                e_max = max(e_max, d1.center_energy + 3 * max(d1.fwhm, 0.1))
                
            if hasattr(region, 'second_distribution') and region.second_distribution:
                d2 = region.second_distribution
                e_min = min(e_min, d2.center_energy - 3 * max(d2.fwhm, 0.1))
                e_max = max(e_max, d2.center_energy + 3 * max(d2.fwhm, 0.1))
                
        except Exception as e:
            print(f"Warning: Error processing surface region {getattr(region, 'region_id', '?')}: {e}")
            # Continue with other regions
    
    return e_min, e_max