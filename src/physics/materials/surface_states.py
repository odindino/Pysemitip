"""
Surface states and surface charge calculations.

This module implements surface state distributions and charge calculations
from SURFRHOMULT in the Fortran code.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, List
from scipy import integrate

from ...utils.constants import PhysicalConstants as PC


@dataclass
class SurfaceStateDistribution:
    """
    Represents a distribution of surface states.
    
    Can be either a constant DOS or Gaussian distribution.
    """
    density: float          # DENS - density of states (cm^-2 eV^-1)
    neutrality_level: float # EN - charge neutrality level (eV)
    fwhm: float = 0.0      # FWHM - full width half maximum for Gaussian (eV)
    center_energy: float = 0.0  # ECENT - center energy for Gaussian (eV)
    
    @property
    def is_gaussian(self) -> bool:
        """Check if this is a Gaussian distribution."""
        return self.fwhm > 0
    
    def density_of_states(self, energy: np.ndarray) -> np.ndarray:
        """
        Calculate density of states at given energies.
        
        Args:
            energy: Energy values (eV)
            
        Returns:
            DOS at each energy (cm^-2 eV^-1)
        """
        if not self.is_gaussian:
            # Constant DOS
            return np.full_like(energy, self.density)
        else:
            # Gaussian distribution
            sigma = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
            return (self.density / (sigma * np.sqrt(2 * np.pi)) * 
                   np.exp(-(energy - self.center_energy)**2 / (2 * sigma**2)))
    
    def integrated_dos(self, energy_min: float, energy_max: float) -> float:
        """
        Integrate DOS over energy range.
        
        Args:
            energy_min: Lower energy bound (eV)
            energy_max: Upper energy bound (eV)
            
        Returns:
            Integrated DOS (cm^-2)
        """
        if not self.is_gaussian:
            # Constant DOS - simple integration
            return self.density * (energy_max - energy_min)
        else:
            # Gaussian - use error function
            sigma = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
            from scipy.special import erf
            
            z1 = (energy_min - self.center_energy) / (np.sqrt(2) * sigma)
            z2 = (energy_max - self.center_energy) / (np.sqrt(2) * sigma)
            
            return self.density * (erf(z2) - erf(z1)) / 2


@dataclass
class SurfaceRegion:
    """
    Represents a surface region with up to two distributions of surface states.
    
    Corresponds to surface state parameters in the /SURF/ common block.
    """
    region_id: int
    position: Tuple[float, float]  # (x, y) position
    
    # Two possible distributions
    distribution1: Optional[SurfaceStateDistribution] = None
    distribution2: Optional[SurfaceStateDistribution] = None
    
    # Temperature
    temperature: float = 300.0  # K
    
    def __post_init__(self):
        """Calculate combined properties."""
        self.kT = PC.thermal_energy(self.temperature)
        self._calculate_combined_neutrality_level()
    
    def _calculate_combined_neutrality_level(self):
        """
        Calculate combined charge neutrality level (EN0).
        
        This implements the ENFIND subroutine functionality.
        """
        if self.distribution1 is None and self.distribution2 is None:
            self.combined_neutrality_level = 0.0
            return
        
        if self.distribution1 is None:
            self.combined_neutrality_level = self.distribution2.neutrality_level
            return
            
        if self.distribution2 is None:
            self.combined_neutrality_level = self.distribution1.neutrality_level
            return
        
        # Both distributions present - need to find combined neutrality level
        # This is where charge from both distributions balance
        d1 = self.distribution1
        d2 = self.distribution2
        
        # For constant DOS, analytical solution exists
        if not d1.is_gaussian and not d2.is_gaussian:
            # Weighted average based on densities
            total_density = d1.density + d2.density
            if total_density > 0:
                self.combined_neutrality_level = (
                    (d1.density * d1.neutrality_level + 
                     d2.density * d2.neutrality_level) / total_density
                )
            else:
                self.combined_neutrality_level = 0.0
        else:
            # For Gaussian or mixed distributions, numerical solution needed
            # This would implement the iterative search in ENFIND
            # For now, use simple average as approximation
            self.combined_neutrality_level = (
                d1.neutrality_level + d2.neutrality_level) / 2
    
    def surface_charge_density(self, fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate surface charge density.
        
        Args:
            fermi_level: Fermi level position (eV)
            potential: Surface potential (V)
            
        Returns:
            Surface charge density (C/cm^2)
        """
        # Adjust Fermi level for surface potential
        ef_surface = fermi_level + potential
        
        charge = 0.0
        
        # Contribution from distribution 1
        if self.distribution1 is not None:
            charge += self._distribution_charge(self.distribution1, ef_surface)
        
        # Contribution from distribution 2
        if self.distribution2 is not None:
            charge += self._distribution_charge(self.distribution2, ef_surface)
        
        # Convert to C/cm^2
        return -PC.E * charge  # Negative because electrons have negative charge
    
    def _distribution_charge(self, dist: SurfaceStateDistribution, 
                           ef_surface: float) -> float:
        """
        Calculate charge from a single distribution.
        
        Args:
            dist: Surface state distribution
            ef_surface: Surface Fermi level (eV)
            
        Returns:
            Charge density (electrons/cm^2)
        """
        if dist.density == 0:
            return 0.0
        
        # Define integration limits
        # Integrate from well below to well above Fermi level
        e_min = ef_surface - 10 * self.kT
        e_max = ef_surface + 10 * self.kT
        
        def integrand(energy):
            """Occupied DOS = DOS * Fermi function."""
            dos = dist.density_of_states(np.array([energy]))[0]
            fermi = 1.0 / (1.0 + np.exp((energy - ef_surface) / self.kT))
            return dos * fermi
        
        # Numerical integration
        charge, _ = integrate.quad(integrand, e_min, e_max)
        
        # Subtract charge at neutrality level
        # (This makes the charge zero when EF = EN)
        def integrand_neutral(energy):
            dos = dist.density_of_states(np.array([energy]))[0]
            fermi = 1.0 / (1.0 + np.exp((energy - dist.neutrality_level) / self.kT))
            return dos * fermi
        
        charge_neutral, _ = integrate.quad(integrand_neutral, e_min, e_max)
        
        return charge - charge_neutral
    
    def tabulate_charge_density(self, energy_array: np.ndarray, 
                              fermi_level: float) -> np.ndarray:
        """
        Tabulate surface charge density vs energy.
        
        This creates the RHOSTAB table used in the Fortran code.
        
        Args:
            energy_array: Array of energy values (eV)
            fermi_level: Bulk Fermi level (eV)
            
        Returns:
            Array of charge densities (C/cm^2)
        """
        charges = np.zeros_like(energy_array)
        
        for i, energy in enumerate(energy_array):
            # Surface potential that would give this energy as surface Fermi level
            potential = energy - fermi_level
            charges[i] = self.surface_charge_density(fermi_level, potential)
        
        return charges


def create_surface_region_from_config(config_region, region_id: int) -> SurfaceRegion:
    """
    Create a SurfaceRegion from configuration data.
    
    Args:
        config_region: Configuration object for a surface region
        region_id: Region identifier
        
    Returns:
        SurfaceRegion object
    """
    # Create distribution 1
    dist1 = None
    if config_region.first_distribution.density > 0:
        dist1 = SurfaceStateDistribution(
            density=config_region.first_distribution.density,
            neutrality_level=config_region.first_distribution.neutrality_level,
            fwhm=config_region.first_distribution.fwhm,
            center_energy=config_region.first_distribution.center_energy
        )
    
    # Create distribution 2
    dist2 = None
    if (config_region.second_distribution and 
        config_region.second_distribution.density > 0):
        dist2 = SurfaceStateDistribution(
            density=config_region.second_distribution.density,
            neutrality_level=config_region.second_distribution.neutrality_level,
            fwhm=config_region.second_distribution.fwhm,
            center_energy=config_region.second_distribution.center_energy
        )
    
    return SurfaceRegion(
        region_id=region_id,
        position=(config_region.position.x, config_region.position.y),
        distribution1=dist1,
        distribution2=dist2,
        temperature=300.0  # Will come from environment config
    )