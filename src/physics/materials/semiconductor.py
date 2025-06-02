"""
Semiconductor material properties and band structure calculations.

This module implements the semiconductor physics from SEMIRHOMULT,
including band structure, carrier densities, and Fermi level calculations.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple
from scipy import optimize
from scipy import special

from ...utils.constants import PhysicalConstants as PC


@dataclass
class SemiconductorRegion:
    """
    Represents a semiconductor region with specific material properties.
    
    Corresponds to the /SEMI/ common block in Fortran.
    """
    # Material identification
    region_id: int
    
    # Doping concentrations (cm^-3)
    donor_concentration: float      # CD in Fortran
    acceptor_concentration: float   # CA in Fortran
    
    # Band structure parameters (eV)
    band_gap: float                # EGAP
    valence_band_offset: float     # DELVB
    electron_affinity: float       # CHI
    
    # Binding energies (eV)
    donor_binding_energy: float    # ED
    acceptor_binding_energy: float # EA
    
    # Effective masses (in units of m0)
    cb_effective_mass: float       # ACB
    vb_effective_mass_heavy: float # AVBH
    vb_effective_mass_light: float # AVBL
    vb_effective_mass_so: float    # AVBSO (split-off)
    
    # Spin-orbit splitting (eV)
    spin_orbit_splitting: float    # ESO
    
    # Material properties
    permittivity: float           # Relative permittivity
    
    # Calculation parameters
    allow_degeneracy: bool = True  # IDEG
    allow_inversion: bool = True   # IINV
    
    # Temperature (K)
    temperature: float = 300.0
    
    def __post_init__(self):
        """Calculate derived quantities."""
        # Average valence band effective mass (Fortran formula)
        self.vb_effective_mass_avg = np.exp(
            2.0 * np.log(np.sqrt(self.vb_effective_mass_heavy**3) + 
                        np.sqrt(self.vb_effective_mass_light**3)) / 3.0
        )
        
        # Thermal energy
        self.kT = PC.thermal_energy(self.temperature)
        
        # Initialize carrier concentration variables
        self.Nc = 0.0
        self.Nv = 0.0
        self.ni = 0.0
        
        # Calculate intrinsic carrier concentration
        self._calculate_intrinsic_concentration()
    
    def _calculate_intrinsic_concentration(self):
        """Calculate intrinsic carrier concentration."""
        # Effective density of states in CB and VB
        self.Nc = 2.0 * (2.0 * np.pi * self.cb_effective_mass * PC.M0 * 
                         PC.KB * self.temperature / PC.H**2)**(3/2) * PC.CM3_TO_M3
        
        self.Nv = 2.0 * (2.0 * np.pi * self.vb_effective_mass_avg * PC.M0 * 
                         PC.KB * self.temperature / PC.H**2)**(3/2) * PC.CM3_TO_M3
        
        # Intrinsic concentration
        self.ni = np.sqrt(self.Nc * self.Nv) * np.exp(-self.band_gap / (2 * self.kT))
    
    @property
    def is_n_type(self) -> bool:
        """Check if the material is n-type."""
        return self.donor_concentration > self.acceptor_concentration
    
    @property
    def is_p_type(self) -> bool:
        """Check if the material is p-type."""
        return self.acceptor_concentration > self.donor_concentration
    
    @property
    def net_doping(self) -> float:
        """Net doping concentration (positive for n-type, negative for p-type)."""
        return self.donor_concentration - self.acceptor_concentration
    
    def conduction_band_edge(self, vacuum_level: float = 0.0) -> float:
        """
        Calculate conduction band edge position.
        
        Args:
            vacuum_level: Reference vacuum level (eV)
            
        Returns:
            Conduction band edge energy (eV)
        """
        return vacuum_level - self.electron_affinity
    
    def valence_band_edge(self, vacuum_level: float = 0.0) -> float:
        """
        Calculate valence band edge position.
        
        Args:
            vacuum_level: Reference vacuum level (eV)
            
        Returns:
            Valence band edge energy (eV)
        """
        return self.conduction_band_edge(vacuum_level) - self.band_gap + self.valence_band_offset
    
    def fermi_level(self, potential: float = 0.0) -> float:
        """
        Calculate Fermi level position (implementation of EFFIND).
        
        Args:
            potential: Electrostatic potential at the point (V)
            
        Returns:
            Fermi level position (eV) relative to vacuum
        """
        # This is a simplified version - full implementation needs iterative solution
        # considering degeneracy and temperature effects
        
        if self.allow_degeneracy:
            # Full Fermi-Dirac statistics
            return self._calculate_fermi_level_fd(potential)
        else:
            # Boltzmann approximation
            return self._calculate_fermi_level_boltzmann(potential)
    
    def _calculate_fermi_level_boltzmann(self, potential: float) -> float:
        """Calculate Fermi level using Boltzmann approximation."""
        Ec = self.conduction_band_edge() - potential
        Ev = self.valence_band_edge() - potential
        
        # Simplified calculation for non-degenerate case
        if self.is_n_type:
            # n-type: EF near conduction band
            n = self.net_doping
            EF = Ec - self.kT * np.log(self.Nc / n)
        elif self.is_p_type:
            # p-type: EF near valence band
            p = -self.net_doping
            EF = Ev + self.kT * np.log(self.Nv / p)
        else:
            # Intrinsic
            EF = (Ec + Ev) / 2 + self.kT / 2 * np.log(self.Nv / self.Nc)
        
        return EF
    
    def _calculate_fermi_level_fd(self, potential: float) -> float:
        """Calculate Fermi level using full Fermi-Dirac statistics."""
        from scipy.optimize import brentq
        
        Ec = self.conduction_band_edge() - potential
        Ev = self.valence_band_edge() - potential
        
        # Define charge neutrality equation
        def charge_neutrality(ef):
            # Electron density
            n = self.carrier_density_cb(ef, potential)
            # Hole density
            p = self.carrier_density_vb(ef, potential)
            # Net charge (should be zero at equilibrium)
            return n - p + self.acceptor_concentration - self.donor_concentration
        
        # Initial guess based on Boltzmann approximation
        ef_guess = self._calculate_fermi_level_boltzmann(potential)
        
        # Search range
        ef_min = Ev - 10 * self.kT
        ef_max = Ec + 10 * self.kT
        
        try:
            # Find root using Brent's method
            ef = brentq(charge_neutrality, ef_min, ef_max, xtol=1e-9)
        except ValueError:
            # If root finding fails, fall back to Boltzmann approximation
            ef = ef_guess
        
        return ef
    
    def carrier_density_cb(self, fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate electron density in conduction band (RHOCB).
        
        Args:
            fermi_level: Fermi level position (eV)
            potential: Electrostatic potential (V)
            
        Returns:
            Electron density (cm^-3)
        """
        Ec = self.conduction_band_edge() - potential
        eta = (fermi_level - Ec) / self.kT
        
        if self.allow_degeneracy:
            # Fermi-Dirac integral F_{1/2}(eta)
            # For now, use approximation for testing
            if eta > 5:
                # Degenerate limit
                n = self.Nc * (2/np.sqrt(np.pi)) * eta**(3/2)
            elif eta < -5:
                # Non-degenerate limit
                n = self.Nc * np.exp(eta)
            else:
                # Intermediate - needs proper Fermi integral
                n = self.Nc * np.exp(eta)  # Placeholder
        else:
            # Boltzmann statistics
            n = self.Nc * np.exp(eta)
        
        return n
    
    def carrier_density_vb(self, fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate hole density in valence band (RHOVB).
        
        Args:
            fermi_level: Fermi level position (eV)
            potential: Electrostatic potential (V)
            
        Returns:
            Hole density (cm^-3)
        """
        Ev = self.valence_band_edge() - potential
        eta = (Ev - fermi_level) / self.kT
        
        if self.allow_degeneracy:
            # Similar to CB calculation but for holes
            if eta > 5:
                p = self.Nv * (2/np.sqrt(np.pi)) * eta**(3/2)
            elif eta < -5:
                p = self.Nv * np.exp(eta)
            else:
                p = self.Nv * np.exp(eta)  # Placeholder
        else:
            p = self.Nv * np.exp(eta)
        
        return p
    
    def total_charge_density(self, fermi_level: float, potential: float = 0.0) -> float:
        """
        Calculate total charge density.
        
        Args:
            fermi_level: Fermi level position (eV)
            potential: Electrostatic potential (V)
            
        Returns:
            Total charge density (C/cm^3)
        """
        n = self.carrier_density_cb(fermi_level, potential)
        p = self.carrier_density_vb(fermi_level, potential)
        
        # Include ionized dopants (complete ionization assumed for now)
        Nd_plus = self.donor_concentration
        Na_minus = self.acceptor_concentration
        
        # Net charge (in units of e/cm^3)
        net_charge = p - n + Nd_plus - Na_minus
        
        # Convert to C/cm^3
        return PC.E * net_charge


def create_semiconductor_from_config(config_region) -> SemiconductorRegion:
    """
    Create a SemiconductorRegion from configuration data.
    
    Args:
        config_region: Configuration object for a semiconductor region
        
    Returns:
        SemiconductorRegion object
    """
    return SemiconductorRegion(
        region_id=config_region.id,
        donor_concentration=config_region.donor_concentration,
        acceptor_concentration=config_region.acceptor_concentration,
        band_gap=config_region.band_gap,
        valence_band_offset=config_region.valence_band_offset,
        electron_affinity=config_region.affinity,
        donor_binding_energy=config_region.donor_binding_energy,
        acceptor_binding_energy=config_region.acceptor_binding_energy,
        cb_effective_mass=config_region.effective_mass.conduction_band,
        vb_effective_mass_heavy=config_region.effective_mass.valence_band_heavy,
        vb_effective_mass_light=config_region.effective_mass.valence_band_light,
        vb_effective_mass_so=config_region.effective_mass.split_off,
        spin_orbit_splitting=config_region.spin_orbit_splitting,
        permittivity=config_region.permittivity,
        allow_degeneracy=config_region.allow_degeneracy,
        allow_inversion=config_region.allow_inversion,
        temperature=300.0  # Will come from environment config
    )