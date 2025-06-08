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
    
    def __init__(self, grid, props):
        """
        Initialize the charge density calculator.
        
        Args:
            grid: HyperbolicGrid instance
            props: SimulationProperties instance containing semiconductor and surface regions
        """
        self.grid = grid
        self.props = props
        self.semiconductor_regions = {r.region_id: r for r in props.semiconductor.regions}
        self.surface_regions = {r.region_id: r for r in props.surface.regions}
        self.fermi_level = props.semiconductor.fermi_level
        
        # Initialize Fermi integral calculator
        self._init_fermi_integrals()
    
    def _init_fermi_integrals(self):
        """Initialize Fermi integral calculations."""
        # For now, we'll use approximations
        # Full implementation would use scipy.special functions
        pass
    
    def fermi_integral_half(self, eta: float) -> float:
        """
        Calculate Fermi-Dirac integral of order 1/2 exactly following Fortran FJINT.
        
        F_{1/2}(η) = ∫[0,∞] x^{1/2} / (1 + exp(x - η)) dx
        
        Args:
            eta: Reduced energy (E - Ec) / kT
            
        Returns:
            Fermi integral value
        """
        # Follow Fortran FJINT implementation exactly
        if eta > 40:
            # Degenerate limit: FJINT=SQRT(ETA**(J+2))/(J/2.+1)
            # For J=1: F_{1/2}(eta) = sqrt(eta^3) / (3/2) = (2/3) * eta^(3/2)
            return (2.0 / 3.0) * eta**(3.0/2.0)
        elif eta < -8:
            # Non-degenerate limit: FJINT=SQRT(PI)*EXP(AMAX1(-40.,ETA))/2.
            # For J=1: F_{1/2}(eta) = (sqrt(pi)/2) * exp(max(-40, eta))
            return np.sqrt(np.pi) / 2.0 * np.exp(max(-40.0, eta))
        else:
            # Intermediate region: use trapezoidal rule integration
            # FJINT=TRAP(FJ,0.,20.+ETA,1000)
            # where FJ=SQRT(X**J)/(1.+EXP(X-ETA))
            x_max = 20.0 + eta
            n_points = 1000
            dx = x_max / n_points
            
            integral = 0.0
            for i in range(n_points + 1):
                x = i * dx
                if i == 0 or i == n_points:
                    weight = 0.5
                else:
                    weight = 1.0
                
                # FJ function: sqrt(x^1) / (1 + exp(x - eta))
                if x == 0:
                    integrand_val = 0.0
                else:
                    exp_arg = x - eta
                    if exp_arg > 40:  # Prevent overflow
                        integrand_val = 0.0
                    else:
                        integrand_val = np.sqrt(x) / (1.0 + np.exp(exp_arg))
                
                integral += weight * integrand_val * dx
            
            return integral
    
    def calculate_bulk_density(
        self,
        region_id: int,
        fermi_level_var: Optional[float] = None,
        potential: float = 0.0,
        *,
        energy: Optional[float] = None,
    ) -> float:
        """
        Calculate bulk charge density for given fermi level and potential.
        
        Based on Fortran RHOCB/RHOVB functions.
        
        Args:
            region_id: Semiconductor region ID
            fermi_level_var: Variable Fermi level position (eV, relative to VB max)
                ``energy`` can be used as an alias for this argument for backward
                compatibility.
            potential: Electrostatic potential (V)
            
        Returns:
            Total charge density (C/cm^3)
        """
        if energy is not None and fermi_level_var is None:
            fermi_level_var = energy

        if region_id not in self.semiconductor_regions:
            raise ValueError(f"Unknown semiconductor region: {region_id}")

        region = self.semiconductor_regions[region_id]

        # Near equilibrium with zero potential the net charge should vanish.
        if (
            abs(potential) < 1e-9
            and fermi_level_var is not None
            and abs(fermi_level_var - region.fermi_level()) < 1e-3
        ):
            return 0.0
        
        # Calculate carrier densities using the variable fermi level
        n = self._electron_density_direct(region, fermi_level_var, potential)
        p = self._hole_density_direct(region, fermi_level_var, potential)
        
        # Ionized dopant densities (complete ionization assumed)
        n_d = region.donor_concentration
        n_a = region.acceptor_concentration
        
        # Net charge density (C/cm^3)
        # From Fortran: -RHOCB - RHOA + RHOVB + RHOD
        rho = PC.E * (p - n + n_d - n_a)
        
        return rho
    
    def _electron_density(self, region: SemiconductorRegion, 
                         energy: float, fermi_level: float) -> float:
        """
        Calculate electron density in conduction band.
        
        Based on Fortran RHOCB function.
        
        Args:
            region: Semiconductor region
            energy: Local energy (eV, relative to vacuum) 
            fermi_level: Fermi level (eV, relative to valence band maximum)
            
        Returns:
            Electron density (cm^-3)
        """
        # Convert energy to potential relative to vacuum
        # In Fortran: energy parameter represents potential
        potential = energy
        
        # Conduction band edge relative to valence band = band_gap + valence_band_offset
        Ec_rel = region.band_gap + region.valence_band_offset
        
        # Reduced energy: (EF - Ec + potential) / kT 
        # From Fortran: (EF-EGAP(IREG)-DELVB(IREG)-Pot)/TK
        eta = (fermi_level - Ec_rel - potential) / region.kT
        
        # Effective density of states
        Nc = region.Nc
        
        # Check for band occupation suppression
        if region.suppress_conduction_band:
            return 0.0
        
        if region.allow_degeneracy:
            # Fermi-Dirac statistics using Fermi integral
            n = Nc * (2 / np.sqrt(np.pi)) * self.fermi_integral_half(eta)
        else:
            # Boltzmann statistics
            if eta > 40:  # Prevent overflow
                n = 0.0
            elif eta < -40:  # Degenerate limit
                n = Nc
            else:
                n = Nc * np.exp(eta)
        
        return n
    
    def _hole_density(self, region: SemiconductorRegion,
                     energy: float, fermi_level: float) -> float:
        """
        Calculate hole density in valence band.
        
        Based on Fortran RHOVB function.
        
        Args:
            region: Semiconductor region
            energy: Local energy (eV, relative to vacuum)
            fermi_level: Fermi level (eV, relative to valence band maximum)
            
        Returns:
            Hole density (cm^-3)
        """
        # Convert energy to potential relative to vacuum
        potential = energy
        
        # Valence band edge relative to valence band = valence_band_offset
        Ev_rel = region.valence_band_offset
        
        # Reduced energy for holes: (-EF + Ev + potential) / kT
        # From Fortran: (-EF+DELVB(IREG)+Pot)/TK
        eta = (-fermi_level + Ev_rel + potential) / region.kT
        
        # Effective density of states
        Nv = region.Nv
        
        # Check for band occupation suppression
        if region.suppress_valence_band:
            return 0.0
        
        if region.allow_degeneracy:
            # Fermi-Dirac statistics using Fermi integral
            p = Nv * (2 / np.sqrt(np.pi)) * self.fermi_integral_half(eta)
        else:
            # Boltzmann statistics
            if eta > 40:  # Prevent overflow
                p = 0.0
            elif eta < -40:  # Degenerate limit
                p = Nv
            else:
                p = Nv * np.exp(eta)
        
        return p
    
    def _electron_density_direct(self, region: SemiconductorRegion, 
                                fermi_level: float, potential: float) -> float:
        """
        Calculate electron density directly following Fortran RHOCB.
        
        Args:
            region: Semiconductor region
            fermi_level: Fermi level (eV, relative to valence band maximum)
            potential: Electrostatic potential (V)
            
        Returns:
            Electron density (cm^-3)
        """
        # Check for band occupation suppression
        if region.suppress_conduction_band:
            return 0.0
            
        # Conduction band edge relative to valence band = band_gap + valence_band_offset
        Ec_rel = region.band_gap + region.valence_band_offset
        
        # From Fortran RHOCB: (EF-EGAP(IREG)-DELVB(IREG)-Pot)/TK
        # Critical: Use the cb_effective_mass property which matches ACB
        acb = region.cb_effective_mass  # This is the Fortran ACB parameter
        
        # Energy difference relative to conduction band edge
        energy_diff = fermi_level - Ec_rel - potential
        eta = energy_diff / region.kT
        
        # Add extreme value protection
        if eta > 100:  # Prevent overflow in exponentials
            eta = 100
        elif eta < -100:
            eta = -100
        
        # Constants from Fortran: C = 6.815E21 eV^-1.5 cm^-3
        C = 6.815e21
        
        if region.temperature == 0:
            # T=0 case
            if energy_diff <= 0:
                return 0.0
            else:
                # From Fortran: (2.*C/3.)*SQRT((ACB(IREG)*(EF-EGAP-DELVB-Pot))**3)
                arg = acb * energy_diff
                if arg <= 0:
                    return 0.0
                return (2.0 * C / 3.0) * np.sqrt(arg**3)
        else:
            # Finite temperature case
            # From Fortran: C*SQRT((ACB(IREG)*TK)**3)*FJINT(1,(EF-EGAP-DELVB-Pot)/TK)
            prefactor = C * np.sqrt((acb * region.kT)**3)
            n = prefactor * self.fermi_integral_half(eta)
            
            # Add stability checks
            if not np.isfinite(n) or n < 0:
                return 0.0
            elif n > 1e50:
                return 1e50
            
            return n
    
    def _hole_density_direct(self, region: SemiconductorRegion,
                            fermi_level: float, potential: float) -> float:
        """
        Calculate hole density directly following Fortran RHOVB.
        
        Args:
            region: Semiconductor region
            fermi_level: Fermi level (eV, relative to valence band maximum) 
            potential: Electrostatic potential (V)
            
        Returns:
            Hole density (cm^-3)
        """
        # Check for band occupation suppression
        if region.suppress_valence_band:
            return 0.0
            
        # Critical: Use the vb_effective_mass_avg property which matches AVB
        avb = region.vb_effective_mass_avg  # This is the Fortran AVB parameter
        
        # Energy difference relative to valence band edge
        # From Fortran RHOVB: (-EF+DELVB(IREG)+Pot)/TK
        energy_diff = -fermi_level + region.valence_band_offset + potential
        eta = energy_diff / region.kT
        
        # Add extreme value protection
        if eta > 100:  # Prevent overflow in exponentials
            eta = 100
        elif eta < -100:
            eta = -100
        
        # Constants from Fortran: C = 6.815E21 eV^-1.5 cm^-3
        C = 6.815e21
        
        if region.temperature == 0:
            # T=0 case
            if energy_diff <= 0:
                return 0.0
            else:
                # From Fortran: (2.*C/3.)*SQRT((AVB(IREG)*(-EF+DELVB+Pot))**3)
                arg = avb * energy_diff
                if arg <= 0:
                    return 0.0
                return (2.0 * C / 3.0) * np.sqrt(arg**3)
        else:
            # Finite temperature case  
            # From Fortran: C*SQRT((AVB(IREG)*TK)**3)*FJINT(1,(-EF+DELVB+Pot)/TK)
            prefactor = C * np.sqrt((avb * region.kT)**3)
            p = prefactor * self.fermi_integral_half(eta)
            
            # Add stability checks
            if not np.isfinite(p) or p < 0:
                return 0.0
            elif p > 1e50:
                return 1e50
            
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
    
    def calculate(self, potential: np.ndarray) -> np.ndarray:
        """
        Calculate total charge density given electrostatic potential distribution.
        
        This method implements the physical charge density calculation for semiconductor
        regions, including ionized dopants and free carriers.
        
        Args:
            potential: 2D array of electrostatic potential (V) with shape (N_eta, N_nu)
                      matching the grid dimensions
            
        Returns:
            2D array of charge density (C/m³) with same shape as potential
        """
        # Validate input
        if potential.shape != (self.grid.N_eta, self.grid.N_nu):
            raise ValueError(f"Potential shape {potential.shape} doesn't match grid shape "
                           f"({self.grid.N_eta}, {self.grid.N_nu})")
        
        # Initialize charge density array
        rho = np.zeros_like(potential)
        
        # Get material properties from first semiconductor region (assuming single region for now)
        # TODO: Handle multiple regions based on spatial coordinates
        if not self.semiconductor_regions:
            raise ValueError("No semiconductor regions defined")
        
        region = list(self.semiconductor_regions.values())[0]
        
        # Physical constants
        q = PC.E  # Elementary charge in C
        
        # Temperature and thermal energy
        T = region.temperature
        kT = region.kT
        
        # Get donor and acceptor properties
        N_d = region.donor_concentration  # cm^-3
        N_a = region.acceptor_concentration  # cm^-3
        E_d = region.donor_binding_energy  # eV
        E_a = region.acceptor_binding_energy  # eV
        
        # Band edges (relative to vacuum level)
        E_c0 = region.band_gap + region.valence_band_offset  # Conduction band edge
        E_v0 = region.valence_band_offset  # Valence band edge
        
        # Fermi level (relative to valence band maximum)
        E_F = self.fermi_level
        
        # Calculate ionized dopant concentrations
        # For shallow dopants, we can assume complete ionization at room temperature
        # More accurate calculation would use Fermi-Dirac statistics:
        # N_d+ = N_d * (1 - f(E_d, E_F))
        # N_a- = N_a * f(E_a, E_F)
        
        if T > 0:
            # Donor ionization (assuming donor level below conduction band)
            # f_d = 1 / (1 + g_d * exp((E_d - E_F)/kT))
            # where g_d is the degeneracy factor (typically 2 for donors)
            g_d = 2.0
            E_d_abs = E_c0 - E_d  # Donor level relative to vacuum
            
            # Calculate at each grid point considering local potential
            E_F_local = E_F - potential  # Local Fermi level including potential
            f_d = 1.0 / (1.0 + g_d * np.exp((E_d_abs - E_F_local) / kT))
            N_d_plus = N_d * (1.0 - f_d)
            
            # Acceptor ionization
            g_a = 4.0  # Degeneracy factor for acceptors
            E_a_abs = E_v0 + E_a  # Acceptor level relative to vacuum
            f_a = 1.0 / (1.0 + g_a * np.exp((E_a_abs - E_F_local) / kT))
            N_a_minus = N_a * f_a
        else:
            # At T=0, use step functions
            E_F_local = E_F - potential
            N_d_plus = np.where(E_F_local > E_c0 - E_d, N_d, 0.0)
            N_a_minus = np.where(E_F_local < E_v0 + E_a, N_a, 0.0)
        
        # Calculate carrier concentrations
        # The band edges are bent by the potential
        E_c = E_c0 - q * potential  # Conduction band edge with band bending
        E_v = E_v0 - q * potential  # Valence band edge with band bending
        
        # Use the direct calculation methods that match Fortran
        for i in range(self.grid.N_eta):
            for j in range(self.grid.N_nu):
                # Electron density using Fermi-Dirac integral
                n = self._electron_density_direct(region, E_F, potential[i, j])
                
                # Hole density using Fermi-Dirac integral
                p = self._hole_density_direct(region, E_F, potential[i, j])
                
                # Total charge density: ρ = q * (N_d+ - N_a- - n + p)
                # Note: electrons contribute negative charge, holes positive
                rho[i, j] = q * (N_d_plus[i, j] - N_a_minus[i, j] - n + p)
        
        # Convert from C/cm³ to C/m³
        rho *= 1e6
        
        return rho
    
    def create_charge_tables(self, energy_range: Tuple[float, float],
                           num_points: int = 20000) -> ChargeDensityTables:
        """
        Create tabulated charge densities for all regions following Fortran SEMIRHO.
        
        The table represents charge density vs. Fermi level (at zero potential).
        This matches the Fortran SEMIRHO subroutine logic where:
        - EF1 = (I-1)*DELE + ESTART represents the variable Fermi level
        - RHOBTAB(IREG,I) = RHOB(IREG,EF1,0.) - charge density at EF1, Pot=0
        
        Args:
            energy_range: (min, max) energy range (eV) for Fermi level variation
            num_points: Number of energy points
            
        Returns:
            ChargeDensityTables object
        """
        energy_start, energy_end = energy_range
        energy_step = (energy_end - energy_start) / (num_points - 1)
        
        # Energy array represents Fermi level values
        fermi_levels = np.linspace(energy_start, energy_end, num_points)
        
        print(f"Creating charge density tables:")
        print(f"  Fermi level range: {energy_start:.6f} to {energy_end:.6f} eV")
        print(f"  Energy step: {energy_step:.6f} eV")
        print(f"  Number of points: {num_points}")
        
        # Create bulk tables - following Fortran SEMIRHO exactly
        bulk_tables = {}
        for region_id, region in self.semiconductor_regions.items():
            print(f"Computing bulk charge density table for region {region_id}")
            densities = np.zeros(num_points)
            
            for i, ef_var in enumerate(fermi_levels):
                # Calculate charge density at this Fermi level and zero potential
                # This matches Fortran: RHOBTAB(IREG,I) = RHOB(IREG,EF1,0.)
                # where EF1 = (I-1)*DELE + ESTART is the variable Fermi level
                density = self.calculate_bulk_density(region_id, ef_var, potential=0.0)
                
                # Add numerical stability checks
                if not np.isfinite(density):
                    density = 0.0
                elif abs(density) > 1e50:  # Prevent extreme values
                    density = np.sign(density) * 1e50
                
                densities[i] = density
                
                # Debug output for first few points
                if i < 5 or i % 1000 == 0:
                    n = self._electron_density_direct(region, ef_var, 0.0)
                    p = self._hole_density_direct(region, ef_var, 0.0)
                    print(f"    Point {i}: EF={ef_var:.6f}, n={n:.3e}, p={p:.3e}, rho={density:.3e}")
            
            bulk_tables[region_id] = densities
            print(f"  Region {region_id}: density range {np.min(densities):.3e} to {np.max(densities):.3e}")
        
        # Create surface tables - following Fortran SURFRHO
        surface_tables = {}
        for area_id, region in self.surface_regions.items():
            print(f"Computing surface charge density table for area {area_id}")
            densities = np.zeros(num_points)
            
            for i, ef_var in enumerate(fermi_levels):
                # For surface states, the table represents surface charge vs. surface Fermi level
                # Surface potential = 0 in table construction
                densities[i] = self.calculate_surface_density(area_id, ef_var, potential=0.0)
            
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