"""
Schrödinger equation solver for tunneling current calculations.

This module implements the quantum mechanical tunneling current calculations
from INTCURR in the Fortran code, using the transfer matrix method and
WKB approximation.
"""

import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
from scipy import integrate, special
from scipy.integrate import solve_ivp
import warnings

from ...utils.constants import PhysicalConstants as PC
from .potential import PotentialProfile


@dataclass
class BandProfile:
    """Container for band edge profiles."""
    z: np.ndarray              # Position array (nm)
    vb_profile: np.ndarray     # Valence band edge profile (eV)
    cb_profile: np.ndarray     # Conduction band edge profile (eV)
    vb_max: float             # Maximum VB energy
    cb_min: float             # Minimum CB energy
    vb_bulk: float            # VB edge in bulk
    cb_bulk: float            # CB edge in bulk


@dataclass
class WaveFunction:
    """Container for wave function data."""
    z: np.ndarray              # Position array
    psi: np.ndarray           # Wave function amplitude
    energy: float             # Energy eigenvalue
    k_parallel: float         # Parallel wave vector
    band: str                 # 'vb' or 'cb'
    localized: bool           # True if bound state


@dataclass
class TunnelCurrent:
    """Container for tunneling current results."""
    total_current: float       # Total current (A)
    vb_current: float         # Valence band contribution
    cb_current: float         # Conduction band contribution
    vb_localized: float       # VB localized state contribution
    cb_localized: float       # CB localized state contribution
    localized_states: List[WaveFunction]  # List of bound states


class SchrodingerSolver:
    """
    Solver for 1D Schrödinger equation in tunneling problems.
    
    Implements transfer matrix method and WKB approximation
    for calculating tunneling current.
    """
    
    def __init__(self, mass_electron: float = PC.M0):
        """
        Initialize Schrödinger solver.
        
        Args:
            mass_electron: Electron mass (kg)
        """
        self.m0 = mass_electron
        self.hbar = PC.HBAR
        
        # Numerical parameters
        self.energy_tolerance = 1e-6  # eV
        self.wavefunction_tolerance = 1e-8
        self.max_iterations = 100
    
    def solve_tunneling_current(self, potential_profile: PotentialProfile,
                              band_params: Dict,
                              bias_voltage: float,
                              temperature: float = 300.0) -> TunnelCurrent:
        """
        Calculate tunneling current through the potential barrier.
        
        Args:
            potential_profile: Potential profile from Poisson solver
            band_params: Band structure parameters
            bias_voltage: Applied bias (V)
            temperature: Temperature (K)
            
        Returns:
            TunnelCurrent object with current components
        """
        # Extract band parameters
        band_gap_semi = band_params['band_gap']
        m_cb = band_params['cb_effective_mass'] * self.m0
        m_vb_heavy = band_params['vb_effective_mass_heavy'] * self.m0
        m_vb_light = band_params['vb_effective_mass_light'] * self.m0
        
        # fermi_level_semi_bulk_rel_VB is EF_semi_bulk - Ev_bulk_semi (from EFFIND)
        fermi_level_semi_bulk_rel_VB = band_params['fermi_level']
        electron_affinity_semi = band_params['electron_affinity'] # CHI_semi

        # Create band profiles relative to EF_semi_bulk_abs = 0 eV
        band_profile = self._create_band_profiles(
            potential_profile,
            band_gap_semi,
            electron_affinity_semi,
            fermi_level_semi_bulk_rel_VB
        )
        
        # Fermi level for semiconductor is now 0.0 eV (our reference)
        fermi_level_semi_ref = 0.0
        # Tip Fermi level relative to semiconductor bulk Fermi level
        tip_fermi_ref = -PC.E_CHARGE * bias_voltage # Assuming sample grounded, V_bias = V_tip

        # Find localized states
        vb_states = self._find_localized_states(band_profile, 'vb', 
                                               [m_vb_light, m_vb_heavy])
        cb_states = self._find_localized_states(band_profile, 'cb', [m_cb])
        
        # Calculate current from extended states
        vb_current_ext = self._calculate_extended_current(
            band_profile, 'vb', [m_vb_light, m_vb_heavy],
            fermi_level_semi_ref, tip_fermi_ref, temperature
        )
        
        cb_current_ext = self._calculate_extended_current(
            band_profile, 'cb', [m_cb],
            fermi_level_semi_ref, tip_fermi_ref, temperature
        )
        
        # Calculate current from localized states
        vb_current_loc = self._calculate_localized_current(
            vb_states, fermi_level_semi_ref, tip_fermi_ref, temperature
        )
        
        cb_current_loc = self._calculate_localized_current(
            cb_states, fermi_level_semi_ref, tip_fermi_ref, temperature
        )
        
        # Total currents
        vb_current = vb_current_ext + vb_current_loc
        cb_current = cb_current_ext + cb_current_loc
        total_current = vb_current + cb_current
        
        return TunnelCurrent(
            total_current=total_current,
            vb_current=vb_current,
            cb_current=cb_current,
            vb_localized=vb_current_loc,
            cb_localized=cb_current_loc,
            localized_states=vb_states + cb_states
        )
    
    def _create_band_profiles(self, potential_profile: PotentialProfile,
                            band_gap_semi: float,
                            electron_affinity_semi: float,
                            fermi_level_semi_bulk_rel_VB: float) -> BandProfile:
        """
        Create band edge profiles relative to semiconductor bulk Fermi level (EF_semi_bulk_abs = 0 eV).

        Args:
            potential_profile: Contains z and phi(z) for vacuum and semiconductor.
                               phi(z) is electrostatic potential with phi(deep_semiconductor)=0.
            band_gap_semi: Semiconductor band gap (eV).
            electron_affinity_semi: Semiconductor electron affinity (CHI_semi) (eV).
            fermi_level_semi_bulk_rel_VB: EF_semi_bulk - Ev_bulk_semi, from EFFIND (eV).
        """
        # Ensure z_semiconductor is increasing for concatenation after flip
        # Typically, z_semiconductor from PotentialProfile is already like [-depth, ..., -dz, 0]
        # If it's [0, -dz, ..., -depth], it needs flipping before use if get_combined_profile isn't used.
        # PotentialProfile.get_combined_profile() handles this.
        z_combined, phi_combined = potential_profile.get_combined_profile() # phi(z) with phi(deep_semi)=0

        # Semiconductor bulk Fermi level (EF_semi_bulk_abs) is our energy reference (0 eV).
        
        # Ev_bulk_offset = Ev_bulk_semi - EF_semi_bulk_abs
        Ev_bulk_offset = -fermi_level_semi_bulk_rel_VB
        # Ec_bulk_offset = Ec_bulk_semi - EF_semi_bulk_abs
        Ec_bulk_offset = Ev_bulk_offset + band_gap_semi

        # Separate semiconductor and vacuum parts of phi_combined based on z_combined
        # Assuming z_combined is sorted and 0.0 marks the interface
        # Find index of z=0 (interface)
        interface_idx = np.argmin(np.abs(z_combined)) # Index where z is closest to 0

        # Handle cases: z_combined might not have an exact 0.
        # If all z < 0 (only semiconductor part in profile, unlikely for tunneling problem)
        # If all z >= 0 (only vacuum part in profile, unlikely)

        # A robust way to split based on original PotentialProfile structure:
        phi_semi = potential_profile.potential_semiconductor # phi(z) for z_semi
        phi_vac = potential_profile.potential_vacuum       # phi(z) for z_vac
        z_semi = potential_profile.z_semiconductor         # z < 0, typically increasing from -depth_max to 0 or near 0
        z_vac = potential_profile.z_vacuum                 # z > 0, typically increasing from near 0 to vac_max

        # Ensure semiconductor z values are negative or zero, and vacuum z are positive or zero
        # And that phi_semi corresponds to z_semi, phi_vac to z_vac.
        # PotentialProfile.get_combined_profile() already gives sorted z_combined and corresponding phi_combined.
        # We need to know which part of phi_combined belongs to semiconductor and which to vacuum.

        # Create masks for semiconductor and vacuum regions based on z_combined
        is_semiconductor = z_combined < 0
        is_vacuum = z_combined >= 0
        # Add a small tolerance for z=0 to be included in vacuum if it's duplicated
        if z_combined[interface_idx] == 0:
            is_vacuum[interface_idx] = True
            if interface_idx > 0 and z_combined[interface_idx-1] == 0: # if z=0 is duplicated due to concat
                 is_semiconductor[interface_idx] = False


        # Band edges in semiconductor (phi_combined is already relative to EF_semi_bulk if PoissonSolver sets bulk potential to 0)
        # PC.E_CHARGE is positive elementary charge. For electron energy, use -e*phi.
        # If phi is potential, energy for electron is -phi (if e=1).
        # If phi is potential in Volts, energy is -e_charge * phi in Joules. Divide by e_charge for eV.
        # So, energy shift in eV due to potential phi (in Volts) is -phi.

        cb_profile_semi_part = Ec_bulk_offset - phi_combined[is_semiconductor]
        vb_profile_semi_part = Ev_bulk_offset - phi_combined[is_semiconductor]

        # Potential energy in vacuum (for electrons)
        # U_vac(z) = (E_vacuum_level_at_interface_no_bending) - phi(z_interface) + CHI_semi - phi(z_vac)
        # U_vac(z) = (Ec_bulk_offset + electron_affinity_semi) - phi_vac_part
        cb_profile_vac_part = (Ec_bulk_offset + electron_affinity_semi) - phi_combined[is_vacuum]

        # Valence band in vacuum: For electron tunneling, this is not directly used.
        # Set it far from the CB to avoid issues, or based on a large vacuum "band gap".
        vb_profile_vac_part = cb_profile_vac_part - 10.0 # Arbitrary large gap (10 eV)

        # Stitch together the profiles. z_combined is already sorted.
        cb_profile = np.zeros_like(z_combined)
        vb_profile = np.zeros_like(z_combined)

        cb_profile[is_semiconductor] = cb_profile_semi_part
        cb_profile[is_vacuum] = cb_profile_vac_part
        vb_profile[is_semiconductor] = vb_profile_semi_part
        vb_profile[is_vacuum] = vb_profile_vac_part
        
        # Find extrema and bulk values (relative to EF_semi_bulk_abs = 0 eV)
        vb_max = np.max(vb_profile) # Max VB energy (could be in accumulation layer)
        cb_min = np.min(cb_profile) # Min CB energy (could be in inversion layer or tip)
        
        # Bulk values are now Ec_bulk_offset and Ev_bulk_offset by definition
        vb_bulk = Ev_bulk_offset
        cb_bulk = Ec_bulk_offset
        
        return BandProfile(
            z=z_combined,
            vb_profile=vb_profile,
            cb_profile=cb_profile,
            vb_max=vb_max,
            cb_min=cb_min,
            vb_bulk=vb_bulk,
            cb_bulk=cb_bulk
        )
    
    def _find_localized_states(self, band_profile: BandProfile,
                             band_type: str,
                             effective_masses: List[float]) -> List[WaveFunction]:
        """
        Find localized (bound) states in the potential well.
        
        Args:
            band_profile: Band edge profiles
            band_type: 'vb' or 'cb'
            effective_masses: List of effective masses for different bands
            
        Returns:
            List of localized wave functions
        """
        localized_states = []
        
        # Get appropriate profile
        if band_type == 'vb':
            profile = band_profile.vb_profile
            e_bulk = band_profile.vb_bulk
            search_range = (e_bulk, band_profile.vb_max)
        else:
            profile = band_profile.cb_profile
            e_bulk = band_profile.cb_bulk
            search_range = (band_profile.cb_min, e_bulk)
        
        # Search for bound states for each effective mass
        for m_eff in effective_masses:
            states = self._find_bound_states_1d(
                band_profile.z, profile, m_eff, search_range
            )
            
            for state in states:
                state.band = band_type
                localized_states.append(state)
        
        return localized_states
    
    def _find_bound_states_1d(self, z: np.ndarray, potential: np.ndarray,
                            m_eff: float,
                            energy_range: Tuple[float, float]) -> List[WaveFunction]:
        """
        Find bound states in 1D potential using shooting method.
        
        Args:
            z: Position array
            potential: Potential profile
            m_eff: Effective mass
            energy_range: (min, max) energy to search
            
        Returns:
            List of bound state wave functions
        """
        bound_states = []
        
        # Energy grid for search
        n_energy = 100
        energies = np.linspace(energy_range[0], energy_range[1], n_energy)
        
        # Look for sign changes in wave function at boundary
        prev_sign = None
        
        for i, energy in enumerate(energies):
            # Skip if classically forbidden everywhere
            if np.all(energy < potential):
                continue
            
            # Solve Schrödinger equation
            psi = self._solve_1d_schrodinger(z, potential, energy, m_eff)
            
            # Check boundary condition (psi -> 0 as z -> ∞)
            boundary_value = psi[-1]
            current_sign = np.sign(boundary_value)
            
            # Sign change indicates eigenvalue
            if prev_sign is not None and current_sign != prev_sign:
                # Refine energy using bisection
                e_refined = self._refine_eigenvalue(
                    z, potential, energies[i-1], energy, m_eff
                )
                
                # Get wave function at refined energy
                psi_refined = self._solve_1d_schrodinger(
                    z, potential, e_refined, m_eff
                )
                
                # Normalize
                psi_normalized = self._normalize_wavefunction(z, psi_refined)
                
                bound_states.append(WaveFunction(
                    z=z,
                    psi=psi_normalized,
                    energy=e_refined,
                    k_parallel=0.0,
                    band='',
                    localized=True
                ))
            
            prev_sign = current_sign
        
        return bound_states
    
    def _solve_1d_schrodinger(self, z: np.ndarray, potential: np.ndarray,
                            energy: float, m_eff: float) -> np.ndarray:
        """
        Solve 1D Schrödinger equation using finite differences.
        
        -ℏ²/2m d²ψ/dz² + V(z)ψ = Eψ
        
        Args:
            z: Position array (nm).
            potential: Potential profile (eV).
            energy: Energy eigenvalue (eV).
            m_eff: Effective mass (kg).
            
        Returns:
            Wave function array (unnormalized).
        """
        n = len(z)
        if n < 3:
            raise ValueError("z array must have at least 3 points for Numerov method.")

        psi = np.zeros(n)
        dz_nm = z[1] - z[0]  # Grid spacing in nm
        dz_m = dz_nm * 1e-9   # Grid spacing in meters

        # k²(z) = 2m*/ħ² * (E - V(z))
        # Ensure (E - V(z)) is converted from eV to Joules
        k_squared_z_values = (2 * m_eff / (self.hbar**2)) * (energy - potential) * PC.E_CHARGE

        # f(z) = 1 + (Δz²/12)k²(z)
        f_z_values = 1.0 + ((dz_m**2) / 12.0) * k_squared_z_values

        # Initial conditions for shooting method (bound state decaying from z_min)
        psi[0] = 0.0
        psi[1] = 1e-5  # Small non-zero value

        # Numerov recurrence relation:
        # ψ(i+1) = ( (12 - 10*f(i))ψ(i) - f(i-1)ψ(i-1) ) / f(i+1)
        for i in range(1, n - 1):
            # Check for potential division by zero if f_z_values[i+1] is too small
            # This might indicate a problem with the potential or energy range
            if abs(f_z_values[i+1]) < 1e-12: # Threshold to avoid division by zero
                # If f(z) is near zero, it's problematic. Wavefunction might be diverging.
                # For simplicity, fill rest with NaNs or large numbers to indicate failure.
                psi[i+1:] = np.nan
                warnings.warn(f"Numerov f(z) near zero at z={z[i+1]:.3f} nm, energy={energy:.3f} eV. Wavefunction may be unreliable.", RuntimeWarning)
                break

            numerator = ((12.0 - 10.0 * f_z_values[i]) * psi[i]) - (f_z_values[i-1] * psi[i-1])
            psi[i+1] = numerator / f_z_values[i+1]

            # Check for NaN/Inf propagation
            if not np.isfinite(psi[i+1]):
                warnings.warn(f"Numerov produced NaN/Inf at z={z[i+1]:.3f} nm, energy={energy:.3f} eV. Previous psi: {psi[i]:.3e}, {psi[i-1]:.3e}", RuntimeWarning)
                # Fill rest with NaNs if propagation occurs
                psi[i+1:] = np.nan
                break
        
        return psi
    
    def _refine_eigenvalue(self, z: np.ndarray, potential: np.ndarray,
                         e1: float, e2: float, m_eff: float) -> float:
        """Refine eigenvalue using bisection method."""
        for _ in range(20):  # Max iterations
            e_mid = (e1 + e2) / 2
            psi_mid = self._solve_1d_schrodinger(z, potential, e_mid, m_eff)
            
            psi1 = self._solve_1d_schrodinger(z, potential, e1, m_eff)
            
            if np.sign(psi_mid[-1]) == np.sign(psi1[-1]):
                e1 = e_mid
            else:
                e2 = e_mid
            
            if abs(e2 - e1) < self.energy_tolerance:
                break
        
        return (e1 + e2) / 2
    
    def _normalize_wavefunction(self, z: np.ndarray, psi: np.ndarray) -> np.ndarray:
        """Normalize wave function."""
        dz = z[1] - z[0]
        norm_sq = np.trapz(np.abs(psi)**2, dx=dz)
        return psi / np.sqrt(norm_sq)
    
    def _calculate_extended_current(self, band_profile: BandProfile,
                                  band_type: str,
                                  effective_masses: List[float],
                                  fermi_level: float,
                                  tip_fermi: float,
                                  temperature: float) -> float:
        """
        Calculate current from extended states using WKB approximation.
        
        Args:
            band_profile: Band edge profiles
            band_type: 'vb' or 'cb'
            effective_masses: List of effective masses
            fermi_level: Semiconductor Fermi level
            tip_fermi: Tip Fermi level
            temperature: Temperature (K)
            
        Returns:
            Current contribution (A)
        """
        kT = PC.thermal_energy(temperature)
        current = 0.0
        
        # Energy integration range
        if band_type == 'vb':
            e_min = band_profile.vb_bulk - 5 * kT
            e_max = min(fermi_level, tip_fermi) + 5 * kT
        else:
            e_min = max(fermi_level, tip_fermi) - 5 * kT
            e_max = band_profile.cb_bulk + 5 * kT
        
        # Skip if no energy window
        if e_max <= e_min:
            return 0.0
        
        # Energy grid
        energies = np.linspace(e_min, e_max, 100)
        
        for m_eff in effective_masses:
            for energy in energies:
                # Transmission coefficient
                T = self._calculate_transmission_wkb(
                    band_profile, band_type, energy, m_eff
                )
                
                # Fermi factors
                f_semi = 1.0 / (1.0 + np.exp((energy - fermi_level) / kT))
                f_tip = 1.0 / (1.0 + np.exp((energy - tip_fermi) / kT))
                
                # Current contribution (simplified - full version needs DOS)
                # I = (2e/h) ∫ T(E) [f_semi(E) - f_tip(E)] dE
                if band_type == 'vb':
                    # Holes tunnel from occupied states
                    dI = T * ((1 - f_semi) - (1 - f_tip))
                else:
                    # Electrons tunnel to empty states
                    dI = T * (f_semi - f_tip)
                
                current += dI
        
        # Convert to proper units and multiply by DOS factors
        # This is simplified - full implementation needs proper DOS
        de = energies[1] - energies[0]
        prefactor = 2 * PC.E / PC.H  # 2e/h for spin
        current *= prefactor * de
        
        return current
    
    def _calculate_transmission_wkb(self, band_profile: BandProfile,
                                  band_type: str, energy: float,
                                  m_eff: float) -> float:
        """
        Calculate transmission coefficient using WKB approximation.
        
        T = exp(-2∫κ(z)dz) where κ = √(2m[V(z)-E])/ℏ
        
        Args:
            band_profile: Band profiles
            band_type: 'vb' or 'cb'
            energy: Tunneling energy
            m_eff: Effective mass
            
        Returns:
            Transmission coefficient
        """
        z = band_profile.z
        
        if band_type == 'vb':
            potential = band_profile.vb_profile
        else:
            potential = band_profile.cb_profile
        
        # Find classically forbidden region
        forbidden = potential > energy
        
        if not np.any(forbidden):
            return 1.0  # No barrier
        
        # Find turning points
        indices = np.where(forbidden)[0]
        if len(indices) == 0:
            return 1.0
        
        z1_idx = indices[0]
        z2_idx = indices[-1]
        
        # WKB integral
        z_barrier = z[z1_idx:z2_idx+1]
        v_barrier = potential[z1_idx:z2_idx+1]
        
        # Calculate local wave vector
        kappa = np.sqrt(2 * m_eff * PC.E * np.maximum(0, v_barrier - energy)) / self.hbar
        
        # Integrate
        integral = np.trapz(kappa, z_barrier * 1e-9)  # Convert nm to m
        
        # Transmission coefficient
        T = np.exp(-2 * integral)
        
        return min(1.0, T)  # Cap at 1
    
    def _calculate_localized_current(self, localized_states: List[WaveFunction],
                                   fermi_level: float, tip_fermi: float,
                                   temperature: float) -> float:
        """
        Calculate current from localized states.
        
        Args:
            localized_states: List of bound states
            fermi_level: Semiconductor Fermi level
            tip_fermi: Tip Fermi level
            temperature: Temperature (K)
            
        Returns:
            Current contribution (A)
        """
        if not localized_states:
            return 0.0
        
        kT = PC.thermal_energy(temperature)
        current = 0.0
        
        for state in localized_states:
            # Matrix element (simplified - needs overlap integral)
            # |M|² ∝ |ψ(0)|² for s-wave tip state
            matrix_element = np.abs(state.psi[0])**2
            
            # Fermi factors
            f_semi = 1.0 / (1.0 + np.exp((state.energy - fermi_level) / kT))
            f_tip = 1.0 / (1.0 + np.exp((state.energy - tip_fermi) / kT))
            
            # Current contribution
            if state.band == 'vb':
                dI = matrix_element * ((1 - f_semi) - (1 - f_tip))
            else:
                dI = matrix_element * (f_semi - f_tip)
            
            current += dI
        
        # Prefactor (simplified - needs proper constants)
        prefactor = 2 * np.pi * PC.E / self.hbar
        current *= prefactor
        
        return current


def create_schrodinger_solver() -> SchrodingerSolver:
    """Create a Schrödinger equation solver."""
    return SchrodingerSolver()