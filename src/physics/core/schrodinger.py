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
        band_gap = band_params['band_gap']
        m_cb = band_params['cb_effective_mass'] * self.m0
        m_vb_heavy = band_params['vb_effective_mass_heavy'] * self.m0
        m_vb_light = band_params['vb_effective_mass_light'] * self.m0
        fermi_level = band_params['fermi_level']
        tip_fermi = band_params['tip_fermi_level']
        
        # Create band profiles
        band_profile = self._create_band_profiles(potential_profile, band_gap)
        
        # Find localized states
        vb_states = self._find_localized_states(band_profile, 'vb', 
                                               [m_vb_light, m_vb_heavy])
        cb_states = self._find_localized_states(band_profile, 'cb', [m_cb])
        
        # Calculate current from extended states
        vb_current_ext = self._calculate_extended_current(
            band_profile, 'vb', [m_vb_light, m_vb_heavy],
            fermi_level, tip_fermi, temperature
        )
        
        cb_current_ext = self._calculate_extended_current(
            band_profile, 'cb', [m_cb],
            fermi_level, tip_fermi, temperature
        )
        
        # Calculate current from localized states
        vb_current_loc = self._calculate_localized_current(
            vb_states, fermi_level, tip_fermi, temperature
        )
        
        cb_current_loc = self._calculate_localized_current(
            cb_states, fermi_level, tip_fermi, temperature
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
                            band_gap: float) -> BandProfile:
        """Create band edge profiles from potential."""
        # Combine z arrays
        z_combined, pot_combined = potential_profile.get_combined_profile()
        
        # Band edges (simplified - assumes electron affinity alignment)
        # In full implementation, would need proper band alignment
        vb_profile = pot_combined - band_gap
        cb_profile = pot_combined
        
        # Find extrema
        vb_max = np.max(vb_profile)
        cb_min = np.min(cb_profile)
        
        # Bulk values (far from surface)
        vb_bulk = vb_profile[-1]
        cb_bulk = cb_profile[-1]
        
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
            z: Position array
            potential: Potential profile
            energy: Energy eigenvalue
            m_eff: Effective mass
            
        Returns:
            Wave function array
        """
        # Convert to dimensionless units for numerical stability
        z_nm = z
        dz = z[1] - z[0]
        
        # Finite difference matrix
        n = len(z)
        H = np.zeros((n, n))
        
        # Kinetic energy term
        ke_factor = -self.hbar**2 / (2 * m_eff * (dz * 1e-9)**2 * PC.E)  # Convert to eV
        
        for i in range(1, n-1):
            H[i, i-1] = ke_factor
            H[i, i] = -2 * ke_factor + potential[i]
            H[i, i+1] = ke_factor
        
        # Boundary conditions
        H[0, 0] = -2 * ke_factor + potential[0]
        H[0, 1] = ke_factor
        H[-1, -2] = ke_factor
        H[-1, -1] = -2 * ke_factor + potential[-1]
        
        # Find wave function using inverse iteration
        # For bound state, use exponential decay boundary conditions
        psi = np.zeros(n)
        psi[0] = 1.0  # Arbitrary normalization
        
        # Simple integration (more sophisticated methods needed for production)
        for i in range(1, n):
            if i == 1:
                psi[i] = psi[0] * np.exp(-np.sqrt(2 * m_eff * PC.E * 
                                                 max(0, potential[0] - energy)) * 
                                        dz * 1e-9 / self.hbar)
            else:
                # Numerov method would be better here
                # Ensure we don't go out of bounds
                pot_idx = min(i, len(potential) - 1)
                k_sq = 2 * m_eff * PC.E * (energy - potential[pot_idx]) / self.hbar**2
                if k_sq > 0:
                    psi[i] = 2 * psi[i-1] - psi[i-2] + k_sq * (dz * 1e-9)**2 * psi[i-1]
                else:
                    kappa = np.sqrt(-k_sq)
                    psi[i] = psi[i-1] * np.exp(-kappa * dz * 1e-9)
        
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