"""
Potential profile extraction and manipulation.

This module implements functionality from POTCUT3 and POTEXPAND for extracting
potential profiles along specific paths and expanding solutions to finer grids.
"""

import numpy as np
from typing import Tuple, Optional, List
from scipy import interpolate
from dataclasses import dataclass


@dataclass
class PotentialProfile:
    """Container for potential profile data."""
    # Coordinate arrays
    z_vacuum: np.ndarray      # Z coordinates in vacuum (positive)
    z_semiconductor: np.ndarray  # Z coordinates in semiconductor (negative)
    
    # Potential arrays
    potential_vacuum: np.ndarray     # Potential in vacuum
    potential_semiconductor: np.ndarray  # Potential in semiconductor
    barrier_height: np.ndarray      # Barrier profile for tunneling
    
    # Metadata
    r_position: float         # Radial position of profile
    phi_position: float       # Angular position of profile
    band_bending: float       # Band bending at interface
    
    def get_combined_profile(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get combined z and potential arrays."""
        z_combined = np.concatenate([np.flip(self.z_semiconductor), self.z_vacuum])
        pot_combined = np.concatenate([np.flip(self.potential_semiconductor), 
                                     self.potential_vacuum])
        return z_combined, pot_combined


class PotentialProcessor:
    """
    Process and extract potential profiles from 3D solutions.
    
    Implements POTCUT3 and POTEXPAND functionality.
    """
    
    def __init__(self, grid: Grid3D):
        """
        Initialize potential processor.
        
        Args:
            grid: 3D computational grid
        """
        self.grid = grid
    
    def extract_profile(self, potential_3d: np.ndarray, 
                       r_position: float = 0.0,
                       phi_position: float = 0.0,
                       include_barrier: bool = True) -> PotentialProfile:
        """
        Extract 1D potential profile at specified position.
        
        Implements POTCUT3 functionality.
        
        Args:
            potential_3d: 3D potential array
            r_position: Radial position for profile (nm)
            phi_position: Angular position for profile (radians)
            include_barrier: Whether to calculate barrier profile
            
        Returns:
            PotentialProfile object
        """
        # Find nearest grid indices
        ir = np.argmin(np.abs(self.grid.r - r_position))
        ip = np.argmin(np.abs(self.grid.phi - phi_position))
        
        # Extract profiles
        ns = self.grid.params.ns
        
        # Semiconductor profile
        z_semi = self.grid.zs
        pot_semi = potential_3d[ir, :ns-1, ip]
        pot_semi = np.flip(pot_semi)  # Make z increasing
        
        # Vacuum profile  
        z_vac = self.grid.zv[1:]  # Skip interface point
        pot_vac = potential_3d[ir, ns:, ip]
        
        # Interface potential (band bending)
        band_bending = potential_3d[ir, ns-1, ip]
        
        # Calculate barrier profile if requested
        if include_barrier:
            barrier = self._calculate_barrier_profile(
                z_vac, pot_vac, z_semi, pot_semi, band_bending
            )
        else:
            barrier = np.zeros_like(np.concatenate([z_semi, z_vac]))
        
        return PotentialProfile(
            z_vacuum=z_vac,
            z_semiconductor=z_semi,
            potential_vacuum=pot_vac,
            potential_semiconductor=pot_semi,
            barrier_height=barrier,
            r_position=self.grid.r[ir],
            phi_position=self.grid.phi[ip],
            band_bending=band_bending
        )
    
    def _calculate_barrier_profile(self, z_vac: np.ndarray, pot_vac: np.ndarray,
                                 z_semi: np.ndarray, pot_semi: np.ndarray,
                                 band_bending: float) -> np.ndarray:
        """
        Calculate tunneling barrier profile.
        
        The barrier includes:
        - Vacuum barrier
        - Band structure
        - Image potential effects (if implemented)
        """
        # For now, simple barrier calculation
        # Full implementation would include image potential
        
        # Combine arrays
        z_combined = np.concatenate([z_semi, z_vac])
        pot_combined = np.concatenate([pot_semi, pot_vac])
        
        # Barrier is potential relative to Fermi level
        # (More sophisticated calculation needed for full implementation)
        barrier = pot_combined.copy()
        
        # In semiconductor: add band edge
        n_semi = len(z_semi)
        # Simplified - assumes conduction band tunneling
        barrier[:n_semi] += 1.42  # GaAs band gap
        
        return barrier
    
    def expand_potential(self, potential_3d: np.ndarray,
                        expansion_factor: int = 2) -> Tuple[np.ndarray, Grid3D]:
        """
        Expand potential to finer grid.
        
        Implements POTEXPAND functionality for more accurate current calculations.
        
        Args:
            potential_3d: Original potential on coarse grid
            expansion_factor: Grid refinement factor
            
        Returns:
            Tuple of (expanded potential, new grid)
        """
        # Create new grid parameters
        old_params = self.grid.params
        
        # New grid dimensions
        nr_new = (old_params.nr - 1) * expansion_factor + 1
        nv_new = (old_params.nv - 1) * expansion_factor + 1
        ns_new = (old_params.ns - 1) * expansion_factor + 1
        np_new = old_params.np  # Keep angular resolution
        
        # New grid spacing
        delr_new = old_params.delr / expansion_factor
        delv_new = old_params.delv / expansion_factor
        dels_new = old_params.dels / expansion_factor
        
        # Create interpolation function
        # Use scipy's RegularGridInterpolator for efficiency
        
        # Original grid points
        r_old = self.grid.r
        z_old = self.grid.z
        phi_old = self.grid.phi
        
        # New grid points
        r_new = np.linspace(0, old_params.rmax, nr_new)
        z_new = np.concatenate([
            -np.linspace(0, old_params.smax, ns_new)[::-1][:-1],
            np.linspace(0, old_params.vmax, nv_new)
        ])
        phi_new = phi_old  # Keep same
        
        # Create interpolator
        interp_func = interpolate.RegularGridInterpolator(
            (r_old, z_old, phi_old),
            potential_3d,
            method='linear',
            bounds_error=False,
            fill_value=0.0
        )
        
        # Generate new grid points
        R_new, Z_new, PHI_new = np.meshgrid(r_new, z_new, phi_new, indexing='ij')
        points = np.column_stack([R_new.ravel(), Z_new.ravel(), PHI_new.ravel()])
        
        # Interpolate
        potential_new = interp_func(points).reshape(nr_new, len(z_new), np_new)
        
        # Create new grid object (simplified - full implementation needed)
        # For now, return the interpolated potential
        return potential_new, None
    
    def extract_multiple_profiles(self, potential_3d: np.ndarray,
                                r_positions: List[float],
                                phi_positions: Optional[List[float]] = None) -> List[PotentialProfile]:
        """
        Extract multiple potential profiles for parallel computation.
        
        Args:
            potential_3d: 3D potential array
            r_positions: List of radial positions
            phi_positions: List of angular positions (default: all 0)
            
        Returns:
            List of PotentialProfile objects
        """
        if phi_positions is None:
            phi_positions = [0.0] * len(r_positions)
        
        profiles = []
        for r, phi in zip(r_positions, phi_positions):
            profile = self.extract_profile(potential_3d, r, phi)
            profiles.append(profile)
        
        return profiles
    
    def get_surface_potential_map(self, potential_3d: np.ndarray) -> np.ndarray:
        """
        Extract 2D surface potential map.
        
        Args:
            potential_3d: 3D potential array
            
        Returns:
            2D array of surface potentials (r, phi)
        """
        ns = self.grid.params.ns
        return potential_3d[:, ns-1, :]
    
    def get_tip_induced_band_bending(self, potential_3d: np.ndarray,
                                    r_cutoff: float = 10.0) -> np.ndarray:
        """
        Calculate tip-induced band bending vs radial position.
        
        Args:
            potential_3d: 3D potential array
            r_cutoff: Maximum radius to consider (nm)
            
        Returns:
            Array of band bending values vs radius
        """
        # Find cutoff index
        ir_max = np.argmin(np.abs(self.grid.r - r_cutoff))
        
        # Extract band bending at interface
        ns = self.grid.params.ns
        bb = potential_3d[:ir_max, ns-1, 0]  # At phi=0
        
        return bb
    
    def calculate_electric_field(self, potential_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate electric field components from potential.
        
        E = -∇φ
        
        Args:
            potential_3d: 3D potential array
            
        Returns:
            Tuple of (Er, Ez, Ephi) field components
        """
        # Gradient in cylindrical coordinates
        dr = self.grid.params.delr
        dz_vac = self.grid.params.delv
        dz_sem = self.grid.params.dels
        dphi = self.grid.params.delp
        
        # Allocate field arrays
        Er = np.zeros_like(potential_3d)
        Ez = np.zeros_like(potential_3d)
        Ephi = np.zeros_like(potential_3d)
        
        # Calculate gradients using central differences
        # Radial component
        Er[1:-1, :, :] = -(potential_3d[2:, :, :] - potential_3d[:-2, :, :]) / (2 * dr)
        
        # Vertical component (different spacing in vacuum/semiconductor)
        ns = self.grid.params.ns
        # Vacuum region
        Ez[:, ns:-1, :] = -(potential_3d[:, ns+1:, :] - potential_3d[:, ns-1:-2, :]) / (2 * dz_vac)
        # Semiconductor region  
        Ez[:, 1:ns-1, :] = -(potential_3d[:, 2:ns, :] - potential_3d[:, :ns-2, :]) / (2 * dz_sem)
        
        # Angular component
        for i in range(self.grid.params.nr):
            if self.grid.r[i] > 0:
                for j in range(potential_3d.shape[1]):
                    for k in range(self.grid.params.np):
                        k_next = (k + 1) % self.grid.params.np
                        k_prev = (k - 1) % self.grid.params.np
                        Ephi[i, j, k] = -(potential_3d[i, j, k_next] - 
                                         potential_3d[i, j, k_prev]) / (2 * dphi * self.grid.r[i])
        
        return Er, Ez, Ephi
    
    def calculate_field_magnitude(self, potential_3d: np.ndarray) -> np.ndarray:
        """
        Calculate electric field magnitude.
        
        Args:
            potential_3d: 3D potential array
            
        Returns:
            Field magnitude array
        """
        Er, Ez, Ephi = self.calculate_electric_field(potential_3d)
        return np.sqrt(Er**2 + Ez**2 + Ephi**2)


def create_potential_processor(grid: Grid3D) -> PotentialProcessor:
    """
    Create a potential processor for the given grid.
    
    Args:
        grid: 3D computational grid
        
    Returns:
        PotentialProcessor object
    """
    return PotentialProcessor(grid)