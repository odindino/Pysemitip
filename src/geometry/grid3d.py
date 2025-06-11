"""
3D Cylindrical Grid Management Module

This module implements the 3D cylindrical coordinate grid system for STM
simulations, including adaptive grid generation, grid point management,
and coordinate transformations.

Based on SEMITIP semitip3-6.1.f implementation.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Union
from dataclasses import dataclass
import warnings
from enum import Enum


class GridRegion(Enum):
    """Grid region identifiers"""
    VACUUM = "vacuum"
    SEMICONDUCTOR = "semiconductor" 
    SURFACE = "surface"


@dataclass
class GridDimensions:
    """Maximum grid dimensions (corresponding to Fortran PARAMETER values)"""
    MAX_R_POINTS: int = 512    # NRDIM: Maximum radial grid points
    MAX_V_POINTS: int = 64     # NVDIM: Maximum vacuum axial grid points  
    MAX_S_POINTS: int = 512    # NSDIM: Maximum semiconductor axial grid points
    MAX_P_POINTS: int = 64     # NPDIM: Maximum angular grid points


@dataclass 
class GridConfig:
    """Configuration for 3D cylindrical grid generation"""
    
    # Grid point numbers (actual usage, must be <= GridDimensions.MAX_*)
    nr_points: int = 64        # Number of radial grid points
    nv_points: int = 16        # Number of vacuum axial grid points
    ns_points: int = 64        # Number of semiconductor axial grid points
    np_points: int = 16        # Number of angular grid points
    
    # Base grid spacing parameters
    delr0: float = 0.1         # Base radial grid spacing [nm]
    dels0: float = 0.1         # Base semiconductor axial grid spacing [nm]
    
    # Domain extents
    max_radius: float = 20.0   # Maximum radial extent [nm]
    vacuum_thickness: float = 5.0      # Vacuum region thickness [nm]
    semiconductor_depth: float = 20.0  # Semiconductor region depth [nm]
    
    # Mirror symmetry flag
    mirror_symmetry: bool = False      # Use mirror symmetry (φ: 0→π vs 0→2π)
    
    # Refinement parameters
    max_refinement_levels: int = 3     # Maximum number of refinement levels
    
    def __post_init__(self):
        """Validate grid configuration"""
        dims = GridDimensions()
        
        if self.nr_points > dims.MAX_R_POINTS:
            raise ValueError(f"nr_points ({self.nr_points}) exceeds MAX_R_POINTS ({dims.MAX_R_POINTS})")
        if self.nv_points > dims.MAX_V_POINTS:
            raise ValueError(f"nv_points ({self.nv_points}) exceeds MAX_V_POINTS ({dims.MAX_V_POINTS})")
        if self.ns_points > dims.MAX_S_POINTS:
            raise ValueError(f"ns_points ({self.ns_points}) exceeds MAX_S_POINTS ({dims.MAX_S_POINTS})")
        if self.np_points > dims.MAX_P_POINTS:
            raise ValueError(f"np_points ({self.np_points}) exceeds MAX_P_POINTS ({dims.MAX_P_POINTS})")
            
        if self.delr0 <= 0 or self.dels0 <= 0:
            raise ValueError("Grid spacing parameters must be positive")


class Grid3D:
    """
    3D Cylindrical Coordinate Grid Manager
    
    This class implements the SEMITIP 3D cylindrical grid system with:
    - Tangent-function radial distribution for high density near axis
    - Separate vacuum and semiconductor axial grids
    - Angular grid with optional mirror symmetry
    - Support for adaptive refinement
    
    Coordinate system:
    - r: radial coordinate [0, r_max]
    - φ: angular coordinate [0, 2π] or [0, π] (mirror symmetry)
    - z: axial coordinate [z_min, z_max] where z=0 is surface
    """
    
    def __init__(self, config: Optional[GridConfig] = None):
        """Initialize 3D cylindrical grid"""
        self.config = config or GridConfig()
        self.dimensions = GridDimensions()
        
        # Current refinement level (0 = coarsest)
        self.refinement_level = 0
        
        # Grid arrays - will be initialized by generate_grid()
        self.r_points = None      # Radial grid points [nm]
        self.z_vacuum_points = None    # Vacuum z grid points [nm] 
        self.z_semiconductor_points = None  # Semiconductor z grid points [nm]
        self.phi_points = None    # Angular grid points [radians]
        
        # Grid spacing arrays
        self.delr = None         # Radial spacing [nm]
        self.delv = None         # Vacuum axial spacing [nm]
        self.dels = None         # Semiconductor axial spacing [nm]
        self.delp = None         # Angular spacing [radians]
        
        # Generate initial grid
        self.generate_grid()
        
    def generate_grid(self):
        """Generate all grid points and spacing arrays"""
        self._generate_radial_grid()
        self._generate_axial_grids()
        self._generate_angular_grid()
        self._calculate_grid_spacing()
        
    def _generate_radial_grid(self):
        """
        Generate radial grid points using tangent mapping
        
        Based on semitip3-6.1.f line 116:
        R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
        """
        nr = self.config.nr_points
        delr0 = self.config.delr0
        
        # Pre-allocate array
        self.r_points = np.zeros(nr, dtype=np.float64)
        
        for i in range(nr):
            # Convert to 1-based indexing for Fortran compatibility
            i_fortran = i + 1
            
            # Tangent mapping for dense distribution near r=0
            self.r_points[i] = (2.0 * nr * delr0 / np.pi) * np.tan(
                np.pi * (i_fortran - 0.5) / (2.0 * nr)
            )
            
        # Note: Do NOT scale radial points to maintain Fortran compatibility
        # The max_radius is just a configuration parameter
        # Actual grid extent is determined by the tangent mapping
            
    def _generate_axial_grids(self):
        """Generate axial grid points for vacuum and semiconductor regions"""
        
        # Vacuum region: z from 0 to -vacuum_thickness
        nv = self.config.nv_points
        vacuum_thickness = self.config.vacuum_thickness
        
        self.z_vacuum_points = np.zeros(nv, dtype=np.float64)
        for j in range(nv):
            # Uniform spacing in vacuum (can be made adaptive later)
            self.z_vacuum_points[j] = -(j + 1) * vacuum_thickness / nv
            
        # Semiconductor region: z from 0 to +semiconductor_depth
        # Using tangent mapping like radial grid for dense distribution near surface
        ns = self.config.ns_points  
        dels0 = self.config.dels0
        
        self.z_semiconductor_points = np.zeros(ns, dtype=np.float64)
        for j in range(ns):
            # Convert to 1-based indexing
            j_fortran = j + 1
            
            # Tangent mapping (semitip3-6.1.f line 164)
            self.z_semiconductor_points[j] = (2.0 * ns * dels0 / np.pi) * np.tan(
                np.pi * (j_fortran - 0.5) / (2.0 * ns)
            )
            
        # Note: Do NOT scale semiconductor points to maintain Fortran compatibility
        # The semiconductor_depth is just a configuration parameter
        # Actual grid extent is determined by the tangent mapping
            
    def _generate_angular_grid(self):
        """Generate angular grid points"""
        np_points = self.config.np_points
        
        if self.config.mirror_symmetry:
            # Mirror symmetry: φ from 0 to π (semitip3-6.1.f line 397-400)
            self.delp = np.pi / np_points
            phi_range = np.pi
        else:
            # Full circle: φ from 0 to 2π
            self.delp = 2.0 * np.pi / np_points
            phi_range = 2.0 * np.pi
            
        self.phi_points = np.linspace(0, phi_range, np_points, endpoint=False, dtype=np.float64)
        
    def _calculate_grid_spacing(self):
        """Calculate grid spacing arrays"""
        
        # Radial spacing
        nr = len(self.r_points)
        self.delr = np.zeros(nr, dtype=np.float64)
        
        for i in range(nr):
            if i == 0:
                # First point spacing
                self.delr[i] = self.r_points[0] if nr == 1 else (self.r_points[1] - self.r_points[0]) / 2.0
            elif i == nr - 1:
                # Last point spacing  
                self.delr[i] = (self.r_points[i] - self.r_points[i-1]) / 2.0
            else:
                # Central difference
                self.delr[i] = (self.r_points[i+1] - self.r_points[i-1]) / 2.0
                
        # Vacuum axial spacing
        nv = len(self.z_vacuum_points)
        self.delv = np.zeros(nv, dtype=np.float64)
        
        for j in range(nv):
            if j == 0:
                # Spacing at surface (can be radially dependent)
                self.delv[j] = abs(self.z_vacuum_points[0]) if nv == 1 else abs(self.z_vacuum_points[1] - self.z_vacuum_points[0]) / 2.0
            elif j == nv - 1:
                self.delv[j] = abs(self.z_vacuum_points[j] - self.z_vacuum_points[j-1]) / 2.0
            else:
                self.delv[j] = abs(self.z_vacuum_points[j+1] - self.z_vacuum_points[j-1]) / 2.0
                
        # Semiconductor axial spacing
        ns = len(self.z_semiconductor_points)
        self.dels = np.zeros(ns, dtype=np.float64)
        
        for j in range(ns):
            if j == 0:
                self.dels[j] = self.z_semiconductor_points[0] if ns == 1 else (self.z_semiconductor_points[1] - self.z_semiconductor_points[0]) / 2.0
            elif j == ns - 1:
                self.dels[j] = (self.z_semiconductor_points[j] - self.z_semiconductor_points[j-1]) / 2.0
            else:
                self.dels[j] = (self.z_semiconductor_points[j+1] - self.z_semiconductor_points[j-1]) / 2.0
                
    def refine_grid(self):
        """
        Refine grid by factor of 2 in all dimensions
        
        Based on semitip3-6.1.f lines 227-240:
        NR=NR*2, NS=NS*2, NV=NV*2, NP=NP*2
        DELR0=DELR0/2., DELS0=DELS0/2.
        """
        if self.refinement_level >= self.config.max_refinement_levels:
            warnings.warn(f"Maximum refinement level ({self.config.max_refinement_levels}) reached")
            return False
            
        # Check if refinement would exceed maximum dimensions
        new_nr = self.config.nr_points * 2
        new_nv = self.config.nv_points * 2
        new_ns = self.config.ns_points * 2
        new_np = self.config.np_points * 2
        
        if (new_nr > self.dimensions.MAX_R_POINTS or 
            new_nv > self.dimensions.MAX_V_POINTS or
            new_ns > self.dimensions.MAX_S_POINTS or
            new_np > self.dimensions.MAX_P_POINTS):
            warnings.warn("Grid refinement would exceed maximum dimensions")
            return False
            
        # Apply refinement
        self.config.nr_points = new_nr
        self.config.nv_points = new_nv
        self.config.ns_points = new_ns
        self.config.np_points = new_np
        
        # Halve base grid spacing
        self.config.delr0 /= 2.0
        self.config.dels0 /= 2.0
        
        # Increment refinement level
        self.refinement_level += 1
        
        # Regenerate grid
        self.generate_grid()
        
        return True
        
    def get_grid_point(self, region: GridRegion, i: int, j: int, k: int = 0) -> Tuple[float, float, float]:
        """
        Get Cartesian coordinates of grid point
        
        Args:
            region: Grid region (VACUUM, SEMICONDUCTOR, SURFACE)
            i: Radial index
            j: Axial index  
            k: Angular index
            
        Returns:
            (x, y, z) coordinates in nm
        """
        if i >= len(self.r_points) or k >= len(self.phi_points):
            raise IndexError("Grid indices out of bounds")
            
        r = self.r_points[i]
        phi = self.phi_points[k]
        
        if region == GridRegion.VACUUM:
            if j >= len(self.z_vacuum_points):
                raise IndexError("Vacuum axial index out of bounds")
            z = self.z_vacuum_points[j]
        elif region == GridRegion.SEMICONDUCTOR:
            if j >= len(self.z_semiconductor_points):
                raise IndexError("Semiconductor axial index out of bounds")
            z = self.z_semiconductor_points[j]
        elif region == GridRegion.SURFACE:
            z = 0.0  # Surface is at z=0
        else:
            raise ValueError(f"Unknown grid region: {region}")
            
        # Convert cylindrical to Cartesian
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        
        return (x, y, z)
        
    def get_cylindrical_coordinates(self, region: GridRegion, i: int, j: int, k: int = 0) -> Tuple[float, float, float]:
        """
        Get cylindrical coordinates of grid point
        
        Returns:
            (r, phi, z) coordinates
        """
        if i >= len(self.r_points) or k >= len(self.phi_points):
            raise IndexError("Grid indices out of bounds")
            
        r = self.r_points[i]
        phi = self.phi_points[k]
        
        if region == GridRegion.VACUUM:
            if j >= len(self.z_vacuum_points):
                raise IndexError("Vacuum axial index out of bounds")
            z = self.z_vacuum_points[j]
        elif region == GridRegion.SEMICONDUCTOR:
            if j >= len(self.z_semiconductor_points):
                raise IndexError("Semiconductor axial index out of bounds")
            z = self.z_semiconductor_points[j]
        elif region == GridRegion.SURFACE:
            z = 0.0
        else:
            raise ValueError(f"Unknown grid region: {region}")
            
        return (r, phi, z)
        
    def find_nearest_grid_point(self, x: float, y: float, z: float) -> Tuple[GridRegion, int, int, int]:
        """
        Find nearest grid point to given Cartesian coordinates
        
        Returns:
            (region, i, j, k): Grid region and indices
        """
        # Convert to cylindrical
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        if phi < 0:
            phi += 2.0 * np.pi
            
        # Find nearest radial index
        i = np.argmin(np.abs(self.r_points - r))
        
        # Find nearest angular index
        k = np.argmin(np.abs(self.phi_points - phi))
        
        # Determine region and find axial index
        if z < 0:
            # Vacuum region
            region = GridRegion.VACUUM
            j = np.argmin(np.abs(self.z_vacuum_points - z))
        else:
            # Semiconductor region  
            region = GridRegion.SEMICONDUCTOR
            j = np.argmin(np.abs(self.z_semiconductor_points - z))
            
        return (region, i, j, k)
        
    def get_grid_volume_element(self, i: int, j: int, k: int, region: GridRegion) -> float:
        """
        Calculate volume element at grid point
        
        Volume element in cylindrical coordinates: dV = r * dr * dφ * dz
        """
        if i >= len(self.r_points):
            raise IndexError("Radial index out of bounds")
            
        r = self.r_points[i]
        dr = self.delr[i]
        dphi = self.delp
        
        if region == GridRegion.VACUUM:
            if j >= len(self.delv):
                raise IndexError("Vacuum axial index out of bounds")
            dz = self.delv[j]
        elif region == GridRegion.SEMICONDUCTOR:
            if j >= len(self.dels):
                raise IndexError("Semiconductor axial index out of bounds") 
            dz = self.dels[j]
        else:
            raise ValueError(f"Cannot calculate volume for region: {region}")
            
        return r * dr * dphi * dz
        
    def get_grid_info(self) -> Dict:
        """Get comprehensive grid information"""
        return {
            'dimensions': {
                'nr_points': self.config.nr_points,
                'nv_points': self.config.nv_points,
                'ns_points': self.config.ns_points,
                'np_points': self.config.np_points
            },
            'domain_extents': {
                'r_min': 0.0,
                'r_max': self.r_points[-1] if self.r_points is not None else 0.0,
                'z_vacuum_min': self.z_vacuum_points[-1] if self.z_vacuum_points is not None else 0.0,
                'z_semiconductor_max': self.z_semiconductor_points[-1] if self.z_semiconductor_points is not None else 0.0,
                'phi_min': 0.0,
                'phi_max': self.phi_points[-1] + self.delp if self.phi_points is not None else 0.0
            },
            'refinement': {
                'current_level': self.refinement_level,
                'max_levels': self.config.max_refinement_levels
            },
            'symmetry': {
                'mirror_symmetry': self.config.mirror_symmetry
            },
            'total_points': self.get_total_grid_points()
        }
        
    def get_total_grid_points(self) -> int:
        """Calculate total number of grid points"""
        total_vacuum = self.config.nr_points * self.config.nv_points * self.config.np_points
        total_semiconductor = self.config.nr_points * self.config.ns_points * self.config.np_points
        total_surface = self.config.nr_points * self.config.np_points
        
        return total_vacuum + total_semiconductor + total_surface
        
    def estimate_memory_usage(self) -> Dict[str, float]:
        """Estimate memory usage for potential arrays"""
        total_points = self.get_total_grid_points()
        
        # Assume double precision (8 bytes per float)
        bytes_per_point = 8
        
        # Multiple arrays (potential, charge density, etc.)
        vacuum_points = self.config.nr_points * self.config.nv_points * self.config.np_points
        semiconductor_points = self.config.nr_points * self.config.ns_points * self.config.np_points
        surface_points = self.config.nr_points * self.config.np_points
        
        return {
            'vacuum_array_mb': vacuum_points * bytes_per_point / 1024**2,
            'semiconductor_array_mb': semiconductor_points * bytes_per_point / 1024**2,
            'surface_array_mb': surface_points * bytes_per_point / 1024**2,
            'total_per_array_mb': total_points * bytes_per_point / 1024**2,
            'estimated_total_mb': total_points * bytes_per_point * 5 / 1024**2  # Multiple arrays
        }
        
    def __str__(self) -> str:
        """String representation"""
        return (f"Grid3D(nr={self.config.nr_points}, nv={self.config.nv_points}, "
               f"ns={self.config.ns_points}, np={self.config.np_points}, "
               f"level={self.refinement_level})")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions for common grid configurations

def create_coarse_grid(mirror_symmetry: bool = False) -> Grid3D:
    """Create coarse grid for initial calculations"""
    config = GridConfig(
        nr_points=32,
        nv_points=8,
        ns_points=32,
        np_points=8,
        mirror_symmetry=mirror_symmetry
    )
    return Grid3D(config)


def create_medium_grid(mirror_symmetry: bool = False) -> Grid3D:
    """Create medium resolution grid"""
    config = GridConfig(
        nr_points=64,
        nv_points=16,
        ns_points=64,
        np_points=16,
        mirror_symmetry=mirror_symmetry
    )
    return Grid3D(config)


def create_fine_grid(mirror_symmetry: bool = False) -> Grid3D:
    """Create high resolution grid"""
    config = GridConfig(
        nr_points=128,
        nv_points=32,
        ns_points=128,
        np_points=32,
        mirror_symmetry=mirror_symmetry
    )
    return Grid3D(config)


if __name__ == "__main__":
    # Demo usage
    print("3D Cylindrical Grid Demo")
    print("=" * 40)
    
    # Create medium grid
    grid = create_medium_grid(mirror_symmetry=True)
    print(f"Created grid: {grid}")
    
    # Get grid info
    info = grid.get_grid_info()
    print(f"\nGrid Information:")
    for key, value in info.items():
        print(f"  {key}: {value}")
        
    # Memory usage
    memory = grid.estimate_memory_usage()
    print(f"\nMemory Usage Estimate:")
    for key, value in memory.items():
        print(f"  {key}: {value:.2f}")
        
    # Test grid point access
    print(f"\nSample grid points:")
    r, phi, z = grid.get_cylindrical_coordinates(GridRegion.VACUUM, 0, 0, 0)
    print(f"  Vacuum (0,0,0): r={r:.3f}, φ={phi:.3f}, z={z:.3f}")
    
    r, phi, z = grid.get_cylindrical_coordinates(GridRegion.SEMICONDUCTOR, 10, 10, 0)
    print(f"  Semiconductor (10,10,0): r={r:.3f}, φ={phi:.3f}, z={z:.3f}")
    
    # Test refinement
    print(f"\nTesting grid refinement...")
    print(f"Before refinement: {grid}")
    success = grid.refine_grid()
    print(f"Refinement successful: {success}")
    print(f"After refinement: {grid}")