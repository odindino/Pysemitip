"""
Grid system and coordinate transformations for SEMITIP simulations.

This module implements the grid structure used in the Fortran code,
including cylindrical coordinates and hyperboloidal transformations.
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from ...utils.constants import PhysicalConstants as PC


@dataclass
class GridParameters:
    """Parameters defining the computational grid."""
    # Grid dimensions
    nr: int          # Number of radial points
    nv: int          # Number of vacuum points
    ns: int          # Number of semiconductor points  
    np: int          # Number of angular points
    
    # Grid spacing
    delr: float      # Radial spacing (nm)
    delv: float      # Vacuum spacing (nm)
    dels: float      # Semiconductor spacing (nm)
    delp: float      # Angular spacing (radians)
    
    # Grid extents
    rmax: float      # Maximum radius (nm)
    vmax: float      # Maximum vacuum extent (nm)
    smax: float      # Maximum semiconductor depth (nm)
    
    # Mirror symmetry
    mirror_symmetry: bool = True  # MIRROR parameter
    
    def __post_init__(self):
        """Calculate derived quantities."""
        # Angular extent
        if self.mirror_symmetry:
            self.pmax = np.pi  # 180 degrees with mirror symmetry
        else:
            self.pmax = 2 * np.pi  # Full 360 degrees
        
        # Ensure consistency
        self.validate()
    
    def validate(self):
        """Validate grid parameters."""
        # Check dimensions
        assert self.nr > 0, "nr must be positive"
        assert self.nv > 0, "nv must be positive"
        assert self.ns > 0, "ns must be positive"
        assert self.np > 0, "np must be positive"
        
        # Check spacing
        assert self.delr > 0, "delr must be positive"
        assert self.delv > 0, "delv must be positive"
        assert self.dels > 0, "dels must be positive"
        assert self.delp > 0, "delp must be positive"
        
        # Check consistency
        assert abs(self.rmax - (self.nr - 1) * self.delr) < 1e-6, "Inconsistent radial grid"
        assert abs(self.vmax - (self.nv - 1) * self.delv) < 1e-6, "Inconsistent vacuum grid"
        assert abs(self.smax - (self.ns - 1) * self.dels) < 1e-6, "Inconsistent semiconductor grid"


class Grid3D:
    """
    3D computational grid for SEMITIP simulations.
    
    Uses cylindrical coordinates (r, z, phi) with special treatment
    for vacuum and semiconductor regions.
    """
    
    def __init__(self, params: GridParameters):
        """
        Initialize the grid.
        
        Args:
            params: Grid parameters
        """
        self.params = params
        self._setup_coordinates()
        self._allocate_arrays()
    
    def _setup_coordinates(self):
        """Set up coordinate arrays."""
        # Radial coordinates
        self.r = np.linspace(0, self.params.rmax, self.params.nr)
        
        # Vacuum z-coordinates (positive z, above surface)
        self.zv = np.linspace(0, self.params.vmax, self.params.nv)
        
        # Semiconductor z-coordinates (negative z, below surface)
        self.zs = -np.linspace(0, self.params.smax, self.params.ns)
        
        # Angular coordinates
        self.phi = np.linspace(0, self.params.pmax, self.params.np)
        
        # Combined z array for convenience
        self.z = np.concatenate([np.flip(self.zs[1:]), self.zv])
        self.nz = len(self.z)
    
    def _allocate_arrays(self):
        """Allocate main data arrays."""
        # Vacuum region array
        self.vac = np.zeros((2, self.params.nr, self.params.nv, self.params.np))
        
        # Semiconductor region array  
        self.sem = np.zeros((2, self.params.nr, self.params.ns, self.params.np))
        
        # Interface array (vacuum-semiconductor boundary)
        self.vsint = np.zeros((2, self.params.nr, self.params.np))
        
        # Tip mask (boolean array indicating tip interior)
        self.tip_mask = np.zeros((self.params.nr, self.params.nv, self.params.np), dtype=bool)
    
    def get_meshgrid(self, region: str = 'all') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get meshgrid for specified region.
        
        Args:
            region: 'vacuum', 'semiconductor', or 'all'
            
        Returns:
            Tuple of (R, Z, PHI) meshgrid arrays
        """
        if region == 'vacuum':
            return np.meshgrid(self.r, self.zv, self.phi, indexing='ij')
        elif region == 'semiconductor':
            return np.meshgrid(self.r, self.zs, self.phi, indexing='ij')
        else:  # 'all'
            return np.meshgrid(self.r, self.z, self.phi, indexing='ij')
    
    def cylindrical_to_cartesian(self, r: np.ndarray, phi: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Convert cylindrical to Cartesian coordinates.
        
        Args:
            r: Radial coordinates
            phi: Angular coordinates
            z: Vertical coordinates
            
        Returns:
            Tuple of (x, y, z) Cartesian coordinates
        """
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        return x, y, z
    
    def hyperboloidal_transform(self, r: np.ndarray, z: np.ndarray, 
                              eta: float, a: float, z0: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transform to hyperboloidal coordinates for tip region.
        
        This implements the coordinate transformation used in SEMITIP3.
        
        Args:
            r: Radial coordinates
            z: Vertical coordinates  
            eta: Hyperboloid parameter
            a: Hyperboloid parameter
            z0: Hyperboloid parameter
            
        Returns:
            Tuple of (xi, zeta) hyperboloidal coordinates
        """
        # Shift z coordinate
        z_shifted = z - z0
        
        # Calculate hyperboloidal coordinates
        rho_sq = r**2 + z_shifted**2
        xi = np.sqrt((rho_sq + a**2 + np.sqrt((rho_sq + a**2)**2 - 4 * a**2 * z_shifted**2)) / 2)
        zeta = np.sqrt((rho_sq - a**2 + np.sqrt((rho_sq - a**2)**2 + 4 * a**2 * r**2)) / 2)
        
        return xi, zeta
    
    def get_grid_indices(self, r: float, z: float) -> Tuple[int, int]:
        """
        Get grid indices for a given (r,z) position.
        
        Args:
            r: Radial position (nm)
            z: Vertical position (nm)
            
        Returns:
            Tuple of (ir, iz) grid indices
        """
        # Radial index
        ir = int(round(r / self.params.delr))
        ir = max(0, min(ir, self.params.nr - 1))
        
        # Vertical index
        if z >= 0:
            # Vacuum region
            iz = int(round(z / self.params.delv))
            iz = max(0, min(iz, self.params.nv - 1))
            is_vacuum = True
        else:
            # Semiconductor region
            iz = int(round(-z / self.params.dels))
            iz = max(0, min(iz, self.params.ns - 1))
            is_vacuum = False
        
        return ir, iz, is_vacuum
    
    def interpolate_to_position(self, r: float, z: float, phi: float, 
                               field: np.ndarray) -> float:
        """
        Interpolate field value to arbitrary position.
        
        Args:
            r: Radial position (nm)
            z: Vertical position (nm)
            phi: Angular position (radians)
            field: Field array to interpolate
            
        Returns:
            Interpolated field value
        """
        # Get grid indices and weights
        ir = r / self.params.delr
        ir0 = int(ir)
        ir1 = min(ir0 + 1, self.params.nr - 1)
        wr = ir - ir0
        
        # Angular interpolation
        ip = phi / self.params.delp
        ip0 = int(ip) % self.params.np
        ip1 = (ip0 + 1) % self.params.np
        wp = ip - int(ip)
        
        # Vertical interpolation depends on region
        if z >= 0:
            # Vacuum
            iz = z / self.params.delv
            iz0 = int(iz)
            iz1 = min(iz0 + 1, self.params.nv - 1)
            wz = iz - iz0
            
            # Trilinear interpolation
            f000 = field[ir0, iz0, ip0]
            f001 = field[ir0, iz0, ip1]
            f010 = field[ir0, iz1, ip0]
            f011 = field[ir0, iz1, ip1]
            f100 = field[ir1, iz0, ip0]
            f101 = field[ir1, iz0, ip1]
            f110 = field[ir1, iz1, ip0]
            f111 = field[ir1, iz1, ip1]
        else:
            # Semiconductor - similar logic
            iz = -z / self.params.dels
            iz0 = int(iz)
            iz1 = min(iz0 + 1, self.params.ns - 1)
            wz = iz - iz0
            
            # Get values from semiconductor array
            # Implementation depends on field structure
            raise NotImplementedError("Semiconductor interpolation to be implemented")
        
        # Trilinear interpolation formula
        f00 = f000 * (1 - wp) + f001 * wp
        f01 = f010 * (1 - wp) + f011 * wp
        f10 = f100 * (1 - wp) + f101 * wp
        f11 = f110 * (1 - wp) + f111 * wp
        
        f0 = f00 * (1 - wz) + f01 * wz
        f1 = f10 * (1 - wz) + f11 * wz
        
        return f0 * (1 - wr) + f1 * wr
    
    def set_tip_mask(self, tip_model):
        """
        Set the tip mask based on tip geometry.
        
        Args:
            tip_model: TipModel object defining tip geometry
        """
        # Get meshgrid for vacuum region
        R, Z, PHI = self.get_meshgrid('vacuum')
        
        # Check each point
        for i in range(self.params.nr):
            for j in range(self.params.nv):
                for k in range(self.params.np):
                    # Convert to Cartesian for tip position check
                    x, y, _ = self.cylindrical_to_cartesian(R[i,j,k], PHI[i,j,k], Z[i,j,k])
                    
                    # Account for tip lateral position
                    r_eff = np.sqrt((x - tip_model.position[0])**2 + 
                                   (y - tip_model.position[1])**2)
                    
                    # Check if inside tip
                    self.tip_mask[i,j,k] = tip_model.is_inside_tip(r_eff, Z[i,j,k])


def create_grid_from_config(config, tip_model) -> Grid3D:
    """
    Create a Grid3D from configuration data.
    
    Args:
        config: Configuration object
        tip_model: TipModel for setting tip mask
        
    Returns:
        Grid3D object
    """
    # Extract grid parameters
    grid_config = config.grid
    
    # Calculate grid spacing
    delr = grid_config.radial_extent / (grid_config.radial_points - 1)
    delv = grid_config.vacuum_extent / (grid_config.vacuum_points - 1)
    dels = grid_config.semiconductor_extent / (grid_config.semiconductor_points - 1)
    delp = np.pi / (grid_config.angular_points - 1) if config.mirror_symmetry else \
           2 * np.pi / (grid_config.angular_points - 1)
    
    # Create parameters
    params = GridParameters(
        nr=grid_config.radial_points,
        nv=grid_config.vacuum_points,
        ns=grid_config.semiconductor_points,
        np=grid_config.angular_points,
        delr=delr,
        delv=delv,
        dels=dels,
        delp=delp,
        rmax=grid_config.radial_extent,
        vmax=grid_config.vacuum_extent,
        smax=grid_config.semiconductor_extent,
        mirror_symmetry=config.mirror_symmetry
    )
    
    # Create grid
    grid = Grid3D(params)
    
    # Set tip mask
    grid.set_tip_mask(tip_model)
    
    return grid