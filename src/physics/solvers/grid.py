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
        
        # Note: For non-uniform grids (tangent-based), delr and dels are input parameters
        # not actual spacing, so the linear consistency check doesn't apply
        # Only check vacuum grid which is still linear
        assert abs(self.vmax - (self.nv - 1) * self.delv) < 1e-6, "Inconsistent vacuum grid"


class Grid3D:
    """
    3D computational grid for SEMITIP simulations.
    
    Uses cylindrical coordinates (r, z, phi) with special treatment
    for vacuum and semiconductor regions.
    """
    
    def __init__(self, params: GridParameters, tip_model=None):
        """
        Initialize the grid.
        
        Args:
            params: Grid parameters
            tip_model: TipModel for hyperboloidal coordinate calculation
        """
        self.params = params
        self._setup_coordinates(tip_model)
        self._allocate_arrays()
    
    def _setup_coordinates(self, tip_model=None):
        """Set up coordinate arrays using Fortran SEMITIP3 formulas."""
        # Radial coordinates using tangent formula from SEMITIP3-6.1.f line 116
        # R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
        self.r = np.zeros(self.params.nr)
        for i in range(self.params.nr):
            # Convert to 1-based indexing like Fortran
            i_fortran = i + 1
            self.r[i] = (2 * self.params.nr * self.params.delr / np.pi) * \
                        np.tan(np.pi * (i_fortran - 0.5) / (2.0 * self.params.nr))
        
        # Vacuum z-coordinates using hyperboloidal coordinates
        # Following Fortran SEMITIP3-6.1.f lines 127-142
        if tip_model is not None:
            eta_tip, a, z0, c = tip_model.hyperboloid_parameters()
            deleta = eta_tip / float(self.params.nv)  # DELETA = ETAT/FLOAT(NV)
            
            # Initialize vacuum coordinates
            self.zv = np.zeros(self.params.nv)
            self.delv_array = np.zeros(self.params.nr)  # DELV varies with radius
            
            for j in range(self.params.nv):
                j_fortran = j + 1
                eta = j_fortran * deleta  # ETA = J*DELETA
                
                # Calculate z coordinates using hyperboloidal transformation
                # For each radial position, we need different z values
                # Store z for the central axis (r=0) as representative
                xsi = 1.0  # XSI = 1 at r=0 (central axis)
                self.zv[j] = a * eta * (xsi + c)  # Z = A*ETA*(XSI+C)
            
            # Calculate DELV array for each radial position (Fortran line 127)
            for i in range(self.params.nr):
                # DELV(I)=(SQRT(A**2+R(I)**2)+C*A)*ETAT/FLOAT(NV)
                self.delv_array[i] = (np.sqrt(a**2 + self.r[i]**2) + c*a) * eta_tip / float(self.params.nv)
            
            # Store DELV for central axis (r=0) for output, matching Fortran
            delv_center = (np.sqrt(a**2) + c*a) * eta_tip / float(self.params.nv)
            self.params.delv = delv_center
        else:
            # Fallback to linear spacing if no tip model
            self.zv = np.linspace(0, self.params.vmax, self.params.nv)
            self.delv_array = np.full(self.params.nr, self.params.delv)
        
        # Semiconductor z-coordinates using tangent formula from SEMITIP3-6.1.f line 164
        # S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
        self.zs_positive = np.zeros(self.params.ns)
        for j in range(self.params.ns):
            # Convert to 1-based indexing like Fortran
            j_fortran = j + 1
            self.zs_positive[j] = (2 * self.params.ns * self.params.dels / np.pi) * \
                                  np.tan(np.pi * (j_fortran - 0.5) / (2.0 * self.params.ns))
        
        # Make semiconductor coordinates negative (below surface)
        self.zs = -self.zs_positive
        
        # Angular coordinates - uniform spacing
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
    Create a Grid3D from configuration data using Fortran SEMITIP3 grid logic.
    
    Args:
        config: Configuration object
        tip_model: TipModel for setting tip mask
        
    Returns:
        Grid3D object
    """
    # Extract grid parameters
    grid_config = config.grid
    
    # Calculate grid spacing following Fortran MultInt3 logic (lines 384-400)
    # This determines the input spacing parameters DELR0, DELS0 for the tangent formulas
    
    # Get depletion width estimate for grid sizing
    # For now, use simple estimate - full implementation would get from simulation
    region = config.semiconductor_regions[0]
    eps = region.permittivity * 8.854e-12  # Approximate EPSILON0
    doping = abs(region.donor_concentration - region.acceptor_concentration) * 1e6  # to m^-3
    bias_est = 2.0  # Typical bias voltage for grid sizing
    if doping > 0:
        w_depl = np.sqrt(2 * eps * bias_est / (1.602e-19 * doping)) * 1e9  # to nm
    else:
        w_depl = 100.0  # Default
    
    # SIZE parameter logic from Fortran - should be 0.5 according to fort.9
    size_factor = getattr(grid_config, 'size_factor', 0.5)  # Default SIZE=0.5
    
    if size_factor > 0:
        # Adaptive sizing based on tip geometry and depletion width
        delr_base = tip_model.radius
        if tip_model.protrusion_radius > 0:
            delr_base = min(tip_model.protrusion_radius, delr_base)
        delr_base = min(delr_base, w_depl / grid_config.radial_points)
        delr = delr_base * size_factor
        
        dels_base = tip_model.radius  
        if tip_model.protrusion_radius > 0:
            dels_base = min(tip_model.protrusion_radius, dels_base)
        dels_base = min(dels_base, w_depl / grid_config.semiconductor_points)
        dels = dels_base * size_factor
    else:
        # Use fixed values from config (if available)
        delr = getattr(grid_config, 'radial_spacing', 0.5)  # DELRIN
        dels = getattr(grid_config, 'semiconductor_spacing', 0.5)  # DELSIN
    
    # Vacuum spacing - simplified for now, should implement hyperboloidal
    delv = grid_config.vacuum_extent / (grid_config.vacuum_points - 1)
    
    # Angular spacing (always uniform)
    if config.mirror_symmetry:
        delp = np.pi / grid_config.angular_points
    else:
        delp = 2 * np.pi / grid_config.angular_points
    
    # Update rmax and smax based on calculated spacing and tangent formula
    # For tangent formula: max = (2*N*del0/π) * tan(π/2) → ∞
    # But in practice, we use the extent from config as upper limit
    rmax = grid_config.radial_extent
    smax = grid_config.semiconductor_extent
    
    # Create parameters
    params = GridParameters(
        nr=grid_config.radial_points,
        nv=grid_config.vacuum_points,
        ns=grid_config.semiconductor_points,
        np=grid_config.angular_points,
        delr=delr,  # Input spacing parameter (DELR0)
        delv=delv,
        dels=dels,  # Input spacing parameter (DELS0)
        delp=delp,
        rmax=rmax,
        vmax=grid_config.vacuum_extent,
        smax=smax,
        mirror_symmetry=config.mirror_symmetry
    )
    
    # Create grid with tip model for hyperboloidal coordinates
    grid = Grid3D(params, tip_model)
    
    # Set tip mask
    grid.set_tip_mask(tip_model)
    
    return grid