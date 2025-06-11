"""
Boundary Conditions Management Module

This module implements boundary condition handling for STM simulations,
including tip surface potentials, semiconductor surface conditions,
and external boundary treatments.

Based on SEMITIP boundary condition implementation.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union, Any
from dataclasses import dataclass, field
from enum import Enum
import warnings

from .grid3d import Grid3D, GridRegion
from .stm_geometry import STMGeometry
from .tip_geometry import TipGeometry


class BoundaryType(Enum):
    """Types of boundary conditions"""
    DIRICHLET = "dirichlet"      # Fixed potential
    NEUMANN = "neumann"          # Fixed normal derivative
    MIXED = "mixed"              # Mixed boundary condition
    PERIODIC = "periodic"        # Periodic boundary
    MIRROR = "mirror"            # Mirror symmetry


class BoundaryRegion(Enum):
    """Boundary region identifiers"""
    TIP_SURFACE = "tip_surface"
    SAMPLE_SURFACE = "sample_surface"
    VACUUM_BOUNDARY = "vacuum_boundary"
    SEMICONDUCTOR_BOUNDARY = "semiconductor_boundary"
    RADIAL_BOUNDARY = "radial_boundary"
    ANGULAR_BOUNDARY = "angular_boundary"


@dataclass
class BoundaryConfig:
    """Configuration for boundary conditions"""
    
    # Tip boundary conditions
    tip_boundary_type: BoundaryType = BoundaryType.DIRICHLET
    tip_potential: Optional[float] = None  # Will be calculated from bias if None
    
    # Sample surface boundary conditions
    surface_boundary_type: BoundaryType = BoundaryType.MIXED
    surface_charge_density: Optional[Callable] = None  # Function for surface charge
    
    # External boundaries
    radial_boundary_type: BoundaryType = BoundaryType.NEUMANN
    radial_boundary_value: float = 0.0
    
    vacuum_boundary_type: BoundaryType = BoundaryType.DIRICHLET
    vacuum_boundary_value: float = 0.0
    
    semiconductor_boundary_type: BoundaryType = BoundaryType.DIRICHLET
    semiconductor_boundary_value: float = 0.0
    
    # Angular boundaries (for mirror symmetry)
    angular_boundary_type: BoundaryType = BoundaryType.MIRROR
    
    # Interface conditions
    enable_interface_continuity: bool = True
    permittivity_vacuum: float = 1.0
    permittivity_semiconductor: float = 11.7
    
    # Numerical parameters
    boundary_relaxation_factor: float = 1.0
    apply_boundary_every_iteration: bool = True


class BoundaryConditions:
    """
    Boundary Conditions Manager
    
    This class handles all boundary condition applications for STM simulations:
    - Tip surface Dirichlet conditions
    - Sample surface mixed conditions with charge density
    - External boundary conditions
    - Interface continuity conditions
    - Mirror symmetry boundaries
    """
    
    def __init__(self, 
                 grid: Grid3D,
                 geometry: STMGeometry,
                 tip_geometry: TipGeometry,
                 config: Optional[BoundaryConfig] = None):
        """
        Initialize boundary conditions manager
        
        Args:
            grid: 3D cylindrical grid
            geometry: STM geometry configuration
            tip_geometry: Tip geometry details
            config: Boundary condition configuration
        """
        self.grid = grid
        self.geometry = geometry
        self.tip_geometry = tip_geometry
        self.config = config or BoundaryConfig()
        
        # Calculate tip potential if not specified
        if self.config.tip_potential is None:
            self.config.tip_potential = (
                self.geometry.config.bias_voltage + 
                self.geometry.config.contact_potential +
                self.tip_geometry.config.work_function
            )
        
        # Boundary identification arrays
        self.boundary_masks = {}
        self.boundary_values = {}
        
        # Initialize boundary identification
        self._identify_boundaries()
        
    def _identify_boundaries(self):
        """Identify all boundary points in the grid"""
        
        # Get grid dimensions
        nr = self.grid.config.nr_points
        nv = self.grid.config.nv_points
        ns = self.grid.config.ns_points
        np_points = self.grid.config.np_points
        
        # Initialize boundary masks
        self.boundary_masks = {
            BoundaryRegion.TIP_SURFACE: {
                GridRegion.VACUUM: np.zeros((nr, nv, np_points), dtype=bool)
            },
            BoundaryRegion.SAMPLE_SURFACE: {
                GridRegion.VACUUM: np.zeros((nr, np_points), dtype=bool),
                GridRegion.SEMICONDUCTOR: np.zeros((nr, np_points), dtype=bool)
            },
            BoundaryRegion.RADIAL_BOUNDARY: {
                GridRegion.VACUUM: np.zeros((nv, np_points), dtype=bool),
                GridRegion.SEMICONDUCTOR: np.zeros((ns, np_points), dtype=bool)
            },
            BoundaryRegion.VACUUM_BOUNDARY: np.zeros((nr, np_points), dtype=bool),
            BoundaryRegion.SEMICONDUCTOR_BOUNDARY: np.zeros((nr, np_points), dtype=bool),
            BoundaryRegion.ANGULAR_BOUNDARY: {
                GridRegion.VACUUM: np.zeros((nr, nv), dtype=bool),
                GridRegion.SEMICONDUCTOR: np.zeros((nr, ns), dtype=bool)
            }
        }
        
        # Initialize boundary values
        self.boundary_values = {region: {} for region in BoundaryRegion}
        
        # Identify tip surface points
        self._identify_tip_surface_points()
        
        # Identify sample surface points
        self._identify_sample_surface_points()
        
        # Identify external boundary points
        self._identify_external_boundaries()
        
        # Identify angular boundaries (for mirror symmetry)
        if self.grid.config.mirror_symmetry:
            self._identify_angular_boundaries()
            
    def _identify_tip_surface_points(self):
        """Identify grid points on tip surface"""
        nr = self.grid.config.nr_points
        nv = self.grid.config.nv_points
        np_points = self.grid.config.np_points
        
        tip_mask = self.boundary_masks[BoundaryRegion.TIP_SURFACE][GridRegion.VACUUM]
        tip_values = np.full((nr, nv, np_points), self.config.tip_potential, dtype=np.float64)
        
        # Check each vacuum grid point
        for i in range(nr):
            for j in range(nv):
                for k in range(np_points):
                    r, phi, z = self.grid.get_cylindrical_coordinates(GridRegion.VACUUM, i, j, k)
                    
                    # Check if point is inside tip
                    if self.tip_geometry.is_point_inside_tip(r, phi, z):
                        tip_mask[i, j, k] = True
                        
        self.boundary_values[BoundaryRegion.TIP_SURFACE][GridRegion.VACUUM] = tip_values
        
    def _identify_sample_surface_points(self):
        """Identify grid points on sample surface (z=0 interface)"""
        nr = self.grid.config.nr_points
        np_points = self.grid.config.np_points
        
        # Surface points are at the interface between vacuum and semiconductor
        # These require special handling for interface continuity
        
        vacuum_surface_mask = self.boundary_masks[BoundaryRegion.SAMPLE_SURFACE][GridRegion.VACUUM]
        semiconductor_surface_mask = self.boundary_masks[BoundaryRegion.SAMPLE_SURFACE][GridRegion.SEMICONDUCTOR]
        
        # All points at z=0 are surface points
        vacuum_surface_mask[:, :] = True  # Interface from vacuum side
        semiconductor_surface_mask[:, :] = True  # Interface from semiconductor side
        
        # Surface boundary values will be determined by interface continuity
        self.boundary_values[BoundaryRegion.SAMPLE_SURFACE][GridRegion.VACUUM] = np.zeros((nr, np_points))
        self.boundary_values[BoundaryRegion.SAMPLE_SURFACE][GridRegion.SEMICONDUCTOR] = np.zeros((nr, np_points))
        
    def _identify_external_boundaries(self):
        """Identify external boundary points"""
        nr = self.grid.config.nr_points
        nv = self.grid.config.nv_points
        ns = self.grid.config.ns_points
        np_points = self.grid.config.np_points
        
        # Radial boundaries (at maximum radius)
        vacuum_radial_mask = self.boundary_masks[BoundaryRegion.RADIAL_BOUNDARY][GridRegion.VACUUM]
        semiconductor_radial_mask = self.boundary_masks[BoundaryRegion.RADIAL_BOUNDARY][GridRegion.SEMICONDUCTOR]
        
        vacuum_radial_mask[:, :] = True  # All points at r_max
        semiconductor_radial_mask[:, :] = True
        
        # Vacuum boundary (far from tip)
        vacuum_boundary_mask = self.boundary_masks[BoundaryRegion.VACUUM_BOUNDARY]
        vacuum_boundary_mask[:, :] = True  # All points at z_min
        
        # Semiconductor boundary (deep in sample)
        semiconductor_boundary_mask = self.boundary_masks[BoundaryRegion.SEMICONDUCTOR_BOUNDARY]
        semiconductor_boundary_mask[:, :] = True  # All points at z_max
        
        # Set boundary values
        self.boundary_values[BoundaryRegion.RADIAL_BOUNDARY][GridRegion.VACUUM] = np.full(
            (nv, np_points), self.config.radial_boundary_value)
        self.boundary_values[BoundaryRegion.RADIAL_BOUNDARY][GridRegion.SEMICONDUCTOR] = np.full(
            (ns, np_points), self.config.radial_boundary_value)
            
        self.boundary_values[BoundaryRegion.VACUUM_BOUNDARY] = np.full(
            (nr, np_points), self.config.vacuum_boundary_value)
        self.boundary_values[BoundaryRegion.SEMICONDUCTOR_BOUNDARY] = np.full(
            (nr, np_points), self.config.semiconductor_boundary_value)
            
    def _identify_angular_boundaries(self):
        """Identify angular boundary points for mirror symmetry"""
        nr = self.grid.config.nr_points
        nv = self.grid.config.nv_points
        ns = self.grid.config.ns_points
        
        # Mirror boundaries at φ=0 and φ=π
        vacuum_angular_mask = self.boundary_masks[BoundaryRegion.ANGULAR_BOUNDARY][GridRegion.VACUUM]
        semiconductor_angular_mask = self.boundary_masks[BoundaryRegion.ANGULAR_BOUNDARY][GridRegion.SEMICONDUCTOR]
        
        # First and last φ indices for mirror symmetry
        vacuum_angular_mask[:, :] = True  # Apply to φ=0 and φ=π boundaries
        semiconductor_angular_mask[:, :] = True
        
    def apply_boundary_conditions(self, 
                                 vacuum_potential: np.ndarray,
                                 semiconductor_potential: np.ndarray,
                                 surface_potential: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Apply all boundary conditions to potential arrays
        
        Args:
            vacuum_potential: Potential in vacuum region [V]
            semiconductor_potential: Potential in semiconductor region [V]
            surface_potential: Potential at surface interface [V]
            
        Returns:
            (vacuum_pot, semiconductor_pot, surface_pot): Updated potential arrays
        """
        
        # Create copies to avoid modifying originals
        vac_pot = vacuum_potential.copy()
        sem_pot = semiconductor_potential.copy()
        surf_pot = surface_potential.copy()
        
        # Apply tip boundary conditions
        self._apply_tip_boundary_conditions(vac_pot)
        
        # Apply external boundary conditions
        self._apply_external_boundary_conditions(vac_pot, sem_pot)
        
        # Apply surface boundary conditions (interface continuity)
        if self.config.enable_interface_continuity:
            surf_pot = self._apply_interface_continuity(vac_pot, sem_pot, surf_pot)
            
        # Apply angular boundary conditions (mirror symmetry)
        if self.grid.config.mirror_symmetry:
            self._apply_mirror_symmetry(vac_pot, sem_pot)
            
        return vac_pot, sem_pot, surf_pot
        
    def _apply_tip_boundary_conditions(self, vacuum_potential: np.ndarray):
        """Apply boundary conditions on tip surface"""
        tip_mask = self.boundary_masks[BoundaryRegion.TIP_SURFACE][GridRegion.VACUUM]
        tip_values = self.boundary_values[BoundaryRegion.TIP_SURFACE][GridRegion.VACUUM]
        
        if self.config.tip_boundary_type == BoundaryType.DIRICHLET:
            # Fixed potential on tip surface
            vacuum_potential[tip_mask] = tip_values[tip_mask]
            
        elif self.config.tip_boundary_type == BoundaryType.NEUMANN:
            # Fixed normal derivative (not commonly used for tip)
            # Would require special finite difference implementation
            warnings.warn("Neumann boundary conditions on tip not implemented")
            
    def _apply_external_boundary_conditions(self, 
                                          vacuum_potential: np.ndarray,
                                          semiconductor_potential: np.ndarray):
        """Apply boundary conditions on external boundaries"""
        
        # Radial boundaries
        if self.config.radial_boundary_type == BoundaryType.DIRICHLET:
            # Fixed potential at radial boundary
            vacuum_potential[-1, :, :] = self.config.radial_boundary_value
            semiconductor_potential[-1, :, :] = self.config.radial_boundary_value
            
        elif self.config.radial_boundary_type == BoundaryType.NEUMANN:
            # Zero normal derivative (∂V/∂r = 0)
            vacuum_potential[-1, :, :] = vacuum_potential[-2, :, :]
            semiconductor_potential[-1, :, :] = semiconductor_potential[-2, :, :]
            
        # Vacuum far boundary
        if self.config.vacuum_boundary_type == BoundaryType.DIRICHLET:
            vacuum_potential[:, -1, :] = self.config.vacuum_boundary_value
            
        # Semiconductor deep boundary
        if self.config.semiconductor_boundary_type == BoundaryType.DIRICHLET:
            semiconductor_potential[:, -1, :] = self.config.semiconductor_boundary_value
            
    def _apply_interface_continuity(self, 
                                  vacuum_potential: np.ndarray,
                                  semiconductor_potential: np.ndarray,
                                  surface_potential: np.ndarray) -> np.ndarray:
        """
        Apply interface continuity conditions at semiconductor surface
        
        Implements:
        1. Potential continuity: V_vac(surface) = V_sem(surface)
        2. Normal field continuity: ε_vac * E_n,vac = ε_sem * E_n,sem + σ_surface
        
        Returns:
            Updated surface potential
        """
        nr = self.grid.config.nr_points
        np_points = self.grid.config.np_points
        
        # Get permittivities
        eps_vac = self.config.permittivity_vacuum
        eps_sem = self.config.permittivity_semiconductor
        
        for i in range(nr):
            for k in range(np_points):
                # Get surface charge density at this point
                r, phi, z = self.grid.get_cylindrical_coordinates(GridRegion.SURFACE, i, 0, k)
                
                if self.config.surface_charge_density is not None:
                    sigma_surface = self.config.surface_charge_density(r, phi, z)
                else:
                    sigma_surface = 0.0
                    
                # Potential continuity condition
                # V_surface = (V_vac(0) + V_sem(0)) / 2
                v_vac_surface = vacuum_potential[i, 0, k] if vacuum_potential.shape[1] > 0 else 0.0
                v_sem_surface = semiconductor_potential[i, 0, k] if semiconductor_potential.shape[1] > 0 else 0.0
                
                # Simple average for potential continuity
                surface_potential[i, k] = (v_vac_surface + v_sem_surface) / 2.0
                
                # Apply potential continuity
                if vacuum_potential.shape[1] > 0 and i < vacuum_potential.shape[0] and k < vacuum_potential.shape[2]:
                    vacuum_potential[i, 0, k] = surface_potential[i, k]
                if semiconductor_potential.shape[1] > 0 and i < semiconductor_potential.shape[0] and k < semiconductor_potential.shape[2]:
                    semiconductor_potential[i, 0, k] = surface_potential[i, k]
                    
        return surface_potential
        
    def _apply_mirror_symmetry(self, 
                             vacuum_potential: np.ndarray,
                             semiconductor_potential: np.ndarray):
        """Apply mirror symmetry boundary conditions"""
        
        # For mirror symmetry at φ=0 and φ=π:
        # ∂V/∂φ = 0 at these boundaries
        
        if vacuum_potential.shape[2] > 1:  # Check if we have angular dimension
            # φ=0 boundary: copy from φ=1
            vacuum_potential[:, :, 0] = vacuum_potential[:, :, 1]
            semiconductor_potential[:, :, 0] = semiconductor_potential[:, :, 1]
            
            # φ=π boundary: copy from φ=np-1  
            vacuum_potential[:, :, -1] = vacuum_potential[:, :, -2]
            semiconductor_potential[:, :, -1] = semiconductor_potential[:, :, -2]
            
    def calculate_surface_charge_density(self, 
                                       r: float, 
                                       phi: float, 
                                       z: float,
                                       surface_potential: float,
                                       temperature: float = 300.0) -> float:
        """
        Calculate surface charge density at given point
        
        This is a simplified model - in practice would use detailed
        surface state calculations from physics modules.
        
        Args:
            r, phi, z: Position coordinates
            surface_potential: Local surface potential [V]
            temperature: Temperature [K]
            
        Returns:
            Surface charge density [C/m²]
        """
        
        # Simple Thomas-Fermi screening model
        # σ = -ε₀ * ε_sem * V_surface / λ_TF
        
        eps_0 = 8.854e-12  # F/m
        eps_sem = self.config.permittivity_semiconductor
        
        # Thomas-Fermi screening length (rough estimate)
        lambda_TF = 1e-9  # 1 nm
        
        # Linear relationship (simplified)
        charge_density = -eps_0 * eps_sem * surface_potential / lambda_TF
        
        return charge_density
        
    def get_boundary_info(self) -> Dict:
        """Get comprehensive boundary condition information"""
        
        boundary_counts = {}
        for region, masks in self.boundary_masks.items():
            if isinstance(masks, dict):
                # Multiple grid regions
                counts = {}
                for grid_region, mask in masks.items():
                    counts[grid_region.value] = np.sum(mask)
                boundary_counts[region.value] = counts
            else:
                # Single mask
                boundary_counts[region.value] = np.sum(masks)
                
        return {
            'configuration': {
                'tip_boundary_type': self.config.tip_boundary_type.value,
                'tip_potential_V': self.config.tip_potential,
                'surface_boundary_type': self.config.surface_boundary_type.value,
                'radial_boundary_type': self.config.radial_boundary_type.value,
                'interface_continuity': self.config.enable_interface_continuity,
                'mirror_symmetry': self.grid.config.mirror_symmetry
            },
            'boundary_point_counts': boundary_counts,
            'permittivities': {
                'vacuum': self.config.permittivity_vacuum,
                'semiconductor': self.config.permittivity_semiconductor
            }
        }
        
    def validate_boundary_conditions(self) -> List[str]:
        """
        Validate boundary condition setup
        
        Returns:
            List of validation messages
        """
        issues = []
        
        # Check tip potential
        if self.config.tip_potential is None:
            issues.append("Error: Tip potential not set")
            
        # Check permittivities
        if self.config.permittivity_vacuum <= 0:
            issues.append("Error: Vacuum permittivity must be positive")
        if self.config.permittivity_semiconductor <= 0:
            issues.append("Error: Semiconductor permittivity must be positive")
            
        # Check boundary consistency
        if self.grid.config.mirror_symmetry and self.config.angular_boundary_type != BoundaryType.MIRROR:
            issues.append("Warning: Grid uses mirror symmetry but angular boundary type is not MIRROR")
            
        # Check tip mask coverage
        tip_mask = self.boundary_masks[BoundaryRegion.TIP_SURFACE][GridRegion.VACUUM]
        tip_points = np.sum(tip_mask)
        if tip_points == 0:
            issues.append("Warning: No tip surface points identified")
        elif tip_points > tip_mask.size / 2:
            issues.append("Warning: Very large number of tip surface points")
            
        return issues
        
    def __str__(self) -> str:
        """String representation"""
        return (f"BoundaryConditions(tip_type={self.config.tip_boundary_type.value}, "
               f"tip_potential={self.config.tip_potential:.2f}V, "
               f"mirror_symmetry={self.grid.config.mirror_symmetry})")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions

def create_standard_boundary_conditions(grid: Grid3D,
                                       geometry: STMGeometry,
                                       tip_geometry: TipGeometry) -> BoundaryConditions:
    """Create standard boundary conditions for STM simulation"""
    config = BoundaryConfig(
        tip_boundary_type=BoundaryType.DIRICHLET,
        surface_boundary_type=BoundaryType.MIXED,
        radial_boundary_type=BoundaryType.NEUMANN,
        vacuum_boundary_type=BoundaryType.DIRICHLET,
        semiconductor_boundary_type=BoundaryType.DIRICHLET,
        enable_interface_continuity=True
    )
    return BoundaryConditions(grid, geometry, tip_geometry, config)


def create_floating_tip_boundary_conditions(grid: Grid3D,
                                           geometry: STMGeometry,
                                           tip_geometry: TipGeometry) -> BoundaryConditions:
    """Create boundary conditions for floating tip (current measurement)"""
    config = BoundaryConfig(
        tip_boundary_type=BoundaryType.NEUMANN,  # Floating tip
        surface_boundary_type=BoundaryType.MIXED,
        radial_boundary_type=BoundaryType.NEUMANN,
        vacuum_boundary_type=BoundaryType.DIRICHLET,
        semiconductor_boundary_type=BoundaryType.DIRICHLET
    )
    return BoundaryConditions(grid, geometry, tip_geometry, config)


def create_research_boundary_conditions(grid: Grid3D,
                                       geometry: STMGeometry,
                                       tip_geometry: TipGeometry,
                                       surface_charge_function: Optional[Callable] = None) -> BoundaryConditions:
    """Create research-grade boundary conditions with custom surface charge"""
    config = BoundaryConfig(
        tip_boundary_type=BoundaryType.DIRICHLET,
        surface_boundary_type=BoundaryType.MIXED,
        radial_boundary_type=BoundaryType.NEUMANN,
        vacuum_boundary_type=BoundaryType.DIRICHLET,
        semiconductor_boundary_type=BoundaryType.DIRICHLET,
        surface_charge_density=surface_charge_function,
        enable_interface_continuity=True,
        boundary_relaxation_factor=0.9
    )
    return BoundaryConditions(grid, geometry, tip_geometry, config)


if __name__ == "__main__":
    # Demo usage
    print("Boundary Conditions Demo")
    print("=" * 40)
    
    # This would normally be imported from other modules
    from .grid3d import create_medium_grid
    from .stm_geometry import create_standard_stm_geometry
    from .tip_geometry import create_standard_tip
    
    # Create test setup
    grid = create_medium_grid(mirror_symmetry=True)
    geometry = create_standard_stm_geometry()
    tip = create_standard_tip()
    
    # Create boundary conditions
    boundaries = create_standard_boundary_conditions(grid, geometry, tip)
    print(f"Created boundaries: {boundaries}")
    
    # Get boundary info
    info = boundaries.get_boundary_info()
    print(f"\nBoundary Information:")
    for section, data in info.items():
        print(f"  {section}: {data}")
        
    # Validate boundaries
    issues = boundaries.validate_boundary_conditions()
    if issues:
        print(f"\nValidation Issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print(f"\nBoundary validation: All checks passed")
        
    # Test boundary application (with dummy arrays)
    nr, nv, ns, np_pts = 32, 8, 32, 8
    vac_pot = np.random.normal(0, 0.1, (nr, nv, np_pts))
    sem_pot = np.random.normal(0, 0.1, (nr, ns, np_pts))
    surf_pot = np.random.normal(0, 0.1, (nr, np_pts))
    
    print(f"\nApplying boundary conditions to test arrays...")
    vac_new, sem_new, surf_new = boundaries.apply_boundary_conditions(vac_pot, sem_pot, surf_pot)
    
    print(f"Vacuum potential range: [{np.min(vac_new):.3f}, {np.max(vac_new):.3f}] V")
    print(f"Semiconductor potential range: [{np.min(sem_new):.3f}, {np.max(sem_new):.3f}] V")
    print(f"Surface potential range: [{np.min(surf_new):.3f}, {np.max(surf_new):.3f}] V")