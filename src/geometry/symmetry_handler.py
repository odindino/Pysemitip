"""
Symmetry Handler Module

This module implements symmetry operations for STM simulations,
including mirror symmetry, rotational symmetry, and optimization
strategies for exploiting geometric symmetries.

Based on SEMITIP mirror symmetry implementation (MIRROR=1).

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union, Any
from dataclasses import dataclass
from enum import Enum
import warnings

from .grid3d import Grid3D, GridRegion


class SymmetryType(Enum):
    """Types of symmetries"""
    NONE = "none"
    MIRROR_XZ = "mirror_xz"       # Mirror symmetry about xz-plane (y=0)
    MIRROR_YZ = "mirror_yz"       # Mirror symmetry about yz-plane (x=0)
    ROTATIONAL = "rotational"     # Rotational symmetry about z-axis
    CYLINDRICAL = "cylindrical"   # Full cylindrical symmetry


class SymmetryOperation(Enum):
    """Symmetry operations"""
    IDENTITY = "identity"
    MIRROR = "mirror"
    ROTATE = "rotate"
    INVERT = "invert"


@dataclass
class SymmetryConfig:
    """Configuration for symmetry handling"""
    
    # Primary symmetry type
    symmetry_type: SymmetryType = SymmetryType.CYLINDRICAL
    
    # Mirror symmetry settings
    enable_mirror_symmetry: bool = True
    mirror_plane_normal: Tuple[float, float, float] = (0.0, 1.0, 0.0)  # y=0 plane
    
    # Rotational symmetry settings  
    rotational_order: int = 1  # N-fold rotational symmetry
    
    # Optimization settings
    use_symmetry_for_computation: bool = True
    symmetry_tolerance: float = 1e-10
    
    # Boundary handling
    apply_symmetry_to_boundaries: bool = True
    symmetry_boundary_relaxation: float = 1.0


class SymmetryHandler:
    """
    Symmetry Operations Manager
    
    This class implements symmetry operations for STM simulations:
    - Mirror symmetry (MIRROR=1 in SEMITIP)
    - Rotational symmetry optimization
    - Symmetry-preserving boundary conditions
    - Computational domain reduction
    """
    
    def __init__(self, 
                 grid: Grid3D,
                 config: Optional[SymmetryConfig] = None):
        """
        Initialize symmetry handler
        
        Args:
            grid: 3D cylindrical grid
            config: Symmetry configuration
        """
        self.grid = grid
        self.config = config or SymmetryConfig()
        
        # Determine effective symmetry based on grid and configuration
        self._determine_effective_symmetry()
        
        # Create symmetry operation matrices and mappings
        self._setup_symmetry_operations()
        
        # Identify symmetric grid points
        self._identify_symmetric_points()
        
    def _determine_effective_symmetry(self):
        """Determine the effective symmetry based on grid and configuration"""
        
        # Check if grid is set up for mirror symmetry
        if self.grid.config.mirror_symmetry:
            self.effective_symmetry = SymmetryType.MIRROR_XZ
            if self.config.symmetry_type != SymmetryType.MIRROR_XZ:
                warnings.warn("Grid has mirror symmetry but config specifies different symmetry")
                
        elif self.config.symmetry_type == SymmetryType.CYLINDRICAL:
            # Check if problem has cylindrical symmetry
            self.effective_symmetry = SymmetryType.CYLINDRICAL
            
        else:
            self.effective_symmetry = self.config.symmetry_type
            
    def _setup_symmetry_operations(self):
        """Setup symmetry operation matrices and transformations"""
        
        self.symmetry_operations = []
        
        if self.effective_symmetry == SymmetryType.MIRROR_XZ:
            # Mirror symmetry about xz-plane (y=0)
            # In cylindrical coordinates: (r, φ, z) → (r, -φ, z)
            self.symmetry_operations.append({
                'type': SymmetryOperation.MIRROR,
                'transform': self._mirror_xz_transform,
                'inverse_transform': self._mirror_xz_transform,  # Self-inverse
                'name': 'Mirror about xz-plane'
            })
            
        elif self.effective_symmetry == SymmetryType.CYLINDRICAL:
            # Full rotational symmetry - any rotation about z-axis
            for n in range(1, self.config.rotational_order):
                angle = 2.0 * np.pi * n / self.config.rotational_order
                self.symmetry_operations.append({
                    'type': SymmetryOperation.ROTATE,
                    'angle': angle,
                    'transform': lambda coords, theta=angle: self._rotate_z_transform(coords, theta),
                    'inverse_transform': lambda coords, theta=angle: self._rotate_z_transform(coords, -theta),
                    'name': f'Rotation by {np.degrees(angle):.1f}°'
                })
                
        elif self.effective_symmetry == SymmetryType.ROTATIONAL:
            # N-fold rotational symmetry
            for n in range(1, self.config.rotational_order):
                angle = 2.0 * np.pi * n / self.config.rotational_order
                self.symmetry_operations.append({
                    'type': SymmetryOperation.ROTATE,
                    'angle': angle,
                    'transform': lambda coords, theta=angle: self._rotate_z_transform(coords, theta),
                    'inverse_transform': lambda coords, theta=angle: self._rotate_z_transform(coords, -theta),
                    'name': f'Rotation by {np.degrees(angle):.1f}°'
                })
                
    def _identify_symmetric_points(self):
        """Identify symmetric grid point pairs"""
        
        self.symmetric_point_pairs = []
        self.independent_points = set()
        self.dependent_points = set()
        
        if self.effective_symmetry == SymmetryType.MIRROR_XZ:
            self._identify_mirror_symmetric_points()
        elif self.effective_symmetry in [SymmetryType.CYLINDRICAL, SymmetryType.ROTATIONAL]:
            self._identify_rotationally_symmetric_points()
            
    def _identify_mirror_symmetric_points(self):
        """Identify mirror symmetric point pairs for φ → -φ symmetry"""
        
        nr = self.grid.config.nr_points
        np_points = self.grid.config.np_points
        
        # For mirror symmetry, we only compute half the φ domain
        # φ ∈ [0, π] instead of [0, 2π]
        
        if not self.grid.config.mirror_symmetry:
            warnings.warn("Grid not configured for mirror symmetry")
            return
            
        # In mirror symmetric grid, all points are independent
        # The symmetry is enforced through boundary conditions
        for i in range(nr):
            for k in range(np_points):
                self.independent_points.add((i, k))
                
    def _identify_rotationally_symmetric_points(self):
        """Identify rotationally symmetric point sets"""
        
        nr = self.grid.config.nr_points
        np_points = self.grid.config.np_points
        
        # For rotational symmetry, identify equivalent φ indices
        phi_equivalence_classes = []
        processed_k = set()
        
        for k in range(np_points):
            if k in processed_k:
                continue
                
            # Find all φ indices equivalent to k under rotations
            equiv_class = [k]
            processed_k.add(k)
            
            for op in self.symmetry_operations:
                if op['type'] == SymmetryOperation.ROTATE:
                    # Transform φ coordinate
                    phi_k = self.grid.phi_points[k]
                    phi_transformed = (phi_k + op['angle']) % (2.0 * np.pi)
                    
                    # Find nearest grid point
                    k_transformed = np.argmin(np.abs(self.grid.phi_points - phi_transformed))
                    
                    if k_transformed not in equiv_class:
                        equiv_class.append(k_transformed)
                        processed_k.add(k_transformed)
                        
            phi_equivalence_classes.append(equiv_class)
            
        # Store equivalence relationships
        self.phi_equivalence_classes = phi_equivalence_classes
        
        # Mark independent and dependent points
        for equiv_class in phi_equivalence_classes:
            # First point in class is independent
            representative = equiv_class[0]
            for i in range(nr):
                self.independent_points.add((i, representative))
                
            # Other points are dependent
            for k in equiv_class[1:]:
                for i in range(nr):
                    self.dependent_points.add((i, k))
                    
    def _mirror_xz_transform(self, coordinates: Tuple[float, float, float]) -> Tuple[float, float, float]:
        """Transform coordinates by mirror reflection about xz-plane"""
        r, phi, z = coordinates
        
        # Mirror reflection: φ → -φ (or equivalently φ → 2π - φ)
        phi_mirrored = -phi
        if phi_mirrored < 0:
            phi_mirrored += 2.0 * np.pi
            
        return (r, phi_mirrored, z)
        
    def _rotate_z_transform(self, coordinates: Tuple[float, float, float], angle: float) -> Tuple[float, float, float]:
        """Transform coordinates by rotation about z-axis"""
        r, phi, z = coordinates
        
        # Rotation: φ → φ + angle
        phi_rotated = (phi + angle) % (2.0 * np.pi)
        
        return (r, phi_rotated, z)
        
    def apply_symmetry_to_potential(self, 
                                  potential_array: np.ndarray,
                                  region: GridRegion) -> np.ndarray:
        """
        Apply symmetry constraints to potential array
        
        Args:
            potential_array: Potential values on grid
            region: Grid region (vacuum, semiconductor, etc.)
            
        Returns:
            Symmetry-constrained potential array
        """
        
        if not self.config.use_symmetry_for_computation:
            return potential_array
            
        result = potential_array.copy()
        
        if self.effective_symmetry == SymmetryType.MIRROR_XZ:
            result = self._apply_mirror_symmetry_to_potential(result)
            
        elif self.effective_symmetry in [SymmetryType.CYLINDRICAL, SymmetryType.ROTATIONAL]:
            result = self._apply_rotational_symmetry_to_potential(result)
            
        return result
        
    def _apply_mirror_symmetry_to_potential(self, potential: np.ndarray) -> np.ndarray:
        """Apply mirror symmetry constraints to potential"""
        
        # For mirror symmetry about xz-plane (y=0):
        # V(r, φ, z) = V(r, -φ, z)
        
        # In a mirror symmetric grid, this is automatically satisfied
        # if boundary conditions are properly applied
        
        if len(potential.shape) >= 3:  # 3D array with φ dimension
            np_points = potential.shape[2]
            
            # Enforce symmetry: V(φ) = V(-φ)
            for k in range(np_points // 2):
                k_mirror = np_points - 1 - k
                if k != k_mirror:
                    # Average symmetric points
                    avg = (potential[:, :, k] + potential[:, :, k_mirror]) / 2.0
                    potential[:, :, k] = avg
                    potential[:, :, k_mirror] = avg
                    
        return potential
        
    def _apply_rotational_symmetry_to_potential(self, potential: np.ndarray) -> np.ndarray:
        """Apply rotational symmetry constraints to potential"""
        
        if not hasattr(self, 'phi_equivalence_classes'):
            return potential
            
        # For each equivalence class, set all values to their average
        for equiv_class in self.phi_equivalence_classes:
            if len(equiv_class) > 1 and len(potential.shape) >= 3:
                # Calculate average over equivalent φ points
                avg_potential = np.mean(potential[:, :, equiv_class], axis=2)
                
                # Set all equivalent points to average
                for k in equiv_class:
                    potential[:, :, k] = avg_potential
                    
        return potential
        
    def reduce_computational_domain(self, full_potential: np.ndarray) -> Tuple[np.ndarray, Dict]:
        """
        Reduce computational domain using symmetry
        
        Args:
            full_potential: Potential on full grid
            
        Returns:
            (reduced_potential, mapping_info): Reduced array and mapping information
        """
        
        if not self.config.use_symmetry_for_computation:
            return full_potential, {'reduction_factor': 1.0, 'symmetry_used': False}
            
        if self.effective_symmetry == SymmetryType.MIRROR_XZ:
            # Domain already reduced in grid generation
            return full_potential, {
                'reduction_factor': 2.0,
                'symmetry_used': True,
                'symmetry_type': 'mirror'
            }
            
        elif self.effective_symmetry == SymmetryType.CYLINDRICAL:
            # For cylindrical symmetry, can reduce to φ-independent
            if len(full_potential.shape) >= 3:
                # Average over all φ and keep only one φ slice
                reduced_potential = np.mean(full_potential, axis=2, keepdims=True)
                reduction_factor = full_potential.shape[2]
                
                return reduced_potential, {
                    'reduction_factor': reduction_factor,
                    'symmetry_used': True,
                    'symmetry_type': 'cylindrical'
                }
                
        elif self.effective_symmetry == SymmetryType.ROTATIONAL:
            # Reduce based on rotational order
            if len(full_potential.shape) >= 3:
                np_points = full_potential.shape[2]
                reduced_np = np_points // self.config.rotational_order
                
                # Keep only fundamental domain
                reduced_potential = full_potential[:, :, :reduced_np]
                
                return reduced_potential, {
                    'reduction_factor': self.config.rotational_order,
                    'symmetry_used': True,
                    'symmetry_type': 'rotational'
                }
                
        return full_potential, {'reduction_factor': 1.0, 'symmetry_used': False}
        
    def expand_from_reduced_domain(self, 
                                 reduced_potential: np.ndarray,
                                 mapping_info: Dict) -> np.ndarray:
        """
        Expand solution from reduced domain to full domain using symmetry
        
        Args:
            reduced_potential: Potential on reduced domain
            mapping_info: Information about domain reduction
            
        Returns:
            Full domain potential
        """
        
        if not mapping_info.get('symmetry_used', False):
            return reduced_potential
            
        symmetry_type = mapping_info.get('symmetry_type')
        
        if symmetry_type == 'mirror':
            # Mirror symmetry - domain already full for computational purposes
            return reduced_potential
            
        elif symmetry_type == 'cylindrical':
            # Expand φ-independent solution to full φ domain
            if len(reduced_potential.shape) >= 3:
                np_points = self.grid.config.np_points
                full_potential = np.repeat(reduced_potential, np_points, axis=2)
                return full_potential
                
        elif symmetry_type == 'rotational':
            # Expand using rotational symmetry
            if len(reduced_potential.shape) >= 3:
                order = self.config.rotational_order
                np_points = self.grid.config.np_points
                
                # Replicate reduced domain 'order' times
                full_potential = np.zeros(reduced_potential.shape[:2] + (np_points,))
                
                reduced_np = reduced_potential.shape[2]
                for n in range(order):
                    start_idx = n * reduced_np
                    end_idx = min((n + 1) * reduced_np, np_points)
                    copy_size = end_idx - start_idx
                    
                    full_potential[:, :, start_idx:end_idx] = reduced_potential[:, :, :copy_size]
                    
                return full_potential
                
        return reduced_potential
        
    def check_symmetry_preservation(self, 
                                  potential: np.ndarray,
                                  tolerance: Optional[float] = None) -> Dict:
        """
        Check if potential preserves expected symmetries
        
        Args:
            potential: Potential array to check
            tolerance: Tolerance for symmetry check
            
        Returns:
            Dictionary with symmetry check results
        """
        
        if tolerance is None:
            tolerance = self.config.symmetry_tolerance
            
        results = {
            'symmetry_preserved': True,
            'max_symmetry_error': 0.0,
            'symmetry_type': self.effective_symmetry.value,
            'details': {}
        }
        
        if self.effective_symmetry == SymmetryType.MIRROR_XZ:
            mirror_error = self._check_mirror_symmetry(potential, tolerance)
            results['max_symmetry_error'] = mirror_error
            results['symmetry_preserved'] = mirror_error < tolerance
            results['details']['mirror_error'] = mirror_error
            
        elif self.effective_symmetry in [SymmetryType.CYLINDRICAL, SymmetryType.ROTATIONAL]:
            rotation_error = self._check_rotational_symmetry(potential, tolerance)
            results['max_symmetry_error'] = rotation_error
            results['symmetry_preserved'] = rotation_error < tolerance
            results['details']['rotation_error'] = rotation_error
            
        return results
        
    def _check_mirror_symmetry(self, potential: np.ndarray, tolerance: float) -> float:
        """Check mirror symmetry preservation"""
        
        if len(potential.shape) < 3:
            return 0.0
            
        np_points = potential.shape[2]
        max_error = 0.0
        
        for k in range(np_points // 2):
            k_mirror = np_points - 1 - k
            if k != k_mirror:
                error = np.max(np.abs(potential[:, :, k] - potential[:, :, k_mirror]))
                max_error = max(max_error, error)
                
        return max_error
        
    def _check_rotational_symmetry(self, potential: np.ndarray, tolerance: float) -> float:
        """Check rotational symmetry preservation"""
        
        if not hasattr(self, 'phi_equivalence_classes') or len(potential.shape) < 3:
            return 0.0
            
        max_error = 0.0
        
        for equiv_class in self.phi_equivalence_classes:
            if len(equiv_class) > 1:
                # Check that all equivalent points have same value
                reference = potential[:, :, equiv_class[0]]
                for k in equiv_class[1:]:
                    error = np.max(np.abs(potential[:, :, k] - reference))
                    max_error = max(max_error, error)
                    
        return max_error
        
    def get_symmetry_info(self) -> Dict:
        """Get comprehensive symmetry information"""
        
        return {
            'configuration': {
                'symmetry_type': self.config.symmetry_type.value,
                'effective_symmetry': self.effective_symmetry.value,
                'enable_mirror_symmetry': self.config.enable_mirror_symmetry,
                'rotational_order': self.config.rotational_order,
                'use_for_computation': self.config.use_symmetry_for_computation
            },
            'operations': [
                {
                    'type': op['type'].value,
                    'name': op['name'],
                    'angle': op.get('angle', 0.0)
                }
                for op in self.symmetry_operations
            ],
            'point_counts': {
                'independent_points': len(self.independent_points),
                'dependent_points': len(self.dependent_points),
                'total_reduction_factor': len(self.independent_points) / 
                                        (len(self.independent_points) + len(self.dependent_points))
                                        if (len(self.independent_points) + len(self.dependent_points)) > 0 else 1.0
            },
            'grid_symmetry': {
                'mirror_symmetry': self.grid.config.mirror_symmetry,
                'phi_range': 'π' if self.grid.config.mirror_symmetry else '2π'
            }
        }
        
    def __str__(self) -> str:
        """String representation"""
        return (f"SymmetryHandler(type={self.effective_symmetry.value}, "
               f"operations={len(self.symmetry_operations)}, "
               f"reduction_factor={len(self.independent_points) / max(1, len(self.independent_points) + len(self.dependent_points)):.2f})")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions

def create_mirror_symmetry_handler(grid: Grid3D) -> SymmetryHandler:
    """Create symmetry handler for mirror symmetry (MIRROR=1)"""
    config = SymmetryConfig(
        symmetry_type=SymmetryType.MIRROR_XZ,
        enable_mirror_symmetry=True,
        use_symmetry_for_computation=True
    )
    return SymmetryHandler(grid, config)


def create_cylindrical_symmetry_handler(grid: Grid3D) -> SymmetryHandler:
    """Create symmetry handler for cylindrical symmetry"""
    config = SymmetryConfig(
        symmetry_type=SymmetryType.CYLINDRICAL,
        rotational_order=8,  # High order for near-cylindrical
        use_symmetry_for_computation=True
    )
    return SymmetryHandler(grid, config)


def create_rotational_symmetry_handler(grid: Grid3D, order: int = 4) -> SymmetryHandler:
    """Create symmetry handler for N-fold rotational symmetry"""
    config = SymmetryConfig(
        symmetry_type=SymmetryType.ROTATIONAL,
        rotational_order=order,
        use_symmetry_for_computation=True
    )
    return SymmetryHandler(grid, config)


def create_no_symmetry_handler(grid: Grid3D) -> SymmetryHandler:
    """Create symmetry handler with no symmetry optimizations"""
    config = SymmetryConfig(
        symmetry_type=SymmetryType.NONE,
        use_symmetry_for_computation=False
    )
    return SymmetryHandler(grid, config)


if __name__ == "__main__":
    # Demo usage
    print("Symmetry Handler Demo")
    print("=" * 40)
    
    # This would normally be imported from other modules
    from .grid3d import create_medium_grid
    
    # Create test setup with mirror symmetry
    grid_mirror = create_medium_grid(mirror_symmetry=True)
    symmetry_mirror = create_mirror_symmetry_handler(grid_mirror)
    print(f"Mirror symmetry: {symmetry_mirror}")
    
    # Create test setup with cylindrical symmetry
    grid_cyl = create_medium_grid(mirror_symmetry=False)
    symmetry_cyl = create_cylindrical_symmetry_handler(grid_cyl)
    print(f"Cylindrical symmetry: {symmetry_cyl}")
    
    # Get symmetry info
    info_mirror = symmetry_mirror.get_symmetry_info()
    print(f"\nMirror Symmetry Information:")
    for section, data in info_mirror.items():
        print(f"  {section}: {data}")
        
    # Test symmetry application on dummy potential
    nr, nv, np_pts = 32, 8, 16
    test_potential = np.random.normal(0, 0.1, (nr, nv, np_pts))
    
    print(f"\nTesting symmetry application...")
    symmetric_potential = symmetry_mirror.apply_symmetry_to_potential(test_potential, GridRegion.VACUUM)
    
    # Check symmetry preservation
    symmetry_check = symmetry_mirror.check_symmetry_preservation(symmetric_potential)
    print(f"Symmetry preservation check: {symmetry_check}")
    
    # Test domain reduction
    reduced_pot, mapping = symmetry_cyl.reduce_computational_domain(test_potential)
    print(f"Domain reduction: {mapping}")
    
    # Test expansion
    expanded_pot = symmetry_cyl.expand_from_reduced_domain(reduced_pot, mapping)
    print(f"Domain expansion shape: {test_potential.shape} → {reduced_pot.shape} → {expanded_pot.shape}")