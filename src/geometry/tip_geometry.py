"""
STM Tip Geometry Module

This module implements the hyperbolic coordinate system for STM tip geometry,
providing precise modeling of tip shape, tip surface identification, and
coordinate transformations.

Based on SEMITIP semitip3-6.1.f hyperbolic coordinate implementation.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union
from dataclasses import dataclass
import warnings
from enum import Enum


class TipRegion(Enum):
    """Tip region identifiers"""
    OUTSIDE = "outside"
    SPHERICAL_END = "spherical_end"
    CONICAL_SHAFT = "conical_shaft"
    HYPERBOLIC_TRANSITION = "hyperbolic_transition"


@dataclass
class TipConfig:
    """Configuration for STM tip geometry"""
    
    # Basic tip parameters
    radius: float = 10.0              # Tip curvature radius [nm]
    radius2: float = 0.5              # Tip end spherical protrusion [nm]
    cone_angle: float = 15.0          # Tip cone half-angle [degrees]
    separation: float = 1.0           # Tip-sample separation [nm]
    
    # Position offsets
    x_offset: float = 0.0             # X position offset [nm]
    y_offset: float = 0.0             # Y position offset [nm]
    
    # Electrical parameters
    work_function: float = 4.5        # Tip work function [eV]
    bias_voltage: float = 1.0         # Applied bias voltage [V]
    contact_potential: float = 0.0    # Contact potential difference [eV]
    
    # Hyperbolic coordinate parameters (calculated)
    slope: Optional[float] = None     # Tip slope (tan of cone angle)
    etat: Optional[float] = None      # ETAT parameter
    A: Optional[float] = None         # Hyperbolic focal parameter
    Z0: Optional[float] = None        # Hyperbolic coordinate origin
    C: Optional[float] = None         # Geometric parameter
    
    def __post_init__(self):
        """Calculate derived parameters and validate configuration"""
        
        # Validate inputs
        if self.radius <= 0:
            raise ValueError("Tip radius must be positive")
        if self.radius2 <= 0:
            raise ValueError("Tip end radius must be positive")
        if self.radius2 >= self.radius:
            warnings.warn("Tip end radius should be smaller than tip radius")
        if not 0 < self.cone_angle < 90:
            raise ValueError("Tip cone angle must be between 0 and 90 degrees")
        if self.separation <= 0:
            raise ValueError("Tip-sample separation must be positive")
            
        # Calculate derived parameters (from semitip3-6.1.f lines 97-101)
        self.slope = np.tan(np.radians(self.cone_angle))
        self.etat = 1.0 / np.sqrt(1.0 + 1.0 / self.slope**2)
        self.A = self.radius * self.slope**2 / self.etat
        sprime = self.A * self.etat
        self.Z0 = self.separation - sprime
        self.C = self.Z0 / sprime if sprime != 0 else 0.0


class TipGeometry:
    """
    STM Tip Geometry Manager
    
    This class implements the hyperbolic coordinate system for precise
    tip geometry modeling as used in SEMITIP. It provides:
    - Tip surface identification
    - Hyperbolic coordinate transformations
    - Electric field enhancement calculations
    - Integration with 3D grid systems
    """
    
    def __init__(self, config: Optional[TipConfig] = None):
        """Initialize tip geometry with configuration"""
        self.config = config or TipConfig()
        
        # Derived geometric parameters
        self._setup_coordinate_system()
        
        # Cache for performance
        self._tip_surface_cache = {}
        
    def _setup_coordinate_system(self):
        """Setup hyperbolic coordinate system parameters"""
        
        # Tip apex position in cylindrical coordinates
        # Tip is at z = -separation from surface
        self.tip_apex_z = -self.config.separation
        self.tip_apex_r = np.sqrt(self.config.x_offset**2 + self.config.y_offset**2)
        
        # Tip potential
        self.tip_potential = (self.config.bias_voltage + 
                             self.config.contact_potential + 
                             self.config.work_function)
                             
    def tip_surface_function(self, r: float) -> float:
        """
        Tip surface function p(r) from semitip3-6.1.f lines 600-605
        
        Args:
            r: Radial coordinate [nm]
            
        Returns:
            Height of tip surface above apex [nm]
        """
        if r < self.config.radius2:
            # Spherical end (semitip3-6.1.f line 602)
            return np.sqrt(self.config.radius2**2 - r**2)
        else:
            # Beyond spherical end - conical section
            return self.config.radius2 + (r - self.config.radius2) * self.config.slope
            
    def is_point_inside_tip(self, r: float, phi: float, z: float) -> bool:
        """
        Check if point (r, φ, z) is inside the tip geometry
        
        Args:
            r: Radial coordinate [nm]
            phi: Angular coordinate [rad]
            z: Axial coordinate [nm]
            
        Returns:
            True if point is inside tip
        """
        # Convert to tip-centered coordinates
        # Account for tip position offsets
        x = r * np.cos(phi) - self.config.x_offset
        y = r * np.sin(phi) - self.config.y_offset
        r_tip = np.sqrt(x**2 + y**2)
        z_tip = z - self.tip_apex_z  # z relative to tip apex
        
        if z_tip <= 0:
            # Below tip apex - inside spherical end
            if r_tip <= self.config.radius2:
                tip_surface_height = self.tip_surface_function(r_tip)
                return z_tip >= -tip_surface_height
                
        elif z_tip > 0:
            # Above tip apex - check conical section
            tip_radius_at_z = self.config.radius2 + z_tip * self.config.slope
            return r_tip <= tip_radius_at_z
            
        return False
        
    def get_tip_region(self, r: float, phi: float, z: float) -> TipRegion:
        """
        Identify which region of tip a point belongs to
        
        Returns:
            TipRegion enum value
        """
        if not self.is_point_inside_tip(r, phi, z):
            return TipRegion.OUTSIDE
            
        # Convert to tip-centered coordinates
        x = r * np.cos(phi) - self.config.x_offset
        y = r * np.sin(phi) - self.config.y_offset
        r_tip = np.sqrt(x**2 + y**2)
        z_tip = z - self.tip_apex_z
        
        if z_tip <= 0 and r_tip <= self.config.radius2:
            return TipRegion.SPHERICAL_END
        elif z_tip > 0:
            return TipRegion.CONICAL_SHAFT
        else:
            return TipRegion.HYPERBOLIC_TRANSITION
            
    def calculate_tip_surface_normal(self, r: float, phi: float, z: float) -> Tuple[float, float, float]:
        """
        Calculate surface normal vector at tip surface point
        
        Returns:
            (n_r, n_phi, n_z): Normal vector components in cylindrical coordinates
        """
        if not self.is_point_inside_tip(r, phi, z):
            return (0.0, 0.0, 0.0)
            
        # Convert to tip-centered coordinates
        x = r * np.cos(phi) - self.config.x_offset
        y = r * np.sin(phi) - self.config.y_offset
        r_tip = np.sqrt(x**2 + y**2)
        z_tip = z - self.tip_apex_z
        
        region = self.get_tip_region(r, phi, z)
        
        if region == TipRegion.SPHERICAL_END:
            # Spherical surface normal points radially outward
            if r_tip > 1e-12:  # Avoid division by zero
                nr = x / r_tip / self.config.radius2
                nphi = y / r_tip / self.config.radius2  
                nz = -z_tip / self.config.radius2
                
                # Convert to cylindrical coordinates
                n_r = nr * np.cos(phi) + nphi * np.sin(phi)
                n_phi = -nr * np.sin(phi) + nphi * np.cos(phi)
                
                return (n_r, n_phi, nz)
            else:
                return (0.0, 0.0, -1.0)  # At apex
                
        elif region == TipRegion.CONICAL_SHAFT:
            # Conical surface normal
            if r_tip > 1e-12:
                # Normal to cone surface
                slope_normal = 1.0 / np.sqrt(1.0 + self.config.slope**2)
                
                n_r = slope_normal
                n_phi = 0.0
                n_z = self.config.slope * slope_normal
                
                return (n_r, n_phi, n_z)
            else:
                return (0.0, 0.0, 1.0)
                
        return (0.0, 0.0, 0.0)
        
    def calculate_field_enhancement(self, r: float, phi: float, z: float) -> float:
        """
        Calculate electric field enhancement factor at given point
        
        Args:
            r, phi, z: Point coordinates
            
        Returns:
            Field enhancement factor (dimensionless)
        """
        if self.is_point_inside_tip(r, phi, z):
            return 0.0  # No field inside conductor
            
        # Distance from tip apex
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        tip_x = self.config.x_offset
        tip_y = self.config.y_offset  
        tip_z = self.tip_apex_z
        
        distance = np.sqrt((x - tip_x)**2 + (y - tip_y)**2 + (z - tip_z)**2)
        
        if distance < self.config.radius2:
            # Very close to tip apex - high enhancement
            return self.config.radius / self.config.radius2
        else:
            # Field enhancement based on tip curvature and distance
            enhancement = self.config.radius / max(distance, self.config.radius2)
            return min(enhancement, 100.0)  # Cap at reasonable value
            
    def hyperbolic_to_cylindrical(self, xi: float, eta: float, phi: float) -> Tuple[float, float, float]:
        """
        Convert hyperbolic coordinates (ξ, η, φ) to cylindrical (r, φ, z)
        
        Based on SEMITIP hyperbolic coordinate system
        
        Args:
            xi: Hyperbolic coordinate ξ
            eta: Hyperbolic coordinate η  
            phi: Angular coordinate φ [rad]
            
        Returns:
            (r, phi, z): Cylindrical coordinates
        """
        # Hyperbolic coordinate transformation
        A = self.config.A
        Z0 = self.config.Z0
        
        # Convert hyperbolic to Cartesian relative to hyperbolic origin
        x_hyp = A * np.sinh(xi) * np.cos(eta)
        z_hyp = A * np.cosh(xi) * np.sin(eta)
        
        # Convert to cylindrical coordinates
        r = np.sqrt(x_hyp**2)  # Assuming axisymmetric for now
        z = Z0 + z_hyp
        
        return (r, phi, z)
        
    def cylindrical_to_hyperbolic(self, r: float, phi: float, z: float) -> Tuple[float, float, float]:
        """
        Convert cylindrical coordinates to hyperbolic coordinates
        
        Args:
            r, phi, z: Cylindrical coordinates
            
        Returns:
            (xi, eta, phi): Hyperbolic coordinates
        """
        # Translate to hyperbolic origin
        A = self.config.A
        Z0 = self.config.Z0
        
        x_hyp = r  # Assuming axisymmetric
        z_hyp = z - Z0
        
        # Convert to hyperbolic coordinates
        if A > 0:
            # Calculate ξ and η
            rho = np.sqrt(x_hyp**2 + z_hyp**2)
            
            if rho > 1e-12:
                cos_sum = z_hyp / rho
                cos_diff = x_hyp / rho
                
                xi = np.arccosh((cos_sum + np.sqrt(cos_sum**2 + cos_diff**2 - 1)) / 2)
                eta = np.arccos((cos_sum - np.sqrt(cos_sum**2 + cos_diff**2 - 1)) / 2)
            else:
                xi = 0.0
                eta = 0.0
        else:
            xi = 0.0
            eta = 0.0
            
        return (xi, eta, phi)
        
    def get_tip_surface_points(self, n_points: int = 100) -> List[Tuple[float, float, float]]:
        """
        Generate points on tip surface for visualization
        
        Args:
            n_points: Number of surface points to generate
            
        Returns:
            List of (r, phi, z) coordinates on tip surface
        """
        surface_points = []
        
        # Generate points on spherical end
        n_sphere = n_points // 3
        for i in range(n_sphere):
            r_tip = self.config.radius2 * i / n_sphere if n_sphere > 0 else 0
            z_surface = self.tip_surface_function(r_tip)
            z_global = self.tip_apex_z + z_surface
            
            # Add point at phi=0 (can be extended for full 3D)
            surface_points.append((r_tip, 0.0, z_global))
            
        # Generate points on conical section
        n_cone = n_points - n_sphere
        max_cone_height = self.config.radius  # Arbitrary limit
        
        for i in range(n_cone):
            z_rel = (i + 1) * max_cone_height / n_cone if n_cone > 0 else 0
            r_tip = self.config.radius2 + z_rel * self.config.slope
            z_global = self.tip_apex_z + z_rel
            
            surface_points.append((r_tip, 0.0, z_global))
            
        return surface_points
        
    def calculate_capacitance(self, grid_spacing: float = 0.1) -> float:
        """
        Estimate tip-sample capacitance using simple geometry
        
        Args:
            grid_spacing: Grid spacing for numerical calculation [nm]
            
        Returns:
            Capacitance estimate [F]
        """
        # Simple parallel plate model with area based on tip end
        epsilon_0 = 8.854e-12  # F/m
        area = np.pi * self.config.radius2**2  # m^2 (convert from nm^2)
        area_m2 = area * 1e-18  # nm^2 to m^2
        separation_m = self.config.separation * 1e-9  # nm to m
        
        # Parallel plate capacitance
        C_parallel = epsilon_0 * area_m2 / separation_m
        
        # Add fringing field correction (rough estimate)
        fringing_factor = 1.0 + self.config.radius2 / self.config.separation
        
        return C_parallel * fringing_factor
        
    def get_tip_info(self) -> Dict:
        """Get comprehensive tip geometry information"""
        return {
            'basic_parameters': {
                'radius_nm': self.config.radius,
                'radius2_nm': self.config.radius2,
                'cone_angle_deg': self.config.cone_angle,
                'separation_nm': self.config.separation
            },
            'position': {
                'x_offset_nm': self.config.x_offset,
                'y_offset_nm': self.config.y_offset,
                'apex_z_nm': self.tip_apex_z
            },
            'hyperbolic_parameters': {
                'slope': self.config.slope,
                'ETAT': self.config.etat,
                'A': self.config.A,
                'Z0': self.config.Z0,
                'C': self.config.C
            },
            'electrical': {
                'work_function_eV': self.config.work_function,
                'bias_voltage_V': self.config.bias_voltage,
                'tip_potential_V': self.tip_potential
            },
            'derived': {
                'capacitance_estimate_F': self.calculate_capacitance()
            }
        }
        
    def validate_tip_geometry(self) -> List[str]:
        """
        Validate tip geometry and return warnings/errors
        
        Returns:
            List of validation messages
        """
        issues = []
        
        # Check physical reasonableness
        if self.config.radius2 > self.config.radius / 2:
            issues.append("Warning: Tip end radius very large compared to tip radius")
            
        if self.config.separation < self.config.radius2:
            issues.append("Warning: Separation smaller than tip end radius")
            
        if self.config.cone_angle > 45:
            issues.append("Warning: Very wide tip cone angle")
            
        # Check hyperbolic parameters
        if self.config.Z0 <= 0:
            issues.append("Error: Hyperbolic coordinate Z0 is not positive")
            
        if abs(self.config.C) > 100:
            issues.append("Warning: Hyperbolic parameter C is very large")
            
        if self.config.A <= 0:
            issues.append("Error: Hyperbolic parameter A is not positive")
            
        # Check tip position
        tip_offset = np.sqrt(self.config.x_offset**2 + self.config.y_offset**2)
        if tip_offset > self.config.separation:
            issues.append("Warning: Tip offset larger than separation distance")
            
        return issues
        
    def __str__(self) -> str:
        """String representation"""
        return (f"TipGeometry(R={self.config.radius:.1f}nm, "
               f"R2={self.config.radius2:.1f}nm, "
               f"angle={self.config.cone_angle:.1f}°, "
               f"sep={self.config.separation:.1f}nm)")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions for common tip configurations

def create_sharp_tip(separation: float = 0.5) -> TipGeometry:
    """Create sharp tip configuration"""
    config = TipConfig(
        radius=5.0,
        radius2=0.2,
        cone_angle=10.0,
        separation=separation
    )
    return TipGeometry(config)


def create_standard_tip(separation: float = 1.0) -> TipGeometry:
    """Create standard tip configuration"""
    config = TipConfig(
        radius=10.0,
        radius2=0.5,
        cone_angle=15.0,
        separation=separation
    )
    return TipGeometry(config)


def create_blunt_tip(separation: float = 2.0) -> TipGeometry:
    """Create blunt tip configuration"""
    config = TipConfig(
        radius=20.0,
        radius2=1.0,
        cone_angle=30.0,
        separation=separation
    )
    return TipGeometry(config)


if __name__ == "__main__":
    # Demo usage
    print("STM Tip Geometry Demo")
    print("=" * 40)
    
    # Create standard tip
    tip = create_standard_tip()
    print(f"Standard tip: {tip}")
    
    # Get tip info
    info = tip.get_tip_info()
    print(f"\nTip Information:")
    for section, params in info.items():
        print(f"  {section}:")
        for key, value in params.items():
            if isinstance(value, float):
                print(f"    {key}: {value:.3e}")
            else:
                print(f"    {key}: {value}")
                
    # Validate geometry
    issues = tip.validate_tip_geometry()
    if issues:
        print(f"\nValidation Issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print(f"\nTip geometry validation: All checks passed")
        
    # Test point inside tip
    print(f"\nPoint tests:")
    test_points = [
        (0.0, 0.0, -1.0),   # Near apex
        (0.3, 0.0, -0.9),   # In spherical end
        (2.0, 0.0, 0.0),    # Outside tip
        (1.0, 0.0, -0.5)    # In conical section
    ]
    
    for r, phi, z in test_points:
        inside = tip.is_point_inside_tip(r, phi, z)
        region = tip.get_tip_region(r, phi, z)
        enhancement = tip.calculate_field_enhancement(r, phi, z)
        print(f"  Point ({r:.1f}, {phi:.1f}, {z:.1f}): inside={inside}, region={region.value}, enhancement={enhancement:.1f}")
        
    # Surface points
    surface_pts = tip.get_tip_surface_points(20)
    print(f"\nGenerated {len(surface_pts)} surface points")
    print(f"Sample surface points:")
    for i in range(0, len(surface_pts), 5):
        r, phi, z = surface_pts[i]
        print(f"  ({r:.3f}, {phi:.3f}, {z:.3f})")