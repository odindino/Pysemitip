"""
STM Geometry Configuration Module

This module implements the basic STM geometry configuration including
tip-sample positioning, coordinate system setup, and fundamental
geometric parameters.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass, field
import warnings


@dataclass
class GeometryConfig:
    """Configuration parameters for STM geometry setup"""
    
    # Basic geometry parameters
    separation: float = 1.0              # Tip-sample separation [nm]
    tip_radius: float = 10.0             # Tip curvature radius [nm]
    tip_cone_angle: float = 15.0         # Tip cone half-angle [degrees]
    tip_radius2: float = 0.5             # Tip end spherical protrusion [nm]
    
    # Tip position offsets
    tip_x_offset: float = 0.0            # X position offset [nm]
    tip_y_offset: float = 0.0            # Y position offset [nm]
    
    # Sample parameters
    sample_work_function: float = 4.1    # Sample work function [eV]
    tip_work_function: float = 4.5       # Tip work function [eV]
    
    # Bias and contact potential
    bias_voltage: float = 1.0            # Applied bias voltage [V]
    contact_potential: float = 0.0       # Contact potential difference [eV]
    
    # Domain size
    max_radius: float = 50.0             # Maximum radial extent [nm]
    vacuum_thickness: float = 10.0       # Vacuum region thickness [nm]
    semiconductor_depth: float = 50.0    # Semiconductor region depth [nm]
    
    # Physical constants
    vacuum_permittivity: float = 8.854e-12   # F/m
    semiconductor_permittivity: float = 11.7  # Relative permittivity (Si)
    
    def __post_init__(self):
        """Validate configuration parameters"""
        if self.separation <= 0:
            raise ValueError("Separation must be positive")
        if self.tip_radius <= 0:
            raise ValueError("Tip radius must be positive")
        if not 0 < self.tip_cone_angle < 90:
            raise ValueError("Tip cone angle must be between 0 and 90 degrees")
        if self.tip_radius2 >= self.tip_radius:
            warnings.warn("Tip protrusion radius should be smaller than tip radius")


class STMGeometry:
    """
    STM Geometry Configuration Manager
    
    This class manages the basic STM geometry setup including tip positioning,
    coordinate system definition, and fundamental geometric calculations.
    
    Based on SEMITIP coordinate system:
    - z = 0: Semiconductor surface
    - z > 0: Semiconductor interior  
    - z < 0: Vacuum region (tip side)
    """
    
    def __init__(self, config: Optional[GeometryConfig] = None):
        """Initialize STM geometry with configuration"""
        self.config = config or GeometryConfig()
        
        # Derived geometric parameters
        self._calculate_derived_parameters()
        
        # Coordinate system bounds
        self._setup_coordinate_bounds()
        
    def _calculate_derived_parameters(self):
        """Calculate derived geometric parameters from basic config"""
        
        # Convert cone angle to slope (tangent)
        self.tip_slope = np.tan(np.radians(self.config.tip_cone_angle))
        
        # Hyperbolic coordinate parameters (from semitip3-6.1.f)
        self.etat = 1.0 / np.sqrt(1.0 + 1.0 / self.tip_slope**2)
        self.A = self.config.tip_radius * self.tip_slope**2 / self.etat
        self.sprime = self.A * self.etat
        self.Z0 = self.config.separation - self.sprime
        self.C = self.Z0 / self.sprime if self.sprime != 0 else 0.0
        
        # Total potential drop across tip-sample junction
        self.total_voltage = (self.config.bias_voltage + 
                             self.config.contact_potential + 
                             (self.config.tip_work_function - self.config.sample_work_function))
        
    def _setup_coordinate_bounds(self):
        """Setup coordinate system boundaries"""
        
        # Z coordinate bounds
        self.z_min = -self.config.vacuum_thickness        # Vacuum side (tip)
        self.z_max = self.config.semiconductor_depth      # Semiconductor side
        self.z_surface = 0.0                              # Surface reference
        
        # Radial bounds
        self.r_min = 0.0
        self.r_max = self.config.max_radius
        
        # Angular bounds (full circle by default)
        self.phi_min = 0.0
        self.phi_max = 2.0 * np.pi
        
    def get_tip_position(self) -> Tuple[float, float, float]:
        """Get tip position in Cartesian coordinates"""
        return (
            self.config.tip_x_offset,
            self.config.tip_y_offset, 
            -self.config.separation  # Tip is at negative z
        )
        
    def get_sample_surface_z(self) -> float:
        """Get z coordinate of sample surface"""
        return self.z_surface
        
    def is_point_in_tip(self, r: float, z: float) -> bool:
        """
        Check if a point (r, z) is inside the tip geometry
        
        Args:
            r: Radial coordinate [nm]
            z: Axial coordinate [nm]
            
        Returns:
            True if point is inside tip
        """
        # Convert to tip-centered coordinates
        z_tip = z + self.config.separation
        
        # Check spherical protrusion at tip end
        if z_tip <= 0:
            if r <= self.config.tip_radius2:
                # Inside spherical end
                tip_surface_z = -np.sqrt(max(0, self.config.tip_radius2**2 - r**2))
                return z_tip >= tip_surface_z
                
        # Check conical section
        if z_tip > 0:
            tip_radius_at_z = self.config.tip_radius2 + z_tip * self.tip_slope
            return r <= tip_radius_at_z
            
        return False
        
    def calculate_local_field_enhancement(self, r: float, z: float) -> float:
        """
        Calculate local field enhancement factor at given position
        
        Simple geometric field enhancement based on tip curvature
        """
        if self.is_point_in_tip(r, z):
            return 0.0  # No field inside tip
            
        # Distance from tip apex
        tip_x, tip_y, tip_z = self.get_tip_position()
        distance = np.sqrt(r**2 + (z - tip_z)**2)
        
        # Simple field enhancement (1/r dependence near tip)
        if distance > 0:
            enhancement = self.config.tip_radius / max(distance, self.config.tip_radius2)
            return min(enhancement, 100.0)  # Cap at reasonable value
        
        return 1.0
        
    def get_electric_field_direction(self, r: float, z: float) -> Tuple[float, float]:
        """
        Get electric field direction (unit vector) at given position
        
        Returns:
            (E_r, E_z): Field components in cylindrical coordinates
        """
        tip_x, tip_y, tip_z = self.get_tip_position()
        
        # Vector from point to tip apex
        dr = -r  # Radial component points inward
        dz = tip_z - z
        
        # Normalize
        magnitude = np.sqrt(dr**2 + dz**2)
        if magnitude > 0:
            return (dr / magnitude, dz / magnitude)
        else:
            return (0.0, 0.0)
            
    def calculate_image_potential(self, r: float, z: float) -> float:
        """
        Calculate image potential contribution at given position
        
        Simple image potential model for tip-sample geometry
        """
        if z >= 0:
            # In semiconductor - image in surface
            distance_to_surface = z
            if distance_to_surface > 1e-3:  # Avoid singularity
                return -0.36 / (4.0 * distance_to_surface)  # eV*nm / (4*distance)
                
        else:
            # In vacuum - consider both surface and tip images
            distance_to_surface = abs(z)
            tip_distance = np.sqrt(r**2 + (z + self.config.separation)**2)
            
            image_surface = 0.0
            image_tip = 0.0
            
            if distance_to_surface > 1e-3:
                image_surface = -0.36 / (4.0 * distance_to_surface)
                
            if tip_distance > 1e-3:
                image_tip = -0.36 / (4.0 * tip_distance)
                
            return image_surface + image_tip
            
        return 0.0
        
    def get_work_function_profile(self, z: float) -> float:
        """
        Get work function as a function of z position
        
        Args:
            z: Axial position [nm]
            
        Returns:
            Work function [eV]
        """
        if z < -self.config.separation + self.config.tip_radius2:
            # Inside or very close to tip
            return self.config.tip_work_function
        elif z >= 0:
            # In semiconductor
            return self.config.sample_work_function
        else:
            # In vacuum - smooth transition
            tip_z = -self.config.separation
            surface_z = 0.0
            
            # Linear interpolation between tip and surface
            fraction = (z - tip_z) / (surface_z - tip_z)
            fraction = max(0.0, min(1.0, fraction))
            
            return (self.config.tip_work_function * (1 - fraction) + 
                   self.config.sample_work_function * fraction)
                   
    def get_geometry_summary(self) -> Dict:
        """Get summary of geometry configuration"""
        return {
            'separation_nm': self.config.separation,
            'tip_radius_nm': self.config.tip_radius,
            'tip_cone_angle_deg': self.config.tip_cone_angle,
            'tip_position': self.get_tip_position(),
            'coordinate_bounds': {
                'z_min': self.z_min,
                'z_max': self.z_max,
                'r_max': self.r_max
            },
            'hyperbolic_parameters': {
                'ETAT': self.etat,
                'A': self.A,
                'Z0': self.Z0,
                'C': self.C
            },
            'bias_voltage_V': self.config.bias_voltage,
            'total_voltage_V': self.total_voltage
        }
        
    def validate_geometry(self) -> List[str]:
        """
        Validate geometry configuration and return list of warnings/errors
        """
        issues = []
        
        # Check physical reasonableness
        if self.config.separation < self.config.tip_radius2:
            issues.append("Warning: Separation smaller than tip protrusion radius")
            
        if self.config.tip_radius2 > self.config.tip_radius / 2:
            issues.append("Warning: Tip protrusion very large compared to tip radius")
            
        if abs(self.config.tip_x_offset) > self.config.max_radius / 2:
            issues.append("Warning: Tip X offset large compared to simulation domain")
            
        if abs(self.config.tip_y_offset) > self.config.max_radius / 2:
            issues.append("Warning: Tip Y offset large compared to simulation domain")
            
        # Check derived parameters
        if self.Z0 <= 0:
            issues.append("Error: Hyperbolic coordinate Z0 is not positive")
            
        if abs(self.C) > 100:
            issues.append("Warning: Hyperbolic parameter C is very large")
            
        return issues
        
    def __str__(self) -> str:
        """String representation of geometry"""
        return (f"STMGeometry(sep={self.config.separation:.2f}nm, "
               f"R_tip={self.config.tip_radius:.1f}nm, "
               f"angle={self.config.tip_cone_angle:.1f}°, "
               f"bias={self.config.bias_voltage:.2f}V)")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions for common configurations

def create_standard_stm_geometry(separation: float = 1.0, 
                               bias_voltage: float = 1.0,
                               tip_radius: float = 10.0) -> STMGeometry:
    """Create standard STM geometry configuration"""
    config = GeometryConfig(
        separation=separation,
        bias_voltage=bias_voltage,
        tip_radius=tip_radius,
        tip_cone_angle=15.0,
        tip_radius2=0.5
    )
    return STMGeometry(config)


def create_sharp_tip_geometry(separation: float = 0.5,
                            bias_voltage: float = 2.0) -> STMGeometry:
    """Create sharp tip geometry for high-resolution imaging"""
    config = GeometryConfig(
        separation=separation,
        bias_voltage=bias_voltage,
        tip_radius=5.0,
        tip_cone_angle=10.0,
        tip_radius2=0.2
    )
    return STMGeometry(config)


def create_blunt_tip_geometry(separation: float = 2.0,
                            bias_voltage: float = 0.5) -> STMGeometry:
    """Create blunt tip geometry for stable tunneling"""
    config = GeometryConfig(
        separation=separation,
        bias_voltage=bias_voltage,
        tip_radius=20.0,
        tip_cone_angle=30.0,
        tip_radius2=1.0
    )
    return STMGeometry(config)


if __name__ == "__main__":
    # Demo usage
    print("STM Geometry Configuration Demo")
    print("=" * 40)
    
    # Create standard geometry
    geom = create_standard_stm_geometry()
    print(f"Standard geometry: {geom}")
    
    # Get summary
    summary = geom.get_geometry_summary()
    print("\nGeometry Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
        
    # Validate
    issues = geom.validate_geometry()
    if issues:
        print("\nValidation Issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("\nGeometry validation: All checks passed")
        
    # Test point checks
    print(f"\nPoint tests:")
    print(f"Point (0, -1) in tip: {geom.is_point_in_tip(0.0, -1.0)}")
    print(f"Point (0.3, -0.9) in tip: {geom.is_point_in_tip(0.3, -0.9)}")
    print(f"Point (2, 0) in tip: {geom.is_point_in_tip(2.0, 0.0)}")