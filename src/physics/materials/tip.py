"""
STM tip model and tip-related calculations.

This module defines the tip geometry and properties for SEMITIP simulations.
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

from ...utils.constants import PhysicalConstants as PC


@dataclass
class TipModel:
    """
    Represents the STM tip geometry and properties.
    
    The tip is modeled as a hyperboloid of revolution with optional
    hemispherical protrusion at the apex.
    """
    # Geometric parameters
    radius: float              # RAD - tip radius of curvature (nm)
    separation: float          # SEP - tip-sample separation (nm)
    slope: float              # SLOPE - tip cone slope parameter
    position: Tuple[float, float]  # (X0, Y0) - lateral position (nm)
    
    # Protrusion parameters
    protrusion_radius: float = 0.0  # RAD2 - hemispherical protrusion radius (nm)
    
    # Electronic properties
    work_function: float = 5.3  # CHI - tip work function (eV)
    fermi_level: float = 0.0   # EFTIP - tip Fermi level (eV)
    
    # Electrical parameters
    contact_potential: float = 0.0  # CPot - contact potential difference (V)
    bias_voltage: float = 0.0      # Applied bias voltage (V)
    
    def __post_init__(self):
        """Calculate derived parameters."""
        # Tip cone angle in degrees
        self.cone_angle = 360.0 * np.arctan(1.0 / self.slope) / (2 * np.pi)
        
        # Effective tip apex position (accounting for protrusion)
        self.apex_height = self.separation + self.protrusion_radius
        
        # Tip potential
        self.tip_potential = self.bias_voltage + self.contact_potential
    
    def hyperboloid_parameters(self) -> Tuple[float, float, float, float]:
        """
        Calculate hyperboloid parameters for the tip shape.
        
        Returns:
            Tuple of (eta, a, z0, c) parameters used in coordinate transformation
        """
        # These formulas match the Fortran code
        eta = 1.0 / np.sqrt(1.0 + self.slope**2)
        a = self.radius / self.slope
        z0 = a / self.slope
        c = np.sqrt(a**2 + z0**2)
        
        return eta, a, z0, c
    
    def is_inside_tip(self, r: float, z: float) -> bool:
        """
        Check if a point is inside the tip volume.
        
        Args:
            r: Radial coordinate (nm)
            z: Vertical coordinate from sample surface (nm)
            
        Returns:
            True if point is inside tip
        """
        # Account for lateral position
        r_eff = np.sqrt((r - self.position[0])**2 + self.position[1]**2)
        
        # Check hemispherical protrusion first
        if self.protrusion_radius > 0:
            # Distance from protrusion center
            dist_to_center = np.sqrt(r_eff**2 + (z - self.apex_height)**2)
            if dist_to_center <= self.protrusion_radius:
                return True
        
        # Check hyperboloid part
        # Using the hyperboloid equation: (z-z0)²/a² - r²/b² = 1
        eta, a, z0, c = self.hyperboloid_parameters()
        
        # Transform to hyperboloid coordinates
        z_tip = z - self.separation
        if z_tip > z0:
            # Above the apex - inside if within cone
            if r_eff <= (z_tip - z0) / self.slope:
                return True
        
        return False
    
    def surface_distance(self, r: float, z: float) -> float:
        """
        Calculate minimum distance from a point to the tip surface.
        
        Args:
            r: Radial coordinate (nm)
            z: Vertical coordinate from sample surface (nm)
            
        Returns:
            Distance to tip surface (nm), negative if inside tip
        """
        # This is a simplified version - full implementation would need
        # to solve for minimum distance to hyperboloid surface
        if self.is_inside_tip(r, z):
            return -1.0  # Inside tip
        
        # Approximate distance calculation
        # For points near apex, use spherical approximation
        r_eff = np.sqrt((r - self.position[0])**2 + self.position[1]**2)
        dist_to_apex = np.sqrt(r_eff**2 + (z - self.apex_height)**2)
        
        return dist_to_apex - self.radius
    
    def potential_contribution(self, r: float, z: float) -> float:
        """
        Calculate the potential contribution from the tip at a point.
        
        This is a placeholder for the full electrostatic calculation
        that would be done in the Poisson solver.
        
        Args:
            r: Radial coordinate (nm)
            z: Vertical coordinate from sample surface (nm)
            
        Returns:
            Potential contribution (V)
        """
        if self.is_inside_tip(r, z):
            return self.tip_potential
        
        # Outside tip - would need full Poisson solution
        # For now return 0 as placeholder
        return 0.0
    
    def adjust_separation(self, delta_sep: float):
        """
        Adjust tip-sample separation (for bias-dependent separation).
        
        Args:
            delta_sep: Change in separation (nm)
        """
        self.separation += delta_sep
        self.apex_height = self.separation + self.protrusion_radius


def create_tip_from_config(config) -> TipModel:
    """
    Create a TipModel from configuration data.
    
    Args:
        config: Configuration object containing tip parameters
        
    Returns:
        TipModel object
    """
    # Extract tip parameters
    tip_config = config.tip
    
    # Calculate slope from radius (if not provided)
    # Default slope = 1.0 corresponds to 90-degree cone angle
    slope = getattr(tip_config, 'slope', 1.0)
    
    # Get position
    x_pos = tip_config.position.x if hasattr(tip_config.position, 'x') else tip_config.x_position
    y_pos = tip_config.position.y if hasattr(tip_config.position, 'y') else tip_config.y_position
    
    return TipModel(
        radius=tip_config.radius,
        separation=tip_config.separation,
        slope=slope,
        position=(x_pos, y_pos),
        protrusion_radius=getattr(tip_config, 'protrusion_radius', 0.0),
        work_function=tip_config.work_function,
        fermi_level=getattr(config, 'tip_fermi_level', 0.0),
        contact_potential=getattr(config, 'contact_potential', 0.0),
        bias_voltage=0.0  # Will be set during simulation
    )