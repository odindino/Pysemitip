"""
Pysemitip Geometry Module

This package implements the 3D geometric modeling and grid management for 
STM simulations, including coordinate system definitions, grid generation,
and boundary condition handling.

Phase 3: Geometry and Grid Layer Implementation

Author: odindino
Date: 2025-06-11
"""

# Core geometry classes
from .stm_geometry import STMGeometry, GeometryConfig
from .grid3d import Grid3D, GridConfig, GridDimensions
from .tip_geometry import TipGeometry, TipConfig

# Advanced grid management
from .adaptive_refinement import AdaptiveGridRefinement, RefinementConfig
from .boundary_conditions import BoundaryConditions, BoundaryConfig

# Symmetry and optimization
from .symmetry_handler import SymmetryHandler, SymmetryConfig

__all__ = [
    # Core geometry
    'STMGeometry', 'GeometryConfig',
    'Grid3D', 'GridConfig', 'GridDimensions', 
    'TipGeometry', 'TipConfig',
    
    # Advanced features
    'AdaptiveGridRefinement', 'RefinementConfig',
    'BoundaryConditions', 'BoundaryConfig',
    'SymmetryHandler', 'SymmetryConfig'
]

__version__ = "0.1.0"
__author__ = "odindino"