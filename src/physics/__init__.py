"""
Pysemitip Physics Models Module

This package implements the core physical models for scanning tunneling microscopy
calculations, including self-consistent Poisson equation solving, charge density
calculations, and tunneling current computations.

Phase 2: Physics Model Layer Implementation

Author: odindino
"""

from .materials import MaterialDatabase, SemiconductorMaterial, MaterialParameters, SurfaceStateParameters
from .charge_density import ChargeDensityCalculator, ChargeDensityConfig
from .poisson import PoissonSolver, PoissonConfig, Grid3D
from .tunneling_current import TunnelingCurrentCalculator, TunnelingConfig

__all__ = [
    'MaterialDatabase', 'SemiconductorMaterial', 'MaterialParameters', 'SurfaceStateParameters',
    'ChargeDensityCalculator', 'ChargeDensityConfig',
    'PoissonSolver', 'PoissonConfig', 'Grid3D',
    'TunnelingCurrentCalculator', 'TunnelingConfig'
]

__version__ = "0.2.0"
__author__ = "odindino"
