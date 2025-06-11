"""
Pysemitip Physics Models Module

This package implements the core physical models for scanning tunneling microscopy
calculations, including self-consistent Poisson equation solving, charge density
calculations, and tunneling current computations.

Phase 2: Physics Model Layer Implementation (Complete)
Phase 3: Integration with Fortran-Equivalent Accuracy (In Progress)

Author: odindino
"""

# Core physics modules
from .materials import MaterialDatabase, SemiconductorMaterial, MaterialParameters, SurfaceStateParameters
from .charge_density import ChargeDensityCalculator, ChargeDensityConfig
from .poisson import PoissonSolver, PoissonConfig, Grid3D

# Tunneling current implementations
from .tunneling_current import TunnelingCurrentCalculator, TunnelingConfig
from .tunneling_current_fortran_equivalent import (
    FortranEquivalentTunnelingCalculator, 
    TunnelingCurrentConfig as FortranTunnelingConfig,
    calculate_fortran_equivalent_current
)

# Unified interface (recommended for new code)
from .tunneling_current_unified import (
    UnifiedTunnelingCurrentCalculator,
    UnifiedTunnelingConfig,
    UnifiedTunnelingResult,
    CalculationMethod,
    AccuracyLevel,
    calculate_tunneling_current_simple,
    calculate_tunneling_current_advanced
)

# Backward compatibility
TunnelingCurrentCalculatorOriginal = TunnelingCurrentCalculator
TunnelingConfigOriginal = TunnelingConfig

__all__ = [
    # Core modules
    'MaterialDatabase', 'SemiconductorMaterial', 'MaterialParameters', 'SurfaceStateParameters',
    'ChargeDensityCalculator', 'ChargeDensityConfig',
    'PoissonSolver', 'PoissonConfig', 'Grid3D',
    
    # Tunneling current implementations
    'TunnelingCurrentCalculator', 'TunnelingConfig',  # Original (backward compatibility)
    'FortranEquivalentTunnelingCalculator', 'FortranTunnelingConfig', 'calculate_fortran_equivalent_current',
    
    # Unified interface (recommended)
    'UnifiedTunnelingCurrentCalculator', 'UnifiedTunnelingConfig', 'UnifiedTunnelingResult',
    'CalculationMethod', 'AccuracyLevel',
    'calculate_tunneling_current_simple', 'calculate_tunneling_current_advanced',
    
    # Backward compatibility aliases
    'TunnelingCurrentCalculatorOriginal', 'TunnelingConfigOriginal'
]

# Version tracking
__version__ = "0.3.0"  # Updated for Phase 3 integration
__author__ = "odindino"

# Module-level convenience function (recommended entry point)
def create_tunneling_calculator(accuracy_level: AccuracyLevel = AccuracyLevel.BALANCED,
                               method: CalculationMethod = CalculationMethod.AUTO) -> UnifiedTunnelingCurrentCalculator:
    """
    Create a tunneling current calculator with specified accuracy and method.
    
    This is the recommended way to create tunneling current calculators in new code.
    
    Args:
        accuracy_level: Required accuracy level
        method: Calculation method (AUTO for automatic selection)
        
    Returns:
        Configured UnifiedTunnelingCurrentCalculator
        
    Example:
        >>> from physics import create_tunneling_calculator, AccuracyLevel
        >>> calculator = create_tunneling_calculator(AccuracyLevel.HIGH_ACCURACY)
    """
    config = UnifiedTunnelingConfig(
        accuracy_level=accuracy_level,
        method=method
    )
    return UnifiedTunnelingCurrentCalculator(config)
