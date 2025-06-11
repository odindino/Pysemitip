"""
Unified Tunneling Current API

This module provides a unified interface for tunneling current calculations,
supporting both the simplified original implementation and the complete
Fortran-equivalent implementation. Users can choose the appropriate level
of accuracy based on their requirements.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Optional, Tuple, Union, List
from dataclasses import dataclass
from enum import Enum
import warnings

from .materials import SemiconductorMaterial, PhysicalConstants
from .tunneling_current import TunnelingCurrentCalculator as OriginalCalculator, TunnelingConfig as OriginalConfig
from .tunneling_current_fortran_equivalent import (
    FortranEquivalentTunnelingCalculator,
    TunnelingCurrentConfig as FortranConfig,
    calculate_fortran_equivalent_current
)


class CalculationMethod(Enum):
    """Available calculation methods"""
    SIMPLIFIED = "simplified"           # Original implementation (fast, approximate)
    FORTRAN_EQUIVALENT = "fortran"     # Complete Fortran-equivalent (accurate, slower)
    AUTO = "auto"                      # Automatically choose based on requirements


class AccuracyLevel(Enum):
    """Accuracy level requirements"""
    FAST = "fast"                      # Speed prioritized over accuracy
    BALANCED = "balanced"              # Balance between speed and accuracy
    HIGH_ACCURACY = "high_accuracy"   # Maximum accuracy (Fortran-equivalent)
    RESEARCH_GRADE = "research_grade"  # Full scientific accuracy with all features


@dataclass
class UnifiedTunnelingConfig:
    """Unified configuration for tunneling current calculations"""
    
    # Method selection
    method: CalculationMethod = CalculationMethod.AUTO
    accuracy_level: AccuracyLevel = AccuracyLevel.BALANCED
    
    # Common parameters
    energy_points: int = 50
    k_points: int = 20
    temperature_tip: float = 4.2      # K
    temperature_sample: float = 4.2   # K
    
    # Advanced options (Fortran-equivalent only)
    include_localized_states: bool = True
    include_image_potential: bool = True
    expansion_factor: float = 20.0
    write_output: int = 0
    
    # Performance tuning
    enable_warnings: bool = True
    
    def to_original_config(self) -> OriginalConfig:
        """Convert to original implementation config"""
        return OriginalConfig(
            energy_points=self.energy_points,
            energy_range=2.0,
            k_points=self.k_points,
            tip_temperature=self.temperature_tip,
            sample_temperature=self.temperature_sample
        )
    
    def to_fortran_config(self) -> FortranConfig:
        """Convert to Fortran-equivalent config"""
        return FortranConfig(
            energy_points=self.energy_points,
            k_points=self.k_points,
            expansion_factor=self.expansion_factor,
            include_image_potential=self.include_image_potential,
            tip_temperature=self.temperature_tip,
            sample_temperature=self.temperature_sample,
            write_output=self.write_output
        )


class UnifiedTunnelingResult:
    """Unified result structure for tunneling current calculations"""
    
    def __init__(self, 
                 method_used: CalculationMethod,
                 total_current: float,
                 extended_current: float = 0.0,
                 localized_current: float = 0.0,
                 band_contributions: Optional[Dict] = None,
                 localized_states_count: Optional[List[int]] = None,
                 calculation_time: float = 0.0,
                 metadata: Optional[Dict] = None):
        
        self.method_used = method_used
        self.total_current = total_current
        self.extended_current = extended_current
        self.localized_current = localized_current
        self.band_contributions = band_contributions or {}
        self.localized_states_count = localized_states_count or []
        self.calculation_time = calculation_time
        self.metadata = metadata or {}
    
    def get_summary(self) -> Dict:
        """Get a summary of the calculation results"""
        return {
            'method_used': self.method_used.value,
            'total_current_A': self.total_current,
            'extended_current_A': self.extended_current,
            'localized_current_A': self.localized_current,
            'has_band_breakdown': bool(self.band_contributions),
            'has_localized_states': bool(self.localized_states_count),
            'calculation_time_s': self.calculation_time,
            'accuracy_level': self.metadata.get('accuracy_level', 'unknown')
        }
    
    def __str__(self) -> str:
        return (f"TunnelingResult(method={self.method_used.value}, "
                f"total_current={self.total_current:.2e} A)")
    
    def __repr__(self) -> str:
        return self.__str__()


class UnifiedTunnelingCurrentCalculator:
    """
    Unified tunneling current calculator with automatic method selection
    
    This class provides a single interface for tunneling current calculations,
    automatically choosing between simplified and Fortran-equivalent methods
    based on the required accuracy level and input parameters.
    """
    
    def __init__(self, config: Optional[UnifiedTunnelingConfig] = None):
        self.config = config or UnifiedTunnelingConfig()
        
        # Initialize both calculators
        self.original_calculator = OriginalCalculator(self.config.to_original_config())
        self.fortran_calculator = FortranEquivalentTunnelingCalculator(self.config.to_fortran_config())
        
        # Performance tracking
        self.calculation_count = 0
        self.total_calculation_time = 0.0
    
    def calculate_current(self,
                         material: SemiconductorMaterial,
                         bias_voltage: float,
                         fermi_level: float,
                         geometry_data: Dict,
                         method: Optional[CalculationMethod] = None) -> UnifiedTunnelingResult:
        """
        Calculate tunneling current with unified interface
        
        Args:
            material: Semiconductor material parameters
            bias_voltage: Applied bias voltage [V]
            fermi_level: Sample Fermi level [eV]
            geometry_data: Geometry and potential data (format depends on method)
            method: Override automatic method selection
            
        Returns:
            UnifiedTunnelingResult with complete calculation results
        """
        import time
        start_time = time.time()
        
        # Determine calculation method
        chosen_method = method or self._choose_method(geometry_data)
        
        if self.config.enable_warnings:
            self._validate_inputs(material, bias_voltage, fermi_level, geometry_data)
        
        # Perform calculation based on chosen method
        if chosen_method == CalculationMethod.SIMPLIFIED:
            result = self._calculate_simplified(material, bias_voltage, fermi_level, geometry_data)
        elif chosen_method == CalculationMethod.FORTRAN_EQUIVALENT:
            result = self._calculate_fortran_equivalent(material, bias_voltage, fermi_level, geometry_data)
        else:
            raise ValueError(f"Unsupported calculation method: {chosen_method}")
        
        # Update performance tracking
        calculation_time = time.time() - start_time
        self.calculation_count += 1
        self.total_calculation_time += calculation_time
        
        # Set metadata
        result.calculation_time = calculation_time
        result.metadata.update({
            'accuracy_level': self.config.accuracy_level.value,
            'calculation_number': self.calculation_count,
            'config_used': self.config
        })
        
        return result
    
    def _choose_method(self, geometry_data: Dict) -> CalculationMethod:
        """Automatically choose calculation method based on requirements"""
        
        if self.config.method != CalculationMethod.AUTO:
            return self.config.method
        
        # Choose based on accuracy level
        if self.config.accuracy_level == AccuracyLevel.FAST:
            return CalculationMethod.SIMPLIFIED
        elif self.config.accuracy_level in [AccuracyLevel.HIGH_ACCURACY, AccuracyLevel.RESEARCH_GRADE]:
            return CalculationMethod.FORTRAN_EQUIVALENT
        elif self.config.accuracy_level == AccuracyLevel.BALANCED:
            # Choose based on problem size and requirements
            has_complex_geometry = len(geometry_data.get('z_positions', [])) > 50
            needs_localized_states = self.config.include_localized_states
            
            if has_complex_geometry or needs_localized_states:
                return CalculationMethod.FORTRAN_EQUIVALENT
            else:
                return CalculationMethod.SIMPLIFIED
        
        # Default to Fortran-equivalent for safety
        return CalculationMethod.FORTRAN_EQUIVALENT
    
    def _calculate_simplified(self,
                             material: SemiconductorMaterial,
                             bias_voltage: float,
                             fermi_level: float,
                             geometry_data: Dict) -> UnifiedTunnelingResult:
        """Calculate using simplified original method"""
        
        # Extract geometry data for original format
        potential_profile = geometry_data.get('potential_profile')
        z_grid = geometry_data.get('z_grid')
        tip_position = geometry_data.get('tip_position', 0.0)
        sample_position = geometry_data.get('sample_position', 1.0)
        
        if potential_profile is None or z_grid is None:
            raise ValueError("Simplified method requires 'potential_profile' and 'z_grid' in geometry_data")
        
        # Calculate using original implementation
        original_result = self.original_calculator.calculate_total_current(
            material, bias_voltage, fermi_level, potential_profile,
            z_grid, tip_position, sample_position
        )
        
        return UnifiedTunnelingResult(
            method_used=CalculationMethod.SIMPLIFIED,
            total_current=original_result['total_current'],
            extended_current=original_result['total_current'],  # All current is "extended" in simplified
            localized_current=0.0,  # No localized states in simplified
            band_contributions=original_result.get('current_components', {}),
            metadata={'original_result': original_result}
        )
    
    def _calculate_fortran_equivalent(self,
                                    material: SemiconductorMaterial,
                                    bias_voltage: float,
                                    fermi_level: float,
                                    geometry_data: Dict) -> UnifiedTunnelingResult:
        """Calculate using Fortran-equivalent method"""
        
        # Extract geometry data for Fortran format
        potential_3d = geometry_data.get('potential_3d')
        z_positions = geometry_data.get('z_positions')
        vacuum_barrier = geometry_data.get('vacuum_barrier')
        separation = geometry_data.get('separation', 1.0)
        tip_fermi_level = geometry_data.get('tip_fermi_level', 4.5)
        
        if any(x is None for x in [potential_3d, z_positions, vacuum_barrier]):
            raise ValueError("Fortran-equivalent method requires 'potential_3d', 'z_positions', "
                           "and 'vacuum_barrier' in geometry_data")
        
        # Calculate using Fortran-equivalent implementation
        fortran_result = self.fortran_calculator.calculate_tunneling_current(
            material, potential_3d, z_positions, vacuum_barrier,
            separation, bias_voltage, fermi_level, tip_fermi_level
        )
        
        return UnifiedTunnelingResult(
            method_used=CalculationMethod.FORTRAN_EQUIVALENT,
            total_current=fortran_result['total_current'],
            extended_current=fortran_result['total_extended'],
            localized_current=fortran_result['total_localized'],
            band_contributions=fortran_result['band_breakdown'],
            localized_states_count=fortran_result['localized_states_count'],
            metadata={'fortran_result': fortran_result}
        )
    
    def _validate_inputs(self,
                        material: SemiconductorMaterial,
                        bias_voltage: float,
                        fermi_level: float,
                        geometry_data: Dict):
        """Validate input parameters and issue warnings if needed"""
        
        # Check for reasonable parameter ranges
        if abs(bias_voltage) > 10.0:
            warnings.warn(f"Large bias voltage ({bias_voltage} V) may lead to numerical issues")
        
        if abs(fermi_level) > 5.0:
            warnings.warn(f"Unusual Fermi level ({fermi_level} eV) - please verify")
        
        # Check geometry data completeness
        if self.config.accuracy_level in [AccuracyLevel.HIGH_ACCURACY, AccuracyLevel.RESEARCH_GRADE]:
            required_fortran_keys = ['potential_3d', 'z_positions', 'vacuum_barrier']
            missing_keys = [key for key in required_fortran_keys if key not in geometry_data]
            if missing_keys:
                warnings.warn(f"High accuracy requested but missing geometry data: {missing_keys}. "
                            "Consider using simplified method or providing complete data.")
    
    def get_performance_stats(self) -> Dict:
        """Get performance statistics"""
        avg_time = self.total_calculation_time / max(self.calculation_count, 1)
        return {
            'total_calculations': self.calculation_count,
            'total_time_s': self.total_calculation_time,
            'average_time_s': avg_time,
            'calculations_per_second': 1.0 / avg_time if avg_time > 0 else 0
        }


# Convenience functions for direct usage

def calculate_tunneling_current_simple(material: SemiconductorMaterial,
                                     bias_voltage: float,
                                     separation: float = 1.0,
                                     accuracy: AccuracyLevel = AccuracyLevel.BALANCED) -> float:
    """
    Simple interface for basic tunneling current calculation
    
    Args:
        material: Semiconductor material
        bias_voltage: Applied bias [V]
        separation: Tip-sample separation [nm]
        accuracy: Required accuracy level
        
    Returns:
        Tunneling current [A]
    """
    from .charge_density import ChargeDensityCalculator
    
    # Set up basic geometry
    z_grid = np.linspace(0, separation, 100)
    potential_profile = np.zeros_like(z_grid)  # Simplified flat potential
    
    # Find Fermi level
    charge_calc = ChargeDensityCalculator()
    fermi_level = charge_calc.find_equilibrium_fermi_level(material)
    
    # Configure calculation
    config = UnifiedTunnelingConfig(accuracy_level=accuracy)
    calculator = UnifiedTunnelingCurrentCalculator(config)
    
    # Prepare geometry data
    geometry_data = {
        'potential_profile': potential_profile,
        'z_grid': z_grid,
        'tip_position': 0.0,
        'sample_position': separation
    }
    
    # Calculate
    result = calculator.calculate_current(material, bias_voltage, fermi_level, geometry_data)
    return result.total_current


def calculate_tunneling_current_advanced(material: SemiconductorMaterial,
                                       potential_3d: np.ndarray,
                                       z_positions: np.ndarray,
                                       vacuum_barrier: np.ndarray,
                                       bias_voltage: float,
                                       fermi_level: float,
                                       separation: float = 1.0,
                                       tip_fermi_level: float = 4.5,
                                       config: Optional[UnifiedTunnelingConfig] = None) -> UnifiedTunnelingResult:
    """
    Advanced interface for complete tunneling current calculation
    
    Args:
        material: Semiconductor material
        potential_3d: 3D potential along central axis [eV]
        z_positions: Z positions in semiconductor [nm]
        vacuum_barrier: Vacuum barrier potential [eV]
        bias_voltage: Applied bias [V]
        fermi_level: Sample Fermi level [eV]
        separation: Tip-sample separation [nm]
        tip_fermi_level: Tip Fermi level [eV]
        config: Optional configuration
        
    Returns:
        Complete UnifiedTunnelingResult
    """
    
    # Use high accuracy by default for advanced interface
    if config is None:
        config = UnifiedTunnelingConfig(
            method=CalculationMethod.FORTRAN_EQUIVALENT,
            accuracy_level=AccuracyLevel.HIGH_ACCURACY
        )
    
    calculator = UnifiedTunnelingCurrentCalculator(config)
    
    # Prepare geometry data
    geometry_data = {
        'potential_3d': potential_3d,
        'z_positions': z_positions,
        'vacuum_barrier': vacuum_barrier,
        'separation': separation,
        'tip_fermi_level': tip_fermi_level
    }
    
    return calculator.calculate_current(material, bias_voltage, fermi_level, geometry_data)


if __name__ == "__main__":
    # Demo usage
    print("Unified Tunneling Current Calculator")
    print("=" * 40)
    print("Author: odindino")
    print("Date: 2025-06-11")
    print()
    print("Features:")
    print("• Unified API for both simplified and Fortran-equivalent methods")
    print("• Automatic method selection based on accuracy requirements")
    print("• Performance tracking and validation")
    print("• Simple and advanced interfaces")
    print("• Complete backward compatibility")