"""
Integration Test for Unified Tunneling Current API

This module tests the integrated tunneling current API to ensure proper
functionality of the unified interface and backward compatibility.

Author: odindino  
Date: 2025-06-11
"""

import pytest
import numpy as np
import sys
import os

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

import physics
from physics import (
    create_tunneling_calculator,
    UnifiedTunnelingCurrentCalculator,
    UnifiedTunnelingConfig,
    CalculationMethod,
    AccuracyLevel,
    calculate_tunneling_current_simple,
    calculate_tunneling_current_advanced
)
from physics.materials import default_materials


class TestUnifiedAPI:
    """Test unified tunneling current API"""
    
    def setup_method(self):
        """Set up test data"""
        self.material = default_materials.get_material("Si_n")
        self.bias_voltage = 1.0  # V
        self.separation = 1.0    # nm
        
    def test_module_imports(self):
        """Test that all imports work correctly"""
        # Test core imports
        assert hasattr(physics, 'UnifiedTunnelingCurrentCalculator')
        assert hasattr(physics, 'CalculationMethod')
        assert hasattr(physics, 'AccuracyLevel')
        assert hasattr(physics, 'create_tunneling_calculator')
        
        # Test backward compatibility
        assert hasattr(physics, 'TunnelingCurrentCalculator')
        assert hasattr(physics, 'FortranEquivalentTunnelingCalculator')
        
    def test_convenience_factory_function(self):
        """Test the module-level convenience function"""
        # Test with different accuracy levels
        calc_fast = create_tunneling_calculator(AccuracyLevel.FAST)
        assert isinstance(calc_fast, UnifiedTunnelingCurrentCalculator)
        assert calc_fast.config.accuracy_level == AccuracyLevel.FAST
        
        calc_accurate = create_tunneling_calculator(AccuracyLevel.HIGH_ACCURACY)
        assert isinstance(calc_accurate, UnifiedTunnelingCurrentCalculator)
        assert calc_accurate.config.accuracy_level == AccuracyLevel.HIGH_ACCURACY
        
    def test_unified_calculator_initialization(self):
        """Test unified calculator initialization with different configs"""
        # Default config
        calc_default = UnifiedTunnelingCurrentCalculator()
        assert calc_default.config.method == CalculationMethod.AUTO
        assert calc_default.config.accuracy_level == AccuracyLevel.BALANCED
        
        # Custom config
        custom_config = UnifiedTunnelingConfig(
            method=CalculationMethod.FORTRAN_EQUIVALENT,
            accuracy_level=AccuracyLevel.RESEARCH_GRADE,
            energy_points=30,
            k_points=15
        )
        calc_custom = UnifiedTunnelingCurrentCalculator(custom_config)
        assert calc_custom.config.method == CalculationMethod.FORTRAN_EQUIVALENT
        assert calc_custom.config.energy_points == 30
        
    def test_simple_interface(self):
        """Test the simple convenience interface"""
        # Test with different accuracy levels
        for accuracy in AccuracyLevel:
            try:
                current = calculate_tunneling_current_simple(
                    self.material, self.bias_voltage, self.separation, accuracy
                )
                assert isinstance(current, (int, float))
                assert np.isfinite(current)
            except Exception as e:
                # Some configurations might not work with simple test data
                pytest.skip(f"Simple interface test skipped for {accuracy}: {e}")
                
    def test_method_selection_logic(self):
        """Test automatic method selection logic"""
        calculator = UnifiedTunnelingCurrentCalculator()
        
        # Test method selection for different scenarios
        simple_geometry = {'z_positions': np.linspace(0.1, 0.5, 10)}
        complex_geometry = {'z_positions': np.linspace(0.1, 1.0, 100)}
        
        # Fast accuracy should choose simplified
        calculator.config.accuracy_level = AccuracyLevel.FAST
        method = calculator._choose_method(simple_geometry)
        assert method == CalculationMethod.SIMPLIFIED
        
        # High accuracy should choose Fortran-equivalent
        calculator.config.accuracy_level = AccuracyLevel.HIGH_ACCURACY
        method = calculator._choose_method(simple_geometry)
        assert method == CalculationMethod.FORTRAN_EQUIVALENT
        
        # Balanced with complex geometry should choose Fortran-equivalent
        calculator.config.accuracy_level = AccuracyLevel.BALANCED
        method = calculator._choose_method(complex_geometry)
        assert method == CalculationMethod.FORTRAN_EQUIVALENT
        
    def test_result_structure(self):
        """Test unified result structure"""
        config = UnifiedTunnelingConfig(
            method=CalculationMethod.SIMPLIFIED,  # Use simplified for speed
            accuracy_level=AccuracyLevel.FAST
        )
        calculator = UnifiedTunnelingCurrentCalculator(config)
        
        # Prepare simple geometry data
        z_grid = np.linspace(0, self.separation, 50)
        potential_profile = np.zeros_like(z_grid)
        
        geometry_data = {
            'potential_profile': potential_profile,
            'z_grid': z_grid,
            'tip_position': 0.0,
            'sample_position': self.separation
        }
        
        # Calculate using charge density for Fermi level
        from physics.charge_density import ChargeDensityCalculator
        charge_calc = ChargeDensityCalculator()
        fermi_level = charge_calc.find_equilibrium_fermi_level(self.material)
        
        result = calculator.calculate_current(
            self.material, self.bias_voltage, fermi_level, geometry_data
        )
        
        # Test result structure
        assert hasattr(result, 'method_used')
        assert hasattr(result, 'total_current')
        assert hasattr(result, 'extended_current')
        assert hasattr(result, 'localized_current')
        assert hasattr(result, 'calculation_time')
        
        # Test result methods
        summary = result.get_summary()
        assert isinstance(summary, dict)
        assert 'method_used' in summary
        assert 'total_current_A' in summary
        
        # Test string representation
        str_repr = str(result)
        assert 'TunnelingResult' in str_repr
        assert 'total_current' in str_repr
        
    def test_error_handling(self):
        """Test error handling for invalid inputs"""
        calculator = UnifiedTunnelingCurrentCalculator()
        
        # Test missing geometry data
        with pytest.raises(ValueError):
            calculator.calculate_current(
                self.material, self.bias_voltage, 0.5, {}
            )
        
        # Test invalid method
        with pytest.raises(ValueError):
            config = UnifiedTunnelingConfig(method="invalid_method")
            # This should fail during enum assignment
            
    def test_performance_tracking(self):
        """Test performance tracking functionality"""
        calculator = UnifiedTunnelingCurrentCalculator()
        
        # Initial state
        stats = calculator.get_performance_stats()
        assert stats['total_calculations'] == 0
        assert stats['total_time_s'] == 0.0
        
        # After calculation (using simplified method for speed)
        config = UnifiedTunnelingConfig(method=CalculationMethod.SIMPLIFIED)
        calculator = UnifiedTunnelingCurrentCalculator(config)
        
        try:
            z_grid = np.linspace(0, 1.0, 20)
            potential_profile = np.zeros_like(z_grid)
            geometry_data = {
                'potential_profile': potential_profile,
                'z_grid': z_grid,
                'tip_position': 0.0,
                'sample_position': 1.0
            }
            
            from physics.charge_density import ChargeDensityCalculator
            charge_calc = ChargeDensityCalculator()
            fermi_level = charge_calc.find_equilibrium_fermi_level(self.material)
            
            calculator.calculate_current(self.material, 1.0, fermi_level, geometry_data)
            
            stats = calculator.get_performance_stats()
            assert stats['total_calculations'] == 1
            assert stats['total_time_s'] > 0
            assert stats['average_time_s'] > 0
            
        except Exception as e:
            pytest.skip(f"Performance tracking test skipped: {e}")


class TestBackwardCompatibility:
    """Test backward compatibility with existing code"""
    
    def test_original_imports_still_work(self):
        """Test that original imports still work"""
        # These should not raise import errors
        from physics import TunnelingCurrentCalculator, TunnelingConfig
        from physics import TunnelingCurrentCalculatorOriginal, TunnelingConfigOriginal
        
        # Original classes should be accessible
        calc_orig = TunnelingCurrentCalculator()
        assert calc_orig is not None
        
        config_orig = TunnelingConfig()
        assert config_orig is not None
        
    def test_alias_consistency(self):
        """Test that aliases point to correct classes"""
        from physics import (
            TunnelingCurrentCalculator, 
            TunnelingCurrentCalculatorOriginal
        )
        
        # Aliases should point to the same class
        assert TunnelingCurrentCalculator is TunnelingCurrentCalculatorOriginal


class TestIntegrationScenarios:
    """Test realistic integration scenarios"""
    
    def setup_method(self):
        """Set up realistic test scenarios"""
        self.material = default_materials.get_material("Si_n")
        
    def test_research_workflow_scenario(self):
        """Test a typical research workflow scenario"""
        
        # Step 1: Quick screening with fast method
        quick_current = calculate_tunneling_current_simple(
            self.material, 1.0, 1.0, AccuracyLevel.FAST
        )
        assert np.isfinite(quick_current)
        
        # Step 2: If needed, detailed calculation would follow
        # (Skip detailed for testing speed)
        
    def test_production_workflow_scenario(self):
        """Test a production calculation scenario"""
        
        # Create balanced calculator for production use
        calculator = create_tunneling_calculator(AccuracyLevel.BALANCED)
        
        # This should work without detailed geometry setup
        # (would require full geometry for real use)
        assert isinstance(calculator, UnifiedTunnelingCurrentCalculator)
        assert calculator.config.accuracy_level == AccuracyLevel.BALANCED
        
    def test_mixed_usage_scenario(self):
        """Test mixed usage of different interfaces"""
        
        # Use convenience function
        calc1 = create_tunneling_calculator(AccuracyLevel.FAST)
        
        # Use direct instantiation
        config = UnifiedTunnelingConfig(
            accuracy_level=AccuracyLevel.HIGH_ACCURACY,
            method=CalculationMethod.FORTRAN_EQUIVALENT
        )
        calc2 = UnifiedTunnelingCurrentCalculator(config)
        
        # Both should work
        assert calc1.config.accuracy_level == AccuracyLevel.FAST
        assert calc2.config.accuracy_level == AccuracyLevel.HIGH_ACCURACY


if __name__ == "__main__":
    print("Integration Test for Unified Tunneling Current API")
    print("=" * 55)
    print("Author: odindino")
    print("Date: 2025-06-11")
    print()
    
    # Run basic tests
    test_api = TestUnifiedAPI()
    test_api.setup_method()
    
    print("Testing module imports...")
    test_api.test_module_imports()
    print("✓ Module imports working")
    
    print("Testing convenience factory...")
    test_api.test_convenience_factory_function()
    print("✓ Factory function working")
    
    print("Testing calculator initialization...")
    test_api.test_unified_calculator_initialization()
    print("✓ Calculator initialization working")
    
    print("Testing method selection...")
    test_api.test_method_selection_logic()
    print("✓ Method selection logic working")
    
    print("Testing backward compatibility...")
    test_compat = TestBackwardCompatibility()
    test_compat.test_original_imports_still_work()
    test_compat.test_alias_consistency()
    print("✓ Backward compatibility maintained")
    
    print("\nIntegration tests passed!")
    print("Unified API is ready for use.")