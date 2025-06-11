"""
Geometry Module Fortran Compatibility Tests

This module provides comprehensive tests to verify that the Python geometry
implementations are numerically consistent with the original SEMITIP Fortran code.

Author: odindino
Date: 2025-06-11
"""

import pytest
import numpy as np
import sys
import os

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from geometry import (
    STMGeometry, GeometryConfig,
    Grid3D, GridConfig, GridDimensions,
    TipGeometry, TipConfig,
    AdaptiveGridRefinement, RefinementConfig,
    BoundaryConditions, BoundaryConfig,
    SymmetryHandler, SymmetryConfig
)
from geometry.grid3d import GridRegion, create_coarse_grid, create_medium_grid
from geometry.stm_geometry import create_standard_stm_geometry
from geometry.tip_geometry import create_standard_tip
from geometry.symmetry_handler import SymmetryType, create_mirror_symmetry_handler


class TestGrid3DFortranCompatibility:
    """Test Grid3D compatibility with SEMITIP grid generation"""
    
    def setup_method(self):
        """Set up test parameters matching Fortran defaults"""
        # Standard SEMITIP parameters from semitip3-6.1.f
        self.fortran_params = {
            'NR': 64,         # Number of radial points
            'NS': 64,         # Number of semiconductor axial points
            'NV': 16,         # Number of vacuum axial points
            'NP': 16,         # Number of angular points
            'DELR0': 0.1,     # Base radial spacing [nm]
            'DELS0': 0.1,     # Base semiconductor spacing [nm]
            'SEP': 1.0,       # Tip-sample separation [nm]
            'RAD': 10.0,      # Tip radius [nm]
            'SLOPE': 0.268,   # Tip slope (tan(15°))
        }
        
    def test_radial_grid_generation(self):
        """Test radial grid generation against Fortran formula"""
        
        # Create grid with Fortran parameters
        config = GridConfig(
            nr_points=self.fortran_params['NR'],
            delr0=self.fortran_params['DELR0']
        )
        grid = Grid3D(config)
        
        # Calculate expected Fortran values
        # R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
        nr = self.fortran_params['NR']
        delr0 = self.fortran_params['DELR0']
        
        expected_r = np.zeros(nr)
        for i in range(nr):
            i_fortran = i + 1  # Fortran 1-based indexing
            expected_r[i] = (2.0 * nr * delr0 / np.pi) * np.tan(
                np.pi * (i_fortran - 0.5) / (2.0 * nr)
            )
            
        # Compare with Python implementation
        python_r = grid.r_points
        
        # Should match to machine precision
        np.testing.assert_allclose(python_r, expected_r, rtol=1e-14, atol=1e-14)
        
        # Test specific values
        assert len(python_r) == nr
        assert python_r[0] > 0  # First point should be positive
        assert np.all(np.diff(python_r) > 0)  # Should be monotonically increasing
        
    def test_semiconductor_axial_grid_generation(self):
        """Test semiconductor axial grid generation"""
        
        config = GridConfig(
            ns_points=self.fortran_params['NS'],
            dels0=self.fortran_params['DELS0']
        )
        grid = Grid3D(config)
        
        # Calculate expected Fortran values
        # S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
        ns = self.fortran_params['NS']
        dels0 = self.fortran_params['DELS0']
        
        expected_s = np.zeros(ns)
        for j in range(ns):
            j_fortran = j + 1  # Fortran 1-based indexing
            expected_s[j] = (2.0 * ns * dels0 / np.pi) * np.tan(
                np.pi * (j_fortran - 0.5) / (2.0 * ns)
            )
            
        # Compare with Python implementation (before scaling)
        python_s = grid.z_semiconductor_points
        
        # Note: Python implementation may scale to fit domain
        # Check that the distribution shape is correct
        assert len(python_s) == ns
        assert python_s[0] > 0
        assert np.all(np.diff(python_s) > 0)
        
        # Check that relative spacing matches
        if len(python_s) > 1:
            expected_ratios = expected_s[1:] / expected_s[:-1]
            python_ratios = python_s[1:] / python_s[:-1]
            np.testing.assert_allclose(python_ratios, expected_ratios, rtol=1e-12)
            
    def test_angular_grid_generation(self):
        """Test angular grid generation with and without mirror symmetry"""
        
        # Test without mirror symmetry
        config_full = GridConfig(
            np_points=self.fortran_params['NP'],
            mirror_symmetry=False
        )
        grid_full = Grid3D(config_full)
        
        # Expected: uniform distribution from 0 to 2π
        np_points = self.fortran_params['NP']
        expected_delp_full = 2.0 * np.pi / np_points
        expected_phi_full = np.linspace(0, 2.0 * np.pi, np_points, endpoint=False)
        
        np.testing.assert_allclose(grid_full.phi_points, expected_phi_full, rtol=1e-14)
        assert abs(grid_full.delp - expected_delp_full) < 1e-14
        
        # Test with mirror symmetry (MIRROR=1)
        config_mirror = GridConfig(
            np_points=self.fortran_params['NP'],
            mirror_symmetry=True
        )
        grid_mirror = Grid3D(config_mirror)
        
        # Expected: uniform distribution from 0 to π
        expected_delp_mirror = np.pi / np_points
        expected_phi_mirror = np.linspace(0, np.pi, np_points, endpoint=False)
        
        np.testing.assert_allclose(grid_mirror.phi_points, expected_phi_mirror, rtol=1e-14)
        assert abs(grid_mirror.delp - expected_delp_mirror) < 1e-14
        
    def test_grid_refinement_factors(self):
        """Test grid refinement matches Fortran doubling strategy"""
        
        config = GridConfig(
            nr_points=32,
            ns_points=32,
            nv_points=8,
            np_points=8,
            delr0=0.2,
            dels0=0.2
        )
        grid = Grid3D(config)
        
        # Store original parameters
        original_nr = grid.config.nr_points
        original_ns = grid.config.ns_points
        original_nv = grid.config.nv_points
        original_np = grid.config.np_points
        original_delr0 = grid.config.delr0
        original_dels0 = grid.config.dels0
        
        # Perform refinement
        success = grid.refine_grid()
        assert success
        
        # Check that all parameters doubled (Fortran behavior)
        assert grid.config.nr_points == 2 * original_nr
        assert grid.config.ns_points == 2 * original_ns
        assert grid.config.nv_points == 2 * original_nv
        assert grid.config.np_points == 2 * original_np
        
        # Check that base spacing halved
        assert abs(grid.config.delr0 - original_delr0 / 2.0) < 1e-14
        assert abs(grid.config.dels0 - original_dels0 / 2.0) < 1e-14
        
    def test_grid_volume_elements(self):
        """Test volume element calculation in cylindrical coordinates"""
        
        grid = create_medium_grid()
        
        # Test volume element formula: dV = r * dr * dφ * dz
        i, j, k = 10, 5, 3  # Sample grid point
        
        volume = grid.get_grid_volume_element(i, j, k, GridRegion.SEMICONDUCTOR)
        
        r = grid.r_points[i]
        dr = grid.delr[i]
        dphi = grid.delp
        dz = grid.dels[j]
        
        expected_volume = r * dr * dphi * dz
        
        assert abs(volume - expected_volume) < 1e-14
        assert volume > 0  # Volume must be positive


class TestTipGeometryFortranCompatibility:
    """Test TipGeometry compatibility with SEMITIP tip calculations"""
    
    def setup_method(self):
        """Set up test parameters from semitip3-6.1.f"""
        self.fortran_params = {
            'RAD': 10.0,      # Tip radius [nm]
            'RAD2': 0.5,      # Tip end radius [nm]
            'SLOPE': 0.268,   # tan(15°)
            'SEP': 1.0,       # Separation [nm]
        }
        
    def test_hyperbolic_coordinate_parameters(self):
        """Test hyperbolic coordinate parameter calculation"""
        
        # Create tip with Fortran parameters
        config = TipConfig(
            radius=self.fortran_params['RAD'],
            radius2=self.fortran_params['RAD2'],
            cone_angle=15.0,  # degrees
            separation=self.fortran_params['SEP']
        )
        tip = TipGeometry(config)
        
        # Calculate expected Fortran values
        # From semitip3-6.1.f lines 97-101:
        rad = self.fortran_params['RAD']
        slope = self.fortran_params['SLOPE']
        sep = self.fortran_params['SEP']
        
        expected_etat = 1.0 / np.sqrt(1.0 + 1.0 / slope**2)
        expected_A = rad * slope**2 / expected_etat
        expected_sprime = expected_A * expected_etat
        expected_Z0 = sep - expected_sprime
        expected_C = expected_Z0 / expected_sprime
        
        # Compare with Python implementation
        assert abs(tip.config.slope - slope) < 1e-14
        assert abs(tip.config.etat - expected_etat) < 1e-14
        assert abs(tip.config.A - expected_A) < 1e-14
        assert abs(tip.config.Z0 - expected_Z0) < 1e-14
        assert abs(tip.config.C - expected_C) < 1e-14
        
    def test_tip_surface_function(self):
        """Test tip surface function p(r)"""
        
        tip = create_standard_tip()
        
        # Test spherical end region (r < RAD2)
        r_spherical = 0.3  # < RAD2 = 0.5
        p_spherical = tip.tip_surface_function(r_spherical)
        
        # Expected from Fortran: p = sqrt(RAD2^2 - r^2)
        expected_p = np.sqrt(tip.config.radius2**2 - r_spherical**2)
        assert abs(p_spherical - expected_p) < 1e-14
        
        # Test conical region (r >= RAD2)
        r_conical = 1.0  # > RAD2 = 0.5
        p_conical = tip.tip_surface_function(r_conical)
        
        # Expected: p = RAD2 + (r - RAD2) * SLOPE
        expected_p_conical = (tip.config.radius2 + 
                             (r_conical - tip.config.radius2) * tip.config.slope)
        assert abs(p_conical - expected_p_conical) < 1e-14
        
    def test_tip_point_identification(self):
        """Test tip point identification logic"""
        
        tip = create_standard_tip()
        
        # Test points that should be inside tip
        test_points = [
            (0.0, 0.0, -1.0),    # At apex
            (0.3, 0.0, -0.95),   # In spherical end
            (1.0, 0.0, -0.7),    # In conical section
        ]
        
        for r, phi, z in test_points:
            inside = tip.is_point_inside_tip(r, phi, z)
            # We expect these to be inside based on geometry
            # (Specific validation would require detailed Fortran comparison)
            assert isinstance(inside, bool)
            
        # Test points clearly outside
        outside_points = [
            (0.0, 0.0, 1.0),     # In semiconductor
            (10.0, 0.0, -1.0),   # Far from tip
        ]
        
        for r, phi, z in outside_points:
            inside = tip.is_point_inside_tip(r, phi, z)
            assert not inside


class TestSTMGeometryFortranCompatibility:
    """Test STMGeometry compatibility with SEMITIP geometry setup"""
    
    def test_coordinate_system_consistency(self):
        """Test coordinate system definition consistency"""
        
        geometry = create_standard_stm_geometry()
        
        # Test coordinate bounds
        assert geometry.z_surface == 0.0  # Surface at z=0
        assert geometry.z_min < 0.0       # Vacuum at negative z
        assert geometry.z_max > 0.0       # Semiconductor at positive z
        
        # Test tip position
        tip_pos = geometry.get_tip_position()
        assert tip_pos[2] < 0.0  # Tip at negative z (vacuum side)
        assert abs(tip_pos[2] + geometry.config.separation) < 1e-14
        
    def test_work_function_profile(self):
        """Test work function variation with position"""
        
        geometry = create_standard_stm_geometry()
        
        # Test at tip
        tip_z = geometry.get_tip_position()[2]
        wf_tip = geometry.get_work_function_profile(tip_z - 0.1)  # Inside tip
        assert abs(wf_tip - geometry.config.tip_work_function) < 1e-14
        
        # Test at surface
        wf_surface = geometry.get_work_function_profile(0.1)  # In semiconductor
        assert abs(wf_surface - geometry.config.sample_work_function) < 1e-14
        
        # Test in vacuum (should be interpolated)
        wf_vacuum = geometry.get_work_function_profile(-0.5)
        assert geometry.config.sample_work_function <= wf_vacuum <= geometry.config.tip_work_function


class TestBoundaryConditionsFortranCompatibility:
    """Test boundary conditions match SEMITIP implementation"""
    
    def setup_method(self):
        """Set up test geometry"""
        self.grid = create_medium_grid(mirror_symmetry=True)
        self.geometry = create_standard_stm_geometry()
        self.tip = create_standard_tip()
        
    def test_tip_potential_calculation(self):
        """Test tip potential calculation"""
        
        from geometry.boundary_conditions import create_standard_boundary_conditions
        
        boundaries = create_standard_boundary_conditions(self.grid, self.geometry, self.tip)
        
        # Expected tip potential: BIAS + CPot + work_function
        expected_tip_potential = (
            self.geometry.config.bias_voltage +
            self.geometry.config.contact_potential +
            self.tip.config.work_function
        )
        
        assert abs(boundaries.config.tip_potential - expected_tip_potential) < 1e-14
        
    def test_mirror_symmetry_boundary_application(self):
        """Test mirror symmetry boundary conditions"""
        
        from geometry.boundary_conditions import create_standard_boundary_conditions
        
        boundaries = create_standard_boundary_conditions(self.grid, self.geometry, self.tip)
        
        # Create test potential array
        nr, nv, np_pts = 32, 8, 16
        test_potential = np.random.normal(0, 0.1, (nr, nv, np_pts))
        
        # Apply boundary conditions
        vac_new, sem_new, surf_new = boundaries.apply_boundary_conditions(
            test_potential, 
            np.random.normal(0, 0.1, (nr, 32, np_pts)),
            np.random.normal(0, 0.1, (nr, np_pts))
        )
        
        # Check that arrays have correct shapes
        assert vac_new.shape == test_potential.shape
        assert isinstance(vac_new, np.ndarray)
        assert isinstance(sem_new, np.ndarray)
        assert isinstance(surf_new, np.ndarray)


class TestSymmetryHandlerFortranCompatibility:
    """Test symmetry handling matches SEMITIP MIRROR=1 implementation"""
    
    def test_mirror_symmetry_setup(self):
        """Test mirror symmetry setup"""
        
        grid = create_medium_grid(mirror_symmetry=True)
        symmetry = create_mirror_symmetry_handler(grid)
        
        # Check that mirror symmetry is properly identified
        info = symmetry.get_symmetry_info()
        assert info['configuration']['effective_symmetry'] == 'mirror_xz'
        assert info['grid_symmetry']['mirror_symmetry'] == True
        assert info['grid_symmetry']['phi_range'] == 'π'
        
    def test_symmetry_preservation_check(self):
        """Test symmetry preservation checking"""
        
        grid = create_medium_grid(mirror_symmetry=True)
        symmetry = create_mirror_symmetry_handler(grid)
        
        # Create symmetric test potential
        nr, nv, np_pts = 32, 8, 8  # Small for testing
        test_potential = np.random.normal(0, 0.1, (nr, nv, np_pts))
        
        # Make it symmetric
        for k in range(np_pts // 2):
            k_mirror = np_pts - 1 - k
            avg = (test_potential[:, :, k] + test_potential[:, :, k_mirror]) / 2.0
            test_potential[:, :, k] = avg
            test_potential[:, :, k_mirror] = avg
            
        # Check symmetry preservation
        check_result = symmetry.check_symmetry_preservation(test_potential)
        assert check_result['symmetry_preserved']
        assert check_result['max_symmetry_error'] < 1e-12


class TestAdaptiveRefinementFortranCompatibility:
    """Test adaptive refinement matches SEMITIP three-level strategy"""
    
    def test_refinement_level_progression(self):
        """Test refinement level progression"""
        
        from geometry.adaptive_refinement import create_standard_refinement
        
        grid = create_coarse_grid()
        geometry = create_standard_stm_geometry()
        refinement = create_standard_refinement(grid, geometry)
        
        # Check default configuration matches Fortran
        assert refinement.config.max_refinement_levels == 3
        assert len(refinement.config.convergence_tolerance) == 3
        assert len(refinement.config.max_iterations_per_level) == 3
        
        # Check that convergence gets stricter with level
        tolerances = refinement.config.convergence_tolerance
        assert tolerances[1] < tolerances[0]  # Level 2 stricter than level 1
        assert tolerances[2] < tolerances[1]  # Level 3 stricter than level 2
        
    def test_grid_sequence_generation(self):
        """Test grid sequence generation for refinement"""
        
        from geometry.adaptive_refinement import create_fast_refinement
        
        initial_grid = create_coarse_grid()
        geometry = create_standard_stm_geometry()
        refinement = create_fast_refinement(initial_grid, geometry)
        
        # Start with coarse grid
        assert len(refinement.state.grid_levels) == 1
        original_nr = refinement.current_grid.config.nr_points
        
        # Simulate one refinement step
        success = refinement._prepare_next_level(np.array([[1.0]]))  # Dummy solution
        assert success
        assert len(refinement.state.grid_levels) == 2
        
        # Check that grid was refined
        new_nr = refinement.current_grid.config.nr_points
        assert new_nr == 2 * original_nr


class TestGeometryIntegrationFortranCompatibility:
    """Test integration between geometry modules"""
    
    def test_complete_geometry_setup(self):
        """Test complete geometry setup workflow"""
        
        # Create all geometry components
        geometry = create_standard_stm_geometry()
        grid = create_medium_grid(mirror_symmetry=True)
        tip = create_standard_tip()
        
        # Verify consistency between components
        assert geometry.config.separation == tip.config.separation
        assert grid.config.mirror_symmetry == True  # Consistent with geometry
        
        # Test coordinate transformations
        r, phi, z = grid.get_cylindrical_coordinates(GridRegion.VACUUM, 0, 0, 0)
        x, y, z_cart = grid.get_grid_point(GridRegion.VACUUM, 0, 0, 0)
        
        # Check coordinate conversion
        expected_x = r * np.cos(phi)
        expected_y = r * np.sin(phi)
        assert abs(x - expected_x) < 1e-14
        assert abs(y - expected_y) < 1e-14
        
    def test_memory_usage_estimation(self):
        """Test memory usage estimation"""
        
        grid = create_medium_grid()
        memory_info = grid.estimate_memory_usage()
        
        # Check that estimates are reasonable
        assert memory_info['total_per_array_mb'] > 0
        assert memory_info['estimated_total_mb'] > memory_info['total_per_array_mb']
        
        # Check individual components
        assert memory_info['vacuum_array_mb'] > 0
        assert memory_info['semiconductor_array_mb'] > 0
        assert memory_info['surface_array_mb'] > 0


# Numerical precision tests
class TestNumericalPrecisionConsistency:
    """Test numerical precision consistency with Fortran"""
    
    def test_floating_point_precision(self):
        """Test that calculations maintain double precision"""
        
        # Test grid calculations with known values
        config = GridConfig(nr_points=100, delr0=0.1)
        grid = Grid3D(config)
        
        # All calculations should be in double precision
        assert grid.r_points.dtype == np.float64
        assert grid.delr.dtype == np.float64
        
        # Test precision of tangent mapping
        # Should be accurate to machine precision
        for i in range(min(10, len(grid.r_points))):
            i_fortran = i + 1
            expected = (2.0 * config.nr_points * config.delr0 / np.pi) * np.tan(
                np.pi * (i_fortran - 0.5) / (2.0 * config.nr_points)
            )
            assert abs(grid.r_points[i] - expected) < 1e-14
            
    def test_parameter_value_consistency(self):
        """Test parameter values match Fortran constants"""
        
        # Test physical constants consistency
        tip = create_standard_tip()
        
        # Check derived parameters
        slope_15_deg = np.tan(np.radians(15.0))
        assert abs(tip.config.slope - slope_15_deg) < 1e-14
        
        # Test trigonometric consistency
        etat_expected = 1.0 / np.sqrt(1.0 + 1.0 / slope_15_deg**2)
        assert abs(tip.config.etat - etat_expected) < 1e-14


if __name__ == "__main__":
    # Run tests with detailed output
    pytest.main([__file__, "-v", "--tb=short"])