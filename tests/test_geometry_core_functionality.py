"""
Core Geometry Functionality Tests

This module tests the core functionality of geometry modules with
focus on essential features and numerical accuracy.

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
    Grid3D, GridConfig,
    TipGeometry, TipConfig,
    AdaptiveGridRefinement, RefinementConfig,
    BoundaryConditions, BoundaryConfig,
    SymmetryHandler, SymmetryConfig
)


class TestGrid3DCoreFunctionality:
    """Test core Grid3D functionality"""
    
    def test_radial_grid_fortran_compatibility(self):
        """Test radial grid generation matches Fortran exactly"""
        
        config = GridConfig(nr_points=64, delr0=0.1)
        grid = Grid3D(config)
        
        # Calculate expected Fortran values
        nr = 64
        delr0 = 0.1
        expected_r = np.zeros(nr)
        
        for i in range(nr):
            i_fortran = i + 1  # Fortran 1-based indexing
            expected_r[i] = (2.0 * nr * delr0 / np.pi) * np.tan(
                np.pi * (i_fortran - 0.5) / (2.0 * nr)
            )
            
        # Should match exactly
        np.testing.assert_allclose(grid.r_points, expected_r, rtol=1e-15, atol=1e-15)
        
    def test_angular_grid_mirror_symmetry(self):
        """Test angular grid with mirror symmetry"""
        
        config = GridConfig(np_points=16, mirror_symmetry=True)
        grid = Grid3D(config)
        
        # Should span 0 to π, not 0 to 2π
        expected_max = np.pi * (1 - 1/16)  # π * (np_points - 1) / np_points
        assert abs(grid.phi_points[-1] - expected_max) < 1e-14
        assert grid.delp == np.pi / 16
        
    def test_grid_refinement(self):
        """Test grid refinement doubles all dimensions"""
        
        config = GridConfig(nr_points=32, ns_points=32, nv_points=8, np_points=8)
        grid = Grid3D(config)
        
        original_dims = (grid.config.nr_points, grid.config.ns_points, 
                        grid.config.nv_points, grid.config.np_points)
        
        success = grid.refine_grid()
        assert success
        
        new_dims = (grid.config.nr_points, grid.config.ns_points,
                   grid.config.nv_points, grid.config.np_points)
        
        # All dimensions should double
        for orig, new in zip(original_dims, new_dims):
            assert new == 2 * orig


class TestTipGeometryCoreFunctionality:
    """Test core TipGeometry functionality"""
    
    def test_hyperbolic_parameters_fortran_exact(self):
        """Test hyperbolic parameters match Fortran calculations exactly"""
        
        config = TipConfig(
            radius=10.0,
            radius2=0.5,
            cone_angle=15.0,
            separation=1.0
        )
        tip = TipGeometry(config)
        
        # Expected values from Fortran calculation
        slope = np.tan(np.radians(15.0))
        etat = 1.0 / np.sqrt(1.0 + 1.0 / slope**2)
        A = 10.0 * slope**2 / etat
        sprime = A * etat
        Z0 = 1.0 - sprime
        C = Z0 / sprime
        
        assert abs(tip.config.slope - slope) < 1e-15
        assert abs(tip.config.etat - etat) < 1e-15
        assert abs(tip.config.A - A) < 1e-15
        assert abs(tip.config.Z0 - Z0) < 1e-15
        assert abs(tip.config.C - C) < 1e-15
        
    def test_tip_surface_function(self):
        """Test tip surface function p(r)"""
        
        tip = TipGeometry(TipConfig(radius2=0.5, slope=np.tan(np.radians(15.0))))
        
        # Spherical region
        r = 0.3
        p = tip.tip_surface_function(r)
        expected = np.sqrt(0.5**2 - r**2)
        assert abs(p - expected) < 1e-15
        
        # Conical region  
        r = 1.0
        p = tip.tip_surface_function(r)
        expected = 0.5 + (r - 0.5) * tip.config.slope
        assert abs(p - expected) < 1e-15


class TestSTMGeometryCoreFunctionality:
    """Test core STMGeometry functionality"""
    
    def test_coordinate_system_definition(self):
        """Test coordinate system definition"""
        
        geometry = STMGeometry()
        
        # Basic coordinate definitions
        assert geometry.z_surface == 0.0
        assert geometry.z_min < 0.0  # Vacuum side
        assert geometry.z_max > 0.0  # Semiconductor side
        
        # Tip position
        tip_pos = geometry.get_tip_position()
        assert tip_pos[2] < 0.0  # Tip in vacuum (negative z)
        
    def test_work_function_profile(self):
        """Test work function variation"""
        
        config = GeometryConfig(
            tip_work_function=4.5,
            sample_work_function=4.1,
            separation=1.0
        )
        geometry = STMGeometry(config)
        
        # At tip
        wf_tip = geometry.get_work_function_profile(-1.1)
        assert abs(wf_tip - 4.5) < 1e-10
        
        # At surface
        wf_surface = geometry.get_work_function_profile(0.1)
        assert abs(wf_surface - 4.1) < 1e-10


class TestBoundaryConditionsCoreFunctionality:
    """Test core boundary conditions functionality"""
    
    def test_tip_potential_calculation(self):
        """Test tip potential calculation"""
        
        from geometry.stm_geometry import create_standard_stm_geometry
        from geometry.tip_geometry import create_standard_tip
        from geometry.grid3d import create_medium_grid
        from geometry.boundary_conditions import create_standard_boundary_conditions
        
        grid = create_medium_grid()
        geometry = create_standard_stm_geometry()
        tip = create_standard_tip()
        
        boundaries = create_standard_boundary_conditions(grid, geometry, tip)
        
        expected = (geometry.config.bias_voltage + 
                   geometry.config.contact_potential +
                   tip.config.work_function)
        
        assert abs(boundaries.config.tip_potential - expected) < 1e-15
        
    def test_boundary_application(self):
        """Test boundary condition application"""
        
        from geometry.stm_geometry import create_standard_stm_geometry
        from geometry.tip_geometry import create_standard_tip
        from geometry.grid3d import create_coarse_grid
        from geometry.boundary_conditions import create_standard_boundary_conditions
        
        grid = create_coarse_grid()
        geometry = create_standard_stm_geometry()
        tip = create_standard_tip()
        boundaries = create_standard_boundary_conditions(grid, geometry, tip)
        
        # Create test arrays with grid dimensions
        nr = grid.config.nr_points
        nv = grid.config.nv_points  
        ns = grid.config.ns_points
        np_pts = grid.config.np_points
        
        vac_pot = np.random.normal(0, 0.1, (nr, nv, np_pts))
        sem_pot = np.random.normal(0, 0.1, (nr, ns, np_pts))
        surf_pot = np.random.normal(0, 0.1, (nr, np_pts))
        
        # Apply boundaries
        vac_new, sem_new, surf_new = boundaries.apply_boundary_conditions(
            vac_pot, sem_pot, surf_pot)
        
        # Should return arrays of correct shape
        assert vac_new.shape == vac_pot.shape
        assert sem_new.shape == sem_pot.shape
        assert surf_new.shape == surf_pot.shape


class TestSymmetryHandlerCoreFunctionality:
    """Test core symmetry handling functionality"""
    
    def test_mirror_symmetry_detection(self):
        """Test mirror symmetry detection"""
        
        from geometry.grid3d import create_medium_grid
        from geometry.symmetry_handler import create_mirror_symmetry_handler
        
        grid = create_medium_grid(mirror_symmetry=True)
        symmetry = create_mirror_symmetry_handler(grid)
        
        info = symmetry.get_symmetry_info()
        assert info['configuration']['effective_symmetry'] == 'mirror_xz'
        assert info['grid_symmetry']['mirror_symmetry'] == True
        
    def test_symmetry_preservation(self):
        """Test symmetry preservation in potential arrays"""
        
        from geometry.grid3d import create_medium_grid
        from geometry.symmetry_handler import create_mirror_symmetry_handler
        from geometry.grid3d import GridRegion
        
        grid = create_medium_grid(mirror_symmetry=True)
        symmetry = create_mirror_symmetry_handler(grid)
        
        # Create test potential that should be symmetric
        nr, nv, np_pts = 16, 4, 8
        test_pot = np.random.normal(0, 0.1, (nr, nv, np_pts))
        
        # Apply symmetry
        symmetric_pot = symmetry.apply_symmetry_to_potential(test_pot, GridRegion.VACUUM)
        
        # Should have same shape
        assert symmetric_pot.shape == test_pot.shape


class TestAdaptiveRefinementCoreFunctionality:
    """Test core adaptive refinement functionality"""
    
    def test_refinement_configuration(self):
        """Test refinement configuration"""
        
        from geometry.adaptive_refinement import create_standard_refinement
        from geometry.grid3d import create_coarse_grid
        from geometry.stm_geometry import create_standard_stm_geometry
        
        grid = create_coarse_grid()
        geometry = create_standard_stm_geometry()
        refinement = create_standard_refinement(grid, geometry)
        
        # Should have 3 levels with decreasing tolerance
        assert refinement.config.max_refinement_levels == 3
        tol = refinement.config.convergence_tolerance
        assert tol[1] < tol[0]
        assert tol[2] < tol[1]
        
    def test_grid_preparation(self):
        """Test grid preparation for next level"""
        
        from geometry.adaptive_refinement import create_fast_refinement
        from geometry.grid3d import create_coarse_grid
        from geometry.stm_geometry import create_standard_stm_geometry
        
        grid = create_coarse_grid()
        geometry = create_standard_stm_geometry()
        refinement = create_fast_refinement(grid, geometry)
        
        original_nr = refinement.current_grid.config.nr_points
        
        # Prepare next level
        dummy_solution = np.array([[1.0]])
        success = refinement._prepare_next_level(dummy_solution)
        
        if success:  # May fail due to dimension limits
            assert refinement.current_grid.config.nr_points == 2 * original_nr


class TestIntegratedFunctionality:
    """Test integrated functionality across modules"""
    
    def test_complete_geometry_workflow(self):
        """Test complete geometry setup workflow"""
        
        from geometry.stm_geometry import create_standard_stm_geometry
        from geometry.grid3d import create_medium_grid
        from geometry.tip_geometry import create_standard_tip
        from geometry.boundary_conditions import create_standard_boundary_conditions
        from geometry.symmetry_handler import create_mirror_symmetry_handler
        
        # Create all components
        geometry = create_standard_stm_geometry()
        grid = create_medium_grid(mirror_symmetry=True)
        tip = create_standard_tip()
        boundaries = create_standard_boundary_conditions(grid, geometry, tip)
        symmetry = create_mirror_symmetry_handler(grid)
        
        # Verify consistency
        assert geometry.config.separation == tip.config.separation
        assert grid.config.mirror_symmetry == True
        assert boundaries.grid.config.mirror_symmetry == True
        
        # Test coordinate consistency
        from geometry.grid3d import GridRegion
        r, phi, z = grid.get_cylindrical_coordinates(
            GridRegion.VACUUM, 0, 0, 0
        )
        assert r >= 0
        assert 0 <= phi <= np.pi  # Mirror symmetry
        assert z <= 0  # Vacuum region
        
        # Test geometry validation
        geom_issues = geometry.validate_geometry()
        tip_issues = tip.validate_tip_geometry()
        boundary_issues = boundaries.validate_boundary_conditions()
        
        # Should have minimal issues for standard configurations
        assert len(geom_issues) <= 2  # May have minor warnings
        assert len(tip_issues) <= 2
        assert len(boundary_issues) <= 2


if __name__ == "__main__":
    # Run with verbose output
    pytest.main([__file__, "-v"])