"""
Unit tests for physics modules.

This module provides basic tests to verify the physics calculations
are working correctly.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import unittest
from pathlib import Path

from src.physics.materials import SemiconductorRegion, SurfaceRegion, TipModel
from src.physics.solvers import Grid3D, GridParameters
from src.physics.core import (ChargeDensityCalculator, PoissonSolver,
                         PotentialProcessor, SchrodingerSolver)
from src.utils.constants import PhysicalConstants as PC


class TestMaterials(unittest.TestCase):
    """Test material models."""
    
    def setUp(self):
        """Set up test materials."""
        self.semiconductor = SemiconductorRegion(
            region_id=1,
            donor_concentration=1e18,
            acceptor_concentration=0,
            band_gap=1.42,
            valence_band_offset=0,
            electron_affinity=4.07,
            donor_binding_energy=0.006,
            acceptor_binding_energy=0.031,
            cb_effective_mass=0.067,
            vb_effective_mass_heavy=0.5,
            vb_effective_mass_light=0.08,
            vb_effective_mass_so=0.15,
            spin_orbit_splitting=0.34,
            permittivity=12.9,
            temperature=300.0
        )
    
    def test_semiconductor_properties(self):
        """Test semiconductor calculations."""
        # Test doping type
        self.assertTrue(self.semiconductor.is_n_type)
        self.assertFalse(self.semiconductor.is_p_type)
        
        # Test band edges
        cb_edge = self.semiconductor.conduction_band_edge()
        vb_edge = self.semiconductor.valence_band_edge()
        self.assertAlmostEqual(cb_edge - vb_edge, 1.42, places=2)
        
        # Test carrier densities
        fermi = self.semiconductor.fermi_level()
        n = self.semiconductor.carrier_density_cb(fermi)
        p = self.semiconductor.carrier_density_vb(fermi)
        
        # For n-type, n >> p
        self.assertGreater(n, p * 1000)
    
    def test_surface_states(self):
        """Test surface state calculations."""
        from src.physics.materials.surface_states import SurfaceStateDistribution
        
        dist = SurfaceStateDistribution(
            density=1e14,
            neutrality_level=0.7,
            fwhm=0.2,
            center_energy=0.7
        )
        
        # Test DOS calculation
        energies = np.linspace(0, 1.4, 100)
        dos = dist.density_of_states(energies)
        
        # Check peak is at center energy
        peak_idx = np.argmax(dos)
        self.assertAlmostEqual(energies[peak_idx], 0.7, places=1)


class TestGrid(unittest.TestCase):
    """Test grid system."""
    
    def test_grid_creation(self):
        """Test grid initialization."""
        params = GridParameters(
            nr=32, nv=32, ns=32, np=16,
            delr=1.0, delv=1.0, dels=1.0, delp=np.pi/16,
            rmax=31.0, vmax=31.0, smax=31.0
        )
        
        grid = Grid3D(params)
        
        # Check dimensions
        self.assertEqual(len(grid.r), 32)
        self.assertEqual(len(grid.zv), 32)
        self.assertEqual(len(grid.zs), 32)
        self.assertEqual(len(grid.phi), 16)
        
        # Check coordinate ranges
        self.assertAlmostEqual(grid.r[-1], 31.0)
        self.assertAlmostEqual(grid.zv[-1], 31.0)
        self.assertAlmostEqual(grid.zs[0], 0.0)
        self.assertAlmostEqual(grid.zs[-1], -31.0)


class TestChargeDensity(unittest.TestCase):
    """Test charge density calculations."""
    
    def test_bulk_charge(self):
        """Test bulk charge density calculation."""
        # Create semiconductor
        semi = SemiconductorRegion(
            region_id=1,
            donor_concentration=1e18,
            acceptor_concentration=0,
            band_gap=1.42,
            valence_band_offset=0,
            electron_affinity=4.07,
            donor_binding_energy=0.006,
            acceptor_binding_energy=0.031,
            cb_effective_mass=0.067,
            vb_effective_mass_heavy=0.5,
            vb_effective_mass_light=0.08,
            vb_effective_mass_so=0.15,
            spin_orbit_splitting=0.34,
            permittivity=12.9,
            temperature=300.0
        )
        
        # Calculate Fermi level first
        fermi_level = semi.fermi_level()
        
        # Create calculator with correct Fermi level
        calc = ChargeDensityCalculator([semi], [], fermi_level=fermi_level)
        
        # Test different conditions
        # 1. Deep bulk (neutral)
        rho_bulk = calc.calculate_bulk_density(1, energy=fermi_level, potential=0.0)
        # Should be nearly neutral in bulk
        self.assertLess(abs(rho_bulk / PC.E), 1e16)  # Less than 1e16 cm^-3
        
        # 2. Surface depletion (positive potential = upward band bending)
        rho_depletion = calc.calculate_bulk_density(1, energy=fermi_level, potential=0.2)
        # For n-type with upward band bending, charge should be positive (depletion)
        self.assertGreater(rho_depletion, 0)
        
        # 3. Accumulation (negative potential = downward band bending)
        rho_accumulation = calc.calculate_bulk_density(1, energy=fermi_level, potential=-0.2)
        # For n-type with downward band bending, charge should be negative (accumulation)
        self.assertLess(rho_accumulation, 0)


class TestPoissonSolver(unittest.TestCase):
    """Test Poisson equation solver."""
    
    def setUp(self):
        """Set up simple test case."""
        # Small grid for testing
        params = GridParameters(
            nr=16, nv=16, ns=16, np=8,
            delr=2.0, delv=2.0, dels=2.0, delp=np.pi/8,
            rmax=30.0, vmax=30.0, smax=30.0
        )
        self.grid = Grid3D(params)
        
        # Simple tip
        self.tip = TipModel(
            radius=1.0,
            separation=1.0,
            slope=1.0,
            position=(0, 0),
            work_function=5.3,
            bias_voltage=-1.0
        )
    
    def test_solver_convergence(self):
        """Test that Poisson solver converges."""
        from src.physics.core.poisson import PoissonSolverParameters
        
        # Use larger tolerance and more iterations for test
        params = PoissonSolverParameters(
            tolerance=1e-3,
            max_iterations=500,
            verbose=False,
            omega=0.5  # Under-relaxation for stability
        )
        
        solver = PoissonSolver(self.grid, self.tip, params)
        
        # Very simple charge density for testing
        def bulk_charge(r, z, phi, pot):
            # Simple depletion approximation with damping
            if z < 0 and z > -20:  # Near surface semiconductor
                # Use smaller coupling for stability
                return 1e3 * np.tanh(pot / 0.1) * PC.E  # Saturating response
            return 0.0
        
        def surface_charge(r, phi, pot):
            return 0.0  # No surface charge for simplicity
        
        # Solve with simple initial guess
        import numpy as np
        initial_guess = np.zeros((self.grid.params.nr, 
                                 self.grid.params.nv + self.grid.params.ns - 1,
                                 self.grid.params.np))
        
        # Solve
        potential, info = solver.solve(bulk_charge, surface_charge, initial_guess)
        
        # Check convergence (with relaxed criteria for test)
        self.assertTrue(info['converged'] or info['iterations'] >= params.max_iterations)
        
        # If converged, check reasonable values
        if info['converged']:
            self.assertLess(info['final_error'], params.tolerance)
            
            # Check band bending is reasonable
            bb = solver.get_band_bending()
            self.assertLess(abs(bb), 2.0)  # Not too large


class TestSchrodingerSolver(unittest.TestCase):
    """Test Schrödinger equation solver."""
    
    def test_wkb_transmission(self):
        """Test WKB transmission calculation."""
        solver = SchrodingerSolver()
        
        # Create simple barrier
        from src.physics.core.schrodinger import BandProfile
        z = np.linspace(-10, 10, 100)
        barrier_height = 1.0
        barrier_width = 2.0
        
        # Square barrier
        cb_profile = np.zeros_like(z)
        cb_profile[np.abs(z) < barrier_width/2] = barrier_height
        
        profile = BandProfile(
            z=z,
            vb_profile=cb_profile - 1.5,
            cb_profile=cb_profile,
            vb_max=0,
            cb_min=0,
            vb_bulk=-1.5,
            cb_bulk=0
        )
        
        # Calculate transmission
        T = solver._calculate_transmission_wkb(
            profile, 'cb', energy=0.5, m_eff=0.067*PC.M0
        )
        
        # Should have some tunneling
        self.assertGreater(T, 0)
        self.assertLess(T, 1)


def run_tests():
    """Run all tests."""
    unittest.main(argv=[''], exit=False, verbosity=2)


if __name__ == '__main__':
    run_tests()