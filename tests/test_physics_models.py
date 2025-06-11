"""
Comprehensive test suite for Phase 2 Physics Models

This test suite validates all the core physics models implemented in Phase 2,
including materials management, charge density calculations, Poisson solving,
and tunneling current calculations.

Author: odindino
"""

import pytest
import numpy as np
import warnings
from typing import Dict

# Import all physics modules
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from physics.materials import (
    SemiconductorMaterial, MaterialDatabase, PhysicalConstants,
    SurfaceStateParameters, MaterialParameters, create_default_surface_states
)
from physics.charge_density import (
    ChargeDensityCalculator, ChargeDensityConfig,
    calculate_intrinsic_fermi_level, calculate_debye_length
)
from physics.poisson import (
    PoissonSolver, PoissonConfig, Grid3D, BoundaryType,
    create_default_grid, setup_basic_stm_solver
)
from physics.tunneling_current import (
    TunnelingCurrentCalculator, TunnelingConfig,
    calculate_simple_stm_current
)


class TestMaterialsModule:
    """Test materials and parameter management"""
    
    def test_semiconductor_material_creation(self):
        """Test semiconductor material dataclass"""
        si = SemiconductorMaterial(
            name="Test_Silicon",
            bandgap=1.12,
            temperature=300.0,
            donor_concentration=1e15
        )
        
        assert si.name == "Test_Silicon"
        assert si.bandgap == 1.12
        assert abs(si.thermal_energy - 0.0259) < 0.001  # kT at 300K
        assert si.net_doping == 1e15
        assert si.intrinsic_concentration > 0
    
    def test_material_database(self):
        """Test material database functionality"""
        db = MaterialDatabase()
        
        # Check available materials
        materials = db.list_materials()
        assert "Si_n" in materials
        assert "Si_p" in materials
        assert "GaAs_n" in materials
        
        # Get specific material
        si_n = db.get_material("Si_n")
        assert si_n.name == "Silicon_n_type"
        assert si_n.bandgap == 1.12
        assert si_n.donor_concentration > 0
        assert si_n.acceptor_concentration == 0
        
        # Test error handling
        with pytest.raises(ValueError):
            db.get_material("NonexistentMaterial")
    
    def test_physical_constants(self):
        """Test physical constants consistency"""
        constants = PhysicalConstants()
        
        # Check fundamental constants are reasonable
        assert 1.6e-19 < constants.ELEMENTARY_CHARGE < 1.61e-19
        assert 8.85e-12 < constants.VACUUM_PERMITTIVITY < 8.86e-12
        assert 6.62e-34 < constants.PLANCK_CONSTANT < 6.63e-34
        
        # Check derived constants
        assert abs(constants.K_B_EV - 8.617333e-5) < 1e-8
        
        # Check MultInt specific constants
        assert constants.EEP > 0
        assert constants.C_KINETIC > 0
        assert constants.RESISTANCE_QUANTUM > 12000
    
    def test_surface_states(self):
        """Test surface state parameters"""
        surface = create_default_surface_states()
        
        assert surface.donor_density > 0
        assert surface.acceptor_density > 0
        assert 0 < surface.donor_energy < 2
        assert 0 < surface.acceptor_energy < 2
        assert surface.temperature_dependent is True


class TestChargeDensityModule:
    """Test charge density calculations"""
    
    def setup_method(self):
        """Setup for charge density tests"""
        self.db = MaterialDatabase()
        self.si_n = self.db.get_material("Si_n")
        self.si_p = self.db.get_material("Si_p")
        self.calculator = ChargeDensityCalculator()
    
    def test_electron_density_calculation(self):
        """Test electron density in conduction band"""
        fermi_level = 0.6  # eV above valence band
        potential = 0.0    # Flat band
        
        n = self.calculator.calculate_electron_density(self.si_n, fermi_level, potential)
        
        assert n > 0  # Should have electrons above bandgap
        assert np.isfinite(n)
        
        # Test vectorized input
        potentials = np.array([0.0, 0.1, 0.2])
        n_array = self.calculator.calculate_electron_density(self.si_n, fermi_level, potentials)
        assert len(n_array) == 3
        assert all(n_val > 0 for n_val in n_array)
    
    def test_hole_density_calculation(self):
        """Test hole density in valence band"""
        fermi_level = 0.5  # eV above valence band
        potential = 0.0
        
        p = self.calculator.calculate_hole_density(self.si_n, fermi_level, potential)
        
        assert p > 0  # Should have holes below Fermi level
        assert np.isfinite(p)
    
    def test_dopant_ionization(self):
        """Test ionized dopant calculations"""
        fermi_level = 0.6  # eV
        potential = 0.0
        
        # Test donor ionization (n-type silicon)
        n_d_plus = self.calculator.calculate_ionized_donor_density(self.si_n, fermi_level, potential)
        assert 0 < n_d_plus <= self.si_n.donor_concentration
        
        # Test acceptor ionization (p-type silicon)
        n_a_minus = self.calculator.calculate_ionized_acceptor_density(self.si_p, fermi_level, potential)
        assert 0 < n_a_minus <= self.si_p.acceptor_concentration
    
    def test_bulk_charge_density(self):
        """Test total bulk charge density calculation"""
        fermi_level = 0.6
        potential = 0.0
        
        rho = self.calculator.calculate_bulk_charge_density(self.si_n, fermi_level, potential)
        assert np.isfinite(rho)
        
        # For n-type, expect positive charge when depleted
        potential_depleted = 0.3  # Upward band bending
        rho_depleted = self.calculator.calculate_bulk_charge_density(
            self.si_n, fermi_level, potential_depleted
        )
        assert rho_depleted > 0  # Positive ionized donors
    
    def test_charge_neutrality(self):
        """Test charge neutrality condition"""
        # Find equilibrium Fermi level
        ef_eq = self.calculator.find_equilibrium_fermi_level(self.si_n)
        
        # Check it's reasonable (within bandgap region)
        assert 0 < ef_eq < self.si_n.bandgap
        
        # Check charge neutrality is satisfied
        balance_func = self.calculator.calculate_charge_neutrality_condition(self.si_n)
        charge_imbalance = balance_func(ef_eq)
        assert abs(charge_imbalance) < 1e10  # Should be nearly zero
    
    def test_intrinsic_fermi_level(self):
        """Test intrinsic Fermi level calculation"""
        si_intrinsic = self.db.get_material("Si_intrinsic")
        ef_intrinsic = calculate_intrinsic_fermi_level(si_intrinsic)
        
        # Should be near mid-gap
        assert abs(ef_intrinsic - si_intrinsic.bandgap/2) < 0.1
    
    def test_debye_length(self):
        """Test Debye screening length calculation"""
        ld = calculate_debye_length(self.si_n, self.si_n.donor_concentration)
        
        # Should be reasonable for typical doping
        assert 1 < ld < 1000  # nm
        assert np.isfinite(ld)
    
    def test_lookup_table_generation(self):
        """Test charge density lookup table"""
        fermi_level = 0.6
        potential_range = (-0.5, 0.5)
        
        v_grid, rho_grid = self.calculator.build_charge_density_lookup_table(
            self.si_n, potential_range, fermi_level
        )
        
        assert len(v_grid) == len(rho_grid)
        assert len(v_grid) > 100  # Should have decent resolution
        assert all(np.isfinite(rho) for rho in rho_grid)
        
        # Test interpolation
        test_potential = 0.1
        rho_interp = self.calculator.interpolate_charge_density(test_potential)
        assert np.isfinite(rho_interp)


class TestPoissonModule:
    """Test Poisson equation solver"""
    
    def setup_method(self):
        """Setup for Poisson tests"""
        self.db = MaterialDatabase()
        self.si_n = self.db.get_material("Si_n")
        self.grid = create_default_grid(tip_radius=5.0, separation=1.0)
        
    def test_grid_creation(self):
        """Test 3D grid setup"""
        assert self.grid.nr > 0
        assert self.grid.nphi > 0
        assert self.grid.nz_vacuum > 0
        assert self.grid.nz_semiconductor > 0
        
        # Check grid arrays
        assert len(self.grid.r) == self.grid.nr
        assert len(self.grid.phi) == self.grid.nphi
        assert len(self.grid.z_vacuum) == self.grid.nz_vacuum
        assert len(self.grid.z_semiconductor) == self.grid.nz_semiconductor
        
        # Check grid spacing
        assert self.grid.dr > 0
        assert self.grid.dphi > 0
        assert self.grid.dz_vacuum > 0
        assert self.grid.dz_semiconductor > 0
    
    def test_poisson_solver_setup(self):
        """Test Poisson solver initialization"""
        solver = setup_basic_stm_solver(self.si_n)
        
        assert solver.grid.nr > 0
        assert 1 in solver.materials
        assert solver.materials[1] == self.si_n
        assert 'tip' in solver.boundary_conditions
        assert 'back' in solver.boundary_conditions
    
    def test_boundary_conditions(self):
        """Test boundary condition setup"""
        solver = PoissonSolver(self.grid)
        
        # Test different boundary types
        solver.set_boundary_condition('test', BoundaryType.DIRICHLET, 1.0)
        assert 'test' in solver.boundary_conditions
        assert solver.boundary_conditions['test']['type'] == BoundaryType.DIRICHLET
        assert solver.boundary_conditions['test']['value'] == 1.0
    
    def test_permittivity_map(self):
        """Test permittivity map generation"""
        solver = PoissonSolver(self.grid)
        solver.set_material(1, self.si_n)
        solver.build_permittivity_map()
        
        assert solver.permittivity_map is not None
        assert solver.permittivity_map.shape == (self.grid.nr, self.grid.nphi, self.grid.nz_total)
        
        # Check vacuum region has ε = 1
        vacuum_region = solver.permittivity_map[:, :, :self.grid.nz_vacuum]
        assert np.allclose(vacuum_region, 1.0)
        
        # Check semiconductor region has ε = εᵣ
        semiconductor_region = solver.permittivity_map[:, :, self.grid.nz_vacuum:]
        assert np.allclose(semiconductor_region, self.si_n.relative_permittivity)
    
    @pytest.mark.slow
    def test_simple_poisson_solve(self):
        """Test simple Poisson equation solve"""
        solver = setup_basic_stm_solver(self.si_n)
        
        # Create simple charge density (uniform)
        nr, nphi, nz = solver.grid.nr, solver.grid.nphi, solver.grid.nz_total
        charge_density = np.zeros((nr, nphi, nz))
        charge_density[:, :, solver.grid.nz_vacuum:] = 1e15  # Uniform charge in semiconductor
        
        # Solve
        potential = solver.solve_poisson_equation(charge_density)
        
        assert potential.shape == (nr, nphi, nz)
        assert np.all(np.isfinite(potential))
        
        # Electric field calculation
        er, ephi, ez = solver.calculate_electric_field(potential)
        assert er.shape == potential.shape
        assert np.all(np.isfinite(er))


class TestTunnelingCurrentModule:
    """Test tunneling current calculations"""
    
    def setup_method(self):
        """Setup for tunneling current tests"""
        self.db = MaterialDatabase()
        self.si_n = self.db.get_material("Si_n")
        self.calculator = TunnelingCurrentCalculator()
        
    def test_tunneling_calculator_creation(self):
        """Test tunneling calculator initialization"""
        config = TunnelingConfig(energy_points=20, k_points=10)
        calc = TunnelingCurrentCalculator(config)
        
        assert calc.config.energy_points == 20
        assert calc.config.k_points == 10
        assert calc.constants is not None
    
    def test_barrier_profile_calculation(self):
        """Test tunnel barrier profile"""
        z_grid = np.linspace(0, 2, 100)  # 2 nm separation
        potential_1d = np.zeros_like(z_grid)  # Flat potential
        
        barrier = self.calculator.calculate_barrier_profile(
            potential_1d, z_grid, tip_position=0.0, sample_position=2.0
        )
        
        assert len(barrier) == len(z_grid)
        assert np.all(np.isfinite(barrier))
        assert barrier[0] > 0  # Work function at tip
        assert barrier[-1] > 0  # Work function at sample
    
    def test_transmission_coefficient(self):
        """Test transmission coefficient calculation"""
        z_grid = np.linspace(0, 1, 50)  # 1 nm barrier
        barrier = np.ones_like(z_grid) * 5.0  # 5 eV barrier
        
        # Low energy electron - should have low transmission
        transmission_low = self.calculator.calculate_transmission_coefficient(
            energy=1.0, k_parallel=0.0, barrier_profile=barrier, z_grid=z_grid
        )
        assert 0 <= transmission_low <= 1
        
        # High energy electron - should have higher transmission
        transmission_high = self.calculator.calculate_transmission_coefficient(
            energy=6.0, k_parallel=0.0, barrier_profile=barrier, z_grid=z_grid
        )
        assert transmission_high >= transmission_low
    
    def test_wavefunction_calculation(self):
        """Test sample wavefunction calculation"""
        energy = 1.5  # eV
        k_parallel = 0.1  # nm^-1
        
        wf_val, wf_deriv = self.calculator.calculate_sample_wavefunction(
            energy, k_parallel, self.si_n, "conduction"
        )
        
        assert wf_val >= 0
        assert np.isfinite(wf_val)
        assert np.isfinite(wf_deriv)
    
    def test_density_of_states(self):
        """Test density of states calculation"""
        energy = 1.5  # eV above valence band
        
        dos_cb = self.calculator.calculate_density_of_states(
            self.si_n, energy, "conduction"
        )
        dos_vb = self.calculator.calculate_density_of_states(
            self.si_n, energy, "valence_heavy"
        )
        
        assert dos_cb >= 0
        assert dos_vb >= 0
        assert np.isfinite(dos_cb)
        assert np.isfinite(dos_vb)
    
    def test_simple_current_calculation(self):
        """Test simple STM current calculation"""
        bias_voltage = 1.0  # V
        
        current = calculate_simple_stm_current(self.si_n, bias_voltage)
        
        assert np.isfinite(current)
        assert current != 0  # Should have some current
        
        # Test different biases
        currents = [calculate_simple_stm_current(self.si_n, v) for v in [-1, 0, 1]]
        
        # Current should change with bias
        assert not all(abs(i - currents[1]) < 1e-15 for i in currents)
    
    @pytest.mark.slow
    def test_extended_state_current(self):
        """Test extended state current calculation"""
        # Create simple 1D potential and barrier
        z_grid = np.linspace(0, 2, 100)
        potential_profile = np.zeros_like(z_grid)
        
        barrier = self.calculator.calculate_barrier_profile(
            potential_profile, z_grid, 0.0, 2.0
        )
        
        fermi_level = 0.6  # eV
        bias_voltage = 1.0  # V
        
        current_components = self.calculator.calculate_extended_state_current(
            self.si_n, bias_voltage, fermi_level, barrier, z_grid
        )
        
        assert 'total' in current_components
        assert 'conduction' in current_components
        assert 'valence_heavy' in current_components
        
        # Check all components are finite
        for component, current in current_components.items():
            assert np.isfinite(current)


class TestPhysicsIntegration:
    """Integration tests for complete physics simulation"""
    
    def setup_method(self):
        """Setup for integration tests"""
        self.db = MaterialDatabase()
        self.si_n = self.db.get_material("Si_n")
        
    @pytest.mark.slow
    def test_complete_stm_simulation(self):
        """Test complete STM simulation workflow"""
        # 1. Setup geometry and materials
        grid = create_default_grid(tip_radius=5.0, separation=1.0)
        solver = setup_basic_stm_solver(self.si_n)
        charge_calc = ChargeDensityCalculator()
        current_calc = TunnelingCurrentCalculator()
        
        # 2. Find equilibrium Fermi level
        fermi_level = charge_calc.find_equilibrium_fermi_level(self.si_n)
        assert 0 < fermi_level < self.si_n.bandgap
        
        # 3. Simplified current calculation (without full self-consistency)
        bias_voltage = 1.0
        current = calculate_simple_stm_current(self.si_n, bias_voltage)
        
        assert np.isfinite(current)
        assert current != 0
        
        print(f"Complete STM simulation successful:")
        print(f"  Material: {self.si_n.name}")
        print(f"  Fermi level: {fermi_level:.3f} eV")
        print(f"  Current at {bias_voltage} V: {current:.2e} A")
    
    def test_material_parameter_consistency(self):
        """Test consistency between different material parameter representations"""
        # Get material from database
        si_from_db = self.db.get_material("Si_n")
        
        # Create equivalent material manually
        si_manual = SemiconductorMaterial(
            name="Silicon_n_type",
            relative_permittivity=11.7,
            bandgap=1.12,
            donor_concentration=1e15,
            temperature=300.0
        )
        
        # Compare key parameters
        assert si_from_db.bandgap == si_manual.bandgap
        assert si_from_db.relative_permittivity == si_manual.relative_permittivity
        assert abs(si_from_db.thermal_energy - si_manual.thermal_energy) < 1e-6
        assert si_from_db.net_doping == si_manual.net_doping
    
    def test_physical_units_consistency(self):
        """Test that physical units are consistent throughout calculations"""
        # Energy units should be eV
        assert self.si_n.bandgap > 1  # Silicon bandgap ~ 1.12 eV
        assert self.si_n.thermal_energy < 0.1  # kT at room temp ~ 0.026 eV
        
        # Concentration units should be cm^-3
        assert 1e14 < self.si_n.donor_concentration < 1e17  # Typical doping
        
        # Length units should be nm in grid calculations
        grid = create_default_grid()
        assert grid.r_max > 10  # Should be tens of nm
        assert grid.dz_vacuum < 10  # Grid spacing should be sub-nm


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
