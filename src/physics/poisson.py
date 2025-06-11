"""
Poisson Equation Solver Module

This module implements finite difference solvers for the 3D Poisson equation
in cylindrical coordinates, as used in STM simulations. The solver handles
complex boundary conditions including tip geometry, semiconductor surfaces,
and material interfaces.

Based on the electrostatic model in semitip3-6.1.f

Author: odindino
"""

import numpy as np
from typing import Dict, Optional, Tuple, Union, List, Callable
from dataclasses import dataclass, field
from enum import Enum
import warnings
from scipy.sparse import csr_matrix, linalg
from scipy.sparse.linalg import spsolve
import time

from .materials import SemiconductorMaterial, PhysicalConstants
from .charge_density import ChargeDensityCalculator

# Import numerical functions with proper path handling
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from utils.numerical import GoldenSectionOptimizer


class BoundaryType(Enum):
    """Types of boundary conditions"""
    DIRICHLET = "dirichlet"  # Fixed potential
    NEUMANN = "neumann"      # Fixed normal derivative
    MIXED = "mixed"          # Mixed boundary condition
    PERIODIC = "periodic"    # Periodic boundary


class SolverMethod(Enum):
    """Poisson solver methods"""
    FINITE_DIFFERENCE = "finite_difference"
    SUCCESSIVE_OVER_RELAXATION = "sor"
    CONJUGATE_GRADIENT = "cg"
    MULTIGRID = "multigrid"


@dataclass
class Grid3D:
    """Three-dimensional grid in cylindrical coordinates (r, φ, z)"""
    
    # Grid dimensions
    nr: int = 64    # Radial points
    nphi: int = 32  # Angular points
    nz_vacuum: int = 32   # Vacuum z points
    nz_semiconductor: int = 64  # Semiconductor z points
    
    # Physical dimensions [nm]
    r_max: float = 100.0
    phi_max: float = 2 * np.pi
    z_vacuum_max: float = 20.0
    z_semiconductor_max: float = 50.0
    
    # Grid spacing
    dr: float = field(init=False)
    dphi: float = field(init=False)  
    dz_vacuum: float = field(init=False)
    dz_semiconductor: float = field(init=False)
    
    # Grid arrays
    r: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)
    z_vacuum: np.ndarray = field(init=False)
    z_semiconductor: np.ndarray = field(init=False)
    
    def __post_init__(self):
        """Initialize grid arrays"""
        # Grid spacing
        self.dr = self.r_max / (self.nr - 1)
        self.dphi = self.phi_max / self.nphi
        self.dz_vacuum = self.z_vacuum_max / (self.nz_vacuum - 1)
        self.dz_semiconductor = self.z_semiconductor_max / (self.nz_semiconductor - 1)
        
        # Grid arrays
        self.r = np.linspace(0, self.r_max, self.nr)
        self.phi = np.linspace(0, self.phi_max, self.nphi, endpoint=False)
        self.z_vacuum = np.linspace(0, self.z_vacuum_max, self.nz_vacuum)
        self.z_semiconductor = np.linspace(0, self.z_semiconductor_max, self.nz_semiconductor)
        
        # Total z grid (vacuum + semiconductor)
        self.nz_total = self.nz_vacuum + self.nz_semiconductor - 1  # Overlap at interface
        self.z_total = np.concatenate([
            self.z_vacuum,
            self.z_vacuum_max + self.z_semiconductor[1:]  # Skip overlapping point
        ])


@dataclass
class PoissonConfig:
    """Configuration for Poisson equation solver"""
    
    # Solver parameters
    method: SolverMethod = SolverMethod.FINITE_DIFFERENCE
    max_iterations: int = 10000
    convergence_tolerance: float = 1e-6
    relaxation_parameter: float = 1.5  # For SOR method
    
    # Grid refinement
    enable_adaptive_mesh: bool = True
    max_refinement_levels: int = 3  # IPMAX in Fortran
    refinement_tolerance: float = 1e-3
    
    # Physical parameters
    include_image_potential: bool = True  # IMPOT flag
    surface_smoothing: bool = True
    
    # Numerical stability
    potential_damping: float = 0.1  # Damping factor for potential updates
    charge_density_threshold: float = 1e-12  # Minimum charge density [cm^-3]


class PoissonSolver:
    """
    Three-dimensional Poisson equation solver for STM simulations.
    
    Solves ∇·[ε(r)∇V(r)] = -ρ(r) using finite difference methods in
    cylindrical coordinates. Handles complex boundary conditions for
    tip-sample geometry and self-consistent iterations.
    
    Based on the implementation in semitip3-6.1.f
    """
    
    def __init__(self, 
                 grid: Grid3D,
                 config: Optional[PoissonConfig] = None):
        self.grid = grid
        self.config = config or PoissonConfig()
        self.constants = PhysicalConstants()
        
        # Solution arrays
        self.potential_vacuum: Optional[np.ndarray] = None
        self.potential_semiconductor: Optional[np.ndarray] = None
        self.charge_density: Optional[np.ndarray] = None
        
        # Material parameters
        self.materials: Dict[int, SemiconductorMaterial] = {}
        self.permittivity_map: Optional[np.ndarray] = None
        
        # Convergence history
        self.convergence_history: List[float] = []
        self.iteration_count: int = 0
        
        # Boundary conditions
        self.boundary_conditions: Dict[str, Dict] = {}
        
        # Golden section optimizer for self-consistent iterations
        self.optimizer = GoldenSectionOptimizer()
        
    def set_material(self, region_id: int, material: SemiconductorMaterial):
        """Set material properties for a region"""
        self.materials[region_id] = material
        
    def set_boundary_condition(self, 
                             region: str,
                             bc_type: BoundaryType, 
                             value: Union[float, Callable, np.ndarray],
                             **kwargs):
        """
        Set boundary condition for a region.
        
        Args:
            region: Boundary region name ('tip', 'surface', 'sides', 'back')
            bc_type: Type of boundary condition
            value: Boundary value (constant, function, or array)
            **kwargs: Additional parameters
        """
        self.boundary_conditions[region] = {
            'type': bc_type,
            'value': value,
            **kwargs
        }
    
    def setup_tip_boundary(self, 
                          tip_bias: float,
                          tip_radius: float,
                          separation: float,
                          work_function_diff: float = 0.0):
        """
        Setup boundary conditions for STM tip.
        
        Args:
            tip_bias: Applied bias voltage [V]
            tip_radius: Tip radius [nm]
            separation: Tip-sample separation [nm]
            work_function_diff: Tip-sample work function difference [eV]
        """
        def tip_potential(r, phi, z):
            """Calculate tip potential including geometry"""
            # Distance from tip apex
            tip_z = separation
            distance = np.sqrt(r**2 + (z - tip_z)**2)
            
            # Potential from spherical tip
            if distance < tip_radius:
                # Inside tip region - constant potential
                return tip_bias + work_function_diff
            else:
                # Outside tip - potential of sphere
                return (tip_bias + work_function_diff) * tip_radius / distance
        
        self.set_boundary_condition('tip', BoundaryType.DIRICHLET, tip_potential)
    
    def setup_sample_boundary(self, 
                            material: SemiconductorMaterial,
                            fermi_level: float):
        """
        Setup boundary conditions for semiconductor sample surface.
        
        Args:
            material: Semiconductor material parameters
            fermi_level: Fermi level relative to valence band [eV]
        """
        # Surface potential includes band bending and work function
        surface_potential = material.electron_affinity + material.bandgap - fermi_level
        
        self.set_boundary_condition('surface', BoundaryType.MIXED, surface_potential,
                                  material=material, fermi_level=fermi_level)
    
    def build_permittivity_map(self):
        """Build 3D permittivity map for heterogeneous materials"""
        # Initialize with vacuum permittivity
        self.permittivity_map = np.ones((self.grid.nr, self.grid.nphi, self.grid.nz_total))
        
        # Set semiconductor regions
        nz_vacuum = self.grid.nz_vacuum
        for region_id, material in self.materials.items():
            if region_id == 1:  # Main semiconductor region
                # Semiconductor part of the grid
                self.permittivity_map[:, :, nz_vacuum:] = material.relative_permittivity
    
    def calculate_finite_difference_coefficients(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate finite difference coefficients for cylindrical coordinates.
        
        Returns:
            Tuple of (coefficient_matrix, rhs_vector)
        """
        nr, nphi, nz = self.grid.nr, self.grid.nphi, self.grid.nz_total
        total_size = nr * nphi * nz
        
        # Create sparse matrix for coefficients
        row_indices = []
        col_indices = []
        data = []
        rhs = np.zeros(total_size)
        
        dr, dphi = self.grid.dr, self.grid.dphi
        
        # Conversion factor: e/ε₀ in SI units, then to convenient units
        # EEP = 1.80943E-20 V⋅cm from Fortran code
        eep_factor = self.constants.EEP
        
        for i in range(nr):
            for j in range(nphi):
                for k in range(nz):
                    # Linear index
                    idx = i * nphi * nz + j * nz + k
                    
                    # Get grid spacing for this region
                    if k < self.grid.nz_vacuum:
                        dz = self.grid.dz_vacuum
                    else:
                        dz = self.grid.dz_semiconductor
                    
                    # Handle boundary conditions
                    if self._is_boundary_point(i, j, k):
                        # Apply boundary condition
                        row_indices.append(idx)
                        col_indices.append(idx)
                        data.append(1.0)
                        rhs[idx] = self._get_boundary_value(i, j, k)
                        continue
                    
                    # Interior point - build finite difference stencil
                    r = self.grid.r[i]
                    
                    # Central coefficient
                    central_coeff = 0.0
                    
                    # Radial derivatives: (1/r) d/dr(r dV/dr)
                    if i > 0 and i < nr - 1:
                        # Central difference
                        r_plus = self.grid.r[i + 1]
                        r_minus = self.grid.r[i - 1]
                        
                        coeff_r_plus = (r + dr/2) / (dr**2)
                        coeff_r_minus = (r - dr/2) / (dr**2)
                        coeff_r_central = -(coeff_r_plus + coeff_r_minus)
                        
                        # Add to matrix
                        if i > 0:
                            row_indices.append(idx)
                            col_indices.append((i-1) * nphi * nz + j * nz + k)
                            data.append(coeff_r_minus)
                        
                        if i < nr - 1:
                            row_indices.append(idx)
                            col_indices.append((i+1) * nphi * nz + j * nz + k)
                            data.append(coeff_r_plus)
                        
                        central_coeff += coeff_r_central
                        
                    elif i == 0:
                        # Special handling at r = 0 (use L'Hôpital's rule)
                        coeff_r = 4.0 / (dr**2)
                        row_indices.append(idx)
                        col_indices.append((i+1) * nphi * nz + j * nz + k)
                        data.append(coeff_r)
                        central_coeff += -coeff_r
                    
                    # Angular derivatives: (1/r²) d²V/dφ²
                    if r > 0 and nphi > 1:
                        coeff_phi = 1.0 / (r**2 * dphi**2)
                        
                        # Previous phi index (periodic)
                        j_prev = (j - 1) % nphi
                        j_next = (j + 1) % nphi
                        
                        row_indices.append(idx)
                        col_indices.append(i * nphi * nz + j_prev * nz + k)
                        data.append(coeff_phi)
                        
                        row_indices.append(idx)
                        col_indices.append(i * nphi * nz + j_next * nz + k)
                        data.append(coeff_phi)
                        
                        central_coeff += -2.0 * coeff_phi
                    
                    # Axial derivatives: d²V/dz²
                    if k > 0 and k < nz - 1:
                        coeff_z = 1.0 / (dz**2)
                        
                        row_indices.append(idx)
                        col_indices.append(i * nphi * nz + j * nz + (k-1))
                        data.append(coeff_z)
                        
                        row_indices.append(idx)
                        col_indices.append(i * nphi * nz + j * nz + (k+1))
                        data.append(coeff_z)
                        
                        central_coeff += -2.0 * coeff_z
                    
                    # Central coefficient
                    row_indices.append(idx)
                    col_indices.append(idx)
                    data.append(central_coeff)
                    
                    # Right-hand side: -ρ/ε
                    if self.charge_density is not None:
                        epsilon = self.permittivity_map[i, j, k] if self.permittivity_map is not None else 1.0
                        rho = self.charge_density[i, j, k]
                        rhs[idx] = -eep_factor * rho / epsilon
        
        # Create sparse matrix
        coefficient_matrix = csr_matrix((data, (row_indices, col_indices)), 
                                      shape=(total_size, total_size))
        
        return coefficient_matrix, rhs
    
    def _is_boundary_point(self, i: int, j: int, k: int) -> bool:
        """Check if point (i,j,k) is on the boundary"""
        nr, nphi, nz = self.grid.nr, self.grid.nphi, self.grid.nz_total
        
        # Boundary conditions at edges
        if (i == 0 or i == nr - 1 or 
            k == 0 or k == nz - 1):
            return True
        
        return False
    
    def _get_boundary_value(self, i: int, j: int, k: int) -> float:
        """Get boundary value at point (i,j,k)"""
        r = self.grid.r[i]
        phi = self.grid.phi[j]
        
        if k < self.grid.nz_vacuum:
            z = self.grid.z_vacuum[k]
        else:
            z_idx = k - self.grid.nz_vacuum + 1
            z = self.grid.z_vacuum_max + self.grid.z_semiconductor[z_idx]
        
        # Apply appropriate boundary condition
        if k == 0:  # Top boundary (tip region)
            if 'tip' in self.boundary_conditions:
                bc = self.boundary_conditions['tip']
                if callable(bc['value']):
                    return bc['value'](r, phi, z)
                else:
                    return bc['value']
        elif k == self.grid.nz_total - 1:  # Bottom boundary
            if 'back' in self.boundary_conditions:
                bc = self.boundary_conditions['back']
                return bc['value']
        elif i == self.grid.nr - 1:  # Outer radial boundary
            if 'sides' in self.boundary_conditions:
                bc = self.boundary_conditions['sides']
                return bc['value']
        
        # Default: zero potential
        return 0.0
    
    def solve_poisson_equation(self, 
                             charge_density: np.ndarray,
                             initial_guess: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Solve the Poisson equation for given charge density.
        
        Args:
            charge_density: 3D charge density array [e·cm^-3]
            initial_guess: Initial potential guess [V]
            
        Returns:
            3D potential array [V]
        """
        self.charge_density = charge_density
        
        # Build permittivity map if not already done
        if self.permittivity_map is None:
            self.build_permittivity_map()
        
        # Build finite difference system
        A, b = self.calculate_finite_difference_coefficients()
        
        # Solve linear system
        if self.config.method == SolverMethod.FINITE_DIFFERENCE:
            solution = spsolve(A, b)
        else:
            raise NotImplementedError(f"Solver method {self.config.method} not implemented")
        
        # Reshape solution to 3D grid
        nr, nphi, nz = self.grid.nr, self.grid.nphi, self.grid.nz_total
        potential_3d = solution.reshape((nr, nphi, nz))
        
        return potential_3d
    
    def solve_self_consistent(self,
                            charge_calculator: ChargeDensityCalculator,
                            materials: Dict[int, SemiconductorMaterial],
                            initial_fermi_level: float,
                            **kwargs) -> Tuple[np.ndarray, float, Dict]:
        """
        Solve self-consistent Poisson-charge system.
        
        Implements the main self-consistent loop from MultInt3-6.4.f
        
        Args:
            charge_calculator: Charge density calculator
            materials: Dictionary of materials by region
            initial_fermi_level: Initial guess for Fermi level [eV]
            **kwargs: Additional parameters
            
        Returns:
            Tuple of (potential_3d, final_fermi_level, convergence_info)
        """
        # Initialize
        fermi_level = initial_fermi_level
        self.materials.update(materials)
        
        # Get main material (region 1)
        main_material = materials.get(1)
        if main_material is None:
            raise ValueError("No material defined for region 1")
        
        # Initialize potential with flat band condition
        nr, nphi, nz = self.grid.nr, self.grid.nphi, self.grid.nz_total
        potential = np.zeros((nr, nphi, nz))
        
        convergence_info = {
            'iterations': 0,
            'convergence_history': [],
            'final_residual': None,
            'success': False
        }
        
        print(f"Starting self-consistent iteration with E_F = {fermi_level:.3f} eV")
        
        for iteration in range(self.config.max_iterations):
            start_time = time.time()
            
            # 1. Calculate charge density for current potential and Fermi level
            charge_density = np.zeros((nr, nphi, nz))
            
            for i in range(nr):
                for j in range(nphi):
                    for k in range(nz):
                        if k >= self.grid.nz_vacuum:  # Semiconductor region
                            pot_val = potential[i, j, k]
                            rho = charge_calculator.calculate_bulk_charge_density(
                                main_material, fermi_level, pot_val
                            )
                            charge_density[i, j, k] = rho
            
            # 2. Solve Poisson equation for new potential
            new_potential = self.solve_poisson_equation(charge_density)
            
            # 3. Check convergence
            potential_change = np.max(np.abs(new_potential - potential))
            convergence_info['convergence_history'].append(potential_change)
            
            # 4. Update potential with damping
            damping = self.config.potential_damping
            potential = (1 - damping) * potential + damping * new_potential
            
            iteration_time = time.time() - start_time
            
            print(f"Iteration {iteration + 1}: ΔV_max = {potential_change:.2e} V "
                  f"({iteration_time:.3f}s)")
            
            # Check convergence
            if potential_change < self.config.convergence_tolerance:
                convergence_info['success'] = True
                convergence_info['final_residual'] = potential_change
                break
        
        convergence_info['iterations'] = iteration + 1
        self.convergence_history = convergence_info['convergence_history']
        
        if not convergence_info['success']:
            warnings.warn(f"Self-consistent iteration did not converge in "
                         f"{self.config.max_iterations} iterations")
        
        return potential, fermi_level, convergence_info
    
    def extract_surface_potential(self, potential_3d: np.ndarray) -> np.ndarray:
        """Extract surface potential (at semiconductor interface)"""
        # Surface is at the vacuum-semiconductor interface
        interface_k = self.grid.nz_vacuum - 1
        return potential_3d[:, :, interface_k]
    
    def calculate_electric_field(self, potential_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate electric field from potential: E = -∇V
        
        Returns:
            Tuple of (E_r, E_phi, E_z) components
        """
        nr, nphi, nz = potential_3d.shape
        
        # Initialize field components
        E_r = np.zeros_like(potential_3d)
        E_phi = np.zeros_like(potential_3d)
        E_z = np.zeros_like(potential_3d)
        
        # Radial component: -dV/dr
        for i in range(1, nr - 1):
            E_r[i, :, :] = -(potential_3d[i+1, :, :] - potential_3d[i-1, :, :]) / (2 * self.grid.dr)
        
        # Angular component: -(1/r) dV/dφ
        for j in range(nphi):
            j_prev = (j - 1) % nphi
            j_next = (j + 1) % nphi
            for i in range(nr):
                if self.grid.r[i] > 0:
                    E_phi[i, j, :] = -(potential_3d[i, j_next, :] - potential_3d[i, j_prev, :]) / \
                                    (2 * self.grid.dphi * self.grid.r[i])
        
        # Axial component: -dV/dz
        for k in range(1, nz - 1):
            if k < self.grid.nz_vacuum:
                dz = self.grid.dz_vacuum
            else:
                dz = self.grid.dz_semiconductor
            E_z[:, :, k] = -(potential_3d[:, :, k+1] - potential_3d[:, :, k-1]) / (2 * dz)
        
        return E_r, E_phi, E_z
    
    def get_solution_summary(self) -> Dict:
        """Get summary of the solution"""
        if self.potential_vacuum is None:
            return {"status": "No solution available"}
        
        # Calculate basic statistics
        summary = {
            "grid_points": self.grid.nr * self.grid.nphi * self.grid.nz_total,
            "convergence_iterations": len(self.convergence_history),
            "final_residual": self.convergence_history[-1] if self.convergence_history else None,
            "potential_range": {
                "min": float(np.min(self.potential_vacuum)) if self.potential_vacuum is not None else None,
                "max": float(np.max(self.potential_vacuum)) if self.potential_vacuum is not None else None
            }
        }
        
        return summary


# Convenience functions

def create_default_grid(tip_radius: float = 5.0, 
                      separation: float = 1.0,
                      sample_thickness: float = 50.0) -> Grid3D:
    """Create default grid for STM geometry"""
    
    # Grid dimensions based on geometry
    r_max = max(20.0, 4 * tip_radius)
    z_vacuum = max(10.0, 3 * separation)
    
    return Grid3D(
        nr=64,
        nphi=32,
        nz_vacuum=32,
        nz_semiconductor=64,
        r_max=r_max,
        z_vacuum_max=z_vacuum,
        z_semiconductor_max=sample_thickness
    )


def setup_basic_stm_solver(material: SemiconductorMaterial,
                         tip_bias: float = 1.0,
                         tip_radius: float = 5.0,
                         separation: float = 1.0) -> PoissonSolver:
    """Setup basic STM Poisson solver with typical parameters"""
    
    grid = create_default_grid(tip_radius, separation)
    solver = PoissonSolver(grid)
    
    # Set material
    solver.set_material(1, material)
    
    # Setup boundaries
    solver.setup_tip_boundary(tip_bias, tip_radius, separation)
    solver.set_boundary_condition('back', BoundaryType.DIRICHLET, 0.0)  # Grounded back
    solver.set_boundary_condition('sides', BoundaryType.NEUMANN, 0.0)   # Zero field at sides
    
    return solver


if __name__ == "__main__":
    # Demo usage
    from .materials import default_materials
    from .charge_density import ChargeDensityCalculator
    
    # Get silicon material
    si_n = default_materials.get_material("Si_n")
    
    # Create solver
    solver = setup_basic_stm_solver(si_n, tip_bias=1.0, separation=1.0)
    
    # Create charge density calculator
    charge_calc = ChargeDensityCalculator()
    
    # Find equilibrium Fermi level
    fermi_level = charge_calc.find_equilibrium_fermi_level(si_n)
    
    print(f"Grid: {solver.grid.nr} × {solver.grid.nphi} × {solver.grid.nz_total}")
    print(f"Equilibrium Fermi level: {fermi_level:.3f} eV")
    
    # Solve self-consistently
    potential, ef_final, info = solver.solve_self_consistent(
        charge_calc, {1: si_n}, fermi_level
    )
    
    print(f"Final solution:")
    print(f"  Iterations: {info['iterations']}")
    print(f"  Converged: {info['success']}")
    print(f"  Final residual: {info['final_residual']:.2e} V")
    print(f"  Potential range: [{np.min(potential):.3f}, {np.max(potential):.3f}] V")
