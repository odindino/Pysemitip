"""
Poisson equation solver for SEMITIP simulations.

This module implements the 3D Poisson equation solver from SEMITIP3,
using finite difference methods with hyperboloidal coordinate transformation
near the tip.
"""

import numpy as np
from typing import Tuple, Optional, Callable
from dataclasses import dataclass
import time

from ...utils.constants import PhysicalConstants as PC
from ..solvers.grid import Grid3D
from ..materials.tip import TipModel
from .charge_density import ChargeDensityTables


@dataclass
class PoissonSolverParameters:
    """Parameters for the Poisson solver."""
    # Convergence criteria
    tolerance: float = 1e-6         # Potential convergence tolerance (V)
    max_iterations: int = 1000      # Maximum iterations (reduced for stability)
    
    # Relaxation parameters
    omega: float = 0.8              # Over-relaxation parameter (balanced)
    adaptive_omega: bool = True     # Use adaptive relaxation
    
    # Solver options
    use_multigrid: bool = False     # Use multigrid acceleration
    verbose: bool = True            # Print progress
    
    # Boundary conditions
    boundary_type: str = 'dirichlet'  # 'dirichlet' or 'neumann'


class PoissonSolver:
    """
    3D Poisson equation solver for semiconductor-vacuum system.
    
    Implements the functionality of SEMITIP3 from the Fortran code.
    """
    
    def __init__(self, grid: Grid3D, tip: TipModel, 
                 params: Optional[PoissonSolverParameters] = None):
        """
        Initialize the Poisson solver.
        
        Args:
            grid: 3D computational grid
            tip: Tip model
            params: Solver parameters
        """
        self.grid = grid
        self.tip = tip
        self.params = params or PoissonSolverParameters()
        
        # Get hyperboloidal transformation parameters
        self.eta, self.a, self.z0, self.c = tip.hyperboloid_parameters()
        
        # Initialize potential arrays
        self._init_potential()
        
        # Set up finite difference coefficients
        self._setup_fd_coefficients()
    
    def _init_potential(self):
        """Initialize potential arrays."""
        # Vacuum potential (2 components for vector potential if needed)
        self.potential_vac = np.zeros_like(self.grid.vac[0])
        
        # Semiconductor potential
        self.potential_sem = np.zeros_like(self.grid.sem[0])
        
        # Interface potential
        self.potential_int = np.zeros_like(self.grid.vsint[0])
        
        # Combined potential for convenience
        nr, nv, ns, np_ = (self.grid.params.nr, self.grid.params.nv,
                          self.grid.params.ns, self.grid.params.np)
        self.potential_full = np.zeros((nr, nv + ns - 1, np_))
    
    def _setup_fd_coefficients(self):
        """Set up finite difference coefficients."""
        # Grid spacing
        dr = self.grid.params.delr
        dv = self.grid.params.delv
        ds = self.grid.params.dels
        dp = self.grid.params.delp
        
        # Coefficients for vacuum region (cylindrical coordinates)
        self.dr2_inv = 1.0 / dr**2
        self.dv2_inv = 1.0 / dv**2
        self.dp2_inv = 1.0 / dp**2
        
        # Coefficients for semiconductor region
        self.ds2_inv = 1.0 / ds**2
        
        # Radial derivative coefficients (1/r terms)
        self.r_coeff = np.zeros(self.grid.params.nr)
        for i in range(1, self.grid.params.nr):
            self.r_coeff[i] = 1.0 / (2.0 * dr * self.grid.r[i])
    
    def solve(self, charge_density_func: Callable, 
              surface_charge_func: Callable,
              initial_guess: Optional[np.ndarray] = None) -> Tuple[np.ndarray, dict]:
        """
        Solve the Poisson equation.
        
        Args:
            charge_density_func: Function to calculate bulk charge density
                                Args: (r, z, phi, potential) -> rho
            surface_charge_func: Function to calculate surface charge density
                                Args: (r, phi, potential) -> sigma
            initial_guess: Initial potential distribution
            
        Returns:
            Tuple of (converged potential, convergence info dict)
        """
        # Initialize potential
        if initial_guess is not None:
            self._set_potential_from_array(initial_guess)
        else:
            self._set_initial_potential()
        
        # Convergence tracking
        iteration = 0
        converged = False
        convergence_history = []
        start_time = time.time()
        
        # Check for zero charge density (common issue)
        test_charge = charge_density_func(1.0, -1.0, 0.0, 0.0)
        if abs(test_charge) < 1e-20:
            print("Warning: Charge density function returns near-zero values")
        
        # Main iteration loop
        while iteration < self.params.max_iterations and not converged:
            # Store old potential
            old_potential = self.potential_full.copy()
            
            # Update potential using successive over-relaxation (SOR)
            self._update_potential_sor(charge_density_func, surface_charge_func)
            
            # Check convergence
            max_change = np.max(np.abs(self.potential_full - old_potential))
            convergence_history.append(max_change)
            
            # Adaptive relaxation parameter
            if self.params.adaptive_omega and iteration > 10:
                self._update_omega(convergence_history)
            
            # Check for convergence
            if max_change < self.params.tolerance:
                converged = True
            
            # Check for divergence or numerical issues
            if not np.isfinite(max_change) or max_change > 1e10:
                print(f"Solution diverged at iteration {iteration}: max_change = {max_change:.2e}")
                break
            
            # Check for immediate convergence (potential issue)
            if iteration == 0 and max_change < self.params.tolerance * 100:
                print(f"Warning: Very fast convergence may indicate insufficient charge density")
            
            # Progress output
            if self.params.verbose and iteration % 100 == 0:
                pot0 = self.potential_int[0, 0]  # Potential at origin
                print(f"Iteration {iteration}: max_change = {max_change:.2e}, "
                      f"Pot(0,0) = {pot0:.6f} V")
            
            iteration += 1
        
        # Convergence info
        convergence_info = {
            'converged': converged,
            'iterations': iteration,
            'final_error': convergence_history[-1] if convergence_history else np.inf,
            'time': time.time() - start_time,
            'convergence_history': np.array(convergence_history)
        }
        
        if self.params.verbose:
            if converged:
                print(f"Converged after {iteration} iterations")
            else:
                print(f"Failed to converge after {iteration} iterations")
            print(f"Band bending at origin: {self.potential_int[0, 0]:.6f} V")
        
        return self.potential_full, convergence_info
    
    def _set_initial_potential(self):
        """Set initial potential distribution."""
        # Simple initial guess: linear drop from tip to ground
        # More sophisticated initialization could use 1D solution
        
        # Set tip potential
        tip_potential = self.tip.tip_potential
        
        # Vacuum region: exponential decay from tip
        for i in range(self.grid.params.nr):
            for j in range(self.grid.params.nv):
                for k in range(self.grid.params.np):
                    if self.grid.tip_mask[i, j, k]:
                        self.potential_vac[i, j, k] = tip_potential
                    else:
                        # Distance from tip apex
                        z = self.grid.zv[j]
                        decay_length = 10.0  # nm
                        self.potential_vac[i, j, k] = tip_potential * np.exp(-z / decay_length)
        
        # Semiconductor region: exponential decay into bulk
        for i in range(self.grid.params.nr):
            for j in range(self.grid.params.ns):
                for k in range(self.grid.params.np):
                    z = -self.grid.zs[j]  # Positive depth
                    decay_length = 50.0  # nm (Debye length estimate)
                    self.potential_sem[i, j, k] = self.potential_int[i, k] * np.exp(-z / decay_length)
        
        # Update full potential array
        self._update_full_potential()
    
    def _update_potential_sor(self, charge_density_func: Callable,
                            surface_charge_func: Callable):
        """
        Update potential using successive over-relaxation.
        
        This implements the finite difference scheme from SEMITIP3.
        """
        omega = self.params.omega
        
        # Update vacuum region
        self._update_vacuum_potential(omega)
        
        # Update semiconductor region
        self._update_semiconductor_potential(charge_density_func, omega)
        
        # Update interface (matching boundary conditions)
        self._update_interface_potential(surface_charge_func, omega)
        
        # Apply boundary conditions
        self._apply_boundary_conditions()
        
        # Update full potential array
        self._update_full_potential()
    
    def _update_vacuum_potential(self, omega: float):
        """Update potential in vacuum region."""
        nr, nv, np_ = self.grid.params.nr, self.grid.params.nv, self.grid.params.np
        
        # Skip boundaries and tip interior
        for i in range(1, nr - 1):
            for j in range(1, nv - 1):
                for k in range(np_):
                    if self.grid.tip_mask[i, j, k]:
                        # Inside tip - fixed potential
                        self.potential_vac[i, j, k] = self.tip.tip_potential
                        continue
                    
                    # Finite difference in cylindrical coordinates
                    # ∇²φ = ∂²φ/∂r² + (1/r)∂φ/∂r + (1/r²)∂²φ/∂φ² + ∂²φ/∂z²
                    
                    # Radial derivatives
                    d2phi_dr2 = (self.potential_vac[i+1, j, k] - 
                                2 * self.potential_vac[i, j, k] + 
                                self.potential_vac[i-1, j, k]) * self.dr2_inv
                    
                    dphi_dr = (self.potential_vac[i+1, j, k] - 
                              self.potential_vac[i-1, j, k]) * self.r_coeff[i]
                    
                    # Angular derivatives (with periodic BC)
                    k_next = (k + 1) % np_
                    k_prev = (k - 1) % np_
                    d2phi_dp2 = (self.potential_vac[i, j, k_next] - 
                                2 * self.potential_vac[i, j, k] + 
                                self.potential_vac[i, j, k_prev]) * self.dp2_inv
                    
                    # Vertical derivatives
                    d2phi_dz2 = (self.potential_vac[i, j+1, k] - 
                                2 * self.potential_vac[i, j, k] + 
                                self.potential_vac[i, j-1, k]) * self.dv2_inv
                    
                    # Laplacian (no charge in vacuum)
                    laplacian = d2phi_dr2 + dphi_dr
                    if i > 0:  # Avoid division by zero at r=0
                        laplacian += d2phi_dp2 / self.grid.r[i]**2
                    laplacian += d2phi_dz2
                    
                    # Update with SOR
                    residual = -laplacian
                    correction = omega * residual / (2 * (self.dr2_inv + self.dv2_inv + 
                                                         self.dp2_inv / max(self.grid.r[i]**2, 1e-10)))
                    self.potential_vac[i, j, k] += correction
    
    def _update_semiconductor_potential(self, charge_density_func: Callable, omega: float):
        """Update potential in semiconductor region."""
        nr, ns, np_ = self.grid.params.nr, self.grid.params.ns, self.grid.params.np
        
        for i in range(1, nr - 1):
            for j in range(1, ns - 1):
                for k in range(np_):
                    # Get charge density at this point
                    r = self.grid.r[i]
                    z = self.grid.zs[j]
                    phi = self.grid.phi[k]
                    pot = self.potential_sem[i, j, k]
                    
                    rho = charge_density_func(r, z, phi, pot)
                    
                    # Poisson equation: ∇²φ = -ρ/ε
                    # Similar finite difference as vacuum
                    d2phi_dr2 = (self.potential_sem[i+1, j, k] - 
                                2 * self.potential_sem[i, j, k] + 
                                self.potential_sem[i-1, j, k]) * self.dr2_inv
                    
                    dphi_dr = (self.potential_sem[i+1, j, k] - 
                              self.potential_sem[i-1, j, k]) * self.r_coeff[i]
                    
                    k_next = (k + 1) % np_
                    k_prev = (k - 1) % np_
                    d2phi_dp2 = (self.potential_sem[i, j, k_next] - 
                                2 * self.potential_sem[i, j, k] + 
                                self.potential_sem[i, j, k_prev]) * self.dp2_inv
                    
                    d2phi_dz2 = (self.potential_sem[i, j+1, k] - 
                                2 * self.potential_sem[i, j, k] + 
                                self.potential_sem[i, j-1, k]) * self.ds2_inv
                    
                    laplacian = d2phi_dr2 + dphi_dr
                    if i > 0:
                        laplacian += d2phi_dp2 / self.grid.r[i]**2
                    laplacian += d2phi_dz2
                    
                    # Include charge density
                    # Assuming uniform permittivity for now
                    eps_r = 12.9  # GaAs
                    source = rho / (PC.EPSILON0 * eps_r)
                    
                    # Update with SOR
                    residual = -laplacian - source
                    correction = omega * residual / (2 * (self.dr2_inv + self.ds2_inv + 
                                                         self.dp2_inv / max(self.grid.r[i]**2, 1e-10)))
                    self.potential_sem[i, j, k] += correction
    
    def _update_interface_potential(self, surface_charge_func: Callable, omega: float):
        """Update potential at vacuum-semiconductor interface."""
        nr, np_ = self.grid.params.nr, self.grid.params.np
        
        for i in range(1, nr - 1):
            for k in range(np_):
                # Surface charge density
                r = self.grid.r[i]
                phi = self.grid.phi[k]
                pot = self.potential_int[i, k]
                sigma = surface_charge_func(r, phi, pot)
                
                # Boundary condition: ε₁∂φ₁/∂z - ε₂∂φ₂/∂z = σ
                # Using finite differences
                
                # Vacuum side derivative
                dphi_dz_vac = (self.potential_vac[i, 1, k] - self.potential_int[i, k]) / self.grid.params.delv
                
                # Semiconductor side derivative
                dphi_dz_sem = (self.potential_int[i, k] - self.potential_sem[i, 1, k]) / self.grid.params.dels
                
                # Permittivities
                eps_vac = 1.0
                eps_sem = 12.9  # GaAs
                
                # Update interface potential
                residual = (eps_vac * dphi_dz_vac - eps_sem * dphi_dz_sem - 
                           sigma / PC.EPSILON0)
                
                # Weighted update based on neighboring points
                weight_vac = eps_vac / self.grid.params.delv
                weight_sem = eps_sem / self.grid.params.dels
                total_weight = weight_vac + weight_sem
                
                new_pot = (weight_vac * self.potential_vac[i, 1, k] + 
                          weight_sem * self.potential_sem[i, 1, k] - 
                          sigma / PC.EPSILON0) / total_weight
                
                self.potential_int[i, k] = (1 - omega) * self.potential_int[i, k] + omega * new_pot
    
    def _apply_boundary_conditions(self):
        """Apply boundary conditions."""
        # Radial boundaries (r = 0 and r = rmax)
        # At r = 0: ∂φ/∂r = 0 (symmetry)
        self.potential_vac[0, :, :] = self.potential_vac[1, :, :]
        self.potential_sem[0, :, :] = self.potential_sem[1, :, :]
        
        # At r = rmax: φ = 0 (ground)
        self.potential_vac[-1, :, :] = 0.0
        self.potential_sem[-1, :, :] = 0.0
        
        # Top boundary (vacuum): φ = 0
        self.potential_vac[:, -1, :] = 0.0
        
        # Bottom boundary (semiconductor): ∂φ/∂z = 0 (deep bulk)
        self.potential_sem[:, -1, :] = self.potential_sem[:, -2, :]
        
        # Angular boundaries (periodic)
        # Handled automatically by periodic indexing
    
    def _update_full_potential(self):
        """Update the combined potential array."""
        # Combine vacuum, interface, and semiconductor potentials
        nv = self.grid.params.nv
        ns = self.grid.params.ns
        
        # Semiconductor part (negative z)
        self.potential_full[:, :ns-1, :] = np.flip(self.potential_sem[:, 1:, :], axis=1)
        
        # Interface
        self.potential_full[:, ns-1, :] = self.potential_int
        
        # Vacuum part (positive z)
        self.potential_full[:, ns:, :] = self.potential_vac[:, 1:, :]
    
    def _update_omega(self, convergence_history: list):
        """Update relaxation parameter based on convergence history."""
        if len(convergence_history) < 20:
            return
        
        # Check if convergence is stalling
        recent = convergence_history[-10:]
        if np.std(recent) / np.mean(recent) < 0.1:
            # Convergence stalling - adjust omega
            if self.params.omega > 1.5:
                self.params.omega *= 0.95
            elif self.params.omega < 1.0:
                self.params.omega *= 1.05
    
    def _set_potential_from_array(self, potential: np.ndarray):
        """Set internal potentials from combined array."""
        ns = self.grid.params.ns
        
        # Split into regions
        self.potential_sem[:, 1:, :] = np.flip(potential[:, :ns-1, :], axis=1)
        self.potential_int[:, :] = potential[:, ns-1, :]
        self.potential_vac[:, 1:, :] = potential[:, ns:, :]
        
        # Set interface values
        self.potential_vac[:, 0, :] = self.potential_int
        self.potential_sem[:, 0, :] = self.potential_int
    
    def get_band_bending(self) -> float:
        """Get band bending at the origin (r=0, z=0)."""
        return self.potential_int[0, 0]
    
    def get_depletion_width(self, threshold: float = 0.01) -> float:
        """
        Estimate depletion width where potential drops to threshold * band_bending.
        
        Args:
            threshold: Fraction of band bending to define depletion edge
            
        Returns:
            Depletion width in nm
        """
        bb = self.get_band_bending()
        target = abs(bb) * threshold
        
        # Search along z-axis at r=0
        for j in range(self.grid.params.ns):
            if abs(self.potential_sem[0, j, 0]) < target:
                return -self.grid.zs[j]
        
        return -self.grid.zs[-1]  # Full depth if not found