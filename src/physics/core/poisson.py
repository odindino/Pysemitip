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


def golden_section_search(func: Callable, xmin: float, xmax: float, 
                         tolerance: float = 1e-6, max_iter: int = 100) -> float:
    """
    Golden section search for finding minimum of a function.
    
    This implements the GSECT subroutine from Fortran SEMITIP.
    
    Args:
        func: Function to minimize
        xmin: Lower bound
        xmax: Upper bound
        tolerance: Convergence tolerance
        max_iter: Maximum iterations
        
    Returns:
        Optimal x value where function is minimized
    """
    if abs(xmax - xmin) < tolerance:
        return (xmin + xmax) / 2.0
    
    # Ensure xmin < xmax
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    
    # Golden ratio constant
    gs = 0.3819660  # (3 - sqrt(5)) / 2
    
    delx = xmax - xmin
    xa = xmin + delx * gs
    fa = func(xa)
    xb = xmax - delx * gs
    fb = func(xb)
    
    for _ in range(max_iter):
        delx_saved = delx
        if delx < tolerance:
            break
            
        if fb < fa:
            # Move to the right interval
            xmin = xa
            delx = xmax - xmin
            if abs(delx - delx_saved) < tolerance * tolerance:
                break
            xa = xb
            fa = fb
            xb = xmax - delx * gs
            fb = func(xb)
        else:
            # Move to the left interval
            xmax = xb
            delx = xmax - xmin
            if abs(delx - delx_saved) < tolerance * tolerance:
                break
            xb = xa
            fb = fa
            xa = xmin + delx * gs
            fa = func(xa)
    
    return (xmin + xmax) / 2.0


@dataclass
class PoissonSolverParameters:
    """Parameters for the Poisson solver matching Fortran SEMITIP3."""
    # Convergence criteria (from Fortran configuration)
    tolerance: float = 1e-5         # EP parameter from Fortran (less strict)
    max_iterations: int = 5000      # ITMAX from Fortran (match Fortran behavior)
    
    # Relaxation parameters
    omega: float = 0.5              # Conservative SOR parameter for stability (Fortran-like)
    adaptive_omega: bool = False    # Keep simple for now
    
    # Nonlinear solver parameters
    nonlinear_tolerance: float = 1e-6  # DELSEM/DELSURF from Fortran
    golden_section_tolerance: float = 1e-6  # For GSECT (less strict)
    
    # Convergence checking (Fortran style)
    convergence_check_interval: int = 100   # Check every 100 iterations
    band_bending_tolerance: float = 1e-4    # Tighter tolerance matching Fortran convergence
    
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
        Solve the Poisson equation using Fortran SEMITIP3 algorithm.
        
        This implements the exact logic from SEMITIP3-6.1.f lines 442-756.
        
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
        
        # Convergence tracking (Fortran style)
        iteration = 0
        converged = False
        convergence_history = []
        start_time = time.time()
        
        # Band bending tracking (Pot0 in Fortran)
        pot0_current = self.get_band_bending()
        pot0_saved = 0.0
        pot0_saved2 = 0.0
        
        # Check for zero charge density (common issue)
        test_charge = charge_density_func(1.0, -1.0, 0.0, 0.0)
        if abs(test_charge) < 1e-20:
            print("Warning: Charge density function returns near-zero values")
        
        # Debug initial state
        if self.params.verbose:
            print(f"Initial Pot0: {pot0_current:.8f}")
            print(f"Tip potential: {self.tip.tip_potential:.6f}")
            print(f"Test charge density: {test_charge:.2e}")
            print(f"Interface potential shape: {self.potential_int.shape}")
            print(f"Interface potential stats: min={np.min(self.potential_int):.6f}, max={np.max(self.potential_int):.6f}, mean={np.mean(self.potential_int):.6f}")
        
        # Constants (from Fortran)
        eep = 1.80943e-20  # e/epsilon_0 in units of V cm, times 1.e-14 cm^2/nm^2
        
        # Main iteration loop (following Fortran DO 500 ITER=1,ITMAX)
        while iteration < self.params.max_iterations and not converged:
            iteration += 1
            
            # 1. Update vacuum potential (Fortran lines 446-567)
            self._update_vacuum_fortran_style()
            
            # 2. Update interface potential with nonlinear solver (Fortran lines 578-622)
            self._update_interface_fortran_style(surface_charge_func, eep)
            
            # 3. Update semiconductor potential with nonlinear solver (Fortran lines 633-731)  
            self._update_semiconductor_fortran_style(charge_density_func, eep)
            
            # 4. Check convergence every 100 iterations (Fortran lines 623-629, 750-751)
            if iteration % self.params.convergence_check_interval == 0:
                pot0_saved2 = pot0_saved
                pot0_saved = pot0_current
                pot0_current = self.get_band_bending()
                
                if self.params.verbose:
                    print(f"ITER,Pot0 = {iteration} {pot0_current:.8f}")
                
                # Fortran convergence check: both current and previous change must be small
                pot_change1 = abs(pot0_current - pot0_saved)
                pot_change2 = abs(pot0_saved - pot0_saved2)
                
                # Exact Fortran convergence condition (from ITER3 line ~500)
                # IF ((MOD(ITER,100).EQ.0.AND.ABS(Pot0-PotSAV).LT.EP.AND.ABS(PotSAV-PotSAV2).LT.2.*EP))
                if (pot_change1 < self.params.band_bending_tolerance and 
                    pot_change2 < 2.0 * self.params.band_bending_tolerance):
                    converged = True
                    if self.params.verbose:
                        print(f"Converged (Fortran-style): change1={pot_change1:.2e}, change2={pot_change2:.2e}")
                
                convergence_history.append(pot_change1)
            
            # Safety check for divergence
            if not np.isfinite(pot0_current) or abs(pot0_current) > 100:
                print(f"Solution diverged at iteration {iteration}: Pot0 = {pot0_current:.2e}")
                break
        
        # Final update of combined potential
        self._update_full_potential()
        
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
            print(f"Final band bending: {pot0_current:.8f} V")
        
        return self.potential_full, convergence_info
    
    def _update_vacuum_fortran_style(self):
        """Update vacuum potential following Fortran SEMITIP3 logic."""
        nr, nv, np_ = self.grid.params.nr, self.grid.params.nv, self.grid.params.np
        
        # Create working arrays (following Fortran VAC(2,...) pattern)
        new_potential = self.potential_vac.copy()
        
        # Update vacuum region (Fortran lines 446-567)
        for k in range(np_):
            for i in range(nr):
                for j in range(1, nv - 1):  # Skip boundaries
                    # Skip if inside tip
                    if self.grid.tip_mask[i, j, k]:
                        continue
                    
                    # Finite difference for Laplace equation in vacuum
                    # d²φ/dr² + (1/r)dφ/dr + (1/r²)d²φ/dφ² + d²φ/dz² = 0
                    
                    # Radial second derivative
                    if i == 0:
                        # At r=0, use symmetry
                        d2phi_dr2 = 2 * (self.potential_vac[1, j, k] - self.potential_vac[0, j, k]) * self.dr2_inv
                        dphi_dr = 0.0  # Symmetry at r=0
                    elif i == nr - 1:
                        # At boundary, use one-sided difference
                        d2phi_dr2 = 2 * (self.potential_vac[i-1, j, k] - self.potential_vac[i, j, k]) * self.dr2_inv
                        dphi_dr = (self.potential_vac[i, j, k] - self.potential_vac[i-1, j, k]) * self.r_coeff[i]
                    else:
                        # Interior points
                        d2phi_dr2 = (self.potential_vac[i+1, j, k] - 2*self.potential_vac[i, j, k] + 
                                    self.potential_vac[i-1, j, k]) * self.dr2_inv
                        dphi_dr = (self.potential_vac[i+1, j, k] - self.potential_vac[i-1, j, k]) * self.r_coeff[i]
                    
                    # Z derivative (vacuum coordinates)
                    d2phi_dz2 = (self.potential_vac[i, j+1, k] - 2*self.potential_vac[i, j, k] + 
                                self.potential_vac[i, j-1, k]) * self.dv2_inv
                    
                    # Angular derivative (with periodic boundary conditions)
                    k_next = (k + 1) % np_
                    k_prev = (k - 1) % np_
                    if i > 0:
                        d2phi_dp2 = (self.potential_vac[i, j, k_next] - 2*self.potential_vac[i, j, k] + 
                                    self.potential_vac[i, j, k_prev]) * self.dp2_inv / self.grid.r[i]**2
                    else:
                        d2phi_dp2 = 0.0  # At r=0, angular terms vanish
                    
                    # Laplacian with numerical stability checks
                    laplacian = d2phi_dr2 + dphi_dr + d2phi_dp2 + d2phi_dz2
                    
                    # Check for numerical issues
                    if not np.isfinite(laplacian):
                        continue  # Skip this point if numerical issues
                    
                    # Update (simple relaxation for vacuum)
                    denominator = 2 * self.dr2_inv + 2 * self.dv2_inv
                    if i > 0:
                        denominator += 2 * self.dp2_inv / self.grid.r[i]**2
                    
                    if abs(denominator) < 1e-15:  # Avoid division by zero
                        continue
                        
                    correction = -laplacian / denominator
                    
                    # Limit correction size for stability
                    max_correction = 0.1  # Volts
                    correction = np.clip(correction, -max_correction, max_correction)
                    
                    new_potential[i, j, k] = self.potential_vac[i, j, k] + self.params.omega * correction
        
        # Copy back to main array
        self.potential_vac = new_potential
        
        # Update interface boundary for vacuum
        self.potential_vac[:, 0, :] = self.potential_int
    
    def _update_interface_fortran_style(self, surface_charge_func: Callable, eep: float):
        """Update interface potential with nonlinear solver (Fortran style)."""
        nr, np_ = self.grid.params.nr, self.grid.params.np
        
        for k in range(np_):
            for i in range(1, nr - 1):  # Skip boundaries
                # Skip if tip is present
                if self.grid.tip_mask[i, 2, k]:  # Check nearby vacuum point
                    continue
                
                # Get coordinates
                r = self.grid.r[i]
                phi = (k - 0.5) * self.grid.params.delp  # Fortran style indexing
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                
                # Store old potential for nonlinear solver
                surf_old = self.potential_int[i, k]
                
                # Calculate finite difference terms (Fortran lines 585-590)
                # 3rd order in vacuum and semiconductor
                stemp = (3.0 * self.potential_vac[i, 1, k] - 
                        (9.0/6.0) * self.potential_vac[i, 2, k] + 
                        (1.0/3.0) * self.potential_vac[i, 3, k]) / self.grid.params.delv
                
                # Add semiconductor contribution
                eps_semi = 12.9  # GaAs permittivity
                stemp += eps_semi * (3.75 * self.potential_sem[i, 1, k] - 
                                   (5.0/6.0) * self.potential_sem[i, 2, k] + 
                                   0.15 * self.potential_sem[i, 3, k]) / self.grid.params.dels
                
                denom = ((11.0/6.0) / self.grid.params.delv + 
                        (46.0/15.0) * eps_semi / self.grid.params.dels)
                
                # Nonlinear solver using golden section search
                def surface_residual(pot_test):
                    """Residual function for surface potential."""
                    rho_surf = surface_charge_func(r, phi, pot_test)
                    temp = stemp - rho_surf * eep * 1e7
                    new_pot = temp / denom
                    return abs(pot_test - new_pot)
                
                # Set search bounds (exact Fortran DELSURF)
                bias = self.tip.bias_voltage
                delta_surf = max(1e-6, abs(bias) / 1e6)  # Exact Fortran DELSURF
                surf_min = surf_old - delta_surf
                surf_max = surf_old + delta_surf
                
                # Use golden section search (GSECT equivalent)
                surf_new = golden_section_search(surface_residual, surf_min, surf_max, 
                                               self.params.golden_section_tolerance)
                
                # Debug: check if potential is changing
                if i == 0 and k == 0:  # Central point
                    print(f"  Surface debug: old={surf_old:.6f}, new={surf_new:.6f}, rho_surf={surface_charge_func(r, phi, surf_old):.2e}")
                
                # Average old and new (Fortran style)
                self.potential_int[i, k] = (surf_old + surf_new) / 2.0
    
    def _update_semiconductor_fortran_style(self, charge_density_func: Callable, eep: float):
        """Update semiconductor potential with nonlinear solver (Fortran style)."""
        nr, ns, np_ = self.grid.params.nr, self.grid.params.ns, self.grid.params.np
        
        # Create working array
        new_potential = self.potential_sem.copy()
        
        for k in range(np_):
            for j in range(1, ns - 1):  # Skip boundaries
                for i in range(1, nr - 1):  # Skip boundaries
                    # Get coordinates
                    r = self.grid.r[i]
                    z = self.grid.zs[j]
                    phi = (k - 0.5) * self.grid.params.delp
                    x = r * np.cos(phi)
                    y = r * np.sin(phi)
                    
                    # Store old potential
                    sem_old = self.potential_sem[i, j, k]
                    
                    # Calculate finite difference terms (Fortran lines 704-708)
                    # Radial derivatives
                    d2phi_dr2 = 2.0 * (self.potential_sem[i+1, j, k] / self.grid.params.delr + 
                                      self.potential_sem[i-1, j, k] / self.grid.params.delr) / (2 * self.grid.params.delr)
                    
                    # Z derivatives
                    d2phi_dz2 = 2.0 * (self.potential_sem[i, j+1, k] / self.grid.params.dels + 
                                      self.potential_sem[i, j-1, k] / self.grid.params.dels) / (2 * self.grid.params.dels)
                    
                    # Radial gradient term
                    dphi_dr = (self.potential_sem[i+1, j, k] - self.potential_sem[i-1, j, k]) / (r * 2 * self.grid.params.delr)
                    
                    # Angular derivatives
                    k_next = (k + 1) % np_
                    k_prev = (k - 1) % np_
                    d2phi_dp2 = (self.potential_sem[i, j, k_next] + self.potential_sem[i, j, k_prev]) / (r**2 * self.grid.params.delp**2)
                    
                    stemp = d2phi_dr2 + d2phi_dz2 + dphi_dr + d2phi_dp2
                    
                    # Check for numerical issues in semiconductor region
                    if not np.isfinite(stemp):
                        continue  # Skip this point if numerical issues
                    
                    # Denominator for finite difference
                    denom = (2.0 * (1.0/self.grid.params.delr + 1.0/self.grid.params.delr) / (2 * self.grid.params.delr) + 
                            2.0 * (1.0/self.grid.params.dels + 1.0/self.grid.params.dels) / (2 * self.grid.params.dels) + 
                            2.0 / (r**2 * self.grid.params.delp**2))
                    
                    if abs(denom) < 1e-15:  # Avoid division by zero
                        continue
                    
                    # Nonlinear solver using golden section search
                    def bulk_residual(pot_test):
                        """Residual function for bulk potential."""
                        rho_bulk = charge_density_func(r, z, phi, pot_test)
                        temp = stemp - rho_bulk * eep / 12.9  # eps_semi = 12.9
                        new_pot = temp / denom
                        return abs(pot_test - new_pot)
                    
                    # Set search bounds (exact Fortran DELSEM)
                    bias = self.tip.bias_voltage
                    delta_sem = max(1e-6, abs(bias) / 1e6)  # Exact Fortran DELSEM
                    sem_min = sem_old - delta_sem
                    sem_max = sem_old + delta_sem
                    
                    # Use golden section search
                    sem_new = golden_section_search(bulk_residual, sem_min, sem_max, 
                                                   self.params.golden_section_tolerance)
                    
                    # Average old and new (Fortran style)
                    new_potential[i, j, k] = (sem_old + sem_new) / 2.0
        
        # Copy back to main array
        self.potential_sem = new_potential
        
        # Update interface boundary for semiconductor
        self.potential_sem[:, 0, :] = self.potential_int
    
    def _set_initial_potential(self):
        """Set initial potential distribution following Fortran SEMITIP3."""
        # Set tip potential in vacuum (Fortran lines 195-206)
        tip_potential = self.tip.tip_potential
        
        # Vacuum region: Initialize with logarithmic potential (Fortran style)
        for k in range(self.grid.params.np):
            for i in range(self.grid.params.nr):
                # Find where tip ends in vacuum
                j_tip_end = 0
                for j in range(self.grid.params.nv - 1, 0, -1):
                    if not self.grid.tip_mask[i, j, k]:
                        j_tip_end = j
                        break
                
                # Set initial guess (Fortran lines 198-204)
                for j in range(j_tip_end, 0, -1):
                    if not self.grid.tip_mask[i, j, k]:
                        # Logarithmic decay from tip (matching Fortran)
                        eta = j * self.eta / float(j_tip_end + 1)
                        if eta < 0.99:  # Avoid singularity
                            cetat = np.log((1. + self.eta) / (1. - self.eta))
                            potential = tip_potential * np.log((1. + eta) / (1. - eta)) / cetat
                        else:
                            potential = tip_potential * 0.1  # Small value near boundaries
                        self.potential_vac[i, j, k] = potential
                    else:
                        # Inside tip
                        self.potential_vac[i, j, k] = tip_potential
        
        # Semiconductor region: Start with zero (Fortran lines 170-178)
        self.potential_sem.fill(0.0)
        
        # Interface: Initialize to zero following Fortran IINIT=1 logic
        # Fortran SEMITIP3-6.1.f lines 155-162: VSINT(1,I,K)=0.
        # This allows the solver to find the correct equilibrium state
        self.potential_int.fill(0.0)
        
        print(f"Initialized interface potential to zero (Fortran IINIT=1 style)")
        
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
        """Get band bending at the origin (r=0, z=0) following Fortran PCENT function."""
        # Fortran PCENT function for JJ=0 (surface potential)
        # SUM=SUM+(9.*VSINT(1,I,K)-VSINT(1,I+1,K))/8.
        np_ = self.grid.params.np
        sum_pot = 0.0
        
        i = 0  # r=0 (first radial point)
        for k in range(np_):
            # Fortran uses I+1 as second point, but arrays are 1-indexed there
            # In Python 0-indexed: I=0, I+1=1
            if self.grid.params.nr > 1:
                weighted_pot = (9.0 * self.potential_int[i, k] - self.potential_int[i+1, k]) / 8.0
            else:
                weighted_pot = self.potential_int[i, k]
            sum_pot += weighted_pot
        
        band_bending = sum_pot / float(np_)
        
        # TODO: Check if sign convention matches Fortran
        # Current Python gives negative values while Fortran gives positive
        # This may require adjusting the reference level or potential initialization
        return band_bending
    
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