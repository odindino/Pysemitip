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
              initial_guess_data: Optional[Tuple[Optional[np.ndarray], Optional['GridParameters']]] = None
              ) -> Tuple[np.ndarray, dict]:
        """
        Solve the Poisson equation using Fortran SEMITIP3 algorithm.
        
        This implements the exact logic from SEMITIP3-6.1.f lines 442-756.
        
        Args:
            charge_density_func: Function to calculate bulk charge density
                                Args: (r, z, phi, potential) -> rho
            surface_charge_func: Function to calculate surface charge density
                                Args: (r, phi, potential) -> sigma
            initial_guess_data: Tuple of (coarse_potential_array, coarse_grid_params).
                                If None, or coarse_potential_array is None, uses default initial potential.
                                If coarse_grid_params is None, assumes coarse_potential_array is for current grid.
            
        Returns:
            Tuple of (converged potential, convergence info dict)
        """
        # Initialize potential
        if initial_guess_data is not None:
            coarse_potential_array, coarse_grid_params = initial_guess_data
            if coarse_potential_array is not None:
                if coarse_grid_params is not None and coarse_grid_params.nr != self.grid.params.nr: # Check if grids differ
                    print(f"Prolonging potential from coarse grid {coarse_grid_params.nr}x{coarse_grid_params.nv+coarse_grid_params.ns-1} to fine grid {self.grid.params.nr}x{self.grid.params.nv+self.grid.params.ns-1}")
                    self._prolong_potential(coarse_potential_array, coarse_grid_params, self.grid.params)
                else: # Coarse_grid_params is None or same as current grid, assume direct use
                    print(f"Using initial guess on the same grid dimensions.")
                    self._set_potential_from_array(coarse_potential_array)
            else: # coarse_potential_array is None
                self._set_initial_potential()
        else: # initial_guess_data is None
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
            # The loop for 'i' should go up to nr-2 (i.e., range(nr-1))
            # to exclude the outer radial boundary i=nr-1 if it's a fixed Dirichlet boundary.
            for i in range(nr - 1):
                for j in range(1, nv - 1):  # Skip j boundaries (interface and far vacuum)
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
                    rho_surf_old = surface_charge_func(r, phi, surf_old)
                    rho_surf_new = surface_charge_func(r, phi, surf_new)
                    residual_old = surface_residual(surf_old)
                    residual_new = surface_residual(surf_new)
                    print(f"  Surface GSECT: old={surf_old:.6f}, new={surf_new:.6f}, change={abs(surf_new-surf_old):.2e}")
                    print(f"    rho_surf: old={rho_surf_old:.2e}, new={rho_surf_new:.2e}")
                    print(f"    residual: old={residual_old:.2e}, new={residual_new:.2e}")
                
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
                    
                    # Calculate linear solution first (Fortran SEMNEW = TEMP/DENOM)
                    # This assumes zero charge density for the linear solution
                    rho_linear = 0.0  # Linear approximation
                    temp_linear = stemp - rho_linear * eep / 12.9  # eps_semi = 12.9
                    sem_new_linear = temp_linear / denom
                    
                    # Nonlinear solver using golden section search
                    def bulk_residual(pot_test):
                        """Residual function for bulk potential."""
                        rho_bulk = charge_density_func(r, z, phi, pot_test)
                        temp = stemp - rho_bulk * eep / 12.9  # eps_semi = 12.9
                        new_pot = temp / denom
                        return abs(pot_test - new_pot)
                    
                    # Set search bounds (correct Fortran approach)
                    # Use sem_old and sem_new_linear as bounds (SEMOLD, SEMNEW)
                    sem_min = min(sem_old, sem_new_linear)
                    sem_max = max(sem_old, sem_new_linear)
                    
                    # delta_sem is used only as convergence tolerance
                    bias = self.tip.bias_voltage
                    delta_sem = max(1e-6, abs(bias) / 1e6)  # Convergence tolerance only
                    
                    # Use golden section search with correct bounds
                    sem_new = golden_section_search(bulk_residual, sem_min, sem_max, 
                                                   delta_sem)  # Use delta_sem as tolerance
                    
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
        ETAT = self.eta  # This is the tip's hyperboloid parameter
        NV = self.grid.params.nv
        NR = self.grid.params.nr
        NP = self.grid.params.np

        if NV == 0: # Should not happen with valid grid params
            # Or handle as error, for now fill with tip_potential if vacuum region exists
            if self.potential_vac is not None:
                 self.potential_vac.fill(tip_potential)
        else:
            deleta = ETAT / float(NV)

            # Pre-calculate denominator for log formula
            log_ETAT_term_numerator = (1. + ETAT)
            log_ETAT_term_denominator = (1. - ETAT)

            if log_ETAT_term_denominator <= 1e-9: # ETAT is very close to 1 or >= 1
                cetat_inv = 0 # Formula is ill-defined or ETAT is too large
            else:
                log_ETAT_term = log_ETAT_term_numerator / log_ETAT_term_denominator
                if log_ETAT_term <= 0:
                    cetat_inv = 0 # Avoid log(<=0)
                else:
                    log_cetat = np.log(log_ETAT_term)
                    if abs(log_cetat) < 1e-9: # Denominator is close to zero
                        cetat_inv = 0 # effectively, potential will be tip_potential
                    else:
                        cetat_inv = 1.0 / log_cetat

            for k_idx in range(NP):
                for i_idx in range(NR):
                    for j_python in range(NV): # j_python is 0 to NV-1 (interface to far vacuum)
                        if self.grid.tip_mask[i_idx, j_python, k_idx]:
                            self.potential_vac[i_idx, j_python, k_idx] = tip_potential
                        else:
                            j_fortran = j_python + 1 # Fortran index for eta calculation (1 to NV)
                            eta_coord = j_fortran * deleta

                            # Ensure eta_coord is strictly less than ETAT for the formula
                            # Using a small factor like 0.99999 to avoid floating point issues at eta_coord == ETAT
                            if eta_coord < ETAT * 0.999999 and cetat_inv != 0:
                                log_eta_coord_term_numerator = (1. + eta_coord)
                                log_eta_coord_term_denominator = (1. - eta_coord)

                                if log_eta_coord_term_denominator <= 1e-9: # eta_coord is very close to 1
                                    # This case implies potential should be tip_potential or very close
                                    self.potential_vac[i_idx, j_python, k_idx] = tip_potential
                                else:
                                    log_eta_coord_term = log_eta_coord_term_numerator / log_eta_coord_term_denominator
                                    if log_eta_coord_term <= 0:
                                        # Fallback if term is problematic (e.g. eta_coord >= 1)
                                        self.potential_vac[i_idx, j_python, k_idx] = tip_potential
                                    else:
                                        potential = tip_potential * np.log(log_eta_coord_term) * cetat_inv
                                        self.potential_vac[i_idx, j_python, k_idx] = potential
                            else:
                                # If eta_coord >= ETAT or cetat_inv is 0 (ETAT is problematic)
                                # point is effectively at or beyond the tip's reference hyperboloid for this guess
                                self.potential_vac[i_idx, j_python, k_idx] = tip_potential
        
        # Semiconductor region: Start with zero (Fortran lines 170-178)
        self.potential_sem.fill(0.0)
        
        # Interface: Initialize to zero following Fortran IINIT=1 logic
        # Fortran SEMITIP3-6.1.f lines 155-162: VSINT(1,I,K)=0.
        self.potential_int.fill(0.0)
        
        # Set far-field boundary conditions explicitly as per typical SEMITIP3 for IINIT=1
        # These are typically Dirichlet boundaries for the initial guess.
        if self.potential_vac is not None and self.potential_vac.size > 0:
            if NV > 0:
                self.potential_vac[:, -1, :] = 0.0  # Top vacuum boundary (j_python = NV-1) to 0
            if NR > 0:
                self.potential_vac[-1, :, :] = 0.0  # Outer radial vacuum boundary (i_idx = NR-1) to 0

        if self.potential_sem is not None and self.potential_sem.size > 0:
            NS = self.grid.params.ns
            if NR > 0:
                self.potential_sem[-1, :, :] = 0.0  # Outer radial semiconductor boundary to 0
            if NS > 0:
                self.potential_sem[:, -1, :] = 0.0  # Deep semiconductor boundary (z -> -inf) to 0

        print(f"Initialized interface potential to zero (Fortran IINIT=1 style)")
        print(f"Applied far-field zero boundary conditions for initial potential.")

        # Update full potential array
        self._update_full_potential()

    def _prolong_potential(self, coarse_potential_full: np.ndarray,
                           coarse_params: 'GridParameters',
                           fine_params: 'GridParameters'):
        """
        Initialize potential on the current fine grid from a coarser grid solution.
        Uses simple injection for now (coarse[i] -> fine[2i]).
        Assumes dimensions are approximately doubled.
        """
        print(f"  Prolongation: Coarse ({coarse_params.nr},{coarse_params.nv},{coarse_params.ns},{coarse_params.np}) to Fine ({fine_params.nr},{fine_params.nv},{fine_params.ns},{fine_params.np})")

        # Deconstruct coarse_potential_full into vac, sem, int parts
        # This logic is reverse of _update_full_potential or similar to _set_potential_from_array

        # coarse_potential_full has shape (NR_c, NV_c + NS_c - 1, NP_c)
        # self.potential_vac has shape (NR_f, NV_f, NP_f)
        # self.potential_sem has shape (NR_f, NS_f, NP_f)
        # self.potential_int has shape (NR_f, NP_f)

        # Approximate scaling factors (should be close to 2)
        scale_nr = fine_params.nr / coarse_params.nr
        scale_nv = fine_params.nv / coarse_params.nv
        scale_ns = fine_params.ns / coarse_params.ns
        scale_np = fine_params.np / coarse_params.np

        # Initialize fine potentials to zero or some default
        self.potential_vac.fill(0.0)
        self.potential_sem.fill(0.0)
        self.potential_int.fill(0.0)

        # Simplified injection:
        # Iterate over coarse grid indices
        for ic in range(coarse_params.nr):
            ifc = int(ic * scale_nr) # Map to fine grid index
            if ifc >= fine_params.nr: continue

            for kpc in range(coarse_params.np):
                kfpc = int(kpc * scale_np)
                if kfpc >= fine_params.np: continue

                # Interface potential
                # coarse_potential_full[:, coarse_ns-1, :] is coarse interface
                self.potential_int[ifc, kfpc] = coarse_potential_full[ic, coarse_params.ns - 1, kpc]

                # Semiconductor part (excluding interface)
                # coarse_potential_full[:, :coarse_ns-1, :] is coarse semiconductor (flipped)
                # self.potential_sem[:, 1:, :] is fine semiconductor bulk (flipped in _set_potential_from_array)
                for jsc in range(coarse_params.ns -1): # Iterate 0 to NS_c-2 (bulk coarse points)
                    jsfc = int((jsc + 1) * scale_ns) # Map to fine grid index (from 1 to NS_f-1 for bulk)
                    if jsfc >= fine_params.ns: continue
                    # coarse_potential_full is indexed [r, z_combined, p]
                    # z_combined for sem runs from 0 (deepest) to ns-2 (near interface)
                    # self.potential_sem is indexed [r, z_sem_idx, p] where z_sem_idx 0 is interface, ns-1 is deep
                    # coarse_potential_full[ic, jsc, kpc] corresponds to coarse_sem_bulk[ic, jsc_flipped, kpc]
                    # coarse_sem_bulk has ns-1 points. zs_idx_coarse_flipped = (coarse_params.ns-2) - jsc
                    idx_coarse_z_in_full_array = jsc # This is index in the coarse_potential_full sem part
                    self.potential_sem[ifc, jsfc, kfpc] = coarse_potential_full[ic, idx_coarse_z_in_full_array, kpc]


                # Vacuum part (excluding interface)
                # coarse_potential_full[:, coarse_ns:, :] is coarse vacuum
                # self.potential_vac[:, 1:, :] is fine vacuum bulk
                for jvc in range(coarse_params.nv -1): # Iterate 0 to NV_c-2 (bulk coarse points in vac part of full)
                    jvfc = int((jvc + 1) * scale_nv) # Map to fine grid index (from 1 to NV_f-1 for bulk)
                    if jvfc >= fine_params.nv: continue
                    idx_coarse_z_in_full_array = coarse_params.ns + jvc
                    self.potential_vac[ifc, jvfc, kfpc] = coarse_potential_full[ic, idx_coarse_z_in_full_array, kpc]

        # Ensure interface values are consistent between sem and vac arrays on fine grid
        self.potential_vac[:, 0, :] = self.potential_int
        self.potential_sem[:, 0, :] = self.potential_int

        # Apply boundary conditions that might have been overwritten or not set by injection
        # Using similar logic from _set_initial_potential for far-field boundaries
        if fine_params.nv > 0:
            self.potential_vac[:, -1, :] = 0.0  # Top vacuum boundary
        if fine_params.nr > 0:
            self.potential_vac[-1, :, :] = 0.0  # Outer radial vacuum boundary

        if fine_params.nr > 0:
            self.potential_sem[-1, :, :] = 0.0  # Outer radial semiconductor boundary
        if fine_params.ns > 0:
            self.potential_sem[:, -1, :] = 0.0  # Deep semiconductor boundary

        # Tip interior should be set to tip_potential
        tip_potential = self.tip.tip_potential
        for k_idx in range(fine_params.np):
            for i_idx in range(fine_params.nr):
                for j_python in range(fine_params.nv):
                    if self.grid.tip_mask[i_idx, j_python, k_idx]:
                        self.potential_vac[i_idx, j_python, k_idx] = tip_potential

        self._update_full_potential() # Ensure self.potential_full is updated
        print(f"  Prolongation completed. Band bending after prolong: {self.get_band_bending():.6f}")


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