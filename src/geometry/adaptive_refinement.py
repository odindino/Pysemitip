"""
Adaptive Grid Refinement Module

This module implements the three-level adaptive grid refinement strategy
used in SEMITIP for achieving convergent solutions. It manages the iterative
refinement process with proper convergence checking.

Based on SEMITIP semitip3-6.1.f refinement implementation.

Author: odindino
Date: 2025-06-11
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union, Any
from dataclasses import dataclass, field
from enum import Enum
import warnings
import time
from copy import deepcopy

from .grid3d import Grid3D, GridConfig
from .stm_geometry import STMGeometry


class RefinementStatus(Enum):
    """Refinement process status"""
    NOT_STARTED = "not_started"
    IN_PROGRESS = "in_progress"
    CONVERGED = "converged"
    MAX_ITERATIONS = "max_iterations"
    FAILED = "failed"


@dataclass
class RefinementConfig:
    """Configuration for adaptive grid refinement"""
    
    # Maximum refinement levels (IPMAX in Fortran)
    max_refinement_levels: int = 3
    
    # Convergence criteria for each level (EP array in Fortran)
    convergence_tolerance: List[float] = field(default_factory=lambda: [1e-4, 1e-5, 1e-6])
    
    # Maximum iterations for each level (ITMAX array in Fortran)
    max_iterations_per_level: List[int] = field(default_factory=lambda: [50, 100, 200])
    
    # Under-relaxation factors for each level
    relaxation_factors: List[float] = field(default_factory=lambda: [0.7, 0.8, 0.9])
    
    # Convergence check interval
    convergence_check_interval: int = 5
    
    # Enable detailed logging
    verbose: bool = True
    
    # Callback functions for monitoring
    level_start_callback: Optional[Callable] = None
    iteration_callback: Optional[Callable] = None
    convergence_callback: Optional[Callable] = None
    
    def __post_init__(self):
        """Validate refinement configuration"""
        if self.max_refinement_levels < 1:
            raise ValueError("Must have at least 1 refinement level")
            
        # Ensure all arrays have correct length
        required_length = self.max_refinement_levels
        
        if len(self.convergence_tolerance) != required_length:
            if len(self.convergence_tolerance) == 1:
                # Extend single value
                self.convergence_tolerance = self.convergence_tolerance * required_length
            else:
                raise ValueError(f"convergence_tolerance must have length {required_length}")
                
        if len(self.max_iterations_per_level) != required_length:
            if len(self.max_iterations_per_level) == 1:
                self.max_iterations_per_level = self.max_iterations_per_level * required_length
            else:
                raise ValueError(f"max_iterations_per_level must have length {required_length}")
                
        if len(self.relaxation_factors) != required_length:
            if len(self.relaxation_factors) == 1:
                self.relaxation_factors = self.relaxation_factors * required_length
            else:
                raise ValueError(f"relaxation_factors must have length {required_length}")


@dataclass
class RefinementState:
    """Current state of refinement process"""
    current_level: int = 0
    current_iteration: int = 0
    status: RefinementStatus = RefinementStatus.NOT_STARTED
    converged_levels: List[bool] = field(default_factory=list)
    convergence_history: List[List[float]] = field(default_factory=list)
    iteration_times: List[List[float]] = field(default_factory=list)
    total_time: float = 0.0
    last_residual: float = float('inf')
    grid_levels: List[Grid3D] = field(default_factory=list)


class AdaptiveGridRefinement:
    """
    Adaptive Grid Refinement Manager
    
    This class implements the SEMITIP three-level adaptive refinement strategy:
    1. Start with coarse grid, iterate to convergence
    2. Refine grid (double resolution), interpolate previous solution as initial guess
    3. Iterate to convergence on refined grid
    4. Repeat for maximum 3 levels
    
    The refinement process is crucial for achieving accurate solutions in SEMITIP.
    """
    
    def __init__(self, 
                 initial_grid: Grid3D,
                 geometry: STMGeometry, 
                 config: Optional[RefinementConfig] = None):
        """
        Initialize adaptive refinement manager
        
        Args:
            initial_grid: Starting coarse grid
            geometry: STM geometry configuration
            config: Refinement configuration
        """
        self.config = config or RefinementConfig()
        self.geometry = geometry
        self.state = RefinementState()
        
        # Initialize with coarse grid
        self.state.grid_levels = [initial_grid]
        self.current_grid = initial_grid
        
        # Initialize convergence tracking
        self.state.converged_levels = [False] * self.config.max_refinement_levels
        self.state.convergence_history = [[] for _ in range(self.config.max_refinement_levels)]
        self.state.iteration_times = [[] for _ in range(self.config.max_refinement_levels)]
        
        # Solution arrays for each refinement level
        self.solutions = []  # Will store potential arrays
        
    def refine_grid_sequence(self, 
                           solver_function: Callable,
                           initial_solution: Optional[np.ndarray] = None) -> Tuple[bool, Dict]:
        """
        Execute complete adaptive refinement sequence
        
        Args:
            solver_function: Function that solves Poisson equation on current grid
                           Should have signature: solver_function(grid, initial_guess) -> solution
            initial_solution: Initial guess for coarsest grid (optional)
            
        Returns:
            (success, results): Success flag and detailed results dictionary
        """
        start_time = time.time()
        self.state.status = RefinementStatus.IN_PROGRESS
        
        try:
            results = {
                'convergence_achieved': [],
                'final_residuals': [],
                'iteration_counts': [],
                'refinement_times': [],
                'total_time': 0.0,
                'grids_used': [],
                'final_solution': None
            }
            
            current_solution = initial_solution
            
            # Loop through refinement levels
            for level in range(self.config.max_refinement_levels):
                self.state.current_level = level
                
                if self.config.verbose:
                    print(f"\n=== Refinement Level {level + 1}/{self.config.max_refinement_levels} ===")
                    grid_info = self.current_grid.get_grid_info()
                    print(f"Grid dimensions: {grid_info['dimensions']}")
                    print(f"Total grid points: {grid_info['total_points']:,}")
                
                # Callback for level start
                if self.config.level_start_callback:
                    self.config.level_start_callback(level, self.current_grid)
                
                # Solve on current refinement level
                level_start_time = time.time()
                
                converged, final_solution, level_results = self._solve_on_level(
                    level, solver_function, current_solution
                )
                
                level_time = time.time() - level_start_time
                
                # Store results
                results['convergence_achieved'].append(converged)
                results['final_residuals'].append(level_results['final_residual'])
                results['iteration_counts'].append(level_results['iterations'])
                results['refinement_times'].append(level_time)
                results['grids_used'].append(deepcopy(self.current_grid.get_grid_info()))
                
                self.state.converged_levels[level] = converged
                self.solutions.append(final_solution)
                
                if self.config.verbose:
                    status = "CONVERGED" if converged else "MAX ITERATIONS"
                    print(f"Level {level + 1} completed: {status}")
                    print(f"Final residual: {level_results['final_residual']:.2e}")
                    print(f"Iterations: {level_results['iterations']}")
                    print(f"Time: {level_time:.2f} seconds")
                
                # Prepare for next level
                if level < self.config.max_refinement_levels - 1:
                    # Refine grid for next level
                    success = self._prepare_next_level(final_solution)
                    if not success:
                        if self.config.verbose:
                            print("Grid refinement failed - stopping refinement sequence")
                        break
                        
                    # Use current solution as initial guess for next level
                    current_solution = self._interpolate_solution_to_refined_grid(
                        final_solution, self.state.grid_levels[level], self.current_grid
                    )
                else:
                    # Final level completed
                    results['final_solution'] = final_solution
                    
            # Overall completion
            self.state.total_time = time.time() - start_time
            results['total_time'] = self.state.total_time
            
            # Determine overall success
            any_converged = any(self.state.converged_levels)
            all_converged = all(self.state.converged_levels)
            
            if all_converged:
                self.state.status = RefinementStatus.CONVERGED
                success = True
            elif any_converged:
                self.state.status = RefinementStatus.CONVERGED  # Partial success
                success = True
            else:
                self.state.status = RefinementStatus.MAX_ITERATIONS
                success = False
                
            if self.config.verbose:
                print(f"\n=== Refinement Sequence Complete ===")
                print(f"Total time: {self.state.total_time:.2f} seconds")
                print(f"Levels converged: {sum(self.state.converged_levels)}/{self.config.max_refinement_levels}")
                print(f"Final status: {self.state.status.value}")
                
            return success, results
            
        except Exception as e:
            self.state.status = RefinementStatus.FAILED
            if self.config.verbose:
                print(f"Refinement sequence failed: {e}")
            return False, {'error': str(e)}
            
    def _solve_on_level(self, 
                       level: int, 
                       solver_function: Callable,
                       initial_guess: Optional[np.ndarray]) -> Tuple[bool, np.ndarray, Dict]:
        """
        Solve Poisson equation on current refinement level
        
        Returns:
            (converged, solution, results): Convergence flag, final solution, and results dict
        """
        tolerance = self.config.convergence_tolerance[level]
        max_iterations = self.config.max_iterations_per_level[level]
        relaxation_factor = self.config.relaxation_factors[level]
        
        self.state.current_iteration = 0
        converged = False
        previous_solution = None
        
        # Initialize solution
        if initial_guess is not None:
            current_solution = initial_guess.copy()
        else:
            # Create zero initial guess based on grid size
            grid_info = self.current_grid.get_grid_info()
            # This would need to be adapted based on actual potential array structure
            current_solution = np.zeros((grid_info['dimensions']['nr_points'],
                                       grid_info['dimensions']['ns_points']), dtype=np.float64)
        
        residual_history = []
        iteration_times = []
        
        # Iterative solution loop
        for iteration in range(max_iterations):
            iter_start_time = time.time()
            self.state.current_iteration = iteration
            
            # Store previous solution for convergence check
            if iteration > 0:
                previous_solution = current_solution.copy()
            
            # Call solver function
            try:
                new_solution = solver_function(self.current_grid, current_solution)
                
                # Apply relaxation
                if previous_solution is not None:
                    current_solution = (relaxation_factor * new_solution + 
                                      (1.0 - relaxation_factor) * previous_solution)
                else:
                    current_solution = new_solution
                    
            except Exception as e:
                if self.config.verbose:
                    print(f"Solver failed at iteration {iteration}: {e}")
                break
                
            iter_time = time.time() - iter_start_time
            iteration_times.append(iter_time)
            
            # Check convergence
            if iteration % self.config.convergence_check_interval == 0 and previous_solution is not None:
                residual = self._calculate_residual(current_solution, previous_solution)
                residual_history.append(residual)
                self.state.last_residual = residual
                
                if self.config.verbose and iteration % (self.config.convergence_check_interval * 2) == 0:
                    print(f"  Iteration {iteration:3d}: residual = {residual:.2e}")
                
                # Convergence check
                if residual < tolerance:
                    converged = True
                    if self.config.verbose:
                        print(f"  Converged at iteration {iteration}")
                    break
                    
                # Callback for iteration monitoring
                if self.config.iteration_callback:
                    self.config.iteration_callback(level, iteration, residual)
                    
        # Store convergence history
        self.state.convergence_history[level] = residual_history
        self.state.iteration_times[level] = iteration_times
        
        # Final results for this level
        results = {
            'converged': converged,
            'iterations': self.state.current_iteration + 1,
            'final_residual': self.state.last_residual,
            'residual_history': residual_history,
            'avg_iteration_time': np.mean(iteration_times) if iteration_times else 0.0
        }
        
        # Convergence callback
        if self.config.convergence_callback:
            self.config.convergence_callback(level, converged, results)
            
        return converged, current_solution, results
        
    def _prepare_next_level(self, current_solution: np.ndarray) -> bool:
        """
        Prepare grid for next refinement level
        
        Args:
            current_solution: Solution from current level
            
        Returns:
            Success flag
        """
        try:
            # Create refined grid
            refined_grid = deepcopy(self.current_grid)
            success = refined_grid.refine_grid()
            
            if success:
                self.state.grid_levels.append(refined_grid)
                self.current_grid = refined_grid
                return True
            else:
                return False
                
        except Exception as e:
            if self.config.verbose:
                print(f"Grid refinement failed: {e}")
            return False
            
    def _interpolate_solution_to_refined_grid(self, 
                                            coarse_solution: np.ndarray,
                                            coarse_grid: Grid3D, 
                                            fine_grid: Grid3D) -> np.ndarray:
        """
        Interpolate solution from coarse grid to refined grid
        
        This is a simplified interpolation - in practice would need to handle
        the full 3D potential arrays properly.
        
        Args:
            coarse_solution: Solution on coarse grid
            coarse_grid: Coarse grid object
            fine_grid: Refined grid object
            
        Returns:
            Interpolated solution on fine grid
        """
        # Simple linear interpolation for demonstration
        # In actual implementation, would need proper 3D interpolation
        
        fine_info = fine_grid.get_grid_info()
        fine_shape = (fine_info['dimensions']['nr_points'],
                     fine_info['dimensions']['ns_points'])
        
        if coarse_solution.shape == fine_shape:
            # Same shape - just copy
            return coarse_solution.copy()
        else:
            # Different shapes - need interpolation
            from scipy.interpolate import RegularGridInterpolator
            
            try:
                # Create interpolator (simplified 2D case)
                coarse_r = coarse_grid.r_points
                coarse_z = coarse_grid.z_semiconductor_points
                
                fine_r = fine_grid.r_points
                fine_z = fine_grid.z_semiconductor_points
                
                # Create meshgrids
                coarse_R, coarse_Z = np.meshgrid(coarse_r, coarse_z, indexing='ij')
                fine_R, fine_Z = np.meshgrid(fine_r, fine_z, indexing='ij')
                
                # Interpolate
                interpolator = RegularGridInterpolator(
                    (coarse_r, coarse_z), coarse_solution, 
                    method='linear', bounds_error=False, fill_value=0.0
                )
                
                points = np.column_stack([fine_R.ravel(), fine_Z.ravel()])
                fine_solution = interpolator(points).reshape(fine_shape)
                
                return fine_solution
                
            except Exception as e:
                if self.config.verbose:
                    print(f"Interpolation failed, using zero initial guess: {e}")
                return np.zeros(fine_shape, dtype=np.float64)
                
    def _calculate_residual(self, new_solution: np.ndarray, old_solution: np.ndarray) -> float:
        """
        Calculate convergence residual between solutions
        
        Args:
            new_solution: Current iteration solution
            old_solution: Previous iteration solution
            
        Returns:
            Relative residual
        """
        if new_solution.shape != old_solution.shape:
            return float('inf')
            
        # Calculate relative change
        diff = np.abs(new_solution - old_solution)
        norm_new = np.abs(new_solution)
        
        # Avoid division by zero
        denominator = np.maximum(norm_new, 1e-12)
        relative_change = diff / denominator
        
        # Return maximum relative change
        return np.max(relative_change)
        
    def get_refinement_summary(self) -> Dict:
        """Get comprehensive summary of refinement process"""
        return {
            'config': {
                'max_levels': self.config.max_refinement_levels,
                'tolerances': self.config.convergence_tolerance,
                'max_iterations': self.config.max_iterations_per_level
            },
            'state': {
                'current_level': self.state.current_level,
                'current_iteration': self.state.current_iteration,
                'status': self.state.status.value,
                'total_time': self.state.total_time
            },
            'convergence': {
                'levels_converged': self.state.converged_levels,
                'final_residuals': [history[-1] if history else float('inf') 
                                  for history in self.state.convergence_history],
                'total_iterations': [len(times) for times in self.state.iteration_times]
            },
            'grids': [grid.get_grid_info() for grid in self.state.grid_levels]
        }
        
    def plot_convergence_history(self, save_path: Optional[str] = None):
        """
        Plot convergence history for all refinement levels
        
        Args:
            save_path: Optional path to save plot
        """
        try:
            import matplotlib.pyplot as plt
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            # Plot residual convergence
            for level, history in enumerate(self.state.convergence_history):
                if history:
                    iterations = np.arange(0, len(history)) * self.config.convergence_check_interval
                    ax1.semilogy(iterations, history, 'o-', label=f'Level {level + 1}')
                    
            ax1.set_xlabel('Iteration')
            ax1.set_ylabel('Residual')
            ax1.set_title('Convergence History by Refinement Level')
            ax1.legend()
            ax1.grid(True)
            
            # Plot iteration times
            for level, times in enumerate(self.state.iteration_times):
                if times:
                    iterations = np.arange(len(times))
                    ax2.plot(iterations, times, 'o-', label=f'Level {level + 1}')
                    
            ax2.set_xlabel('Iteration')
            ax2.set_ylabel('Time per Iteration (s)')
            ax2.set_title('Iteration Times by Refinement Level')
            ax2.legend()
            ax2.grid(True)
            
            plt.tight_layout()
            
            if save_path:
                plt.savefig(save_path, dpi=150, bbox_inches='tight')
                print(f"Convergence plot saved to: {save_path}")
                
            plt.show()
            
        except ImportError:
            print("Matplotlib not available - cannot plot convergence history")
        except Exception as e:
            print(f"Plotting failed: {e}")
            
    def __str__(self) -> str:
        """String representation"""
        return (f"AdaptiveGridRefinement(levels={self.config.max_refinement_levels}, "
               f"current_level={self.state.current_level}, "
               f"status={self.state.status.value})")
               
    def __repr__(self) -> str:
        return self.__str__()


# Factory functions

def create_standard_refinement(initial_grid: Grid3D, 
                             geometry: STMGeometry,
                             verbose: bool = True) -> AdaptiveGridRefinement:
    """Create standard three-level refinement configuration"""
    config = RefinementConfig(
        max_refinement_levels=3,
        convergence_tolerance=[1e-4, 1e-5, 1e-6],
        max_iterations_per_level=[50, 100, 200],
        verbose=verbose
    )
    return AdaptiveGridRefinement(initial_grid, geometry, config)


def create_fast_refinement(initial_grid: Grid3D, 
                         geometry: STMGeometry) -> AdaptiveGridRefinement:
    """Create fast two-level refinement for testing"""
    config = RefinementConfig(
        max_refinement_levels=2,
        convergence_tolerance=[1e-3, 1e-4],
        max_iterations_per_level=[30, 50],
        relaxation_factors=[0.7, 0.8],
        verbose=False
    )
    return AdaptiveGridRefinement(initial_grid, geometry, config)


def create_high_precision_refinement(initial_grid: Grid3D, 
                                   geometry: STMGeometry) -> AdaptiveGridRefinement:
    """Create high-precision refinement for research applications"""
    config = RefinementConfig(
        max_refinement_levels=3,
        convergence_tolerance=[1e-5, 1e-6, 1e-7],
        max_iterations_per_level=[100, 200, 400],
        relaxation_factors=[0.6, 0.7, 0.8],
        verbose=True
    )
    return AdaptiveGridRefinement(initial_grid, geometry, config)


if __name__ == "__main__":
    # Demo usage
    print("Adaptive Grid Refinement Demo")
    print("=" * 40)
    
    # This would normally be imported from other modules
    from .grid3d import create_coarse_grid
    from .stm_geometry import create_standard_stm_geometry
    
    # Create test setup
    grid = create_coarse_grid()
    geometry = create_standard_stm_geometry()
    
    # Create refinement manager
    refinement = create_standard_refinement(grid, geometry)
    print(f"Created refinement: {refinement}")
    
    # Get summary
    summary = refinement.get_refinement_summary()
    print(f"\nRefinement Summary:")
    for section, data in summary.items():
        print(f"  {section}: {data}")
        
    # Demo solver function (placeholder)
    def dummy_solver(grid, initial_guess):
        """Dummy solver that just returns slightly modified input"""
        if initial_guess is not None:
            return initial_guess * 0.95 + np.random.normal(0, 0.01, initial_guess.shape)
        else:
            info = grid.get_grid_info()
            shape = (info['dimensions']['nr_points'], info['dimensions']['ns_points'])
            return np.random.normal(0, 0.1, shape)
    
    print(f"\nNote: In actual usage, refinement.refine_grid_sequence(actual_solver)")
    print(f"would be called with a real Poisson solver function.")