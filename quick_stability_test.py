#!/usr/bin/env python3
import sys
import time
import numpy as np
from src.physics.core.poisson import PoissonSOREquation
from src.physics.solvers.grid import HyperbolicGrid

def test_stability():
    """Quick test of numerical stability"""
    print("Testing Poisson solver stability...")
    
    try:
        # Create a small grid for testing
        grid = HyperbolicGrid(N_eta=8, N_nu=4, R=1.0, Z_TS=1.0, shank_slope=1.0)
        
        # Create minimal props object
        class MinimalProps:
            def __init__(self):
                class SemiconductorProps:
                    epsilon_r = 12.9
                    Ev_offset_eV = 0.0
                self.semiconductor_props = SemiconductorProps()
        
        props = MinimalProps()
        
        # Create solver
        solver = PoissonSOREquation(grid, props)
        
        print(f"Grid size: {grid.N_eta} x {grid.N_nu}")
        print(f"SOR coefficients range: A_P max={np.max(solver.A_P):.3e}, min={np.min(solver.A_P):.3e}")
        
        # Test Laplace equation with reasonable parameters
        V_tip = -2.0
        V_sample = 0.0
        
        start_time = time.time()
        potential, iterations, max_error = solver.solve_laplace(
            V_tip, V_sample, max_iterations=100, tolerance=1e-4
        )
        end_time = time.time()
        
        print(f"Laplace test completed in {end_time - start_time:.2f} seconds")
        print(f"Iterations: {iterations}, Max error: {max_error:.3e}")
        print(f"Potential range: [{np.min(potential):.3f}, {np.max(potential):.3f}] V")
        
        # Check for numerical stability
        if np.any(np.isnan(potential)) or np.any(np.isinf(potential)):
            print("❌ FAILED: NaN or Inf values detected")
            return False
        
        if np.max(np.abs(potential)) > 100:
            print("❌ FAILED: Potential values too large (>100V)")
            return False
            
        pot0 = solver._calculate_pot0_fortran_style(potential)
        print(f"Pot0 (band bending): {pot0:.6f} V")
        
        if abs(pot0) < 1e-10:
            print("⚠️  WARNING: Pot0 is essentially zero")
        elif abs(pot0) > 10:
            print("⚠️  WARNING: Pot0 magnitude is very large")
        else:
            print("✅ PASSED: Numerical stability test")
            
        return True
        
    except Exception as e:
        print(f"❌ FAILED: Exception occurred: {e}")
        return False

if __name__ == "__main__":
    success = test_stability()
    sys.exit(0 if success else 1)