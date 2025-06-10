#!/usr/bin/env python3
"""
Test the interface potential fix with a simple iteration
"""

import numpy as np
import sys
import os
import logging

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from simulation.multint import MultIntSimulation
from core.filereader import YamlConfigReader

def test_interface_iteration():
    """Test multiple SOR iterations to see potential development"""
    
    # Load the quick test config
    config_path = "data/input/examples/test/quick_test.yaml"
    reader = YamlConfigReader()
    config = reader.load_config(config_path)
    
    # Create simulation
    output_dir = "data/output/test_interface_fix"
    os.makedirs(output_dir, exist_ok=True)
    sim = MultIntSimulation(config, output_dir)
    
    grid = sim.poisson_solver.grid
    print(f"Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    
    # Test case: -1V bias
    V_tip = -1.0
    V_sample = 0.0
    
    # Initialize potential
    potential = np.zeros((grid.N_eta, grid.N_nu))
    
    print(f"\n=== Testing Interface Development with Multiple Iterations ===")
    print(f"Bias: {V_tip:.1f} V")
    
    # Test multiple iterations to see potential development
    for iter_num in range(10):
        # Apply boundary conditions (tip only, not sample surface)
        potential = sim.poisson_solver._apply_boundary_conditions(potential, V_tip, V_sample)
        
        # Update interface potential
        potential = sim.poisson_solver._update_interface_potential(potential)
        
        # Simple SOR iteration for interior points
        if iter_num < 5:  # Only do a few SOR steps
            potential_old = np.copy(potential)
            omega = 1.0
            for i in range(1, grid.N_eta - 1):
                for j in range(1, grid.N_nu - 1):
                    # Simple Laplace equation update
                    laplacian = (potential_old[i+1, j] + potential_old[i-1, j] + 
                               potential_old[i, j+1] + potential_old[i, j-1] - 4*potential_old[i, j])
                    potential[i, j] = potential_old[i, j] + omega * 0.25 * laplacian
        
        # Calculate Pot0
        pot0 = sim.poisson_solver._calculate_pot0_fortran_style(potential)
        
        print(f"Iteration {iter_num}: Pot0 = {pot0:.8f} V")
        if iter_num < 3:
            print(f"  Sample surface: {potential[:5, -1]}")
    
    print(f"\nFinal potential at sample surface (nu={grid.N_nu-1}):")
    for i in range(min(8, grid.N_eta)):
        print(f"  eta[{i}]: {potential[i, -1]:.6f} V")

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)  # Reduce log noise
    test_interface_iteration()