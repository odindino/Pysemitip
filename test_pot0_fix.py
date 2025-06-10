#!/usr/bin/env python3
"""
Test the corrected Pot0 calculation with a simple Poisson solve
"""

import numpy as np
import sys
import os
import logging

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from simulation.multint import MultIntSimulation
from core.filereader import YamlConfigReader

def test_corrected_pot0():
    """Test the corrected Pot0 calculation with actual Poisson solving"""
    
    # Load the quick test config
    config_path = "data/input/examples/test/quick_test.yaml"
    reader = YamlConfigReader()
    config = reader.load_config(config_path)
    
    # Create simulation
    output_dir = "data/output/test_pot0_fix"
    os.makedirs(output_dir, exist_ok=True)
    sim = MultIntSimulation(config, output_dir)
    
    print("=== Testing Corrected Pot0 Calculation ===")
    
    # Test different bias voltages
    bias_voltages = [-1.0, 0.0, 1.0]
    
    for bias in bias_voltages:
        print(f"\n--- Testing bias = {bias:.1f} V ---")
        
        # Update tip potential
        contact_potential = 0.0  # From config
        tip_potential = bias + contact_potential
        print(f"Tip potential: {tip_potential:.6f} V")
        
        # Initialize with simple potential distribution for testing
        grid = sim.poisson_solver.grid
        potential = np.zeros((grid.N_eta, grid.N_nu))
        
        # Apply boundary conditions
        potential = sim.poisson_solver._apply_boundary_conditions(potential, tip_potential, 0.0)
        
        print(f"After boundary conditions:")
        print(f"  Tip surface potential[0, :]: {potential[0, :5]}")
        print(f"  Sample surface potential[:, -1]: {potential[:5, -1]}")
        
        # Calculate Pot0 with corrected method
        pot0_corrected = sim.poisson_solver._calculate_pot0_fortran_style(potential)
        print(f"Corrected Pot0: {pot0_corrected:.8f} V")
        
        # Test with a simple SOR iteration to see behavior
        print("\nTesting with one SOR iteration...")
        N_eta, N_nu = grid.N_eta, grid.N_nu
        omega = 1.0
        
        # Simple SOR update (no charge density)
        potential_old = np.copy(potential)
        for i in range(1, N_eta - 1):
            for j in range(1, N_nu - 1):
                # Simple Laplace equation update
                laplacian = (potential_old[i+1, j] + potential_old[i-1, j] + 
                           potential_old[i, j+1] + potential_old[i, j-1] - 4*potential_old[i, j])
                potential[i, j] = potential_old[i, j] + omega * 0.25 * laplacian
        
        # Reapply boundary conditions
        potential = sim.poisson_solver._apply_boundary_conditions(potential, tip_potential, 0.0)
        
        pot0_after_sor = sim.poisson_solver._calculate_pot0_fortran_style(potential)
        print(f"Pot0 after 1 SOR iteration: {pot0_after_sor:.8f} V")
        print(f"Change in Pot0: {pot0_after_sor - pot0_corrected:.8f} V")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_corrected_pot0()