#!/usr/bin/env python3
"""
Detailed debug of the Pot0 calculation to understand what's happening
"""

import numpy as np
import sys
import os
import logging

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from simulation.multint import MultIntSimulation
from core.filereader import YamlConfigReader

def debug_pot0_detailed():
    """Detailed debug of the Pot0 calculation"""
    
    # Load the quick test config
    config_path = "data/input/examples/test/quick_test.yaml"
    reader = YamlConfigReader()
    config = reader.load_config(config_path)
    
    # Create simulation
    output_dir = "data/output/debug_pot0_detailed"
    os.makedirs(output_dir, exist_ok=True)
    sim = MultIntSimulation(config, output_dir)
    
    grid = sim.poisson_solver.grid
    print(f"Grid: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    
    # Create a test potential with realistic values
    potential = np.zeros((grid.N_eta, grid.N_nu))
    
    # Test case: -1V bias
    V_tip = -1.0
    V_sample = 0.0
    
    # Apply boundary conditions
    potential = sim.poisson_solver._apply_boundary_conditions(potential, V_tip, V_sample)
    
    print("\n=== Boundary Conditions Applied ===")
    print(f"Tip surface (eta=0): {potential[0, :]}")
    print(f"Sample surface (nu={grid.N_nu-1}): {potential[:, -1]}")
    print(f"Axis (nu=0): {potential[:5, 0]}")
    print(f"Far field (eta={grid.N_eta-1}): {potential[-1, :5]}")
    
    # Create a realistic potential distribution by linear interpolation
    print("\n=== Creating Linear Potential Distribution ===")
    for i in range(grid.N_eta):
        for j in range(grid.N_nu):
            # Linear interpolation from tip to sample
            if j == grid.N_nu - 1:  # Sample surface
                potential[i, j] = V_sample
            elif j == 0:  # Axis - interpolate between tip center and far field
                eta_frac = i / (grid.N_eta - 1)
                potential[i, j] = V_tip * (1 - eta_frac) + V_sample * eta_frac
            else:  # Interior points
                nu_frac = j / (grid.N_nu - 1)
                eta_frac = i / (grid.N_eta - 1)
                # Interpolate in both directions
                potential[i, j] = V_tip * (1 - nu_frac) * (1 - eta_frac) + V_sample * nu_frac
    
    # Reapply boundary conditions to ensure consistency
    potential = sim.poisson_solver._apply_boundary_conditions(potential, V_tip, V_sample)
    
    print("After linear interpolation:")
    print(f"Sample surface values: {potential[:8, -1]}")
    print(f"Axis values: {potential[:8, 0]}")
    
    # Now test the PCENT calculation in detail
    print("\n=== Detailed PCENT Calculation ===")
    interface_nu_idx = grid.N_nu - 1
    
    print(f"Interface nu index: {interface_nu_idx}")
    print(f"Potential at interface (sample surface):")
    for i in range(min(8, grid.N_eta)):
        print(f"  eta[{i}]: {potential[i, interface_nu_idx]:.6f} V")
    
    # Manual PCENT calculation
    sum_val = 0.0
    valid_points = 0
    
    print("\nManual PCENT calculation:")
    for I in range(min(grid.N_eta - 1, 8)):
        if I + 1 < grid.N_eta:
            v1 = potential[I, interface_nu_idx]
            v2 = potential[I + 1, interface_nu_idx]
            weighted_value = (9.0 * v1 - v2) / 8.0
            sum_val += weighted_value
            valid_points += 1
            print(f"  I={I}: V1={v1:.6f}, V2={v2:.6f}, weighted={(9.0*v1-v2)/8.0:.6f}")
    
    pot0_manual = sum_val / valid_points if valid_points > 0 else 0.0
    print(f"Manual Pot0: {pot0_manual:.8f} V (sum={sum_val:.6f}, points={valid_points})")
    
    # Compare with function
    pot0_function = sim.poisson_solver._calculate_pot0_fortran_style(potential)
    print(f"Function Pot0: {pot0_function:.8f} V")
    
    # Test with different potential profiles
    print("\n=== Testing Different Potential Profiles ===")
    
    # Profile 1: Strong gradient from tip to sample
    potential1 = np.copy(potential)
    for i in range(grid.N_eta):
        for j in range(grid.N_nu):
            if j == grid.N_nu - 1:  # Sample surface - create variation
                potential1[i, j] = V_sample + 0.1 * i / grid.N_eta  # Small gradient
            elif j == 0:  # Axis
                potential1[i, j] = V_tip + 0.5 * i / grid.N_eta  # Larger gradient
    
    potential1 = sim.poisson_solver._apply_boundary_conditions(potential1, V_tip, V_sample)
    pot0_test1 = sim.poisson_solver._calculate_pot0_fortran_style(potential1)
    print(f"Test profile 1 - Sample surface with gradient: {pot0_test1:.8f} V")
    print(f"  Sample surface values: {potential1[:5, -1]}")
    
    # Profile 2: Different approach
    potential2 = np.ones((grid.N_eta, grid.N_nu)) * V_sample
    potential2[0, :] = V_tip  # Tip surface
    potential2[:, -1] = np.linspace(V_sample, V_sample + 0.2, grid.N_eta)  # Sample surface gradient
    
    pot0_test2 = sim.poisson_solver._calculate_pot0_fortran_style(potential2)
    print(f"Test profile 2 - Linear sample surface gradient: {pot0_test2:.8f} V")
    print(f"  Sample surface values: {potential2[:5, -1]}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)  # Reduce log noise
    debug_pot0_detailed()