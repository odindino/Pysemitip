#!/usr/bin/env python3
"""
Debug script to analyze the Pot0 calculation differences between Python and Fortran
"""

import numpy as np
import sys
import os
import logging

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from simulation.multint import MultIntSimulation
from core.filereader import YamlConfigReader

def debug_pot0_calculation():
    """Debug the Pot0 calculation to understand coordinate system and array access"""
    
    # Load the quick test config
    config_path = "data/input/examples/test/quick_test.yaml"
    reader = YamlConfigReader()
    config = reader.load_config(config_path)
    
    # Create simulation
    output_dir = "data/output/debug_pot0"
    os.makedirs(output_dir, exist_ok=True)
    sim = MultIntSimulation(config, output_dir)
    
    # Get grid information
    grid = sim.poisson_solver.grid
    print(f"Grid dimensions: N_eta={grid.N_eta}, N_nu={grid.N_nu}")
    print(f"Grid physical parameters: R={grid.R:.3f} nm, Z_TS={grid.Z_TS:.3f} nm")
    print(f"Grid coordinate system: f={grid.f:.6f}, eta_tip={grid.eta_tip:.6f}")
    
    # Show coordinate arrays near the interface
    print("\n=== Grid Coordinates Analysis ===")
    print("Sample surface points (nu = N_nu-1 = last nu index):")
    nu_sample_idx = grid.N_nu - 1
    for i in range(min(5, grid.N_eta)):
        r_val = grid.r[i, nu_sample_idx]
        z_val = grid.z[i, nu_sample_idx]
        print(f"  eta_idx={i}: r={r_val:.6f} nm, z={z_val:.6f} nm")
    
    print("\nTip surface points (eta_idx = 0):")
    for j in range(min(5, grid.N_nu)):
        r_val = grid.r[0, j]
        z_val = grid.z[0, j]
        nu_val = grid.nu[0, j]
        print(f"  nu_idx={j}, nu={nu_val:.6f}: r={r_val:.6f} nm, z={z_val:.6f} nm")
    
    # Initialize potential with some test values to understand indexing
    print("\n=== Testing Potential Array Structure ===")
    potential = np.zeros((grid.N_eta, grid.N_nu))
    
    # Set tip potential
    V_tip = -1.0  # Test bias
    potential[0, :] = V_tip
    print(f"Tip potential (eta_idx=0): {potential[0, :5]}")
    
    # Set sample potential  
    V_sample = 0.0
    potential[-1, :] = V_sample  # Last eta index
    print(f"Sample potential (eta_idx={grid.N_eta-1}): {potential[-1, :5]}")
    
    # Test our PCENT calculation
    print("\n=== PCENT Calculation Analysis ===")
    
    # Set up a simple linear potential for testing
    for i in range(grid.N_eta):
        for j in range(grid.N_nu):
            # Linear interpolation between tip and sample
            potential[i, j] = V_tip + (V_sample - V_tip) * (i / (grid.N_eta - 1))
    
    # Our current Python PCENT implementation
    I = 0  # First eta index (tip surface)
    if I + 1 < grid.N_eta:
        sum_val = 0.0
        for k in range(grid.N_nu):
            weighted_value = (9.0 * potential[I, k] - potential[I + 1, k]) / 8.0
            sum_val += weighted_value
        pot0_python = sum_val / grid.N_nu
        print(f"Python PCENT calculation: {pot0_python:.8f} V")
        print(f"  Using I={I} (tip surface), I+1={I+1}")
        print(f"  potential[{I}, :] = {potential[I, :5]}")
        print(f"  potential[{I+1}, :] = {potential[I+1, :5]}")
    
    # Test what happens if we use interface indices differently
    print("\n=== Testing Different Interface Interpretations ===")
    
    # Maybe the interface is not at eta_idx=0 but somewhere else?
    # Let's check if there's a specific interface region
    
    # Try interface at sample surface (last nu index)
    print("Testing interface at sample surface (nu = N_nu-1):")
    interface_nu = grid.N_nu - 1
    sum_val_alt = 0.0
    for i in range(grid.N_eta - 1):  # Avoid boundary
        weighted_value = (9.0 * potential[i, interface_nu] - potential[i + 1, interface_nu]) / 8.0
        sum_val_alt += weighted_value
    pot0_alt = sum_val_alt / (grid.N_eta - 1)
    print(f"Alternative PCENT (at sample surface): {pot0_alt:.8f} V")
    
    # Check Fortran coordinate mapping
    print("\n=== Fortran Coordinate System Analysis ===")
    print("Fortran PCENT uses VSINT(JJ, I, K) with JJ=0 for interface")
    print("In Fortran SEMITIP:")
    print("  - JJ=0 is the interface (sample surface)")
    print("  - I runs from 1 to NR (radial direction)")
    print("  - K runs from 1 to NP (angular direction)")
    print("  - I=1 in Fortran corresponds to I=0 in Python")
    
    print(f"\nIn our Python grid:")
    print(f"  - eta direction (I): 0 to {grid.N_eta-1} (0=tip, {grid.N_eta-1}=far field)")
    print(f"  - nu direction (J): 0 to {grid.N_nu-1} (0=axis, {grid.N_nu-1}=sample surface)")
    
    print("\nProblem identified: We may be using wrong array indices!")
    print("Fortran VSINT(JJ=0, I, K) likely refers to potential at sample surface")
    print("Our current implementation uses potential[I=0, K] which is tip surface")
    
    # Test the correct interpretation: sample surface interface
    print("\n=== Corrected PCENT Calculation ===")
    # Interface should be at sample surface: potential[:, -1]
    interface_potential = potential[:, -1]  # Sample surface
    
    if grid.N_eta >= 2:
        sum_val_corrected = 0.0
        # Use first two points in eta direction at sample surface
        I_interface = 0  # Still use I=0 as per Fortran, but at sample surface
        weighted_value = (9.0 * interface_potential[I_interface] - interface_potential[I_interface + 1]) / 8.0
        # For PCENT, we need to average over angular points, but since we're at sample surface,
        # we have only one "angular" point. Let's check what happens with more points.
        
        # Actually, let's reinterpret: maybe K loop is over eta points at interface
        sum_val_corrected = 0.0
        valid_points = 0
        for i in range(min(grid.N_eta - 1, 10)):  # Limit to avoid boundary issues
            weighted_value = (9.0 * potential[i, -1] - potential[i + 1, -1]) / 8.0
            sum_val_corrected += weighted_value
            valid_points += 1
        
        if valid_points > 0:
            pot0_corrected = sum_val_corrected / valid_points
            print(f"Corrected PCENT (sample surface interface): {pot0_corrected:.8f} V")
            print(f"  Used {valid_points} points along eta at sample surface")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    debug_pot0_calculation()