#!/usr/bin/env python3
"""
Test the fixed grid generation to compare with Fortran values
"""

import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from src.core.filereader import load_yaml_config
from src.simulation.multint import MultIntSimulation

def test_grid_generation():
    """Test the new grid generation logic"""
    
    print("="*60)
    print("測試修正後的網格生成")
    print("="*60)
    
    # Load configuration
    config_path = "data/input/examples/quick_test.yaml"
    config = load_yaml_config(config_path)
    
    # Create simulation
    sim = MultIntSimulation(config)
    
    # Check new grid parameters
    grid = sim.grid
    print(f"\n修正後的網格參數:")
    print(f"   Grid dimensions: nr={grid.params.nr}, nv={grid.params.nv}, ns={grid.params.ns}, np={grid.params.np}")
    print(f"   Input spacing parameters:")
    print(f"     delr (DELR0): {grid.params.delr:.6f}")
    print(f"     dels (DELS0): {grid.params.dels:.6f}")
    print(f"     delv: {grid.params.delv:.6f}")
    print(f"     delp: {grid.params.delp:.6f}")
    
    # Check actual coordinate spacing
    print(f"\n實際坐標分佈:")
    print(f"   R values (first 5): {grid.r[:5]}")
    print(f"   R spacing (first 4): {np.diff(grid.r[:5])}")
    print(f"   ZS values (first 5): {grid.zs_positive[:5]}")
    print(f"   ZS spacing (first 4): {np.diff(grid.zs_positive[:5])}")
    
    # Compare with expected Fortran values
    print(f"\n與Fortran比較:")
    print(f"   delr: {grid.params.delr:.6f} (target: ~0.5)")
    print(f"   dels: {grid.params.dels:.6f} (target: ~0.5)")
    print(f"   delv: {grid.params.delv:.6f} (target: ~0.25)")
    print(f"   delp: {grid.params.delp:.6f} (target: ~0.393)")
    
    # Check if grid coordinates make sense
    print(f"\n網格坐標範圍:")
    print(f"   R range: 0 to {grid.r[-1]:.3f} nm")
    print(f"   ZS range: 0 to {grid.zs_positive[-1]:.3f} nm")
    print(f"   ZV range: 0 to {grid.zv[-1]:.3f} nm")
    
    # Test depletion width calculation with new grid
    test_bias = -2.0
    depl_width = sim._estimate_depletion_width_1d(test_bias)
    print(f"\nDepletion width測試:")
    print(f"   Python (修正後): {depl_width:.6f} nm")
    print(f"   Fortran目標: 54.337425 nm")
    print(f"   比例: {depl_width / 54.337425:.3f}")

if __name__ == "__main__":
    test_grid_generation()