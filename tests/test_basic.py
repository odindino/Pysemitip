#!/usr/bin/env python3
"""
Basic test to verify the installation and imports are working.
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

print("Testing basic imports...")

try:
    # Test utils
    from src.utils.constants import PhysicalConstants as PC
    print("[OK] Utils imported successfully")
    print(f"  Elementary charge: {PC.E:.6e} C")
    
    # Test materials
    from src.physics.materials import SemiconductorRegion, TipModel
    print("[OK] Materials imported successfully")
    
    # Test grid
    from src.physics.solvers import Grid3D, GridParameters
    print("[OK] Solvers imported successfully")
    
    # Test core physics
    from src.physics.core import PoissonSolver, SchrodingerSolver
    print("[OK] Core physics imported successfully")
    
    # Test simulation
    from src.simulation.multint import MultIntSimulation
    print("[OK] Simulation module imported successfully")
    
    # Test visualization
    from src.physics.visualization import SEMITIPPlotter
    print("[OK] Visualization imported successfully")
    
    print("\nAll imports successful! The package is ready to use.")
    
    # Quick functionality test
    print("\nTesting basic functionality...")
    
    # Create a semiconductor
    semi = SemiconductorRegion(
        region_id=1,
        donor_concentration=1e18,
        acceptor_concentration=0,
        band_gap=1.42,
        valence_band_offset=0,
        electron_affinity=4.07,
        donor_binding_energy=0.006,
        acceptor_binding_energy=0.031,
        cb_effective_mass=0.067,
        vb_effective_mass_heavy=0.5,
        vb_effective_mass_light=0.08,
        vb_effective_mass_so=0.15,
        spin_orbit_splitting=0.34,
        permittivity=12.9,
        temperature=300.0
    )
    
    print(f"[OK] Created semiconductor: n-type = {semi.is_n_type}")
    print(f"  Fermi level ~ {semi.fermi_level():.3f} eV")
    
    print("\nBasic tests passed!")
    
except ImportError as e:
    print(f"[FAIL] Import error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"[FAIL] Error: {e}")
    sys.exit(1)