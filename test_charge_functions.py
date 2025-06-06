#!/usr/bin/env python3
"""
Simple test to verify charge density functions are working and return non-zero values.
"""

import sys
import os
import numpy as np

# Add both src and current directory to path 
src_path = os.path.join(os.path.dirname(__file__), 'src')
sys.path.insert(0, src_path)
sys.path.insert(0, os.path.dirname(__file__))

try:
    from src.physics.materials.semiconductor import SemiconductorRegion
    from src.physics.core.charge_density import ChargeDensityCalculator
except ImportError:
    # Try alternative import paths
    import src.physics.materials.semiconductor as semiconductor_module
    import src.physics.core.charge_density as charge_module
    SemiconductorRegion = semiconductor_module.SemiconductorRegion
    ChargeDensityCalculator = charge_module.ChargeDensityCalculator

def test_semiconductor_creation():
    """Test that we can create a SemiconductorRegion with the exact Fortran parameters."""
    print("Creating SemiconductorRegion with exact Fortran parameters...")
    
    # Use exact parameters from Fortran simulation
    region = SemiconductorRegion(
        region_id=1,
        donor_concentration=1e18,  # 1e24 m^-3 = 1e18 cm^-3 in Fortran
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
    
    print(f"Successfully created region with ID: {region.region_id}")
    print(f"  Band gap: {region.band_gap} eV")
    print(f"  CB effective mass: {region.cb_effective_mass}")
    print(f"  VB effective mass (avg): {region.vb_effective_mass_avg:.6f}")
    print(f"  kT: {region.kT:.6f} eV")
    print(f"  Nc: {region.Nc:.3e} cm^-3")
    print(f"  Nv: {region.Nv:.3e} cm^-3")
    print(f"  ni: {region.ni:.3e} cm^-3")
    
    return region

def test_charge_density_calculation(region):
    """Test charge density calculation functions."""
    print("\nTesting charge density calculations...")
    
    # Create calculator
    fermi_level = 1.4186435  # From Fortran output
    calculator = ChargeDensityCalculator([region], [], fermi_level)
    
    # Test at zero potential 
    potential = 0.0
    
    print(f"  Using Fermi level: {fermi_level:.6f} eV")
    print(f"  Using potential: {potential} V")
    
    # Test electron density
    try:
        n = calculator._electron_density_direct(region, fermi_level, potential)
        print(f"  Electron density: {n:.6e} cm^-3")
        
        # Check if non-zero and reasonable
        if n > 0 and n < 1e30:
            print("  ✓ Electron density is non-zero and reasonable")
        else:
            print("  ✗ Electron density is zero or unreasonable")
            
    except Exception as e:
        print(f"  ✗ Error calculating electron density: {e}")
    
    # Test hole density
    try:
        p = calculator._hole_density_direct(region, fermi_level, potential)
        print(f"  Hole density: {p:.6e} cm^-3")
        
        # Check if non-zero and reasonable
        if p > 0 and p < 1e30:
            print("  ✓ Hole density is non-zero and reasonable")
        else:
            print("  ✗ Hole density is zero or unreasonable")
            
    except Exception as e:
        print(f"  ✗ Error calculating hole density: {e}")
    
    # Test bulk charge density
    try:
        rho = calculator.calculate_bulk_density(region.region_id, fermi_level, potential)
        print(f"  Bulk charge density: {rho:.6e} C/cm^3")
        
        # Check if finite
        if np.isfinite(rho):
            print("  ✓ Bulk charge density is finite")
        else:
            print("  ✗ Bulk charge density is not finite")
            
    except Exception as e:
        print(f"  ✗ Error calculating bulk charge density: {e}")
    
    return n, p, rho

def test_fermi_level_range(region):
    """Test charge densities across a range of Fermi levels."""
    print("\nTesting charge densities across Fermi level range...")
    
    # Test range around the expected Fortran value
    fermi_levels = np.linspace(1.4, 1.45, 11)
    potential = 0.0
    
    print("EF (eV)   | n (cm^-3)    | p (cm^-3)    | rho (C/cm^3)")
    print("-" * 60)
    
    all_zero_n = True
    all_zero_p = True
    
    for ef in fermi_levels:
        calculator = ChargeDensityCalculator([region], [], ef)
        
        try:
            n = calculator._electron_density_direct(region, ef, potential)
            p = calculator._hole_density_direct(region, ef, potential)
            rho = calculator.calculate_bulk_density(region.region_id, ef, potential)
            
            print(f"{ef:.6f} | {n:.3e} | {p:.3e} | {rho:.3e}")
            
            if n > 0:
                all_zero_n = False
            if p > 0:
                all_zero_p = False
                
        except Exception as e:
            print(f"{ef:.6f} | ERROR: {e}")
    
    # Summary
    print("\nSummary:")
    if not all_zero_n:
        print("  ✓ Non-zero electron densities found")
    else:
        print("  ✗ All electron densities are zero")
        
    if not all_zero_p:
        print("  ✓ Non-zero hole densities found")
    else:
        print("  ✗ All hole densities are zero")

def main():
    """Main test function."""
    print("=" * 60)
    print("CHARGE DENSITY FUNCTION TEST")
    print("=" * 60)
    
    try:
        # Test 1: Create semiconductor region
        region = test_semiconductor_creation()
        
        # Test 2: Basic charge density calculation
        n, p, rho = test_charge_density_calculation(region)
        
        # Test 3: Fermi level range test
        test_fermi_level_range(region)
        
        print("\n" + "=" * 60)
        print("TEST COMPLETED")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
