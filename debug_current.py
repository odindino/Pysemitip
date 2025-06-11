#!/usr/bin/env python3
"""
Debug tunneling current calculation
"""

import numpy as np
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from physics.materials import MaterialDatabase
from physics.charge_density import ChargeDensityCalculator
from physics.tunneling_current import TunnelingCurrentCalculator

def debug_tunneling_current():
    """Debug the tunneling current calculation step by step"""
    
    # Get material
    db = MaterialDatabase()
    material = db.get_material("Si_n")
    
    print(f"Material: {material.name}")
    print(f"Bandgap: {material.bandgap:.3f} eV")
    print(f"Net doping: {material.net_doping:.1e} cm⁻³")
    print()
    
    # Create calculator
    calculator = TunnelingCurrentCalculator()
    charge_calc = ChargeDensityCalculator()
    
    # Find Fermi level
    fermi_level = charge_calc.find_equilibrium_fermi_level(material)
    print(f"Fermi level: {fermi_level:.3f} eV")
    
    # Test simple barrier
    separation = 1.0  # nm
    z_grid = np.linspace(0, separation, 50)
    bias_voltage = 1.0  # V
    
    # Create simple triangular barrier
    work_function_tip = 4.5  # eV
    work_function_sample = 4.0  # eV
    barrier_height = 2.0  # eV additional barrier
    
    potential_profile = np.zeros_like(z_grid)
    for i, z in enumerate(z_grid):
        # Linear potential drop
        linear_drop = work_function_tip + (work_function_sample + bias_voltage - work_function_tip) * z / separation
        potential_profile[i] = linear_drop + barrier_height
    
    print(f"Barrier at tip (z=0): {potential_profile[0]:.3f} eV")
    print(f"Barrier at sample (z={separation}): {potential_profile[-1]:.3f} eV")
    print()
    
    # Test transmission coefficient
    energy = fermi_level + bias_voltage/2  # Mid-gap energy
    k_parallel = 0.1  # nm^-1
    
    transmission = calculator.calculate_transmission_coefficient(
        energy, k_parallel, potential_profile, z_grid
    )
    print(f"Test transmission coefficient:")
    print(f"  Energy: {energy:.3f} eV")
    print(f"  k_parallel: {k_parallel:.3f} nm⁻¹")
    print(f"  Transmission: {transmission:.2e}")
    print()
    
    # Test density of states
    dos = calculator.calculate_density_of_states(material, energy, "conduction")
    print(f"Density of states (conduction): {dos:.2e} eV⁻¹⋅cm⁻³")
    
    dos_valence = calculator.calculate_density_of_states(material, energy, "valence_heavy")
    print(f"Density of states (valence): {dos_valence:.2e} eV⁻¹⋅cm⁻³")
    print()
    
    # Test wavefunction
    wf_val, wf_deriv = calculator.calculate_sample_wavefunction(
        energy, k_parallel, material, "conduction"
    )
    print(f"Sample wavefunction:")
    print(f"  Value: {wf_val:.2e}")
    print(f"  Derivative: {wf_deriv:.2e}")
    print()
    
    # Calculate full current
    result = calculator.calculate_total_current(
        material, bias_voltage, fermi_level, potential_profile,
        z_grid, 0.0, separation
    )
    
    print(f"Current calculation results:")
    print(f"  Total current: {result['total_current']:.2e} A")
    print(f"  Extended states: {result['extended_states']['total']:.2e} A")
    print(f"  Localized states: {result['localized_states']:.2e} A")
    print()
    
    for band, current in result['current_components'].items():
        if abs(current) > 1e-20:
            print(f"  {band}: {current:.2e} A")

if __name__ == "__main__":
    debug_tunneling_current()
