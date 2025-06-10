#!/usr/bin/env python3
"""
Debug numerical differences between Python and Fortran versions.
"""

import numpy as np
from src.physics.materials.semiconductor import SemiconductorRegion as SemiconductorRegionPhysics
from src.utils.constants import PhysicalConstants as PC

def debug_carrier_density_calculation():
    """Debug carrier density calculation to find differences with Fortran."""
    print("=== Debugging Carrier Density Calculation ===\n")
    
    # Create semiconductor region with exact parameters
    # Using donor concentration close to Fortran value
    region = SemiconductorRegionPhysics(
        region_id=1,
        donor_concentration=9.99999984e17,  # Exact Fortran value
        acceptor_concentration=0.0,
        band_gap=1.42,
        valence_band_offset=0.0,
        donor_binding_energy=0.006,
        acceptor_binding_energy=0.028,
        cb_effective_mass=0.0635,
        vb_effective_mass_heavy=0.643,
        vb_effective_mass_light=0.081,
        vb_effective_mass_so=0.145,
        spin_orbit_splitting=0.341,
        permittivity=12.9,
        temperature=300.0
    )
    
    print(f"Region parameters:")
    print(f"  Donor concentration: {region.donor_concentration:.8e} cm^-3")
    print(f"  Band gap: {region.band_gap} eV")
    print(f"  CB effective mass: {region.cb_effective_mass}")
    print(f"  VB effective mass (avg): {region.vb_effective_mass_avg:.6f}")
    print(f"  Temperature: {region.temperature} K")
    print()
    
    # Temperature conversion
    fortran_tk = 300.0 * 8.617e-5  # Fortran: TK=TEM*8.617E-5
    python_tk = PC.KB_EV * region.temperature  # Python version
    print(f"Temperature conversion:")
    print(f"  Fortran kT: {fortran_tk:.8e} eV")
    print(f"  Python kT:  {python_tk:.8e} eV")
    print(f"  Difference: {abs(fortran_tk - python_tk):.2e} eV")
    print()
    
    # Effective density of states
    C_fortran = 6.815e21  # Fortran constant
    region.kT = fortran_tk  # Use Fortran temperature for comparison
    Nc_fortran_style = C_fortran * np.sqrt((region.cb_effective_mass * fortran_tk)**3)
    Nv_fortran_style = C_fortran * np.sqrt((region.vb_effective_mass_avg * fortran_tk)**3)
    
    print(f"Effective density of states (using Fortran kT):")
    print(f"  Nc: {Nc_fortran_style:.8e} cm^-3")
    print(f"  Nv: {Nv_fortran_style:.8e} cm^-3")
    print()
    
    # Calculate Fermi level using our method
    try:
        fermi_level = region.fermi_level()
        print(f"Fermi level calculation:")
        print(f"  Python result: {fermi_level:.7f} eV")
        print(f"  Fortran result: 1.4186435 eV")
        print(f"  Difference: {abs(fermi_level - 1.4186435):.7f} eV")
        print()
        
        # Calculate carrier densities
        n_cb = region.carrier_density_cb(fermi_level)
        n_vb = region.carrier_density_vb(fermi_level)
        
        print(f"Carrier densities:")
        print(f"  Python CB: {n_cb:.8e} cm^-3")
        print(f"  Python VB: {n_vb:.6f} cm^-3")
        print(f"  Fortran CB: 2.94679424e+17 cm^-3")
        print(f"  Fortran VB: 57.446033 cm^-3")
        print(f"  CB difference: {abs(n_cb - 2.94679424e17):.2e} cm^-3")
        print(f"  VB difference: {abs(n_vb - 57.446033):.6f} cm^-3")
        print()
        
    except Exception as e:
        print(f"Error in fermi level calculation: {e}")
        print()

def debug_fortran_constants():
    """Check if our constants match Fortran exactly."""
    print("=== Debugging Physical Constants ===\n")
    
    # Temperature conversion
    fortran_factor = 8.617e-5
    python_factor = PC.KB_EV
    print(f"Temperature conversion factor:")
    print(f"  Fortran: {fortran_factor:.6e} eV/K")
    print(f"  Python:  {python_factor:.6e} eV/K") 
    print(f"  Relative difference: {abs(fortran_factor - python_factor) / fortran_factor:.2e}")
    print()
    
    # Carrier density constant
    fortran_c = 6.815e21
    print(f"Carrier density constant:")
    print(f"  Fortran: {fortran_c:.3e} eV^-1.5 cm^-3")
    print(f"  Formula: (2/√π) × 2 × (m/(2π*ℏ²))^1.5")
    
    # Calculate the theoretical value
    m0 = PC.M0  # kg
    hbar = PC.HBAR  # J·s
    hbar_ev = PC.HBAR_EV  # eV·s
    
    # Theoretical constant in SI then convert to eV^-1.5 cm^-3
    theoretical_si = (2/np.sqrt(np.pi)) * 2 * (m0/(2*np.pi*hbar**2))**(1.5)
    # Convert: multiply by (J to eV)^1.5 and (m^-3 to cm^-3)
    theoretical_ev_cm3 = theoretical_si * (PC.E**1.5) * 1e-6
    
    print(f"  Theoretical: {theoretical_ev_cm3:.3e} eV^-1.5 cm^-3")
    print(f"  Difference from Fortran: {abs(theoretical_ev_cm3 - fortran_c):.2e}")
    print()

if __name__ == "__main__":
    debug_fortran_constants()
    debug_carrier_density_calculation()