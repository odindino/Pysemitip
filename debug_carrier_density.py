#!/usr/bin/env python3
"""
Debug script to test carrier density calculations.
"""

import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

from src.physics.materials.semiconductor import SemiconductorRegion
from src.physics.core.charge_density import ChargeDensityCalculator

# Create a simple semiconductor region matching EXACT Fortran parameters
region = SemiconductorRegion(
    region_id=1,
    donor_concentration=1e18,  # CD=1e18 cm^-3
    acceptor_concentration=0.0,  # CA=0
    band_gap=1.42,  # EGAP=1.42 eV
    electron_affinity=4.07,
    valence_band_offset=0.0,  # DELVB=0
    donor_binding_energy=0.006,  # ED=0.006 eV
    acceptor_binding_energy=0.028,  # EA=0.028 eV (EXACT from Fortran line 15)
    cb_effective_mass=0.0635,  # ACB=0.0635 (EXACT from Fortran line 16)
    vb_effective_mass_heavy=0.643,  # AVBH=0.643 (EXACT from Fortran line 17)
    vb_effective_mass_light=0.081,  # AVBL=0.081 (EXACT from Fortran line 18)
    vb_effective_mass_so=0.172,  # AVBSO=0.172 (EXACT from Fortran line 19)
    spin_orbit_splitting=0.341,  # ESO=0.341 (EXACT from Fortran line 20)
    permittivity=12.9,
    allow_degeneracy=False,  # Use Boltzmann for simplicity
    temperature=300.0
)

print(f"Region properties:")
print(f"  Net doping: {region.net_doping:.3e} cm^-3")
print(f"  Is n-type: {region.is_n_type}")
print(f"  CB effective mass: {region.conduction_band_effective_mass}")
print(f"  VB effective mass: {region.valence_band_effective_mass}")
print(f"  kT: {region.kT:.6f} eV")

# Calculate Fermi level
fermi_level = region.fermi_level()
print(f"  Fermi level: {fermi_level:.6f} eV")

# Debug: check effective masses used
print(f"  AVB (valence band effective mass): {region.vb_effective_mass_avg:.6f}")
print(f"  ACB (conduction band effective mass): {region.cb_effective_mass:.6f}")

# Debug: manually check Fortran EF calculation for intrinsic case 
# EF = EGAP/2 + 0.75*TK*ALOG(AVB/ACB)
ef_manual = region.band_gap/2.0 + 0.75*region.kT*np.log(region.vb_effective_mass_avg/region.cb_effective_mass)
print(f"  Manual EF (intrinsic formula): {ef_manual:.6f} eV")

# For doped case, use simple approximation
# n-type: EF ≈ EGAP - kT*ln(Nc/Nd)
# But we need to use the Fortran constant properly
C = 6.815e21
Nc_fortran = C * np.sqrt((region.cb_effective_mass * region.kT)**3)
ef_doped = region.band_gap - region.kT * np.log(Nc_fortran / region.donor_concentration)
print(f"  Manual EF (doped approximation): {ef_doped:.6f} eV")

# Create charge density calculator  
calculator = ChargeDensityCalculator([region], [], fermi_level)

# Test electron and hole densities at zero potential
potential = 0.0
n = calculator._electron_density_direct(region, fermi_level, potential)
p = calculator._hole_density_direct(region, fermi_level, potential)

print(f"\nCarrier densities at EF={fermi_level:.6f} eV, Pot={potential} V:")
print(f"  Electron density: {n:.6e} cm^-3")
print(f"  Hole density: {p:.6e} cm^-3")

# Compare with expected values from Fortran output:
# CARRIER DENSITY IN CB, VB = 2.94679424E+17   57.446033
print(f"\nExpected from Fortran:")
print(f"  Fermi level: 1.4186435 eV")
print(f"  Electron density: 2.947e+17 cm^-3")
print(f"  Hole density: 57.4 cm^-3")

# Test with exact Fortran Fermi level
fortran_ef = 1.4186435
n_fortran = calculator._electron_density_direct(region, fortran_ef, potential)
p_fortran = calculator._hole_density_direct(region, fortran_ef, potential)

print(f"\nUsing exact Fortran EF={fortran_ef:.6f} eV:")
print(f"  Electron density: {n_fortran:.6e} cm^-3")
print(f"  Hole density: {p_fortran:.6e} cm^-3")

# Test different Fermi levels to understand the sensitivity
test_efs = [1.410, 1.412, 1.414, 1.4150, 1.4155, 1.4160, 1.4165]
print(f"\nFermi level sensitivity test:")
for test_ef in test_efs:
    n_test = calculator._electron_density_direct(region, test_ef, potential)
    p_test = calculator._hole_density_direct(region, test_ef, potential)
    print(f"  EF={test_ef:.6f}: n={n_test:.3e}, p={p_test:.3e}")

# Find the EF that gives n = 2.947e17
target_n = 2.947e17
best_ef = None
min_diff = float('inf')
for ef_test in np.linspace(1.410, 1.420, 100):
    n_test = calculator._electron_density_direct(region, ef_test, potential)
    diff = abs(n_test - target_n)
    if diff < min_diff:
        min_diff = diff
        best_ef = ef_test

print(f"\nBest EF for target n={target_n:.3e}: {best_ef:.6f}")
n_best = calculator._electron_density_direct(region, best_ef, potential)
p_best = calculator._hole_density_direct(region, best_ef, potential)
print(f"  n={n_best:.3e}, p={p_best:.3e}")

# Check the eta values for Fortran EF
fortran_ef = 1.4186435
eta_cb = (fortran_ef - region.band_gap - region.valence_band_offset - potential) / region.kT
eta_vb = (-fortran_ef + region.valence_band_offset + potential) / region.kT
print(f"\nFor Fortran EF={fortran_ef:.6f}:")
print(f"  eta_CB = {eta_cb:.6f}")
print(f"  eta_VB = {eta_vb:.6f}")

# Manual hole density calculation following Fortran RHOVB exactly
C = 6.815e21
AVB = region.vb_effective_mass_avg
TK = region.kT

# Since eta_vb = -54.9 < -8, use non-degenerate formula:
# RHOVB = C*SQRT((AVB*TK)**3) * SQRT(PI)/2 * EXP(max(-40, eta_vb))
manual_p = C * np.sqrt((AVB * TK)**3) * np.sqrt(np.pi)/2.0 * np.exp(max(-40, eta_vb))
print(f"\nManual hole density calculation:")
print(f"  C = {C:.3e}")
print(f"  AVB = {AVB:.6f}")
print(f"  TK = {TK:.6f}")
print(f"  sqrt((AVB*TK)^3) = {np.sqrt((AVB * TK)**3):.6e}")
print(f"  sqrt(pi)/2 = {np.sqrt(np.pi)/2.0:.6f}")
print(f"  exp(max(-40, eta_vb)) = {np.exp(max(-40, eta_vb)):.6e}")
print(f"  Manual p = {manual_p:.6e}")