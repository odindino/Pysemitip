#!/usr/bin/env python3
"""
Detailed debugging script for Poisson solver band bending issue.
This script will create detailed logs to help identify why band bending remains 0.
"""

import sys
import os
import logging
from pathlib import Path
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from src.core.filereader import load_yaml_config
from src.simulation.multint import MultIntSimulation

def setup_detailed_logging():
    """Set up detailed logging to file and console."""
    # Create logs directory
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    
    # Set up logging
    log_file = log_dir / "poisson_debug_detailed.log"
    
    # Clear previous log
    if log_file.exists():
        log_file.unlink()
    
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return log_file

def debug_poisson_solver():
    """Run detailed debugging of Poisson solver."""
    
    log_file = setup_detailed_logging()
    logger = logging.getLogger(__name__)
    
    logger.info("="*60)
    logger.info("Starting detailed Poisson solver debugging")
    logger.info("="*60)
    
    try:
        # Load configuration
        config_path = "data/input/examples/quick_test.yaml"
        logger.info(f"Loading configuration: {config_path}")
        config = load_yaml_config(config_path)
        
        # Create simulation
        logger.info("Creating MultInt simulation...")
        sim = MultIntSimulation(config)
        
        # Test with first bias voltage
        bias_voltages = [-2.0, 2.0]  # From config
        test_bias = bias_voltages[0]
        logger.info(f"Testing with bias voltage: {test_bias} V")
        
        # Set up simulation for this bias
        sim.tip.bias_voltage = test_bias
        logger.info(f"Tip potential: {sim.tip.tip_potential} V")
        
        # Create charge density tables (simplified)
        logger.info("Setting up charge density...")
        sim.fermi_level = 1.418687  # From config
        
        # Try to access components directly
        logger.info("Accessing grid and solvers...")
        grid = sim.grid
        poisson_solver = sim.poisson_solver
        
        # Debug grid details
        logger.info(f"Grid parameters:")
        logger.info(f"  nr={grid.params.nr}, nv={grid.params.nv}, ns={grid.params.ns}, np={grid.params.np}")
        logger.info(f"  delr={grid.params.delr}, delv={grid.params.delv}, dels={grid.params.dels}, delp={grid.params.delp}")
        logger.info(f"  Grid r values: {grid.r[:5]} ... {grid.r[-5:]}")
        
        # Debug tip mask
        tip_count = np.sum(sim.grid.tip_mask)
        logger.info(f"Tip mask: {tip_count} points inside tip out of {np.prod(sim.grid.tip_mask.shape)} total")
        
        # Debug Poisson solver initial state
        logger.info("Debugging Poisson solver initial state...")
        poisson = poisson_solver
        
        # Initialize potential
        logger.info("Setting initial potential...")
        poisson._set_initial_potential()
        
        # Check initial potentials
        logger.info(f"Initial potential stats:")
        logger.info(f"  Vacuum: min={np.min(poisson.potential_vac):.6f}, max={np.max(poisson.potential_vac):.6f}, mean={np.mean(poisson.potential_vac):.6f}")
        logger.info(f"  Interface: min={np.min(poisson.potential_int):.6f}, max={np.max(poisson.potential_int):.6f}, mean={np.mean(poisson.potential_int):.6f}")
        logger.info(f"  Semiconductor: min={np.min(poisson.potential_sem):.6f}, max={np.max(poisson.potential_sem):.6f}, mean={np.mean(poisson.potential_sem):.6f}")
        
        # Test band bending calculation
        initial_bb = poisson.get_band_bending()
        logger.info(f"Initial band bending: {initial_bb:.8f} V")
        
        # Simple dummy charge functions for testing
        def dummy_bulk_charge(r, z, phi, pot):
            return 0.1  # Small non-zero value
        
        def dummy_surface_charge(r, phi, pot):
            return 0.01  # Small non-zero value
        
        # Test charge functions
        logger.info("Testing dummy charge functions...")
        bulk_test = dummy_bulk_charge(1.0, -1.0, 0.0, 0.0)
        surface_test = dummy_surface_charge(0.0, 0.0, 0.0)
        logger.info(f"Test bulk charge density: {bulk_test:.2e}")
        logger.info(f"Test surface charge density: {surface_test:.2e}")
        
        # Run a few iterations manually to see what changes
        logger.info("Running 5 manual iterations...")
        eep = 1.80943e-20
        
        for iteration in range(1, 6):
            logger.info(f"--- Iteration {iteration} ---")
            
            # Store old potentials
            old_vac = poisson.potential_vac.copy()
            old_int = poisson.potential_int.copy()
            old_sem = poisson.potential_sem.copy()
            
            # Update vacuum
            poisson._update_vacuum_fortran_style()
            vac_change = np.max(np.abs(poisson.potential_vac - old_vac))
            logger.info(f"  Vacuum potential max change: {vac_change:.2e}")
            
            # Update interface  
            poisson._update_interface_fortran_style(dummy_surface_charge, eep)
            int_change = np.max(np.abs(poisson.potential_int - old_int))
            logger.info(f"  Interface potential max change: {int_change:.2e}")
            
            # Update semiconductor
            poisson._update_semiconductor_fortran_style(dummy_bulk_charge, eep)
            sem_change = np.max(np.abs(poisson.potential_sem - old_sem))
            logger.info(f"  Semiconductor potential max change: {sem_change:.2e}")
            
            # Check band bending
            bb = poisson.get_band_bending()
            logger.info(f"  Band bending: {bb:.8f} V")
            
            # Update full potential
            poisson._update_full_potential()
        
        # Test what happens if we artificially set interface potential
        logger.info("Testing artificial interface potential...")
        poisson.potential_int[0, 0] = 0.1  # Set 0.1V at origin
        bb_artificial = poisson.get_band_bending()
        logger.info(f"Artificial band bending (set to 0.1V): {bb_artificial:.8f} V")
        
        # Test interface calculation details
        logger.info("Testing interface potential calculation details...")
        nr, np_ = poisson.grid.params.nr, poisson.grid.params.np
        i, k = 0, 0  # Origin point
        r = poisson.grid.r[i]
        phi = (k - 0.5) * poisson.grid.params.delp
        
        # Calculate stemp manually
        stemp = (3.0 * poisson.potential_vac[i, 1, k] - 
                (9.0/6.0) * poisson.potential_vac[i, 2, k] + 
                (1.0/3.0) * poisson.potential_vac[i, 3, k]) / poisson.grid.params.delv
        
        eps_semi = 12.9
        stemp += eps_semi * (3.75 * poisson.potential_sem[i, 1, k] - 
                           (5.0/6.0) * poisson.potential_sem[i, 2, k] + 
                           0.15 * poisson.potential_sem[i, 3, k]) / poisson.grid.params.dels
        
        denom = ((11.0/6.0) / poisson.grid.params.delv + 
                (46.0/15.0) * eps_semi / poisson.grid.params.dels)
        
        logger.info(f"Manual interface calculation at origin:")
        logger.info(f"  stemp = {stemp:.6f}")
        logger.info(f"  denom = {denom:.6f}")
        logger.info(f"  dummy surface charge at origin = {dummy_surface_charge(r, phi, 0.0):.2e}")
        
        # Check vacuum potentials near interface
        logger.info(f"Vacuum potentials near interface (i=0):")
        for j in range(min(4, poisson.grid.params.nv)):
            logger.info(f"  VAC[0,{j},0] = {poisson.potential_vac[0, j, 0]:.6f}")
        
        # Check semiconductor potentials near interface
        logger.info(f"Semiconductor potentials near interface (i=0):")
        for j in range(min(4, poisson.grid.params.ns)):
            logger.info(f"  SEM[0,{j},0] = {poisson.potential_sem[0, j, 0]:.6f}")
        
    except Exception as e:
        logger.error(f"Error during debugging: {e}", exc_info=True)
    
    logger.info("="*60)
    logger.info(f"Debugging complete. Detailed log saved to: {log_file}")
    logger.info("="*60)

if __name__ == "__main__":
    debug_poisson_solver()