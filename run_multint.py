#!/usr/bin/env python3
"""
Main entry point for running MultInt simulations.

Usage:
    python run_multint.py config.yaml [options]
    
Options:
    --plot          Generate plots after simulation
    --output DIR    Output directory (default: output/)
    --biases LIST   Override bias voltages (comma-separated)
"""

import sys
import os
import argparse
import logging
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

from src.simulation.multint import run_multint_simulation
from src.physics.visualization import plot_simulation_results, create_contour_plots


def main():
    """Main function."""
    # Parse arguments
    parser = argparse.ArgumentParser(description='Run MultInt SEMITIP simulation')
    parser.add_argument('config', help='Path to YAML configuration file')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--output', help='Output directory (default: auto-generated)')
    parser.add_argument('--name', help='Simulation name (default: from config filename)')
    parser.add_argument('--biases', help='Override bias voltages (comma-separated)')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, 
                       format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Determine output directory
    if args.output:
        output_dir = Path(args.output)
    else:
        # Auto-generate output directory based on timestamp and config name
        from datetime import datetime
        config_name = Path(args.config).stem
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        sim_name = args.name or config_name
        output_dir = Path('data/output/results') / f'{sim_name}_{timestamp}'
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse bias voltages if provided
    bias_voltages = None
    if args.biases:
        try:
            bias_voltages = [float(b.strip()) for b in args.biases.split(',')]
        except ValueError:
            print(f"Error: Invalid bias voltages: {args.biases}")
            sys.exit(1)
    
    # Change to output directory
    original_dir = os.getcwd()
    print(f"Running MultInt simulation with config: {args.config}")
    print(f"Output directory: {output_dir}")
    os.chdir(output_dir)
    
    try:
        
        config_path = Path(original_dir) / args.config
        results = run_multint_simulation(str(config_path))
        
        print(f"\nSimulation completed successfully!")
        print(f"Results saved to: {output_dir / 'multint_results.pkl'}")
        
        # Generate plots if requested
        if args.plot:
            print("\nGenerating plots...")
            
            # Create subdirectories - use absolute paths
            figures_dir = Path('figures')
            contours_dir = Path('contours')
            figures_dir.mkdir(exist_ok=True)
            contours_dir.mkdir(exist_ok=True)
            
            # Generate standard plots
            plot_simulation_results('multint_results.pkl', str(figures_dir))
            
            # Generate contour plots
            create_contour_plots('multint_results.pkl', str(contours_dir))
            
            print(f"Plots saved to: {figures_dir} and {contours_dir}")
    
    finally:
        # Return to original directory
        os.chdir(original_dir)


if __name__ == '__main__':
    main()