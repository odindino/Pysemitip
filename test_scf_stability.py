#!/usr/bin/env python3
"""Test SCF loop stability with new fixes"""
import sys
import os
sys.path.append(os.path.dirname(__file__))

from src.simulation.multint import MultInt
from src.core.config_schema import SemitipConfig
import logging
import yaml

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_scf_stability():
    """Test a few SCF iterations with the new stable solver"""
    try:
        # Load config
        config_path = "data/input/examples/test/quick_test.yaml"
        logger.info(f"Loading config from {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)
        config = SemitipConfig(**config_data)
        
        # Create output directory
        output_dir = "data/output/results/scf_stability_test"
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize MultInt
        multint = MultInt(config, output_dir)
        
        logger.info("Testing reduced SCF loop (5 iterations max)...")
        
        # Modify computation settings for quick test
        if hasattr(config, 'computation'):
            config.computation.max_iterations = [5, 5, 5, 5]  # Very short test
            config.computation.convergence_parameters = [1e-2, 1e-2, 1e-2, 1e-2]  # Loose tolerance
        
        # Run a short SCF test
        multint.run_self_consistent_loop()
        
        logger.info("✅ SCF stability test completed successfully!")
        return True
        
    except Exception as e:
        logger.error(f"❌ SCF stability test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_scf_stability()
    sys.exit(0 if success else 1)