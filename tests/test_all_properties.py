#!/usr/bin/env python3
"""
Test all root-level property access.
"""

import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

def test_all_properties():
    """Test all the backward compatibility properties."""
    try:
        from src.core.filereader import YamlConfigReader
        
        # Test loading a config file
        config_path = "data/input/examples/quick_test.yaml"
        if Path(config_path).exists():
            reader = YamlConfigReader()
            config = reader.load_config(config_path)
            
            print(f"[OK] Configuration loaded successfully")
            print(f"  Contact potential: {config.contact_potential}")
            print(f"  Mirror symmetry: {config.mirror_symmetry}")
            print(f"  Charge table points: {config.charge_table_points}")
            print(f"  Max iterations: {config.max_iterations}")
            print(f"  Convergence tolerance: {config.convergence_tolerance}")
            print(f"  Temperature: {config.temperature}")
            print(f"  Dielectric constant: {config.dielectric_constant}")
            print(f"  Number of semiconductor regions: {len(config.semiconductor_regions)}")
            print(f"  Number of surface regions: {len(config.surface_regions)}")
            
            return True
        else:
            print(f"[FAIL] Config file not found: {config_path}")
            return False
            
    except Exception as e:
        print(f"[FAIL] Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_all_properties()
    sys.exit(0 if success else 1)