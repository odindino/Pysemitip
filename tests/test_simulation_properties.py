#!/usr/bin/env python3
"""
Comprehensive test for all simulation properties required by MultInt.
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_simulation_properties():
    """Test all properties required by the MultInt simulation."""
    try:
        from src.core.filereader import YamlConfigReader
        
        # Test loading a config file
        config_path = project_root / "data/input/examples/test/quick_test.yaml"
        if config_path.exists():
            reader = YamlConfigReader()
            config = reader.load_config(config_path)
            
            print(f"[OK] Configuration loaded successfully")
            
            # Test all properties used by MultInt simulation
            properties_to_test = [
                # From _log_parameters()
                ('config.tip.radius', lambda: config.tip.radius),
                ('config.tip.slope', lambda: config.tip.slope),
                ('config.contact_potential', lambda: config.contact_potential),
                ('config.tip.position.x', lambda: config.tip.position.x),
                ('config.tip.position.y', lambda: config.tip.position.y),
                ('config.semiconductor_regions', lambda: config.semiconductor_regions),
                ('config.surface_regions', lambda: config.surface_regions),
                ('config.mirror_symmetry', lambda: config.mirror_symmetry),
                ('config.temperature', lambda: config.temperature),
                
                # From voltage scan
                ('config.voltage_scan.voltages', lambda: config.voltage_scan.voltages),
                
                # From charge density calculation
                ('config.charge_table_points', lambda: config.charge_table_points),
            ]
            
            all_passed = True
            for prop_name, prop_getter in properties_to_test:
                try:
                    value = prop_getter()
                    print(f"  ✓ {prop_name}: {value}")
                except Exception as e:
                    print(f"  ✗ {prop_name}: ERROR - {e}")
                    all_passed = False
            
            if all_passed:
                print(f"\n[OK] All simulation properties are accessible!")
                return True
            else:
                print(f"\n[FAIL] Some properties failed")
                return False
        else:
            print(f"[FAIL] Config file not found: {config_path}")
            return False
            
    except Exception as e:
        print(f"[FAIL] Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_simulation_properties()
    sys.exit(0 if success else 1)