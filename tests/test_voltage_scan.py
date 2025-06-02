#!/usr/bin/env python3
"""
Test VoltageScanConfig voltages property.
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def test_voltage_scan():
    """Test the voltages property in VoltageScanConfig."""
    try:
        from src.core.filereader import YamlConfigReader
        
        # Test loading a config file
        config_path = project_root / "data/input/examples/quick_test.yaml"
        if config_path.exists():
            reader = YamlConfigReader()
            config = reader.load_config(config_path)
            
            print(f"[OK] Configuration loaded successfully")
            print(f"  Voltage scan start: {config.voltage_scan.start}")
            print(f"  Voltage scan end: {config.voltage_scan.end}")
            print(f"  Voltage scan points: {config.voltage_scan.points}")
            print(f"  Generated voltages: {config.voltage_scan.voltages}")
            print(f"  Voltage count: {len(config.voltage_scan.voltages)}")
            
            # Test voltages property access (this was failing before)
            voltages = config.voltage_scan.voltages
            print(f"[OK] Voltages property accessible: {len(voltages)} voltages")
            
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
    success = test_voltage_scan()
    sys.exit(0 if success else 1)