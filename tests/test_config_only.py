#!/usr/bin/env python3
"""
Test configuration loading only (without numpy dependencies).
"""

import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

def test_config_loading():
    """Test only the configuration loading part."""
    try:
        # Test the filereader import
        from src.core.filereader import YamlConfigReader
        print("[OK] YamlConfigReader imported successfully")
        
        # Test configuration schema import
        from src.core.config_schema import SemitipConfig
        print("[OK] SemitipConfig imported successfully")
        
        # Test loading a config file
        config_path = "data/input/examples/quick_test.yaml"
        if Path(config_path).exists():
            reader = YamlConfigReader()
            
            # Load YAML manually first to check types
            import yaml
            with open(config_path, 'r', encoding='utf-8') as f:
                yaml_data = yaml.safe_load(f)
            print(f"[DEBUG] First region donor_concentration type: {type(yaml_data['semiconductor_regions'][0]['donor_concentration'])}")
            print(f"[DEBUG] First region donor_concentration value: {yaml_data['semiconductor_regions'][0]['donor_concentration']}")
            
            config = reader.load_config(config_path)
            print(f"[OK] Configuration loaded successfully")
            print(f"  Simulation type: {config.simulation_type}")
            print(f"  Temperature: {config.temperature}")
            print(f"  Tip radius: {config.tip.radius}")
            print(f"  Tip work function: {config.tip.work_function}")
            print(f"  Number of semiconductor regions: {len(config.semiconductor_regions)}")
            print(f"  First region band gap: {config.semiconductor_regions[0].band_gap}")
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
    success = test_config_loading()
    sys.exit(0 if success else 1)