#!/usr/bin/env python3
"""
Test contact_potential property access.
"""

import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, os.path.dirname(__file__))

def test_contact_potential():
    """Test the contact_potential property access."""
    try:
        from src.core.filereader import YamlConfigReader
        
        # Test loading a config file
        config_path = "data/input/examples/quick_test.yaml"
        if Path(config_path).exists():
            reader = YamlConfigReader()
            config = reader.load_config(config_path)
            
            print(f"[OK] Configuration loaded successfully")
            print(f"  Contact potential (tip): {config.tip.contact_potential}")
            print(f"  Contact potential (root): {config.contact_potential}")
            print(f"  Both should be equal: {config.tip.contact_potential == config.contact_potential}")
            
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
    success = test_contact_potential()
    sys.exit(0 if success else 1)