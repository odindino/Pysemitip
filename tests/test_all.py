#!/usr/bin/env python3
"""
Comprehensive test script for Pysemitip.

This script runs all tests and basic simulations to verify the installation.
"""

import sys
import os
import subprocess
from pathlib import Path

# Colors for output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def run_command(cmd, description):
    """Run a command and report results."""
    print(f"\n{BLUE}Testing: {description}{RESET}")
    print(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"{GREEN}[OK] Success{RESET}")
            if result.stdout:
                print(result.stdout)
            return True
        else:
            print(f"{RED}[FAIL] Failed{RESET}")
            print(result.stderr)
            if result.stdout:
                print(result.stdout)
            return False
    except Exception as e:
        print(f"{RED}[FAIL] Error: {e}{RESET}")
        return False

def main():
    """Run all tests."""
    print(f"{BLUE}=== Pysemitip Comprehensive Test ==={RESET}")
    
    # Change to project directory
    project_dir = Path(__file__).parent
    os.chdir(project_dir)
    
    all_passed = True
    
    # 1. Setup directories
    if run_command("python setup_directories.py", "Setting up directory structure"):
        pass
    else:
        print(f"{YELLOW}Warning: Could not set up directories{RESET}")
    
    # 2. Test basic imports
    if not run_command("python test_basic.py", "Basic import test"):
        all_passed = False
        print(f"{RED}Basic imports failed. Please check your installation.{RESET}")
        return 1
    
    # 3. Run unit tests
    if not run_command("python run_tests.py", "Unit tests"):
        all_passed = False
        print(f"{YELLOW}Some unit tests failed, but continuing...{RESET}")
    
    # 4. Test quick simulation
    print(f"\n{BLUE}Running quick simulation test...{RESET}")
    quick_test_config = "data/input/examples/quick_test.yaml"
    
    if Path(quick_test_config).exists():
        cmd = f"python run_multint.py {quick_test_config} --name quick_test"
        if run_command(cmd, "Quick simulation"):
            print(f"{GREEN}[OK] Simulation completed successfully{RESET}")
            
            # Check if output files were created
            output_dir = Path("data/output/results")
            if output_dir.exists():
                results = list(output_dir.glob("quick_test_*/multint_results.pkl"))
                if results:
                    print(f"{GREEN}[OK] Results file created: {results[0]}{RESET}")
                else:
                    print(f"{YELLOW}Warning: No results file found{RESET}")
        else:
            all_passed = False
            print(f"{RED}Simulation failed{RESET}")
    else:
        print(f"{YELLOW}Warning: Quick test config not found{RESET}")
    
    # 5. Summary
    print(f"\n{BLUE}=== Test Summary ==={RESET}")
    if all_passed:
        print(f"{GREEN}[OK] All tests passed!{RESET}")
        print("\nYou can now run simulations with:")
        print("  python run_multint.py data/input/examples/test_multint.yaml --plot")
        return 0
    else:
        print(f"{YELLOW}[WARN] Some tests failed, but the package may still be usable.{RESET}")
        print("Check the error messages above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())