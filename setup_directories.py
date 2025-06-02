#!/usr/bin/env python3
"""
Setup directory structure for Pysemitip project.

This script creates the necessary directories for inputs, outputs, and tests.
"""

import os
from pathlib import Path

def setup_directories():
    """Create directory structure."""
    # Base directories
    base_dir = Path(__file__).parent
    
    # Directory structure
    directories = [
        # Input/Output directories
        "data/input/configs",      # User configuration files
        "data/input/examples",     # Example configurations
        "data/output/results",     # Simulation results
        "data/output/logs",        # Log files
        "data/output/figures",     # Generated plots
        "data/output/contours",    # Contour plots
        
        # Test directories
        "tests/data/input",        # Test input files
        "tests/data/output",       # Test output files
        "tests/data/reference",    # Reference results for validation
        
        # Documentation
        "docs/examples",           # Example notebooks/scripts
        "docs/api",               # API documentation
        "docs/tutorials",         # User tutorials
    ]
    
    # Create directories
    for dir_path in directories:
        full_path = base_dir / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"Created: {dir_path}")
    
    # Move example configs to proper location
    example_configs = [
        ("config/examples/test_multint.yaml", "data/input/examples/test_multint.yaml"),
        ("config/examples/example_multint.yaml", "data/input/examples/example_multint.yaml"),
        ("config/examples/example_multint_fixed.yaml", "data/input/examples/example_multint_fixed.yaml"),
        ("config/examples/example_multplane.yaml", "data/input/examples/example_multplane.yaml"),
        ("config/examples/example_multplane_fixed.yaml", "data/input/examples/example_multplane_fixed.yaml"),
    ]
    
    for src, dst in example_configs:
        src_path = base_dir / src
        dst_path = base_dir / dst
        if src_path.exists() and not dst_path.exists():
            import shutil
            shutil.copy2(src_path, dst_path)
            print(f"Copied: {src} -> {dst}")
    
    # Create README files
    readme_content = {
        "data/input/configs/README.md": """# User Configuration Files

Place your SEMITIP configuration files here.

Example usage:
```bash
python run_multint.py data/input/configs/my_config.yaml
```
""",
        "data/output/results/README.md": """# Simulation Results

This directory contains simulation output files:
- `*.pkl` - Pickled simulation results
- `*.log` - Simulation log files

Results are organized by timestamp or configuration name.
""",
        "tests/data/reference/README.md": """# Reference Results

This directory contains reference results for validation testing.
These files are used to verify that the Python implementation
produces results consistent with the Fortran version.
"""
    }
    
    for readme_path, content in readme_content.items():
        full_path = base_dir / readme_path
        with open(full_path, 'w') as f:
            f.write(content)
        print(f"Created: {readme_path}")
    
    print("\nDirectory structure setup complete!")
    
    # Print directory tree
    print("\nProject structure:")
    print_tree(base_dir)

def print_tree(directory, prefix="", max_depth=3, current_depth=0):
    """Print directory tree."""
    if current_depth >= max_depth:
        return
    
    # Skip some directories
    skip_dirs = {'.git', '__pycache__', '.pytest_cache', 'src'}
    
    items = sorted(directory.iterdir())
    for i, item in enumerate(items):
        if item.name.startswith('.') and item.name not in {'.gitignore'}:
            continue
        if item.name in skip_dirs:
            continue
            
        is_last = i == len(items) - 1
        current_prefix = "└── " if is_last else "├── "
        print(f"{prefix}{current_prefix}{item.name}")
        
        if item.is_dir():
            extension = "    " if is_last else "│   "
            print_tree(item, prefix + extension, max_depth, current_depth + 1)

if __name__ == "__main__":
    setup_directories()