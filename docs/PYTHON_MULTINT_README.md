# Python MultInt - SEMITIP Simulation Package

## Overview

This is a Python implementation of the SEMITIP MultInt program for simulating scanning tunneling microscopy (STM) of semiconductors. It solves the 3D Poisson equation self-consistently with quantum mechanical tunneling calculations.

## Installation

1. Ensure you have Python 3.8+ installed
2. Install required packages:
```bash
pip install numpy scipy matplotlib pyyaml
```

## Quick Start

### 1. Test Installation
```bash
python test_basic.py
```

### 2. Run Tests
```bash
python run_tests.py
```

### 3. Run a Simulation
```bash
python run_multint.py config/examples/test_multint.yaml --plot
```

## Usage

### Basic Simulation

```python
from src.simulation.multint import run_multint_simulation

# Run simulation from config file
results = run_multint_simulation('config.yaml')

# Access results
for result in results:
    print(f"Bias: {result.bias_voltage} V")
    print(f"Current: {result.current} A")
    print(f"Band bending: {result.band_bending} eV")
```

### Custom Analysis

```python
from src.physics.visualization import SEMITIPPlotter

# Create plotter
plotter = SEMITIPPlotter()

# Plot I-V curve
fig = plotter.plot_current_voltage(results)
fig.savefig('iv_curve.png')

# Plot potential profile
fig = plotter.plot_potential_profile(results[0].potential_profile)
fig.savefig('potential.png')
```

## Configuration

Configuration files use YAML format. See `config/examples/` for examples.

Key parameters:
- `semiconductor_regions`: Material properties
- `surface_regions`: Surface state distributions  
- `tip`: STM tip parameters
- `voltage_scan`: Bias voltage range
- `grid`: Computational grid settings

## Project Structure

```
Pysemitip/
├── src/
│   ├── physics/
│   │   ├── core/          # Core physics calculations
│   │   ├── materials/     # Material models
│   │   ├── solvers/       # Numerical solvers
│   │   └── visualization/ # Plotting tools
│   ├── simulation/        # Main simulation programs
│   └── utils/            # Constants and utilities
├── config/               # Configuration files
│   └── examples/         # Example configs
├── tests/               # Unit tests
└── docs/                # Documentation
```

## Physics Models

- **Poisson Solver**: 3D finite difference with cylindrical coordinates
- **Charge Density**: Fermi-Dirac statistics with degeneracy
- **Tunneling Current**: WKB approximation and transfer matrix method
- **Band Structure**: Effective mass approximation

## Output Files

- `multint_results.pkl`: Pickled simulation results
- `output_MultInt.log`: Detailed calculation log
- `figures/`: Generated plots
- `contours/`: Potential contour plots

## Troubleshooting

### Import Errors
Make sure you're running from the project root directory:
```bash
cd /path/to/Pysemitip
python run_multint.py config/examples/test_multint.yaml
```

### Memory Issues
Reduce grid size in configuration:
```yaml
grid:
  radial_points: 32  # Reduce from default
  vacuum_points: 32
  semiconductor_points: 32
```

### Convergence Problems
Adjust solver parameters:
```yaml
max_iterations: 1000  # Increase iterations
convergence_tolerance: 1.0e-5  # Relax tolerance
```

## Differences from Fortran Version

- Uses Python's object-oriented design
- Modular architecture for easy extension
- Built-in visualization tools
- YAML configuration instead of fort.9
- Some numerical differences due to different implementations

## Contributing

To add new features:
1. Create new modules in appropriate directories
2. Follow existing code style
3. Add unit tests
4. Update documentation

## License

This is a research code. Please cite the original SEMITIP papers when using.