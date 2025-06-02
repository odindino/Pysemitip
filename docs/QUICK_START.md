# Pysemitip Quick Start Guide

## Installation

1. **Install Python 3.8+** and required packages:
```bash
pip install numpy scipy matplotlib pyyaml
```

2. **Set up directories**:
```bash
python setup_directories.py
```

3. **Run tests** to verify installation:
```bash
python test_all.py
```

## Directory Structure

```
Pysemitip/
├── data/
│   ├── input/
│   │   ├── configs/      # Your configuration files
│   │   └── examples/     # Example configurations
│   └── output/
│       ├── results/      # Simulation results (auto-organized by date/time)
│       ├── logs/         # Log files
│       ├── figures/      # Generated plots
│       └── contours/     # Contour plots
├── src/                  # Source code
├── tests/               # Test files
└── docs/                # Documentation
```

## Running Simulations

### Quick Test
```bash
python run_multint.py data/input/examples/quick_test.yaml
```

### Full Simulation with Plots
```bash
python run_multint.py data/input/examples/test_multint.yaml --plot
```

### Custom Configuration
```bash
python run_multint.py data/input/configs/my_config.yaml --name my_simulation --plot
```

## Output Files

Results are saved in `data/output/results/[simulation_name]_[timestamp]/`:
- `multint_results.pkl` - Main results (Python pickle format)
- `output_MultInt.log` - Detailed calculation log
- `figures/` - If --plot is used
- `contours/` - If --plot is used

## Configuration Files

Configuration files use YAML format. Key sections:

```yaml
# Material properties
semiconductor_regions:
  - id: 1
    donor_concentration: 1.0e18  # cm^-3
    band_gap: 1.42               # eV

# Surface states
surface_regions:
  - id: 1
    first_distribution:
      density: 4.4e14            # cm^-2 eV^-1
      neutrality_level: 0.125    # eV

# Tip parameters
tip:
  radius: 1.0                    # nm
  separation: 1.0                # nm

# Voltage scan
voltage_scan:
  start: -2.0                    # V
  end: 2.0                       # V
  points: 21                     # Number of bias points
```

## Analyzing Results

### Load Results in Python
```python
import pickle
import matplotlib.pyplot as plt

# Load results
with open('data/output/results/my_simulation_20241206_120000/multint_results.pkl', 'rb') as f:
    data = pickle.load(f)

results = data['results']

# Plot I-V curve
biases = [r.bias_voltage for r in results]
currents = [r.current for r in results]

plt.figure()
plt.plot(biases, currents, 'o-')
plt.xlabel('Bias (V)')
plt.ylabel('Current (A)')
plt.show()
```

### Use Built-in Plotting
```python
from src.physics.visualization import plot_simulation_results

plot_simulation_results('path/to/multint_results.pkl')
```

## Tips for Fast Testing

1. Use `quick_test.yaml` for debugging
2. Reduce grid points for faster computation:
   ```yaml
   grid:
     radial_points: 16
     vacuum_points: 16
   ```
3. Use fewer bias points:
   ```yaml
   voltage_scan:
     points: 5
   ```
4. Disable degeneracy for speed:
   ```yaml
   allow_degeneracy: false
   ```

## Troubleshooting

### Import Errors
Run from the project root directory:
```bash
cd /path/to/Pysemitip
python run_multint.py [config_file]
```

### Convergence Issues
- Increase `max_iterations` in config
- Relax `convergence_tolerance`
- Reduce grid size

### Memory Issues
- Reduce grid points
- Use fewer energy points in `charge_table_points`

## Common Tasks

### Compare with Fortran Results
Place Fortran output in `tests/data/reference/` and use comparison scripts.

### Batch Processing
```bash
for config in data/input/configs/*.yaml; do
    python run_multint.py "$config" --plot
done
```

### Custom Analysis
See `docs/examples/` for Jupyter notebooks with analysis examples.

## Getting Help

1. Check error messages in log files
2. Run with `--verbose` for detailed output
3. See `PYTHON_MULTINT_README.md` for detailed documentation
4. Check unit tests for usage examples