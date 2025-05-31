# Fort.9 to YAML Converter

This tool converts old fort.9 files from SEMITIP Fortran programs to the new YAML configuration format.

## Usage

```bash
python3 fileconverter.py <fort.9 file> [output.yaml]
```

### Examples

Convert MultInt fort.9 file:
```bash
python3 fileconverter.py Fortran/MultInt/fort_MultInt.9
# Output: converted_MultInt_config.yaml
```

Convert MultPlane fort.9 file:
```bash
python3 fileconverter.py Fortran/MultPlane/fort_MultPlane.9
# Output: converted_MultPlane_config.yaml
```

Specify custom output filename:
```bash
python3 fileconverter.py fort_MultInt.9 my_config.yaml
```

## Features

- Automatically detects simulation type (MultInt or MultPlane) from filename
- Preserves all parameter values from fort.9 format
- Generates properly structured YAML with appropriate nesting
- Handles comments and mixed format lines in fort.9 files
- Converts voltage arrays to start/end voltage format
- Maps old parameter structure to new hierarchical YAML structure

## Supported Conversions

### MultInt
- Basic tip parameters (slope, separation, radius, etc.)
- Multiple semiconductor regions with all material properties
- Surface state distributions (first and second)
- Grid parameters specific to MultInt
- Computation parameters (scaling steps, iterations, convergence)
- Voltage scan settings
- MultInt-specific parameters (parallel wavevectors, energy points, etc.)
- Output settings (contours, angles, etc.)

### MultPlane
- All MultInt parameters plus:
- Vacuum width and spacing
- Maximum energies for different bands
- Compute all bands option

## Output Format

The generated YAML files follow the new hierarchical structure:
- `environment`: temperature, dielectric constant
- `tip`: all tip-related parameters including position
- `semiconductor`: regions with material properties
- `surface`: surface state regions
- `grid`: computational grid settings
- `computation`: solver parameters
- `voltage_scan`: voltage sweep settings
- `multint_specific` or `multplane_specific`: simulation-specific parameters
- `output`: output control parameters

## Notes

- The converter automatically handles the different parameter orders between MultInt and MultPlane fort.9 files
- Numeric values are preserved exactly as in the original file
- Boolean flags (0/1) are converted to true/false in YAML
- The generated YAML includes Chinese comments for clarity