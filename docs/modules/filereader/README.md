# FileReader Module Documentation

## Overview

The `filereader.py` module provides a robust YAML configuration reader for SEMITIP simulations. It handles loading, parsing, and validating configuration files with comprehensive error handling and logging support.

## Features

- **YAML Parsing**: Load and parse YAML configuration files
- **Schema Validation**: Automatic validation against the SEMITIP configuration schema
- **Error Handling**: Detailed error messages for debugging
- **Type Safety**: Ensures all configuration values have correct types
- **Backward Compatibility**: Supports legacy configuration formats
- **Logging**: Comprehensive logging for debugging

## Main Components

### YamlConfigReader Class

The main class for reading YAML configuration files.

```python
from src.core.filereader import YamlConfigReader

reader = YamlConfigReader()
config = reader.load_config('path/to/config.yaml')
```

#### Methods

- `load_config(file_path: Union[str, Path]) -> SemitipConfig`
  - Loads and validates a YAML configuration file
  - Returns a validated SemitipConfig object
  - Raises FileNotFoundError if file doesn't exist
  - Raises yaml.YAMLError for parsing errors
  - Raises ValueError for validation errors

### Convenience Functions

#### load_yaml_config

Quick function to load a configuration file:

```python
from src.core.filereader import load_yaml_config

config = load_yaml_config('path/to/config.yaml')
```

#### save_yaml_config

Save a configuration object back to YAML:

```python
from src.core.filereader import save_yaml_config

save_yaml_config(config, 'path/to/output.yaml')
```

## Usage Examples

### Basic Usage

```python
from src.core.filereader import YamlConfigReader

# Create reader instance
reader = YamlConfigReader()

# Load configuration
config = reader.load_config('data/input/MultInt_config.yaml')

# Access configuration values
print(f"Temperature: {config.temperature} K")
print(f"Simulation type: {config.simulation_type}")
```

### Error Handling

```python
from src.core.filereader import YamlConfigReader
import yaml

reader = YamlConfigReader()

try:
    config = reader.load_config('config.yaml')
except FileNotFoundError:
    print("Configuration file not found")
except yaml.YAMLError as e:
    print(f"YAML parsing error: {e}")
except ValueError as e:
    print(f"Validation error: {e}")
```

### Modifying and Saving Configuration

```python
from src.core.filereader import load_yaml_config, save_yaml_config

# Load configuration
config = load_yaml_config('input.yaml')

# Modify values
config.temperature = 350.0
config.tip.radius = 2.0

# Save modified configuration
save_yaml_config(config, 'output.yaml')
```

## Configuration Validation

The filereader automatically validates configurations during loading:

1. **Type Validation**: Ensures all fields have correct data types
2. **Range Validation**: Checks that numeric values are within valid ranges
3. **Required Fields**: Verifies all required fields are present
4. **Consistency Checks**: Ensures related fields are consistent

### Validation Examples

```python
# This will raise ValueError during loading:
# - temperature: -100  # Must be positive
# - tip.radius: 0      # Must be positive
# - grid.radial_points: -10  # Must be positive
```

## Supported Configuration Types

- **MultInt**: Multi-region integration simulations
- **MultPlane**: Multi-plane simulations

Each type has specific parameters and validation rules.

## File Format

Configuration files must be in YAML format with the following structure:

```yaml
version: "1.0"
simulation_type: "MultInt"  # or "MultPlane"
environment:
  temperature: 300.0
  dielectric_constant: 12.9
# ... additional configuration
```

## Error Messages

The module provides detailed error messages:

- **File not found**: Clear indication of missing file path
- **YAML syntax errors**: Line number and description of syntax issues
- **Validation errors**: Specific field and constraint violations

## Logging

Enable detailed logging for debugging:

```python
import logging
logging.basicConfig(level=logging.INFO)

# Now filereader will output detailed logs
reader = YamlConfigReader()
config = reader.load_config('config.yaml')
```

## Best Practices

1. **Always use absolute paths** when possible
2. **Handle exceptions** appropriately in production code
3. **Validate configurations** before running simulations
4. **Use logging** for debugging configuration issues
5. **Keep backups** of working configurations

## See Also

- [Config Schema Documentation](../config_schema/README.md) - Details on configuration structure
- [FileConverter Documentation](../fileconverter/README.md) - Converting between formats