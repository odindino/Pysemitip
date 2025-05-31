# Core Modules

This directory contains the core Python modules for the Pysemitip project.

## Modules

- **config_schema.py**: Defines the configuration schema and validation for SEMITIP simulations
- **fileconverter.py**: Converts between different file formats (see [detailed documentation](../../docs/modules/fileconverter/README.md))
- **filereader.py**: YAML configuration file reader with validation support

## Quick Start

```python
from src.core.filereader import YamlConfigReader
from src.core.config_schema import SemitipConfig

# Load a configuration file
reader = YamlConfigReader()
config = reader.load_config('path/to/config.yaml')
```

## Documentation

For detailed documentation on each module, see:
- [FileConverter Documentation](../../docs/modules/fileconverter/README.md)
- [FileReader Documentation](../../docs/modules/filereader/README.md) 
- [Config Schema Documentation](../../docs/modules/config_schema/README.md)