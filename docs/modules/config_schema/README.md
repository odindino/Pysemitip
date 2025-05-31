# Config Schema Module Documentation

## Overview

The `config_schema.py` module defines the complete configuration schema for SEMITIP simulations. It provides structured data classes with validation, type safety, and backward compatibility for both MultInt and MultPlane simulation types.

## Architecture

The module uses Python's `dataclasses` with custom validation to ensure configuration integrity. The schema is hierarchical, with specialized classes for different aspects of the simulation.

## Main Classes

### SemitipConfig

The root configuration class that contains all simulation parameters.

```python
@dataclass
class SemitipConfig:
    version: str
    simulation_type: str
    environment: Environment
    grid: Grid
    tip: Tip
    voltage_scan: VoltageScan
    semiconductor_regions: List[SemiconductorRegion]
    surface_regions: List[SurfaceRegion]
    multint_specific: Optional[MultIntSpecific] = None
    multplane_specific: Optional[MultPlaneSpecific] = None
    # ... additional fields
```

### Environment

Environmental parameters for the simulation:

```python
@dataclass
class Environment:
    temperature: float  # Kelvin
    dielectric_constant: float
```

### Grid

Grid parameters for numerical calculations:

```python
@dataclass
class Grid:
    radial_points: int
    vacuum_points: int
    angular_points: int
    energy_points: int
    radial_extent: float  # nm
    vacuum_extent: float  # nm
    angular_extent: float  # degrees
```

### Tip

STM tip configuration:

```python
@dataclass
class Tip:
    radius: float  # nm
    separation: float  # nm
    work_function: float  # eV
    position: TipPosition
    # Backward compatibility properties
    @property
    def x_position(self) -> float:
        return self.position.x
```

### SemiconductorRegion

Defines semiconductor material properties:

```python
@dataclass
class SemiconductorRegion:
    id: int
    donor_concentration: float  # cm^-3
    acceptor_concentration: float  # cm^-3
    band_gap: float  # eV
    affinity: float  # eV
    permittivity: float
    effective_mass: EffectiveMass
    work_function: float  # eV
    # ... additional properties
```

### SurfaceRegion

Surface state distributions:

```python
@dataclass
class SurfaceRegion:
    id: int
    position: Position
    first_distribution: SurfaceDistribution
    second_distribution: Optional[SurfaceDistribution] = None
```

## Validation

Each class includes validation methods to ensure data integrity:

```python
def validate(self) -> bool:
    """Validate configuration parameters."""
    errors = []
    
    if self.temperature <= 0:
        errors.append("溫度必須大於 0")
    
    if self.grid.radial_points <= 0:
        errors.append("徑向網格點數必須大於 0")
        
    # ... additional validations
    
    if errors:
        raise ValueError(f"配置驗證失敗:\\n" + "\\n".join(f"- {e}" for e in errors))
    
    return True
```

## Usage Examples

### Creating a Configuration

```python
from src.core.config_schema import (
    SemitipConfig, Environment, Grid, Tip, TipPosition,
    VoltageScan, SemiconductorRegion, EffectiveMass
)

# Create configuration components
env = Environment(temperature=300.0, dielectric_constant=12.9)
grid = Grid(
    radial_points=201,
    vacuum_points=251,
    angular_points=61,
    energy_points=40,
    radial_extent=60.0,
    vacuum_extent=30.0,
    angular_extent=180.0
)
tip = Tip(
    radius=1.0,
    separation=1.0,
    work_function=5.3,
    position=TipPosition(x=0.0, y=0.0)
)

# Create main configuration
config = SemitipConfig(
    version="1.0",
    simulation_type="MultInt",
    environment=env,
    grid=grid,
    tip=tip,
    # ... additional parameters
)
```

### Accessing Properties

```python
# Direct access
temp = config.temperature  # Via property
temp = config.environment.temperature  # Via nested structure

# Backward compatibility
x_pos = config.tip.x_position  # Legacy property
x_pos = config.tip.position.x  # New structure
```

### Validation

```python
try:
    config.validate()
    print("Configuration is valid")
except ValueError as e:
    print(f"Validation failed: {e}")
```

## Simulation-Specific Parameters

### MultInt Specific

For multi-region integration simulations:

```python
@dataclass
class MultIntSpecific:
    parallel_wavevectors: int
    energy_points: int
    expansion_factor: int
    integration_cutoff: float
```

### MultPlane Specific

For multi-plane simulations:

```python
@dataclass
class MultPlaneSpecific:
    vacuum_width: float
    vacuum_spacing: float
    max_energies: int
    periodic_copies: int
```

## Properties and Backward Compatibility

The schema includes many convenience properties for backward compatibility:

```python
# Top-level properties
config.temperature  # -> config.environment.temperature
config.dielectric_constant  # -> config.environment.dielectric_constant
config.tip_radius  # -> config.tip.radius
config.radial_grid_points  # -> config.grid.radial_points

# Nested properties
config.tip.x_position  # -> config.tip.position.x
config.tip.y_position  # -> config.tip.position.y
```

## Conversion Methods

### to_dict()

Convert configuration to dictionary:

```python
config_dict = config.to_dict()
# Save to YAML or JSON
```

### from_dict()

Create configuration from dictionary:

```python
config = SemitipConfig.from_dict(config_dict)
```

## Validation Rules

### Global Rules
- Version must be a valid version string
- Simulation type must be "MultInt" or "MultPlane"
- Temperature must be positive
- All physical dimensions must be positive

### Grid Rules
- All point counts must be positive integers
- Extents must be positive floats
- Angular extent typically 0-360 degrees

### Material Rules
- Band gap must be positive
- Concentrations must be non-negative
- Work functions typically 3-6 eV
- Effective masses must be positive

### Surface State Rules
- Densities must be non-negative
- Energy levels relative to band edges
- Distribution types must be valid

## Best Practices

1. **Always validate** configurations before use
2. **Use type hints** when creating configurations programmatically
3. **Prefer new property names** over backward-compatible ones
4. **Document custom validations** in your code
5. **Keep configurations versioned** for reproducibility

## Extending the Schema

To add new parameters:

1. Add field to appropriate dataclass
2. Add validation logic
3. Update to_dict/from_dict methods
4. Add backward compatibility if needed
5. Update documentation

## Error Handling

The module provides detailed validation errors:

```python
try:
    config.validate()
except ValueError as e:
    # Error message will list all validation failures
    # e.g., "配置驗證失敗:\n- 溫度必須大於 0\n- 探針半徑必須大於 0"
```

## See Also

- [FileReader Documentation](../filereader/README.md) - Loading configurations from YAML
- [FileConverter Documentation](../fileconverter/README.md) - Converting configurations
- [SEMITIP V6 Technical Manual](../../Fortran-semitip/SEMITIP%20V6,%20Technical%20Manual.pdf) - Physical parameter details