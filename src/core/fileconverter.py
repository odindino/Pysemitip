#!/usr/bin/env python3
"""
Fort.9 to YAML Converter for SEMITIP

This module converts old fort.9 files (fort_MultInt.9 or fort_MultPlane.9)
to the new YAML configuration format.

Author: Odindino
Date: 2025-05-31
"""

import sys
import yaml
from typing import Dict, List, Any, Optional
from pathlib import Path


class Fort9ToYamlConverter:
    """Converter for fort.9 files to YAML format"""
    
    def __init__(self):
        self.data = {}
        self.simulation_type = None
        
    def read_fort9(self, filename: str) -> Dict[str, Any]:
        """Read fort.9 file and parse its content"""
        lines = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                lines.append(line)
        
        # Determine simulation type from filename
        if 'MultInt' in filename:
            self.simulation_type = 'MultInt'
            return self._parse_multint(lines)
        elif 'MultPlane' in filename:
            self.simulation_type = 'MultPlane'
            return self._parse_multplane(lines)
        else:
            raise ValueError("Cannot determine simulation type from filename")
    
    def _parse_value(self, line: str) -> Any:
        """Parse a value from a line, handling comments"""
        # First, try to find a number at the beginning of the line
        parts = line.split()
        if not parts:
            return 0
        
        # Try each part until we find a valid number
        for part in parts:
            try:
                # Check if it's a float
                if '.' in part or 'e' in part.lower():
                    return float(part)
                else:
                    # Try integer
                    return int(part)
            except ValueError:
                # If first part is not a number, it might be after it
                continue
        
        # If no number found, return 0
        return 0
    
    def _parse_voltage_array(self, line: str) -> List[float]:
        """Parse voltage array from comma-separated values"""
        # Remove any comments
        for delimiter in ['(', '#', '!']:
            if delimiter in line:
                line = line.split(delimiter)[0]
        
        # Split by comma and convert to floats
        values = []
        for v in line.split(','):
            v = v.strip()
            if v:
                try:
                    values.append(float(v))
                except ValueError:
                    pass
        return values
    
    def _parse_multint(self, lines: List[str]) -> Dict[str, Any]:
        """Parse MultInt fort.9 format"""
        idx = 0
        
        # Basic parameters
        nparm = int(self._parse_value(lines[idx])); idx += 1
        
        # Tip parameters
        shank_slope = self._parse_value(lines[idx]); idx += 1
        separation = self._parse_value(lines[idx]); idx += 1
        radius = self._parse_value(lines[idx]); idx += 1
        protrusion_radius = self._parse_value(lines[idx]); idx += 1
        contact_potential = self._parse_value(lines[idx]); idx += 1
        x_position = self._parse_value(lines[idx]); idx += 1
        y_position = self._parse_value(lines[idx]); idx += 1
        
        # Semiconductor regions
        num_regions = int(self._parse_value(lines[idx])); idx += 1
        regions = []
        
        for i in range(num_regions):
            region = {
                'id': i + 1,
                'donor_concentration': self._parse_value(lines[idx]),
                'acceptor_concentration': self._parse_value(lines[idx + 1]),
                'band_gap': self._parse_value(lines[idx + 2]),
                'valence_band_offset': self._parse_value(lines[idx + 3]),
                'donor_binding_energy': self._parse_value(lines[idx + 4]),
                'acceptor_binding_energy': self._parse_value(lines[idx + 5]),
                'effective_mass': {
                    'conduction_band': self._parse_value(lines[idx + 6]),
                    'heavy_hole': self._parse_value(lines[idx + 7]),
                    'light_hole': self._parse_value(lines[idx + 8]),
                    'split_off_hole': self._parse_value(lines[idx + 9]),
                },
                'spin_orbit_splitting': self._parse_value(lines[idx + 10]),
                'degeneracy_indicator': self._parse_value(lines[idx + 11]),
                'inversion_indicator': self._parse_value(lines[idx + 12]),
            }
            idx += 13
            regions.append(region)
        
        # Environment
        dielectric_constant = self._parse_value(lines[idx]); idx += 1
        temperature = self._parse_value(lines[idx]); idx += 1
        
        # Surface regions
        num_surface_regions = int(self._parse_value(lines[idx])); idx += 1
        surface_regions = []
        
        for i in range(num_surface_regions):
            surface_region = {
                'id': i + 1,
                'first_distribution': {
                    'density': self._parse_value(lines[idx]),
                    'neutrality_level': self._parse_value(lines[idx + 1]),
                    'fwhm': self._parse_value(lines[idx + 2]),
                    'centroid_energy': self._parse_value(lines[idx + 3]),
                },
                'second_distribution': {
                    'density': self._parse_value(lines[idx + 4]),
                    'neutrality_level': self._parse_value(lines[idx + 5]),
                    'fwhm': self._parse_value(lines[idx + 6]),
                    'centroid_energy': self._parse_value(lines[idx + 7]),
                }
            }
            idx += 8
            surface_regions.append(surface_region)
        
        temp_dependence = self._parse_value(lines[idx]); idx += 1
        
        # Grid parameters
        mirror_plane = self._parse_value(lines[idx]); idx += 1
        radial_points = self._parse_value(lines[idx]); idx += 1
        vacuum_points = self._parse_value(lines[idx]); idx += 1
        semiconductor_points = self._parse_value(lines[idx]); idx += 1
        angular_points = self._parse_value(lines[idx]); idx += 1
        initial_grid_size = self._parse_value(lines[idx]); idx += 1
        
        # Computation parameters
        scaling_steps = int(self._parse_value(lines[idx])); idx += 1
        max_iterations = [self._parse_value(x) for x in lines[idx].split()]; idx += 1
        convergence_params = [self._parse_value(x) for x in lines[idx].split()]; idx += 1
        charge_density_table_size = self._parse_value(lines[idx]); idx += 1
        
        # Output parameters
        output_param = self._parse_value(lines[idx]); idx += 1
        
        # Voltage scan
        num_voltage_points = int(self._parse_value(lines[idx])); idx += 1
        voltage_array = self._parse_voltage_array(lines[idx]); idx += 1
        
        # Contour parameters
        num_contours = self._parse_value(lines[idx]); idx += 1
        contour_spacing = self._parse_value(lines[idx]); idx += 1
        contour_angle = self._parse_value(lines[idx]); idx += 1
        
        # Additional parameters
        electron_affinity = self._parse_value(lines[idx]); idx += 1
        fermi_energy = self._parse_value(lines[idx]); idx += 1
        
        # MultInt specific parameters
        parallel_wavevectors = self._parse_value(lines[idx]); idx += 1
        energy_points = self._parse_value(lines[idx]); idx += 1
        expansion_factor = self._parse_value(lines[idx]); idx += 1
        semiconductor_depth_fraction = self._parse_value(lines[idx]); idx += 1
        
        # Voltage parameters
        modulation_voltage = self._parse_value(lines[idx]); idx += 1
        negative_ramp = self._parse_value(lines[idx]); idx += 1
        positive_ramp = self._parse_value(lines[idx]); idx += 1
        start_voltage = self._parse_value(lines[idx]); idx += 1
        
        # Build YAML structure
        yaml_data = {
            'version': '1.0',
            'simulation_type': 'MultInt',
            'environment': {
                'temperature': temperature,
                'dielectric_constant': dielectric_constant
            },
            'tip': {
                'shank_slope': shank_slope,
                'separation': separation,
                'radius': radius,
                'protrusion_radius': protrusion_radius,
                'contact_potential': contact_potential,
                'position': {
                    'x': x_position,
                    'y': y_position
                },
                'fermi_energy': fermi_energy
            },
            'semiconductor': {
                'regions': regions,
                'electron_affinity': electron_affinity
            },
            'surface': {
                'regions': surface_regions,
                'temperature_dependence': bool(temp_dependence)
            },
            'grid': {
                'mirror_plane': bool(mirror_plane),
                'radial_points': radial_points,
                'vacuum_points': vacuum_points,
                'semiconductor_points': semiconductor_points,
                'angular_points': angular_points,
                'initial_grid_size': initial_grid_size
            },
            'computation': {
                'scaling_steps': scaling_steps,
                'max_iterations': max_iterations,
                'convergence_parameters': convergence_params,
                'charge_density_table_size': charge_density_table_size
            },
            'voltage_scan': {
                'points': num_voltage_points,
                'start_voltage': voltage_array[0] if voltage_array else start_voltage,
                'end_voltage': voltage_array[-1] if voltage_array else 2.0,
                'modulation_voltage': modulation_voltage,
                'negative_ramp': negative_ramp,
                'positive_ramp': positive_ramp
            },
            'multint_specific': {
                'parallel_wavevectors': parallel_wavevectors,
                'energy_points': energy_points,
                'expansion_factor': expansion_factor,
                'semiconductor_depth_fraction': semiconductor_depth_fraction
            },
            'output': {
                'basic_output': output_param == 1,
                'equipotential_contours': output_param == 2,
                'full_potential': output_param == 3,
                'num_contours': num_contours,
                'contour_spacing': contour_spacing,
                'contour_angle': contour_angle
            }
        }
        
        return yaml_data
    
    def _parse_multplane(self, lines: List[str]) -> Dict[str, Any]:
        """Parse MultPlane fort.9 format"""
        idx = 0
        
        # Basic parameters
        nparm = int(self._parse_value(lines[idx])); idx += 1
        
        # Tip parameters
        shank_slope = self._parse_value(lines[idx]); idx += 1
        separation = self._parse_value(lines[idx]); idx += 1
        radius = self._parse_value(lines[idx]); idx += 1
        protrusion_radius = self._parse_value(lines[idx]); idx += 1
        contact_potential = self._parse_value(lines[idx]); idx += 1
        x_position = self._parse_value(lines[idx]); idx += 1
        y_position = self._parse_value(lines[idx]); idx += 1
        
        # Semiconductor regions
        num_regions = int(self._parse_value(lines[idx])); idx += 1
        regions = []
        
        for i in range(num_regions):
            region = {
                'id': i + 1,
                'donor_concentration': self._parse_value(lines[idx]),
                'acceptor_concentration': self._parse_value(lines[idx + 1]),
                'band_gap': self._parse_value(lines[idx + 2]),
                'valence_band_offset': self._parse_value(lines[idx + 3]),
                'donor_binding_energy': self._parse_value(lines[idx + 4]),
                'acceptor_binding_energy': self._parse_value(lines[idx + 5]),
                'effective_mass': {
                    'conduction_band': self._parse_value(lines[idx + 6]),
                    'heavy_hole': self._parse_value(lines[idx + 7]),
                    'light_hole': self._parse_value(lines[idx + 8]),
                    'split_off_hole': self._parse_value(lines[idx + 9]),
                },
                'spin_orbit_splitting': self._parse_value(lines[idx + 10]),
                'degeneracy_indicator': self._parse_value(lines[idx + 11]),
                'inversion_indicator': self._parse_value(lines[idx + 12]),
            }
            idx += 13
            regions.append(region)
        
        # Environment
        dielectric_constant = self._parse_value(lines[idx]); idx += 1
        temperature = self._parse_value(lines[idx]); idx += 1
        
        # Surface regions
        num_surface_regions = int(self._parse_value(lines[idx])); idx += 1
        surface_regions = []
        
        for i in range(num_surface_regions):
            surface_region = {
                'id': i + 1,
                'first_distribution': {
                    'density': self._parse_value(lines[idx]),
                    'neutrality_level': self._parse_value(lines[idx + 1]),
                    'fwhm': self._parse_value(lines[idx + 2]),
                    'centroid_energy': self._parse_value(lines[idx + 3]),
                },
                'second_distribution': {
                    'density': self._parse_value(lines[idx + 4]),
                    'neutrality_level': self._parse_value(lines[idx + 5]),
                    'fwhm': self._parse_value(lines[idx + 6]),
                    'centroid_energy': self._parse_value(lines[idx + 7]),
                }
            }
            idx += 8
            surface_regions.append(surface_region)
        
        temp_dependence = self._parse_value(lines[idx]); idx += 1
        
        # Grid parameters
        mirror_plane = self._parse_value(lines[idx]); idx += 1
        radial_points = self._parse_value(lines[idx]); idx += 1
        vacuum_points = self._parse_value(lines[idx]); idx += 1
        semiconductor_points = self._parse_value(lines[idx]); idx += 1
        angular_points = self._parse_value(lines[idx]); idx += 1
        initial_grid_size = self._parse_value(lines[idx]); idx += 1
        
        # Computation parameters
        scaling_steps = int(self._parse_value(lines[idx])); idx += 1
        max_iterations = [self._parse_value(x) for x in lines[idx].split()]; idx += 1
        convergence_params = [self._parse_value(x) for x in lines[idx].split()]; idx += 1
        charge_density_table_size = self._parse_value(lines[idx]); idx += 1
        
        # Output parameters
        output_param = self._parse_value(lines[idx]); idx += 1
        
        # Voltage scan
        num_voltage_points = int(self._parse_value(lines[idx])); idx += 1
        voltage_array = self._parse_voltage_array(lines[idx]); idx += 1
        
        # Contour parameters
        num_contours = self._parse_value(lines[idx]); idx += 1
        contour_spacing = self._parse_value(lines[idx]); idx += 1
        contour_angle = self._parse_value(lines[idx]); idx += 1
        
        # Additional parameters
        electron_affinity = self._parse_value(lines[idx]); idx += 1
        fermi_energy = self._parse_value(lines[idx]); idx += 1
        
        # MultPlane specific parameters
        vacuum_width = self._parse_value(lines[idx]); idx += 1
        vacuum_spacing = self._parse_value(lines[idx]); idx += 1
        max_energies = [self._parse_value(x) for x in lines[idx].split()]; idx += 1
        compute_all_bands = self._parse_value(lines[idx]); idx += 1
        
        # Voltage parameters
        modulation_voltage = self._parse_value(lines[idx]); idx += 1
        negative_ramp = self._parse_value(lines[idx]); idx += 1
        positive_ramp = self._parse_value(lines[idx]); idx += 1
        start_voltage = self._parse_value(lines[idx]); idx += 1
        
        # Build YAML structure
        yaml_data = {
            'version': '1.0',
            'simulation_type': 'MultPlane',
            'environment': {
                'temperature': temperature,
                'dielectric_constant': dielectric_constant
            },
            'tip': {
                'shank_slope': shank_slope,
                'separation': separation,
                'radius': radius,
                'protrusion_radius': protrusion_radius,
                'contact_potential': contact_potential,
                'position': {
                    'x': x_position,
                    'y': y_position
                },
                'fermi_energy': fermi_energy
            },
            'semiconductor': {
                'regions': regions,
                'electron_affinity': electron_affinity
            },
            'surface': {
                'regions': surface_regions,
                'temperature_dependence': bool(temp_dependence)
            },
            'grid': {
                'mirror_plane': bool(mirror_plane),
                'radial_points': radial_points,
                'vacuum_points': vacuum_points,
                'semiconductor_points': semiconductor_points,
                'angular_points': angular_points,
                'initial_grid_size': initial_grid_size
            },
            'computation': {
                'scaling_steps': scaling_steps,
                'max_iterations': max_iterations,
                'convergence_parameters': convergence_params,
                'charge_density_table_size': charge_density_table_size
            },
            'voltage_scan': {
                'points': num_voltage_points,
                'start_voltage': voltage_array[0] if voltage_array else start_voltage,
                'end_voltage': voltage_array[-1] if voltage_array else 2.5,
                'modulation_voltage': modulation_voltage,
                'negative_ramp': negative_ramp,
                'positive_ramp': positive_ramp
            },
            'multplane_specific': {
                'vacuum_width': vacuum_width,
                'vacuum_spacing': vacuum_spacing,
                'max_energies': {
                    'light_hole': max_energies[0] if len(max_energies) > 0 else 6.0,
                    'heavy_hole': max_energies[1] if len(max_energies) > 1 else 1.8,
                    'split_off': max_energies[2] if len(max_energies) > 2 else 4.0,
                    'conduction_band': max_energies[3] if len(max_energies) > 3 else 7.5
                },
                'compute_all_bands': bool(compute_all_bands)
            },
            'output': {
                'basic_output': output_param == 1,
                'equipotential_contours': output_param == 2,
                'full_potential': output_param == 3,
                'num_contours': num_contours,
                'contour_spacing': contour_spacing,
                'contour_angle': contour_angle
            }
        }
        
        return yaml_data
    
    def save_yaml(self, data: Dict[str, Any], output_filename: str):
        """Save data to YAML file with proper formatting"""
        # Add header comments
        header = f"""# SEMITIP {data['simulation_type']} 模擬參數配置檔
# 自動從 fort.9 檔案轉換
# 版本: {data['version']}
# 模擬類型: {data['simulation_type']}
"""
        
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(header)
            f.write('\n')
            yaml.dump(data, f, 
                     default_flow_style=False, 
                     allow_unicode=True,
                     sort_keys=False,
                     width=1000)  # Prevent line wrapping
        
        print(f"Successfully converted to: {output_filename}")


def main():
    """Main function to handle command line conversion"""
    if len(sys.argv) < 2:
        print("Usage: python fileconverter.py <fort.9 file> [output.yaml]")
        print("Example: python fileconverter.py fort_MultInt.9")
        print("         python fileconverter.py fort_MultPlane.9 my_config.yaml")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Determine output filename
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Generate output filename based on input
        base_name = Path(input_file).stem
        if 'MultInt' in base_name:
            output_file = 'converted_MultInt_config.yaml'
        elif 'MultPlane' in base_name:
            output_file = 'converted_MultPlane_config.yaml'
        else:
            output_file = 'converted_config.yaml'
    
    # Create converter and process file
    converter = Fort9ToYamlConverter()
    
    try:
        yaml_data = converter.read_fort9(input_file)
        converter.save_yaml(yaml_data, output_file)
        
        print(f"\nConversion completed successfully!")
        print(f"Input: {input_file}")
        print(f"Output: {output_file}")
        print(f"Simulation type: {converter.simulation_type}")
        
    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()