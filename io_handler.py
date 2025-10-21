import yaml
import os

class IOHandler:
    """Handles file input/output operations."""
    
    def load_parameters(self, file_path):
        """Load parameters from YAML file."""
        self._check_file_extension(file_path)
        # TODO: Validate YAML structure if needed

        try:
            with open(file_path, 'r') as file:
                return yaml.safe_load(file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Parameter file not found: {file_path}")
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML file: {e}")
    
    def save_parameters(self, params, file_path):
        """Save parameters to YAML file."""
        try:
            with open(file_path, 'w') as file:
                yaml.safe_dump(params, file, default_flow_style=False)
        except IOError as e:
            raise IOError(f"Error writing to file {file_path}: {e}")
    
    def setup_output_file(self, file_path, mode='w'):
        """Setup output file for simulation results."""
        try:
            return open(file_path, mode)
        except IOError as e:
            raise IOError(f"Cannot open output file {file_path}: {e}")
    
    def write_simulation_header(self, file_handle, headers):
        """Write header line to simulation output file."""
        if isinstance(headers, list):
            file_handle.write('\t'.join(headers) + '\n')
        else:
            file_handle.write(str(headers) + '\n')
    
    def write_simulation_step(self, file_handle, data):
        """Write a single simulation step to output file."""
        if isinstance(data, list):
            file_handle.write('\t'.join(map(str, data)) + '\n')
        else:
            file_handle.write(str(data) + '\n')
    
    def _check_file_extension(self, file_path):
        """Check that file has .yaml or .yml extension."""
        _, ext = os.path.splitext(file_path)
        if ext.lower() not in ['.yaml', '.yml']:
            raise ValueError(f"Only YAML files (.yaml/.yml) are supported. Got: {ext}")