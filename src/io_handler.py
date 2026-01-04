import os
import csv
import json
import yaml
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

class IOHandler:
    """Handles file input/output operations."""
    
    def __init__(self, detail_level: int = 0, debug_level: int = 0, output_dir: str = 'output'):
        """Initialize IO handler.
        
        Args:
            detail_level: 0=basic (summary+layer), 1=detailed (+constituents), 2=full
            debug_level: 0=none, 1=logs, 2=debug files
            output_dir: Directory for output files
        """
        self.detail_level = detail_level
        self.debug_level = debug_level
        self.output_dir = Path(output_dir)
        self.output_writers = {}  # Track open file handles and writers
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # PARAMETER LOADING (Input)
    # =========================================================================
    
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

    
    # =========================================================================
    # SETUP OUTPUT FILES
    # =========================================================================
    
    def setup_simulation_summary(self):
        """Setup summary CSV file (always created)."""
        filepath = self.output_dir / 'simulation_summary.csv'
        file_handle = open(filepath, 'w', newline='')
        
        headers = ['time', 'layer_name', 'a', 'h', 'stress_theta', 'stress_z', 'rho', 'pressure', 'flow_rate', 'wss']
        writer = csv.DictWriter(file_handle, fieldnames=headers)
        writer.writeheader()
        
        self.output_writers['summary'] = {
            'handle': file_handle,
            'writer': writer,
            'headers': headers
        }
    
    def setup_layer_output(self, layer_name: str):
        """Setup per-layer CSV file (always created)."""
        safe_name = layer_name.replace(' ', '_').lower()
        filepath = self.output_dir / f'layer_{safe_name}.csv'
        file_handle = open(filepath, 'w', newline='')
        
        # TODO: create registry for outputs
        headers = ['time', 'a', 'h', 'lambda_z', 'rho', 'J', 
                'F_radial', 'F_circ', 'F_axial',
                'stress_radial', 'stress_circ', 'stress_axial', 
                'stress', 'pressure', 'flow_rate', 'wss']
        writer = csv.DictWriter(file_handle, fieldnames=headers)
        writer.writeheader()
        
        self.output_writers[f'layer_{safe_name}'] = {
            'handle': file_handle,
            'writer': writer,
            'headers': headers
        }
    
    def setup_constituent_output(self, constituent_name: str, layer_name: str = ''):
        """Setup per-constituent CSV file (only if detail_level >= 1)."""
        if self.detail_level < 1:
            return
        
        safe_name = constituent_name.replace(' ', '_').lower()
        filepath = self.output_dir / f'constituent_{safe_name}.csv'
        file_handle = open(filepath, 'w', newline='')
        
        headers = ['time', 'rho_alpha', 'k_alpha', 'mR_alpha', 'stimulus_function', 'stress_r', 'stress_theta', 'stress_z', 
                'sigma_hat_r', 'sigma_hat_theta', 'sigma_hat_z', 'sigma_hat_act', 'a_act']
        writer = csv.DictWriter(file_handle, fieldnames=headers)
        writer.writeheader()
        
        self.output_writers[f'constituent_{safe_name}'] = {
            'handle': file_handle,
            'writer': writer,
            'headers': headers
        }
    
    def setup_metadata_output(self):
        """Setup metadata JSON file."""
        self.output_writers['metadata'] = {
            'filepath': self.output_dir / 'metadata.json',
            'data': {
                'simulation_info': {
                    'start_time': datetime.now().isoformat(),
                    'software': 'svGrowth'
                },
                'configuration': {},
                'statistics': {}
            }
        }

    # =========================================================================
    # WRITE DATA
    # =========================================================================
    
    def write_summary_row(self, data: Dict[str, Any]):
        """Write one row to summary CSV."""
        if 'summary' not in self.output_writers:
            raise RuntimeError("Summary output not setup")
        
        writer = self.output_writers['summary']['writer']
        writer.writerow(data)
        self.output_writers['summary']['handle'].flush()
    
    def write_layer_row(self, layer_name: str, data: Dict[str, Any]):
        """Write one row to layer CSV."""
        safe_name = layer_name.replace(' ', '_').lower()
        key = f'layer_{safe_name}'
        
        if key not in self.output_writers:
            raise RuntimeError(f"Layer output for '{layer_name}' not setup")
        
        writer = self.output_writers[key]['writer']
        writer.writerow(data)
        self.output_writers[key]['handle'].flush()
    
    def write_constituent_row(self, constituent_name: str, data: Dict[str, Any], layer_name: str = ''):
        """Write one row to constituent CSV."""
        if self.detail_level < 1:
            return
        
        safe_name = constituent_name.replace(' ', '_').lower()
        key = f'constituent_{safe_name}'
        
        if key not in self.output_writers:
            raise RuntimeError(f"Constituent output for '{constituent_name}' not setup")
        
        writer = self.output_writers[key]['writer']
        writer.writerow(data)
        self.output_writers[key]['handle'].flush()
    
    def update_metadata(self, section: str, data: Dict[str, Any]):
        """Update metadata section."""
        if 'metadata' not in self.output_writers:
            return
        
        self.output_writers['metadata']['data'][section].update(data)
    
    def write_metadata(self):
        """Write metadata to file."""
        if 'metadata' not in self.output_writers:
            return
        
        metadata = self.output_writers['metadata']
        metadata['data']['simulation_info']['end_time'] = datetime.now().isoformat()
        
        with open(metadata['filepath'], 'w') as f:
            json.dump(metadata['data'], f, indent=2)

    # =========================================================================
    # CLEANUP
    # =========================================================================
    
    def close_all(self):
        """Close all open file handles."""
        # Write final metadata
        if 'metadata' in self.output_writers:
            self.write_metadata()
        
        # Close all file handles
        for key, writer_info in self.output_writers.items():
            if 'handle' in writer_info:
                writer_info['handle'].close()
        
        self.output_writers.clear()


    # Other methods...
    def save_parameters(self, params, file_path):
        """Save parameters to YAML file."""
        try:
            with open(file_path, 'w') as file:
                yaml.safe_dump(params, file, default_flow_style=False)
        except IOError as e:
            raise IOError(f"Error writing to file {file_path}: {e}")
    
    def _check_file_extension(self, file_path):
        """Check that file has .yaml or .yml extension."""
        _, ext = os.path.splitext(file_path)
        if ext.lower() not in ['.yaml', '.yml']:
            raise ValueError(f"Only YAML files (.yaml/.yml) are supported. Got: {ext}")