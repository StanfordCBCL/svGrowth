import os
import csv
import json
import yaml
import copy
from pathlib import Path
from typing import Dict, Any, List, Union 
from datetime import datetime

class IOHandler:
    """Handles file input/output operations."""
    
    def __init__(self, 
                 detail_level: int = 1, 
                 log_level: str = 'INFO', 
                 output_dir: str = None):
        """Initialize IO handler.
        
        Args:
            detail_level: 0=basic (summary+layer), 1=detailed (+constituents), 2=full
            log_level: Console logging level (DEBUG, INFO, WARNING, ERROR) 
            output_dir: Directory for output files
        """
        self.detail_level = detail_level
        self.log_level = log_level
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.output_writers = {}  # Track open file handles and writers
        
        # Create output directory if specified
        if self.output_dir is not None:
            self.output_dir.mkdir(parents=True, exist_ok=True)

    def set_output_dir(self, output_dir: str):
        """Set output directory after initialization.
        
        Args:
            output_dir: Path to output directory
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # INPUT: PARAMETER LOADING
    # =========================================================================
    
    def load_parameters(self, file_path: Union[str, Path]) -> Dict[str, Any]:
        """Load parameters from YAML file.
        
        Args:
            file_path: Path to YAML configuration file (can be relative or absolute)
            
        Returns:
            Dictionary of parameters loaded from YAML
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file is not .yaml/.yml or has invalid YAML syntax
        """
        # Resolve path to absolute and verify existence
        resolved_path = self._resolve_file_path(file_path)
        self._check_file_extension(resolved_path)
        
        try:
            with open(resolved_path, 'r') as file:
                params = yaml.safe_load(file)
            
            # TODO: Validate YAML structure if needed
    
            return params
            
        except yaml.YAMLError as e:
            raise ValueError(
                f"Error parsing YAML file: {resolved_path}\n"
                f"YAML error: {e}"
            )

    def apply_cli_overrides(self, params: Dict[str, Any], output_dir: str = None, log_level: str = None) -> Dict[str, Any]:
        """Apply CLI overrides to parameters.
        
        Args:
            params: Original parameters from YAML
            output_dir: CLI override for output directory
            log_level: CLI override for log level
            
        Returns:
            Updated parameters with CLI overrides applied
        """
        params_copy = copy.deepcopy(params)
        
        # Apply output directory override
        if output_dir is not None:
            params_copy['simulation']['output_directory'] = output_dir
        
        # Apply log level override
        if log_level is not None:
            params_copy['simulation']['log_level'] = log_level
        
        return params_copy

    # =========================================================================
    # OUTPUT: PARAMETER SAVING
    # =========================================================================
         
    def save_parameters(self, params: Dict[str, Any]) -> None:
        """Save simulation parameters with safe names to output directory.
        
        Currently not including any command-line overrides.
        
        Args:
            params: Original parameters from YAML
        """
        final_params = copy.deepcopy(params)

        # Add safe names for layer
        # TODO: can improve YAML readability by adding safe_name in the lines after layer_name in the YAML file
        if 'layers' in final_params:
            for layer_config in final_params['layers']:
                layer_config['safe_name'] = layer_config['layer_name'].replace(' ', '_').lower()
                                
        # Save to output directory
        params_file = self.output_dir / 'sim_params.yaml'
        with open(params_file, 'w') as f:
            yaml.dump(final_params, f, default_flow_style=False, sort_keys=False)

    # =========================================================================
    # OUTPUT: FILE SETUP
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
    # OUTPUT: DATA WRITING
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
    # OUTPUT: FINALIZATION
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

    # =========================================================================
    # UTILITIES: VALIDATION HELPERS
    # =========================================================================

    def _resolve_file_path(self, file_path: Union[str, Path]) -> Path:
        """Resolve file path to absolute path and verify it exists.
        
        Args:
            file_path: Path to file (relative or absolute, string or Path object)
            
        Returns:
            Absolute Path object
            
        Raises:
            FileNotFoundError: If file doesn't exist
        """
        path = Path(file_path)
        
        # Resolve to absolute path (handles .., ., symlinks)
        if not path.is_absolute():
            path = path.resolve()
        
        # Verify file exists
        if not path.exists():
            raise FileNotFoundError(
                f"File not found: {path}\n"
                f"Original path provided: {file_path}\n"
                f"Please check the path and try again."
            )
        
        return path
             
    def _check_file_extension(self, file_path: Path) -> None:
        """Check that file has .yaml or .yml extension.
        
        Args:
            file_path: Path object to check
            
        Raises:
            ValueError: If file extension is not .yaml or .yml
        """
        ext = file_path.suffix.lower()
        if ext not in ['.yaml', '.yml']:
            raise ValueError(
                f"Only YAML files (.yaml/.yml) are supported.\n"
                f"Got: {file_path.name} with extension '{ext}'"
            )

