"""
Postprocessing pipeline for svGrowth simulation results.

Reads simulation outputs and generates:
- Ratio CSV files (values normalized by homeostatic/initial state)
- Plots (future)
"""

import argparse
import sys
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional


class SimulationPostprocessor:
    """Postprocess svGrowth simulation results.
    
    Reads output files and generates:
    - Normalized/ratio CSV files
    - Statistical summaries
    """
    
    def __init__(
        self, 
        sim_results_dir: str, 
        params_filename: str = 'simulation_params.yaml',
        postproc_output_dir: str = 'postprocessing'
    ):
        """Initialize postprocessor.
        
        Args:
            sim_results_dir: Path to simulation results directory
            params_filename: Name of parameter file in results directory
            postproc_output_dir: Subdirectory name for postprocessing outputs
        """
        self.sim_results_dir = Path(sim_results_dir)
        self.params_filename = params_filename
        
        if not self.sim_results_dir.exists():
            raise FileNotFoundError(f"Results directory not found: {sim_results_dir}")
        
        # Create postprocessing output directory
        self.postproc_output_dir = self.sim_results_dir / postproc_output_dir
        self.postproc_output_dir.mkdir(exist_ok=True)
        
        # Load simulation parameters
        self.params = self._load_params()
        
        # Extract layer and constituent names
        self.layer_names = self._get_layer_names()
        self.constituent_names = self._get_constituent_names()
    
    def _load_params(self) -> Dict:
        """Load simulation parameters from results directory."""
        
        params_file = self.sim_results_dir / self.params_filename
        
        if not params_file.exists():
            raise FileNotFoundError(
                f"Parameter file not found: {params_file}\n"
                f"Expected '{self.params_filename}' in {self.sim_results_dir}"
            )
        
        with open(params_file, 'r') as f:
            return yaml.safe_load(f)
    
    def _get_layer_names(self) -> List[str]:
        """Extract layer names from parameters."""
        layers = []
        if 'layers' in self.params:
            for layer_config in self.params['layers']:
                layers.append(layer_config['safe_name'])
        return layers
    
    def _get_constituent_names(self) -> List[str]:
        """Extract constituent names from parameters."""
        constituents = []
        if 'layers' in self.params:
            for layer_config in self.params['layers']:
                if 'constituents' in layer_config:
                    for const_name, const_props in layer_config['constituents'].items():
                        # Handle multi-fiber families
                        if const_props.get('constituent_type') == 'multi_fiber_family':
                            for family_name in const_props.get('fiber_families', {}).keys():
                                full_name = f"{const_name}_{family_name}"
                                constituents.append(full_name)
                        else:
                            constituents.append(const_name)
        return constituents
    
    def _compute_ratios_from_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute ratios of all columns relative to initial state (t=0).
        
        Args:
            df: DataFrame with 'time' column and data columns
            
        Returns:
            DataFrame with time and ratio columns
        """
        # Get homeostatic values (t=0)
        homeostatic_rows = df[df['time'] == 0.0]
        
        if homeostatic_rows.empty:
            raise ValueError(
                "No homeostatic state (time=0) found in data. "
                "Cannot compute ratios without reference homeostatic state."
            )
        
        homeostatic = homeostatic_rows.iloc[0]

        # Columns to compute ratios for
        ratio_columns = [col for col in df.columns if col != 'time']
        
        ratios = pd.DataFrame()
        ratios['time'] = df['time']
        for col in ratio_columns:
            # Check for valid homeostatic value
            if pd.notna(homeostatic[col]) and homeostatic[col] != 0:
                ratios[f'{col}_ratio'] = df[col] / homeostatic[col]
            else:
                # Can't compute ratio (division by zero or NaN)
                ratios[f'{col}_ratio'] = np.nan
        
        return ratios
    
    def compute_layer_ratios(self, layer_name: str) -> Optional[pd.DataFrame]:
        """Compute ratios of layer quantities relative to initial state (t=0)."""
        
        layer_file = self.sim_results_dir / f"layer_{layer_name}.csv"
        if not layer_file.exists():
            return None
        
        df = pd.read_csv(layer_file)
        return self._compute_ratios_from_dataframe(df)
    
    def compute_constituent_ratios(self, constituent_name: str) -> Optional[pd.DataFrame]:
        """Compute ratios of constituent quantities relative to initial state (t=0)."""
        
        const_file = self.sim_results_dir / f"constituent_{constituent_name}.csv"
        if not const_file.exists():
            return None
        
        df = pd.read_csv(const_file)
        return self._compute_ratios_from_dataframe(df)
    
    def save_ratio_csvs(self):
        """Generate and save all ratio CSV files."""
        # Process layers
        for layer_name in self.layer_names:
            ratios = self.compute_layer_ratios(layer_name)
            if ratios is not None:
                output_file = self.postproc_output_dir / f"layer_{layer_name}_ratios.csv"
                ratios.to_csv(output_file, index=False)
        
        # Process constituents
        for const_name in self.constituent_names:
            ratios = self.compute_constituent_ratios(const_name)
            if ratios is not None:
                output_file = self.postproc_output_dir / f"constituent_{const_name}_ratios.csv"
                ratios.to_csv(output_file, index=False)
    
    def run_all(self):
        """Run complete postprocessing pipeline."""
        self.save_ratio_csvs()


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='svGrowth Postprocessing - Analyze simulation results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
            python postprocessing.py outputdir
        """
    )
    
    parser.add_argument(
        'sim_results_dir',
        metavar='RESULTS_DIR',
        help='Path to simulation results directory'
    )
    
    parser.add_argument(
        '-p', '--params',
        dest='params_filename',
        default='simulation_params.yaml',
        metavar='FILE',
        help='Parameter file name in results directory (default: simulation_params.yaml)'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='postproc_output_dir',
        default='postprocessing',
        metavar='DIR',
        help='Subdirectory for postprocessing outputs (default: postprocessing)'
    )
    
    return parser.parse_args()


def main():
    """Main entry point for postprocessing."""
    args = parse_args()
    
    postprocessor = SimulationPostprocessor(
        sim_results_dir=args.sim_results_dir,
        params_filename=args.params_filename,
        postproc_output_dir=args.postproc_output_dir
    )
    postprocessor.run_all()

    return 0


if __name__ == "__main__":
    sys.exit(main())