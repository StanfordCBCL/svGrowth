"""
E2E test for Latorre 2018 paper validation.

Reference: Latorre & Humphrey (2018). "Modeling biological growth and remodeling..."
Journal of Biomechanics, DOI: 10.1016/j.jbiomech.2018.02.017
"""

import pytest
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import shutil
import subprocess
import sys

# Get path to src directory
SRC_DIR = Path(__file__).parent.parent.parent / 'src'

class TestLatorre2018:
    """E2E validation against Latorre 2018 paper."""
    
    @pytest.fixture
    def fixture_dir(self):
        """Get fixture directory for this test."""
        return Path(__file__).parent / 'fixtures' / 'latorre2018'
    
    @pytest.fixture
    def comparison_config(self, fixture_dir):
        """Load comparison configuration."""
        with open(fixture_dir / 'comparison_config.yaml', 'r') as f:
            return yaml.safe_load(f)
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create temporary output directory."""
        temp_dir = tempfile.mkdtemp(prefix='e2e_test_latorre2018')
        yield temp_dir
        shutil.rmtree(temp_dir)  # Cleanup after test
    
    def test_latorre2018_hypertension(self, fixture_dir, comparison_config, temp_output_dir):
        """Run simulation and compare to Latorre 2018 reference."""
        
        print(f"\n{'='*70}")
        print("E2E Test: Latorre 2018")
        print(f"{'='*70}")
        print(f"Fixture dir: {fixture_dir}")
        print(f"Output dir:  {temp_output_dir}")
        
        # Step 1: Run simulation using CLI
        input_yaml = fixture_dir / 'sim_params.yaml'
        
        cmd = [
            sys.executable,
            str(SRC_DIR / 'main.py'),
            '-i', str(input_yaml),
            '-o', str(temp_output_dir),
        ]
        
        print(f"\nRunning simulation...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"\nSimulation FAILED:")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            pytest.fail(f"Simulation failed with return code {result.returncode}")
        
        print(f"✓ Simulation completed successfully")
        
        # Step 2: Run postprocessing
        print(f"\nRunning postprocessing...")
        
        sys.path.insert(0, str(SRC_DIR))
        from postprocessing import SimulationPostprocessor
        
        postprocessor = SimulationPostprocessor(
            sim_results_dir=temp_output_dir,
            params_filename='simulation_params.yaml',
            postproc_output_dir='postprocessing'
        )
        postprocessor.run_all()
        
        print(f"✓ Postprocessing completed")
        
        # Step 3: Compare results
        print(f"\n{'='*70}")
        print("Comparing results to reference data")
        print(f"{'='*70}")
        
        postproc_dir = Path(temp_output_dir) / 'postprocessing'
        
        # Get time filtering options from config
        compare_time_points = comparison_config.get('comparen_time_points')
        compare_time_range = comparison_config.get('compare_time_range')
        
        # Print what we're comparing
        if compare_time_points:
            print(f"Comparing at specific time points: {compare_time_points}")
        elif compare_time_range:
            print(f"Comparing time range: {compare_time_range['start']} to {compare_time_range['end']} (step {compare_time_range['step']})")
        else:
            print("Comparing all overlapping time points")
        
        for comparison in comparison_config['comparisons']:
            output_file = postproc_dir / comparison['output_file']
            reference_file = fixture_dir / comparison['reference_file']
            
            if not output_file.exists():
                pytest.fail(f"Output file not found: {output_file}")
            
            if not reference_file.exists():
                pytest.fail(f"Reference file not found: {reference_file}")
            
            # Load data
            output_df = pd.read_csv(output_file)
            reference_df = pd.read_csv(reference_file)
            
            print(f"\nComparing: {comparison['output_file']}")
            
            for var_config in comparison['variables']:
                var_name = var_config['name']
                ref_name = var_config.get('reference_name', var_name)
                rel_tol = var_config.get('tolerance', {}).get('relative', 0.05)
                
                self._compare_variable(
                    output_df, 
                    reference_df, 
                    var_name, 
                    ref_name,
                    rel_tol,
                    comparison_config.get('comparison_time_points'),
                    comparison_config.get('comparison_time_range')
                )
        
        print(f"\n{'='*70}")
        print("✅ PASSED E2E Test: Latorre 2018")
        print(f"{'='*70}\n")
    
    def _compare_variable(
        self, 
        output_df: pd.DataFrame, 
        reference_df: pd.DataFrame,
        output_var: str,
        reference_var: str,
        rel_tol: float,
        compare_time_points: list = None,
        compare_time_range: dict = None
    ):
        """Compare a single variable between output and reference.
        
        Args:
            output_df: Output DataFrame
            reference_df: Reference DataFrame
            output_var: Variable name in output
            reference_var: Variable name in reference
            rel_tol: Relative tolerance
            comparison_time_points: Specific time points to compare (list)
            comparison_time_range: Time range to compare (dict with 'start', 'end', 'step')
        """
        
        # Normalize time column types (avoid int/float merge warning)
        output_df = output_df.copy()
        reference_df = reference_df.copy()
        output_df['time'] = output_df['time'].astype(float)
        reference_df['time'] = reference_df['time'].astype(float)
        
        # Merge dataframes on time column
        merged = pd.merge(
            output_df[['time', output_var]],
            reference_df[['time', reference_var]],
            on='time',
            how='inner'
        )
        
        # Filter time points based on what's provided
        if compare_time_points is not None:
            # Use specific time points
            merged = merged[merged['time'].isin(compare_time_points)]
        elif compare_time_range is not None:
            # Generate time points from range
            start = compare_time_range.get('start', 0)
            end = compare_time_range.get('end', merged['time'].max())
            step = compare_time_range.get('step', 1)
            
            # Create array of times to compare
            comparison_times = np.arange(start, end + step, step)
            merged = merged[merged['time'].isin(comparison_times)]
    
        # Check that we have overlapping data
        if merged.empty:
            pytest.fail(
                f"No overlapping time points for {output_var}. "
                f"Output times: {output_df['time'].values[:5]}..., "
                f"Reference times: {reference_df['time'].values[:5]}..."
            )
        
        # Extract values
        output_vals = merged[output_var].values
        reference_vals = merged[reference_var].values
        
        # Compute relative error (handle division by zero)
        rel_error = np.where(
            np.abs(reference_vals) > 1e-10,
            np.abs(output_vals - reference_vals) / np.abs(reference_vals),
            0.0
        )
        
        # Check if all errors are within tolerance
        passed = rel_error <= rel_tol
        
        if not np.all(passed):
            failed_indices = np.where(~passed)[0]
            failure_info = []
            
            for idx in failed_indices[:5]:  # Show first 5 failures
                t = merged.iloc[idx]['time']
                out_val = output_vals[idx]
                ref_val = reference_vals[idx]
                rel_err = rel_error[idx]
                
                failure_info.append(
                    f"    t={t:.1f}: output={out_val:.6f}, "
                    f"reference={ref_val:.6f}, "
                    f"rel_err={rel_err:.4f} (tolerance={rel_tol})"
                )
            
            pytest.fail(
                f"\n  {output_var} comparison FAILED at {len(failed_indices)}/{len(merged)} time points:\n" +
                "\n".join(failure_info)
            )
        
        # Print success summary
        max_rel_err = rel_error.max()
        mean_rel_err = rel_error.mean()
        n_points = len(merged)
        print(f"  ✓ {output_var}: {n_points} points, max_rel_err={max_rel_err:.4f}, mean_rel_err={mean_rel_err:.4f} (tol={rel_tol})")