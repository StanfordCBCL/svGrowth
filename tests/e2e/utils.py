"""Shared utilities for E2E tests."""

import pandas as pd
import numpy as np
from typing import Optional, List


def compare_timeseries(
    output_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    variable_map: dict,
    relative_tolerance: float = 0.05,
    absolute_tolerance: float = 0.01,
    time_points: Optional[List[float]] = None
) -> dict:
    """Compare time series data.
    
    Args:
        output_df: Simulation output
        reference_df: Reference data
        variable_map: Dict mapping output var → reference var
        relative_tolerance: Relative error tolerance
        absolute_tolerance: Absolute error tolerance
        time_points: Specific times to compare (None = all)
        
    Returns:
        Dict with comparison results for each variable
    """
    results = {}
    
    for output_var, ref_var in variable_map.items():
        # Merge on time
        merged = pd.merge(
            output_df[['time', output_var]],
            reference_df[['time', ref_var]],
            on='time',
            how='inner'
        )
        
        # Filter time points
        if time_points:
            merged = merged[merged['time'].isin(time_points)]
        
        # Compute errors
        abs_error = np.abs(merged[output_var] - merged[ref_var])
        rel_error = abs_error / np.abs(merged[ref_var])
        
        results[output_var] = {
            'passed': np.all((abs_error <= absolute_tolerance) | 
                           (rel_error <= relative_tolerance)),
            'max_abs_error': abs_error.max(),
            'max_rel_error': rel_error.max(),
            'mean_abs_error': abs_error.mean(),
            'mean_rel_error': rel_error.mean()
        }
    
    return results