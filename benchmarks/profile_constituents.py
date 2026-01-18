#!/usr/bin/env python3
"""
Profile constituent-level calculations in detail.

This script focuses on profiling the most computationally intensive
constituent operations: stress computation, mass density calculation,
and survival function evaluation.

Usage:
    python benchmarks/profile_constituents.py
"""

import sys
from pathlib import Path
from datetime import datetime

# Add src to Python path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

try:
    from pyinstrument import Profiler
except ImportError:
    print("❌ PyInstrument not installed!")
    print("   Install with: pip install pyinstrument")
    sys.exit(1)

from io_handler import IOHandler
from configuration import Configuration


def profile_constituent_operations():
    """Profile constituent stress and mass density calculations."""
    
    print("🔍 Profiling constituent-level operations")
    print("-" * 60)
    
    # Load configuration
    script_dir = Path(__file__).parent.parent / 'src'
    config_file = script_dir / "latorre2018.yaml"
    
    io_handler = IOHandler()
    params = io_handler.load_parameters(str(config_file))
    config = Configuration.from_parameters(params)
    
    # Get a constituent to profile
    layer = config.layers[0]
    constituent = layer.constituents[0]
    
    print(f"📦 Profiling constituent: {constituent.name}")
    print(f"📍 In layer: {layer.name}")
    print("-" * 60)
    
    # Create profiler
    profiler = Profiler()
    
    # Profile constituent operations
    profiler.start()
    
    # Simulate multiple timestep calculations
    print("\n⏱️  Running calculations for 100 timesteps...")
    for timestep in range(100):
        # These are typically the most expensive operations
        if hasattr(constituent, 'compute_sigma_alpha'):
            constituent.compute_sigma_alpha(current_timestep=timestep)
        
        if hasattr(constituent, 'compute_rhoR_alpha'):
            constituent.compute_rhoR_alpha(current_timestep=timestep)
    
    profiler.stop()
    
    print("✅ Calculations completed")
    print("-" * 60)
    
    # Print results
    print("\n📊 Profiling Results:\n")
    print(profiler.output_text(unicode=True, color=True))
    
    # Save HTML report
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    html_path = results_dir / f"profile_constituents_{timestamp}.html"
    
    with open(html_path, 'w') as f:
        f.write(profiler.output_html())
    
    print(f"\n💾 HTML report saved to: {html_path}")
    print(f"   Open with: open {html_path}")


if __name__ == "__main__":
    try:
        profile_constituent_operations()
    except Exception as e:
        print(f"\n❌ Profiling failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
