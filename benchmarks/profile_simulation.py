#!/usr/bin/env python3
"""
Profile a complete svGrowth simulation using PyInstrument.

Usage:
    python benchmarks/profile_simulation.py -i latorre2018.yaml
    python benchmarks/profile_simulation.py --input latorre2018.yaml --browser
    python benchmarks/profile_simulation.py --help

Output:
    - Console: Text summary of profiling results
    - HTML (optional): Interactive flame graph in benchmarks/results/
"""

import sys
import argparse
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

from main import main


def profile_simulation(params_file: str = "latorre2018.yaml", 
                       output_html: bool = True,
                       output_text: bool = True,
                       open_browser: bool = False):
    """Profile a complete simulation run.
    
    Args:
        params_file: Path to parameter file (YAML format)
        output_html: Save interactive HTML report
        output_text: Save text report
        open_browser: Open HTML report in browser automatically
    
    Returns:
        Path to HTML report if generated, None otherwise
    """
    print(f"🔍 Profiling simulation")
    print(f"📁 Parameters: {params_file}")
    print(f"⏱️  Starting profiling at {datetime.now().strftime('%H:%M:%S')}")
    print("-" * 60)
    
    # Create profiler
    profiler = Profiler()
    
    # Run simulation with profiling
    profiler.start()
    try:
        main(params_file)  # Run your simulation
        profiler.stop()
        print(f"✅ Simulation completed at {datetime.now().strftime('%H:%M:%S')}")
    except Exception as e:
        profiler.stop()
        print(f"❌ Simulation failed: {e}")
        raise
    
    print("-" * 60)
    
    # Create results directory
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    # Generate timestamp for filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    params_name = Path(params_file).stem
    
    # Print to console
    print("\n📊 Profiling Results:\n")
    print(profiler.output_text(unicode=True, color=True))
    
    # Save HTML report
    html_path = None
    if output_html:
        html_path = results_dir / f"profile_{params_name}_{timestamp}.html"
        with open(html_path, 'w') as f:
            f.write(profiler.output_html())
        print(f"\n💾 HTML report saved to: {html_path}")
        print(f"   Open with: open {html_path}")
    
    # Save text report
    if output_text:
        text_path = results_dir / f"profile_{params_name}_{timestamp}.txt"
        with open(text_path, 'w') as f:
            f.write(profiler.output_text(unicode=False, color=False))
        print(f"💾 Text report saved to: {text_path}")
    
    # Open in browser
    if open_browser and html_path:
        import webbrowser
        webbrowser.open(f'file://{html_path.absolute()}')
        print(f"🌐 Opened in browser")
    
    return html_path


def main_cli():
    """Command-line interface for profiling."""
    parser = argparse.ArgumentParser(
        description="Profile svGrowth simulation performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Profile with default parameter file
    python benchmarks/profile_simulation.py -i latorre2018.yaml
    
    # Profile and open in browser
    python benchmarks/profile_simulation.py -i latorre2018.yaml --browser
    
    # Profile without saving HTML (console only)
    python benchmarks/profile_simulation.py -i latorre2018.yaml --no-html
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        dest='params_file',
        default='latorre2018.yaml',
        metavar='FILE',
        help='Parameter file to profile (default: latorre2018.yaml)'
    )
    
    parser.add_argument(
        '--browser', '-b',
        action='store_true',
        help='Open HTML report in browser after profiling'
    )
    
    parser.add_argument(
        '--no-html',
        action='store_true',
        help='Do not generate HTML report (console output only)'
    )
    
    parser.add_argument(
        '--no-text',
        action='store_true',
        help='Do not save text report'
    )
    
    args = parser.parse_args()
    
    try:
        profile_simulation(
            params_file=args.params_file,
            output_html=not args.no_html,
            output_text=not args.no_text,
            open_browser=args.browser
        )
    except Exception as e:
        print(f"\n❌ Profiling failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main_cli()
