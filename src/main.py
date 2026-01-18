#!/usr/bin/env python3
"""
svGrowth - Growth and Remodeling Simulation Framework
"""

import sys
import argparse
from io_handler import IOHandler
from configuration import Configuration
from simulation import Simulation
from custom_logging import init_logging

def main(params_file: str, 
         log_level: str = None, 
         output_dir: str = None):
        
    # Read input parameter file
    io_handler = IOHandler()
    params = io_handler.load_parameters(params_file)
    params = io_handler.apply_cli_overrides(params, output_dir=output_dir, log_level=log_level)

    sim_params = params['simulation']    
    init_logging(sim_params['log_level'])

    # Save input parameter file to output directory
    io_handler.set_output_dir(sim_params['output_directory']) 
    io_handler.save_parameters(params)

    # Create and initialize vessel configuration 
    config = Configuration.from_parameters(params)

    # Initialize and run simulation
    sim = Simulation(
        configuration=config,
        io_handler=io_handler,
        simulation_name=sim_params['simulation_name'],
        output_directory=sim_params['output_directory'],
        dt=sim_params['dt'],
        n_steps=int(sim_params['n_days'] / sim_params['dt']),
        tolerance=float(sim_params.get('tolerance', 1e-12)),
        max_iterations=int(sim_params.get('max_iterations', 50)),
        integration_method=sim_params.get('integration_method', 'trapezoidal'),
        survival_function_computation=sim_params.get('survival_function_computation', 'backward'), 
        verbose=sim_params.get('verbose', False),
        detail_level=sim_params.get('detail_level', 1),
        log_level=sim_params['log_level']  
    )
    
    sim.run()

def parse_args():
    parser = argparse.ArgumentParser(
        description='svGrowth - Growth and Remodeling Simulation Framework',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Example:
        python main.py -i latorre2018.yaml -o output_dir
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        dest='params_file',  
        required=True,
        metavar='FILE',
        help='Input parameter file (YAML format)'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_dir',
        default=None,
        metavar='DIR',
        help='Output directory (overrides input YAML parameter file)'
    )
    
    parser.add_argument(
        '--log-level',
        dest='log_level',
        default=None,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Set logging level (overrides input YAML parameter file)'
    )

    parser.add_argument(
        '--version',
        action='version',
        version='svGrowth 0.1.0'
    )
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    main(
        params_file=args.params_file, 
        log_level=args.log_level,
        output_dir=args.output_dir
    )
    sys.exit(0)