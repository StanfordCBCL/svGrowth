from io_handler import IOHandler
from configuration import Configuration
from simulation import Simulation

def main():
    """Main entry point for svGrowth simulation."""
    # Read input parameter file
    io_handler = IOHandler()
    params = io_handler.load_parameters("latorre2018_updated.yaml")

    # Create and initialize vessel configuration 
    config = Configuration.from_parameters(params)

    # Extract simulation parameters
    sim_params = params['simulation']
    
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
        debug_level=sim_params.get('debug_level', 0)    
    )
    
    sim.run()

if __name__ == "__main__":
    main()