from io_handler import IOHandler
from configuration import Configuration
from simulation import Simulation

def main():
    """Main entry point for pyGrowth simulation."""
    # Read input parameter file
    io_handler = IOHandler()
    params = io_handler.load_parameters("latorre2018_updated.yaml")

    # Create and initialize vessel configuration 
    config = Configuration.from_parameters(params)

    # Extract simulation control parameters
    sim_params = params['simulation']
    
    # Initialize and run simulation
    sim = Simulation(
        configuration=config,
        io_handler=io_handler,
        simulation_name=sim_params['simulation_name'],
        output_directory=sim_params['output_directory'],
        dt=sim_params['dt'],
        n_steps=int(sim_params['n_days'] / sim_params['dt']),
        tolerance=sim_params.get('tolerance', 1e-6),
        max_iterations=sim_params.get('max_iterations', 50),
        verbose=sim_params.get('verbose', False)
    )
    
    sim.run()
    sim.save_final_configuration()

if __name__ == "__main__":
    main()