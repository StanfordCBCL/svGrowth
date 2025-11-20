import os
from configuration import Configuration
from io_handler import IOHandler

class Simulation:
    def __init__(self, 
                 configuration: Configuration, 
                 io_handler=None, 
                 simulation_name=None, 
                 output_directory=None, 
                 dt=0.1, 
                 n_steps=100, 
                 tolerance=1e-6, 
                 max_iterations=50, 
                 integration_method='simpson',
                 survival_function_computation='backward',  
                 verbose=False):
        
        self.configuration = configuration
        self.dt = dt
        self.n_steps = n_steps
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.integration_method = integration_method
        self.survival_function_computation = survival_function_computation
        self.verbose = verbose
        
        self.current_timestep = 0  # Current timestep index (0 = initial/homeostatic state)
        
        # Simulation metadata (owned by Simulation)
        self.simulation_name = simulation_name or getattr(configuration, 'simulation_name', 'default_simulation')
        self.output_directory = output_directory or getattr(configuration, 'output_directory', 'output')
        
        # IO service (used by Simulation)
        self.io_handler = io_handler or IOHandler()
        
        # Output file handle (will be set up during run)
        self.output_file = None

    def run(self):
        """Run the complete G&R simulation."""        
        for step in range(self.n_steps):
            print(f"\n--- Simulation step {step+1} (timestep {self.current_timestep} → {self.current_timestep+1}) ---")

            self._initialize_timestep()
            while self._convergence_not_reached:
                self.configuration.compute_all_rhoR_alpha(
                    next_timestep, 
                    self.dt, 
                    self.integration_method,
                    self.survival_function_computation
                )
                self.configuration.compute_all_stress(next_timestep, self.dt, self.integration_method)
                self.configuration.compute_radius(next_timestep)
                self._check_convergence(self.configuration.get_rhoR)             


    def _initialize_timestep(self):
        next_timestep = self.current_timestep + 1
        print(f"Advancing from timestep {self.current_timestep} to {next_timestep}")
        
        # STEP 1: Initial guesses for next timestep using configured method
        self.configuration.guess_all_rhoR_alpha(next_timestep, guess_method="from_previous_timestep")
        # TODO: Make sure guess geometry includes active radius for active stress
        self.configuration.guess_geometry(next_timestep, guess_method="from_previous_timestep")
        # TODO: Not needed once we pre-allocate arrays with 0. Guess stress for mass production.
        self.configuration.guess_stress_and_wss(next_timestep)
        
        # STEP 2: Compute rhoR_alpha for next timestep
        self.configuration.compute_all_rhoR_alpha(
            next_timestep, 
            self.dt, 
            self.integration_method,
            self.survival_function_computation 
        )

        # STEP 3: Refine geometry based on incompressibility constraint
        #self.configuration.update_geometry_from_density(next_timestep, guess_variable="mid_radius")

        # STEP 4: Compute geometry for next timestep
        #self.configuration.compute_theoretical_stress(next_timestep)

        # STEP 5: Compute sigma for next timestep
        self.configuration.compute_all_stress(
            next_timestep,
            self.dt,
            self.integration_method,
            self.survival_function_computation  # Reuse same survival values!
        )

        # STEP 5: Compute radius for next timestep
        self.configuration.compute_radius(next_timestep)

    def _setup_output(self):
        """Setup output files at simulation start."""
        # Create output directory if it doesn't exist
        os.makedirs(self.output_directory, exist_ok=True)
        
        # Setup main results file
        output_filename = f"{self.simulation_name}_results.txt"
        output_path = os.path.join(self.output_directory, output_filename)
        
        self.output_file = self.io_handler.setup_output_file(output_path)
        
        # Write header
        headers = ["step", "time_days", "inner_radius_mm", "thickness_mm", "cauchy_stress_kPa", "total_mass"]
        self.io_handler.write_simulation_header(self.output_file, headers)
        
        print(f"Output will be written to: {output_path}")

    def _write_step_data(self, step):
        """Write comprehensive simulation data."""
        vessel_geometry = self.configuration.get_vessel_geometry()
        
        if vessel_geometry:
            inner_radius_mm = vessel_geometry['inner_radius'] * 1000
            thickness_mm = vessel_geometry['thickness'] * 1000
        else:
            inner_radius_mm, thickness_mm = 1.40, 0.12
        
        # Get mechanical state
        vessel_layer = None
        for layer in self.configuration.layers:
            if layer.layer_type == "biological":
                vessel_layer = layer
                break
        
        cauchy_stress_kPa = vessel_layer.cauchy_stress / 1000 if vessel_layer else 0.0
        total_mass = self.configuration.total_mass
        
        step_data = [
            step + 1,
            (step + 1) * self.dt,
            inner_radius_mm,
            thickness_mm,
            cauchy_stress_kPa,
            total_mass
        ]
        
        self.io_handler.write_simulation_step(self.output_file, step_data)

    def _cleanup_output(self):
        """Clean up output files at simulation end."""
        if self.output_file:
            self.output_file.close()
            print(f"Simulation complete. Results saved to {self.output_directory}")

    def save_final_configuration(self, filename=None):
        """Save the final configuration parameters to YAML file."""
        if filename is None:
            filename = f"{self.simulation_name}_final_config.yaml"
        
        output_path = os.path.join(self.output_directory, filename)
        
        # Get final geometry from vessel layer
        vessel_geometry = self.configuration.get_vessel_geometry()
        
        # Create a dictionary with final configuration state
        final_config = {
            'simulation_info': {
                'name': self.simulation_name,
                'total_steps': self.n_steps,
                'dt': self.dt,
                'final_time': self.n_steps * self.dt
            },
            'final_geometry': {
                'inner_radius_mm': vessel_geometry['inner_radius'] * 1000 if vessel_geometry else 1.40,
                'thickness_mm': vessel_geometry['thickness'] * 1000 if vessel_geometry else 0.12,
                'axial_stretch': vessel_geometry['axial_stretch'] if vessel_geometry else 1.0
            },
            'layers': []
        }
        
        # Add layer information
        for layer in self.configuration.layers:
            layer_info = {
                'name': layer.name,
                'type': layer.layer_type,
                'constituents': [c.name for c in layer.constituents]
            }
            final_config['layers'].append(layer_info)
        
        self.io_handler.save_parameters(final_config, output_path)
        print(f"Final configuration saved to: {output_path}")