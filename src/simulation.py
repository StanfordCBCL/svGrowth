import os
from typing import Dict
from configuration import Configuration
from io_handler import IOHandler
from custom_logging import get_logger
from fixed_point_solver import FixedPointSolver

class Simulation:
    def __init__(self, 
                 configuration: Configuration, 
                 io_handler=None, 
                 simulation_name=None, 
                 output_directory=None, 
                 dt=None, 
                 n_steps=None,
                 
                 # Fixed-point solver parameters
                 fixed_point_tolerance=1e-6,
                 fixed_point_max_iterations=50,
                 
                 # Equilibrium solver parameters
                 equilibrium_solver_method='brentq',
                 equilibrium_tolerance=1e-5,
                 equilibrium_verbose=False,
                 
                 # Numerical methods
                 integration_method='trapezoidal',
                 survival_function_computation='backward',
                 
                 # Verbosity
                 verbose=False,
                 detail_level=None,
                 log_level='INFO'):

        self.logger = get_logger(__name__)

        self.configuration = configuration
        self.dt = dt
        self.n_steps = n_steps
        
        # Fixed-point solver parameters
        self.fixed_point_tolerance = float(fixed_point_tolerance)
        self.fixed_point_max_iterations = int(fixed_point_max_iterations)
        
        # Equilibrium solver parameters
        self.equilibrium_solver_method = equilibrium_solver_method
        self.equilibrium_tolerance = float(equilibrium_tolerance)
        self.equilibrium_verbose = equilibrium_verbose
        
        # Numerical methods
        self.integration_method = integration_method
        self.survival_function_computation = survival_function_computation
        
        # Verbosity
        self.verbose = verbose
        self.log_level = log_level
        self.current_timestep = 0
        
        # Simulation metadata
        self.simulation_name = simulation_name 
        self.output_directory = output_directory 

        # Fixed-point solver
        self.fixed_point_solver = FixedPointSolver(
            tolerance=self.fixed_point_tolerance,
            max_iterations=self.fixed_point_max_iterations,
            perform_initial_solve=True,  
            verbose=self.verbose
        )

        # IO service
        self.io_handler = IOHandler(
            detail_level=detail_level,
            log_level=log_level,
            output_dir=output_directory
        )


    def run(self):
        """Run the complete G&R simulation."""

        self._setup_outputs()

        self.logger.section("Starting G&R simulation")
        self.logger.info(f"Steps: {self.n_steps}, dt: {self.dt}")

        # Write initial state (t=0, homeostatic)
        self.logger.section("TIMESTEP 0 (Homeostatic State)")
        time_0 = 0
        self._write_timestep_data(self.current_timestep, time_0)

        #TODO: refactor with current_timestep variable
        # March through time     
        for step in range(self.n_steps):

            next_timestep = self.current_timestep + 1
            time = next_timestep * self.dt
            self.logger.section(f"TIMESTEP {next_timestep}/{self.n_steps}")
            
            # Step 1: Make initial guesses for fixed point iteration
            self._make_initial_guesses(next_timestep)

            # Step 2: Apply perturbations
            # TODO: Streamline this in combination with guesses. Currenty we need to recompute wss and F after perturbations, if we perturb axial stretch or flow.
            # TODO: In 3D, every time we compute F we might need to update J as well, see when J = rhoR_alpha/rhoR_h is valid.
            self.configuration.apply_perturbations(next_timestep, time)
            self.configuration.update_wss_and_F(next_timestep)
                        
            # Step 3: Solve Coupled Mass-Geometry System
            result = self.fixed_point_solver.solve(
                configuration=self.configuration,
                timestep=next_timestep,
                dt=self.dt,
                integration_method=self.integration_method,
                survival_function_computation=self.survival_function_computation,
                equilibrium_solver_params={
                    'solver_method': self.equilibrium_solver_method,
                    'tolerance': self.equilibrium_tolerance,
                    'verbose': self.equilibrium_verbose
                }
            )
            
            if not result['converged']:
                self.logger.warning(
                    f"⚠️  Timestep {next_timestep}: {result['message']}"
                )

            self._write_timestep_data(next_timestep, time)
            self.current_timestep = next_timestep  

        self._finalize_outputs()    

    def _make_initial_guesses(
        self, 
        timestep: int, 
        guess_method: str = "from_previous_timestep"
    ) -> None:
        """Make initial guesses for all variables at target timestep.
        
        Initializes:
        - Mass densities (ρ_α)
        - Geometry (a, h, λ_z)
        - Loading variables (P, Q)
        - Stress and WSS
        
        Args:
            timestep: Target timestep to make guesses for
            guess_method: Guessing strategy (default: from previous timestep)
        """
        # TODO: refactor guess_method options with registry and expose to input params file
        # TODO: refactor which guesses are needed based on perturbations and kinetics
        self.configuration.guess_all_rhoR_alpha(timestep, guess_method)
        self.configuration.guess_geometry(timestep, guess_method)
        self.configuration.guess_loading_variables(timestep)
        self.configuration.guess_stress_and_wss(timestep)
                
    def _setup_outputs(self):
        """Setup all output files at start."""
        self.io_handler.setup_simulation_summary()
        self.io_handler.setup_metadata_output()
        self.configuration.setup_outputs(self.io_handler)
        
    
    def _write_timestep_data(self, timestep: int, time: float):
        """Write data for current timestep."""
        # Write summary (one row per layer)
        for summary_row in self.configuration.to_dict_summary(timestep, time):
            self.io_handler.write_summary_row(summary_row)
        
        # Write detailed data (layer + constituent files)
        self.configuration.write_timestep_data(self.io_handler, timestep, time)
    
    def _finalize_outputs(self):
        """Finalize and close all outputs."""
        self.io_handler.update_metadata('statistics', {
            'total_steps': self.current_timestep,
            'final_time': self.current_timestep * self.dt
        })
        self.io_handler.close_all()
        
        print(f"✓ Simulation complete. Results saved to {self.io_handler.output_dir}/")

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
