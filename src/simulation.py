import os
from configuration import Configuration
from io_handler import IOHandler

class Simulation:
    def __init__(self, 
                 configuration: Configuration, 
                 io_handler=None, 
                 simulation_name=None, 
                 output_directory=None, 
                 dt=1, 
                 n_steps=100, 
                 tolerance=1e-12, 
                 max_iterations=50, 
                 integration_method='simpson',
                 survival_function_computation='naive',  
                 verbose=False,
                 detail_level=None,
                 debug_level=None):
        
        self.configuration = configuration
        self.dt = dt
        self.n_steps = n_steps
        self.tolerance = float(tolerance)
        self.max_iterations = int(max_iterations)
        self.integration_method = integration_method
        self.survival_function_computation = survival_function_computation
        self.verbose = verbose
        
        self.current_timestep = 0  # Current timestep index (0 = initial/homeostatic state)
        
        # Simulation metadata (owned by Simulation)
        self.simulation_name = simulation_name or getattr(configuration, 'simulation_name', 'default_simulation')
        self.output_directory = output_directory or getattr(configuration, 'output_directory', 'output')
        
        # IO service (used by Simulation)
        self.io_handler = IOHandler(
            detail_level=detail_level,
            debug_level=debug_level,
            output_dir=output_directory
        )
        
        # Output file handle (will be set up during run)
        self.output_file = None

    def run(self):
        """Run the complete G&R simulation."""

        self._setup_outputs()
        
        #TODO: refactor step vs current_timestep        
        for step in range(self.n_steps):

            next_timestep = self.current_timestep + 1
            time = next_timestep * self.dt
            print(f"Advancing from timestep {self.current_timestep} to {next_timestep}")
            
            self.configuration.apply_perturbations(next_timestep, time)           
            # ========================================================================
            # STEP 1: Initial Guesses
            # ========================================================================
            self.configuration.guess_all_rhoR_alpha(
                next_timestep, guess_method="from_previous_timestep"
            )

            # TODO: Make sure guess geometry includes active radius for active stress
            self.configuration.guess_geometry(
                next_timestep, guess_method="from_previous_timestep"
            )

            self.configuration.guess_loading_variables(next_timestep)
            # TODO: Not needed if we pre-allocate arrays with 0. Guess stress for mass production.
            self.configuration.guess_stress_and_wss(next_timestep)

            self.configuration.apply_perturbations(next_timestep, time)

            # TODO: Streamline this in combination with guesses. Currenty we need to recompute wss and F after perturbations, if we perturb axial stretch or flow.
            # TODO: In 3D, every time we compute F we might need to update J as well, see when J = rhoR_alpha/rhoR_h is valid.
            self.configuration.update_wss_and_F(next_timestep)

            # Compute kinetics with initial guesses
            self.configuration.compute_all_rhoR(
                next_timestep, 
                self.dt, 
                self.integration_method,
                self.survival_function_computation 
            )
            
            # Solve for equilibrium geometry with initial kinetics
            result = self.configuration.solve_equilibrium_geometry(
                timestep=next_timestep,
                dt=self.dt,
                integration_method=self.integration_method,
                survival_function_computation=self.survival_function_computation,
                solver_method='brentq',
                tolerance=1e-5,
                verbose=self.verbose
            )

            if not result['all_converged']:
                raise RuntimeError(
                    f"Initial equilibrium solver failed at timestep {next_timestep}"
                )
                        
           # ====================================================================
            # STEP 4: Fixed-Point Iteration (Mass-Geometry Coupling)
            # ====================================================================           
            iteration = 0
            converged = False

            while iteration < self.max_iterations:
                iteration += 1

                # Store old mass densities BEFORE updating them
                rho_old = {}
                for layer in self.configuration.layers:
                    rho_old[layer.name] = layer.get_density(next_timestep)

                # Compute Mass Densities (Heredity Integrals)
                self.configuration.compute_all_rhoR(
                    next_timestep, 
                    self.dt, 
                    self.integration_method,
                    self.survival_function_computation 
                )
                    
                # Solve for Equilibrium Geometry
                # 
                # This encapsulates an iterative solve that repeatedly executes:
                #   TODO: check order of step 3, or if can be excluded from loop
                #   3. Update geometry from trial value (uses incompressibility J=ρ_h/ρ)
                #   4. Update WSS (depends on inner radius)
                #   5. Compute mixture stress (via constituent heredity integrals)
                #   6. Check residual: σ_θθ(mixture) - σ_θθ(theoretical)
                #
                # The solver finds geometry where residual ≈ 0.
                result = self.configuration.solve_equilibrium_geometry(
                    timestep=next_timestep,
                    dt=self.dt,
                    integration_method=self.integration_method,
                    survival_function_computation=self.survival_function_computation,
                    solver_method='brentq', # exposure to yaml file
                    tolerance=1e-5, # exposure to yaml file
                    verbose=self.verbose
                )
                
                if not result['all_converged']:
                    raise RuntimeError(
                        f"Equilibrium solver failed at timestep {next_timestep}. "
                        f"Check layer results for details."
                    )
                
                # (d) Check convergence on mass density for each layer
                all_layers_converged = True
                max_relative_change = 0.0
                
                for layer in self.configuration.layers:
                    rho_new = layer.get_density(next_timestep)
                    rho_prev = rho_old[layer.name]
                    
                    # Compute relative change
                    if rho_prev > 0:
                        relative_change = abs(rho_new - rho_prev) / rho_prev
                    else:
                        relative_change = 0.0
                    
                    max_relative_change = max(max_relative_change, relative_change)
                    
                    # Check if this layer has converged
                    layer_converged = relative_change <= self.tolerance
                    
                    if not layer_converged:
                        all_layers_converged = False

                # (e) Check stopping criteria
                if all_layers_converged:
                    print(f"\n✓ Converged in {iteration} iteration(s)")
                    print(f"  Max Δρ/ρ = {max_relative_change:.2e} < {self.tolerance:.2e}")
                    converged = True
                    break

            # Check if we exceeded max iterations
            if not converged:
                print(f"\n⚠️  Warning: Max iterations ({self.max_iterations}) reached")
                print(f"   Max Δρ/ρ = {max_relative_change:.2e} > {self.tolerance:.2e}")
                print(f"   Continuing to next timestep...")
            
            # (f) Write results for this timestep
            self._write_timestep_data(next_timestep, time)

            # Timestep complete
            self.current_timestep = next_timestep  

        self._finalize_outputs()    
        
    def _setup_outputs(self):
        """Setup all output files at start."""
        self.io_handler.setup_simulation_summary()
        self.io_handler.setup_metadata_output()
        self.configuration.setup_outputs(self.io_handler)
        
        # Store simulation config in metadata
        self.io_handler.update_metadata('configuration', {
            'dt': self.dt,
            'n_steps': self.n_steps,
            'integration_method': self.integration_method
        })
    
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
        self.io_handler.update_metadata('statistics', {
            'total_steps': self.current_timestep,
            'final_time': self.current_timestep * self.dt
        })
        self.io_handler.close_all()
        
        print(f"✓ Simulation complete. Results saved to {self.io_handler.output_dir}/")

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