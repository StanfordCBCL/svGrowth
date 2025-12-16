from typing import List, Dict, Any, Optional
from io_handler import IOHandler
from layer import Layer
from constituent import Constituent

class Configuration:
    def __init__(self, params=None):
        self.layers: List[Layer] = []
        self.params = params or {}
        
        # Global configuration state (computed from layers)
        self.total_mass = 0.0
        self.global_stress_state = None
    
    @classmethod
    def from_parameters(cls, params: Dict[str, Any]) -> 'Configuration':
        """Create and fully initialize Configuration from parameter dictionary."""
        #print("Initializing vessel configuration...")
        
        # Validate parameters first
        cls._validate_parameters(params)
        
        config = cls(params)
        
        # Initialize each layer (layer initializes its own constituents)
        layers_data = params['layers']
        for layer_data in layers_data:
            layer = Layer.from_parameters(layer_data)
            config.add_layer(layer)
        
        # Compute homeostatic stress (direct method - no simulation params needed)
        # TODO: integrate this into Layer.from_parameters.
        #print("\nComputing homeostatic stresses...")
        for layer in config.layers:
            layer.compute_homeostatic_stress_direct()
    
        #print(f"\nConfiguration initialized with {len(config.layers)} layers")
        return config

    @staticmethod
    def _validate_parameters(params):
        # TODO: Refactor at the end.
        """Validate parameter structure."""
        if 'layers' not in params:
            raise ValueError("Missing required parameter section: layers")
        if 'simulation' not in params :
            raise ValueError("Missing required section: simulation")
    
    def _enforce_layer_interactions(self):
        """Enforce constraints between layers (contact, compatibility, etc.)."""
        # [FUTURE] Implement layer interaction constraints here.

    def add_layer(self, layer):
        self.layers.append(layer)

    def guess_all_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess mass densities for all constituents in all layers using specified method."""
        #print(f"Guessing rhoR_alpha for all constituents at timestep {target_timestep} using method '{guess_method}' method")
        
        for layer in self.layers:
            layer.guess_all_rhoR_alpha(target_timestep, guess_method)

    def guess_loading_variables(self, target_timestep: int) -> None:
        """Guess loading variables for all layers at target timestep."""
        for layer in self.layers:
            layer.guess_loading_variables(target_timestep)

    def compute_all_rhoR(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute mass densities for all constituents in all layers.
        
        Returns:
            List[float]: Total mass density for each layer (kg/m³)
        """
        #print(f"Computing rhoR for all layers at timestep {target_timestep}")
        
        rhoR_for_all_layers = []
        for layer in self.layers:
            rhoR = layer.compute_all_rhoR_alpha(
                target_timestep, 
                dt, 
                integration_method, 
                survival_function_computation
            )
            rhoR_for_all_layers.append(rhoR)
        
        return rhoR_for_all_layers

    def guess_stress_and_wss(self, target_timestep: int) -> None:
        """Guess stress and wall shear stress for all layers at target timestep."""
        #print(f"  Guessing stress and WSS for all layers at timestep {target_timestep}")

        for layer in self.layers:
            layer.guess_stress_and_wss(target_timestep)

    def guess_geometry(self, timestep: int, 
                  guess_method: str = "from_previous_timestep") -> None:
        """Guess geometry for all layers at target timestep."""
        for layer in self.layers:
            layer.guess_geometry(timestep, guess_method)

    def compute_all_stress(self, target_timestep: int, dt: float, 
                      integration_method: str,
                      survival_function_computation: str) -> None:
        """Compute stresses for all layers."""
        #print(f"Computing stress for all layers at target timestep {target_timestep}")
        
        for layer in self.layers:
            layer.compute_all_stress(
                target_timestep,
                dt,
                integration_method,
                survival_function_computation
            )

    def update_geometry_from_volume_ratio(self, timestep: int, guess_variable: str = "mid_radius") -> None:
        """Update geometry for all layers based on incompressibility constraint.
        
        Uses the volume ratio J = ρ_h / ρ(t) to refine geometry guess.
        
        Args:
            timestep: Target timestep to update geometry for
            guess_variable: Which variable to hold constant during refinement
                        - "mid_radius": Hold mid-radius constant, update inner_radius and thickness
                        - "inner_radius": [FUTURE] Hold inner radius constant
        """
        #print(f"Updating geometry from density at timestep {timestep} (method: {guess_variable})")
        
        for layer in self.layers:
            layer.update_geometry_from_volume_ratio(timestep, guess_variable)

    def update_wss(self, timestep: int) -> None:
        """Update wall shear stress for all layers at given timestep.
        
        Called after geometry updates (e.g., from incompressibility constraint) 
        to recompute WSS based on new inner radius.
        """
        #print(f"Updating WSS for all layers at timestep {timestep}")
        
        for layer in self.layers:
            layer.update_wss(timestep)


    def solve_equilibrium_geometry(
        self,
        timestep: int,
        dt: float,
        integration_method: str,
        survival_function_computation: str,
        solver_method: str = 'brentq',
        tolerance: float = 1e-5,
        verbose: bool = False
    ) -> dict:
        """Solve for equilibrium geometry for all layers.
        
        Encapsulates the iterative equilibrium solution process.
        For each trial geometry value, the solver executes:
            1. Update geometry from incompressibility (J = ρ_h/ρ)
            2. Update WSS (depends on inner radius)
            3. Compute mixture stress (via heredity integrals)
            4. Compute theoretical stress (from Laplace law)
            5. Check residual: σ_θθ(mixture) - σ_θθ(theoretical)
        
        Finds geometry where residual ≈ 0.
        
        Each layer's kinematics defines its own equilibrium condition.
        
        Args:
            timestep: Target timestep
            dt: Time step size
            integration_method: Integration method for heredity integrals
            survival_function_computation: Survival computation method
            solver_method: Root-finding method ('brentq', 'toms748', etc.)
            tolerance: Convergence tolerance
            verbose: Print solver iteration details
            
        Returns:
            Dictionary with:
                - 'all_converged': Boolean (True if all layers converged)
                - 'layer_results': List of result dicts from each layer
        """
        #print(f"Solving equilibrium geometry at timestep {timestep}")
        
        layer_results = []
        for layer in self.layers:
            result = layer.solve_equilibrium_geometry(
                timestep=timestep,
                dt=dt,
                integration_method=integration_method,
                survival_function_computation=survival_function_computation,
                solver_method=solver_method,
                tolerance=tolerance,
                verbose=verbose
            )
            
            if result['converged']:
                print(f"  Layer '{layer.name}': ✓ Converged "
                    f"(geometry={result['solution']*1000:.4f} mm, "
                    f"iterations={result['iterations']})")
            
            layer_results.append(result)
        
        return {
            'all_converged': all(r['converged'] for r in layer_results),
            'layer_results': layer_results
        }
    

    def setup_outputs(self, io_handler: 'IOHandler'):
        """Setup outputs for all layers and constituents.
        
        Args:
            io_handler: IOHandler instance to setup
        """
        for layer in self.layers:
            io_handler.setup_layer_output(layer.name)
            
            if io_handler.detail_level >= 1:
                for constituent in layer.constituents:
                    io_handler.setup_constituent_output(
                        constituent.name,
                        layer_name=layer.name
                    )
    
    def to_dict_summary(self, timestep: int, time: float) -> List[Dict[str, Any]]:
        """Get summary data for all layers.
        
        Returns:
            List of dicts, one per layer
        """
        summary_rows = []
        for layer in self.layers:
            summary_rows.append({
                'time': time,
                'layer_name': layer.name,
                'a': layer.get_inner_radius(timestep),
                'h': layer.get_thickness(timestep),
                'stress_theta': layer.get_stress(timestep)[1,1],
                'stress_z': layer.get_stress(timestep)[2,2],
                'rho': layer.get_density(timestep),
                'pressure': layer.get_pressure(timestep),
                'flow_rate': layer.get_flow_rate(timestep),
                'wss': layer.get_wss(timestep)
            })
        return summary_rows
    
    def write_timestep_data(self, io_handler: 'IOHandler', timestep: int, time: float):
        """Write detailed timestep data for all layers.
        
        Args:
            io_handler: IOHandler to write to
            timestep: Current timestep index
            time: Current simulation time
        """
        for layer in self.layers:
            # Write layer data
            layer_data = layer.to_dict(timestep, time)
            io_handler.write_layer_row(layer.name, layer_data)
            
            # Write constituent data (if detail level permits)
            if io_handler.detail_level >= 1:
                for constituent in layer.constituents:
                    const_data = constituent.to_dict(timestep, time)
                    io_handler.write_constituent_row(
                        constituent.name,
                        const_data,
                        layer_name=layer.name
                    )
    def apply_perturbations(self, timestep: int, current_time: float):
        """Apply perturbations to all layers.
        
        This is called from Simulation at the start of each timestep to update
        loading conditions based on active perturbations.
        
        Args:
        current_time: Current physical time (days)
        
        Side Effects:
            Updates layer histories with perturbed values
        """
        for layer in self.layers:
            layer.apply_perturbations(timestep, current_time)

    def update_wss_and_F(self, timestep: int) -> None:
        """Update wall shear stress and deformation gradient for all layers at given timestep.
        
        Called after perturbations to recompute WSS and F based on new flow rate and axial stretch.
        """
        for layer in self.layers:
            layer.update_wss(timestep)
            layer.update_deformation_gradient(timestep)