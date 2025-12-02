from typing import List, Dict, Any, Optional
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
        print("Initializing vessel configuration...")
        
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
        print("\nComputing homeostatic stresses...")
        for layer in config.layers:
            layer.compute_homeostatic_stress_direct()
    
        print(f"\nConfiguration initialized with {len(config.layers)} layers")
        return config

    @staticmethod
    def _validate_parameters(params):
        # TODO: Refactor at the end.
        """Validate parameter structure."""
        if 'layers' not in params:
            raise ValueError("Missing required parameter section: layers")
        if 'simulation' not in params :
            raise ValueError("Missing required section: simulation")
        
    def _initialize_configuration_from_parameters(self, params):
        """Initialize configuration: validate and create layers."""
        # Validate parameters
        self._validate_parameters(params)
        
        # Initialize each layer (layer initializes its own constituents)
        layer_data = params['layer']
        layer = Layer.from_parameters(layer_data, params)
        self.add_layer(layer)
    
    def _enforce_layer_interactions(self):
        """Enforce constraints between layers (contact, compatibility, etc.)."""
        # Example: Ensure layers remain in contact
        for i in range(len(self.layers) - 1):
            inner_layer = self.layers[i]
            outer_layer = self.layers[i + 1]
            
            # Enforce geometric compatibility
            outer_layer.enforce_contact_with(inner_layer) 

    def add_layer(self, layer):
        self.layers.append(layer)

    def guess_all_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess mass densities for all constituents in all layers using specified method."""
        print(f"Guessing rhoR_alpha for all constituents at timestep {target_timestep} using method '{guess_method}' method")
        
        for layer in self.layers:
            layer.guess_all_rhoR_alpha(target_timestep, guess_method)

    def compute_all_rhoR_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute mass densities for all constituents in all layers."""
        print(f"Computing rhoR_alpha for all constituents at timestep {target_timestep}")
        
        for layer in self.layers:
            layer.compute_all_rhoR_alpha(
                target_timestep, 
                dt, 
                integration_method, 
                survival_function_computation
            )

    def guess_stress_and_wss(self, target_timestep: int) -> None:
        """Guess stress and wall shear stress for all layers at target timestep."""
        print(f"  Guessing stress and WSS for all layers at timestep {target_timestep}")

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
        print(f"Computing stress for all layers at target timestep {target_timestep}")
        
        for layer in self.layers:
            layer.compute_all_stress(
                target_timestep,
                dt,
                integration_method,
                survival_function_computation
            )