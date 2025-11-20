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
        
        # Add external support layers if present
        # if 'external_supports' in params:
        #    config._initialize_support_layers(params['external_supports'])
        
        print(f"Configuration initialized with {len(config.layers)} layers")
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
        
        # Add external support layers if present
        if 'external_supports' in params:
            self._initialize_support_layers(params['external_supports'])
    
    def _initialize_support_layers(self, support_params):
        """Initialize external support layers."""
        for support_name, properties in support_params.items():
            support_layer = Layer(support_name, "synthetic", properties)
            
            if 'constituents' in properties:
                for constituent_name, constituent_props in properties['constituents'].items():
                    constituent = Constituent.from_parameters(constituent_name, constituent_props)
                    support_layer.add_constituent(constituent)
            
            self.add_layer(support_layer)

    def advance_time(self, current_time, dt):
        """Advance configuration one timestep."""
        print(f"Advancing configuration by dt={dt}")
        current_time += dt

        # 1. Update all constituents (bottom-up kinetics)
        self._update_all_constituents(dt)
        
        # 2. Update all layers (mechanical state from constituents)
        self._update_all_layers(dt)
        
        # 3. Handle layer interactions and constraints
        self._enforce_layer_interactions()
        
        # 4. Update global configuration properties
        self._update_global_state()
    
    def _update_all_constituents(self, dt):
        """Update kinetics for all constituents in all layers."""
        for layer in self.layers:
            layer.update_constituent_kinetics(dt)
    
    def _update_all_layers(self, dt):
        """Update mechanical state for all layers."""
        for layer in self.layers:
            layer.update_mechanical_state(dt)
    
    def _enforce_layer_interactions(self):
        """Enforce constraints between layers (contact, compatibility, etc.)."""
        # Example: Ensure layers remain in contact
        for i in range(len(self.layers) - 1):
            inner_layer = self.layers[i]
            outer_layer = self.layers[i + 1]
            
            # Enforce geometric compatibility
            outer_layer.enforce_contact_with(inner_layer)
    
    def _update_global_state(self):
        """Update global configuration properties from layer states."""
        self.total_mass = sum(layer.get_total_mass() for layer in self.layers)
        # Update other global properties as needed    

    def add_layer(self, layer):
        self.layers.append(layer)

    def get_vessel_geometry(self):
        """Get geometry from the vessel layer."""
        for layer in self.layers:
            if layer.layer_type == "biological":
                return layer.get_geometry()

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

    def compute_all_stress(self, timestep: int, dt: float, 
                      integration_method: str,
                      survival_function_computation: str) -> None:
        """Compute stresses for all layers."""
        print(f"Computing stress for all layers at timestep {timestep}")
        
        for layer in self.layers:
            layer.compute_all_stress(
                timestep,
                dt,
                integration_method,
                survival_function_computation
            )