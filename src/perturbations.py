import math
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any, List


class Perturbation(ABC):
    """Base class for perturbation strategies.
    
    A perturbation modifies a homeostatic value over time based on a specific
    temporal profile (step, linear ramp, exponential, etc.).
    """
    
    def __init__(self, perturbation_time: float, perturbation_percentage: float):
        """Initialize perturbation.
        
        Args:
            perturbation_time: Time at which perturbation begins (days)
            perturbation_percentage: Magnitude as percentage of homeostatic value
                                    (positive = increase, negative = decrease)
        """
        self.perturbation_time = perturbation_time
        self.perturbation_percentage = perturbation_percentage
    
    @abstractmethod
    def apply(self, current_time: float, homeostatic_value: float) -> float:
        """Compute perturbed value at current time.
        
        Args:
            current_time: Current simulation time (days)
            homeostatic_value: Baseline homeostatic value
            
        Returns:
            Perturbed value
        """
        pass
    
    def is_active(self, current_time: float) -> bool:
        """Check if perturbation is currently active.
        
        Args:
            current_time: Current simulation time (days)
            
        Returns:
            True if perturbation has started, False otherwise
        """
        return current_time >= self.perturbation_time
    
    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"t={self.perturbation_time} days, "
                f"magnitude={self.perturbation_percentage:+.1f}%)")


class StepPerturbation(Perturbation):
    """Instantaneous step change in value.
    
    At t = t_pert, value jumps immediately:
        V(t) = V_h                           for t < t_pert
        V(t) = V_h · (1 + percentage/100)    for t >= t_pert
    """
    
    def apply(self, current_time: float, homeostatic_value: float) -> float:
        """Apply step perturbation."""
        if not self.is_active(current_time):
            return homeostatic_value
        
        perturbation_magnitude = homeostatic_value * (self.perturbation_percentage / 100.0)
        return homeostatic_value + perturbation_magnitude


class LinearPerturbation(Perturbation):
    """Linear ramp to target value over specified duration.
    
    Value changes linearly from V_h to target over duration:
        V(t) = V_h                                      for t < t_pert
        V(t) = V_h + slope · (t - t_pert)              for t_pert <= t < t_pert + duration
        V(t) = V_h · (1 + percentage/100)              for t >= t_pert + duration
    """
    
    def __init__(self, perturbation_time: float, perturbation_percentage: float, 
                 duration: float):
        """Initialize linear perturbation.
        
        Args:
            perturbation_time: Time at which ramp begins (days)
            perturbation_percentage: Target magnitude as percentage
            duration: Duration of linear ramp (days)
        """
        super().__init__(perturbation_time, perturbation_percentage)
        
        if duration <= 0:
            raise ValueError(f"Duration must be positive, got {duration}")
        
        self.duration = duration
    
    def apply(self, current_time: float, homeostatic_value: float) -> float:
        """Apply linear ramp perturbation."""
        if not self.is_active(current_time):
            return homeostatic_value
        
        elapsed_time = current_time - self.perturbation_time
        perturbation_magnitude = homeostatic_value * (self.perturbation_percentage / 100.0)
        target_value = homeostatic_value + perturbation_magnitude
        
        # Linear ramp
        if elapsed_time < self.duration:
            slope = perturbation_magnitude / self.duration
            return homeostatic_value + slope * elapsed_time
        else:
            # Ramp complete - return final value
            return target_value
    
    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"t={self.perturbation_time} days, "
                f"magnitude={self.perturbation_percentage:+.1f}%, "
                f"duration={self.duration} days)")


class Latorre2018Perturbation(Perturbation):
    """Exponential approach to 50% increase (Latorre 2018 paper).
    
    Specific profile from Latorre et al. (2018):
        V(t) = V_h                                           for t < t_pert
        V(t) = V_h · [1 + 0.5 · (1 - exp(-(t-t_pert)/10))]  for t >= t_pert
    
    Asymptotically approaches 50% increase with time constant τ = 10 days.
    """
    
    def __init__(self, perturbation_time: float):
        """Initialize Latorre 2018 perturbation.
        
        Args:
            perturbation_time: Time at which perturbation begins (days)
        
        Note:
            Magnitude (50%) and time constant (10 days) are fixed by the model.
        """
        # Fixed 50% increase for Latorre 2018 model
        super().__init__(perturbation_time, perturbation_percentage=50.0)
        self.time_constant = 10.0  # days
    
    def apply(self, current_time: float, homeostatic_value: float) -> float:
        """Apply Latorre 2018 exponential perturbation."""
        if not self.is_active(current_time):
            return homeostatic_value
        
        elapsed_time = current_time - self.perturbation_time
        factor = 1.0 - math.exp(-elapsed_time / self.time_constant)
        
        return homeostatic_value * (1.0 + 0.5 * factor)
    
    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"t={self.perturbation_time} days, "
                f"50% increase, τ={self.time_constant} days)")


class PerturbationFactory:
    """Factory for creating perturbations using registry pattern.
    
    Benefits:
    - Extensible: register new perturbation types easily
    - Centralized: single place to manage all perturbation types
    - Clean: no scattered if-statements
    """
    
    # Registry of available perturbation types
    _registry = {
        'step': {
            'class': StepPerturbation,
            'required_params': ['perturbation_time', 'perturbation_percentage'],
            'optional_params': {}
        },
        'linear': {
            'class': LinearPerturbation,
            'required_params': ['perturbation_time', 'perturbation_percentage', 'duration'],
            'optional_params': {}
        },
        'latorre2018': {
            'class': Latorre2018Perturbation,
            'required_params': ['perturbation_time'],
            'optional_params': {}
        }
    }
    
    @classmethod
    def create(cls, perturbation_config: Dict[str, Any]) -> Optional[Perturbation]:
        """Create perturbation from configuration dictionary.
        
        Args:
            perturbation_config: Dict with keys:
                - 'type': Type of perturbation ('step', 'linear', etc.)
                - 'perturbation_time': Time when perturbation starts (days)
                - 'perturbation_percentage': Magnitude (for step/linear)
                - 'duration': Duration (for linear only)
                - Other type-specific parameters
        
        Returns:
            Perturbation instance, or None if no perturbation specified
            
        Raises:
            ValueError: If invalid type or missing required parameters
        """
        if perturbation_config is None or perturbation_config.get('type') is None:
            return None
        
        perturbation_type = perturbation_config['type'].lower()
        
        if perturbation_type not in cls._registry:
            available = ', '.join(cls._registry.keys())
            raise ValueError(
                f"Unknown perturbation type: '{perturbation_type}'. "
                f"Available types: {available}"
            )
        
        # Extract parameters and create instance
        spec = cls._registry[perturbation_type]
        perturbation_class = spec['class']

        params = {}
        # Add required parameters
        for param_name in spec['required_params']:
            params[param_name] = perturbation_config[param_name]

        # Add optional parameters with defaults
        for param_name, default_value in spec['optional_params'].items():
            if param_name in perturbation_config:
                params[param_name] = perturbation_config[param_name]
            else:
                params[param_name] = default_value
        
        return perturbation_class(**params)
    
    @classmethod
    def register(cls, name: str, perturbation_class):
        """Register a new perturbation type.
        
        Example:
            >>> class CustomPerturbation(Perturbation):
            ...     # implementation
            >>> PerturbationFactory.register('custom', CustomPerturbation)
        """
        cls._registry[name.lower()] = perturbation_class
    
    @classmethod
    def available_types(cls) -> List[str]:
        """Get list of available perturbation types."""
        return list(cls._registry.keys())


class PerturbableVariableRegistry:
    """Central registry of all variables that can be perturbed.
    
    Organized by scope (layer, constituent, kinetics, etc.).
    This is the single source of truth for what can be perturbed.
    """
    
    # Layer-level variables
    LAYER_VARIABLES = {
        'pressure': 'homeostatic_pressure',
        'flow_rate': 'homeostatic_flow_rate',
        'axial_stretch': 'homeostatic_axial_stretch',
    }
    
    # Constituent-level variables
    CONSTITUENT_VARIABLES = {
        'deposition_stretch': 'deposition_stretch',
    }
    
    # Kinetics-level variables (for future use)
    # TODO: currently kinetics perturbations are not implemented 
    # TODO: potentially code kinetics variable changes over time directly to kinetics
    KINETICS_VARIABLES = {
        'production_gain_stress': 'production_gain_params.intramural_stress',
        'degradation_gain_stress': 'degradation_gain_params.intramural_stress',
    }
    
    @classmethod
    def validate(cls, scope: str, variable_name: str) -> None:
        """Validate that a variable is perturbable in the given scope.
        
        Args:
            scope: Scope name ('layer', 'constituent', 'kinetics')
            variable_name: Variable name from YAML
        
        Raises:
            ValueError: If scope or variable name is invalid
        """
        valid_scopes = {
            'layer': cls.LAYER_VARIABLES,
            'constituent': cls.CONSTITUENT_VARIABLES,
            'kinetics': cls.KINETICS_VARIABLES,
        }
        
        if scope not in valid_scopes:
            raise ValueError(
                f"Invalid scope: '{scope}'. "
                f"Valid scopes: {', '.join(valid_scopes.keys())}"
            )
        
        scope_vars = valid_scopes[scope]
        
        if variable_name not in scope_vars:
            valid_names = ', '.join(scope_vars.keys())
            raise ValueError(
                f"Invalid variable '{variable_name}' for scope '{scope}'. "
                f"Valid variables: {valid_names}"
            )

    @classmethod
    def get_value(cls, scope: str, variable_name: str, obj) -> float:
        #TODO: to rethink, scope vs obj?
        """Extract homeostatic value from object using explicit attribute access.
        
        This method encapsulates all the mapping logic. It knows how to get
        values from objects based on scope and variable name.
        
        Args:
            scope: Scope ('layer', 'constituent', etc.)
            variable_name: Variable name (e.g., 'pressure')
            obj: Object to extract value from (Layer, Constituent, etc.)
        
        Returns:
            Homeostatic value
        
        Example:
            >>> layer = Layer("test")
            >>> value = PerturbableVariableRegistry.get_value('layer', 'pressure', layer)
            >>> # Returns layer.homeostatic_pressure
        """
        # Validate first
        cls.validate(scope, variable_name)
        
        # Get attribute name from registry
        attr_name = cls.get_attribute_name(scope, variable_name)
        
        # Extract value using explicit dispatch based on scope
        if scope == 'layer':
            return cls._get_layer_value(obj, attr_name)
        elif scope == 'constituent':
            return cls._get_constituent_value(obj, attr_name)
        else:
            raise ValueError(f"Unknown scope: {scope}")

    @classmethod
    def get_attribute_name(cls, scope: str, variable_name: str) -> str:
        """Get the attribute name for a variable in a given scope.
        
        Args:
            scope: Scope name ('layer', 'constituent', etc.)
            variable_name: Variable name from YAML
        
        Returns:
            Attribute name on the object
        """
        cls.validate(scope, variable_name)
        
        if scope == 'layer':
            return cls.LAYER_VARIABLES[variable_name]
        elif scope == 'constituent':
            return cls.CONSTITUENT_VARIABLES[variable_name]
        elif scope == 'kinetics':
            return cls.KINETICS_VARIABLES[variable_name]
    
    @classmethod
    def get_valid_variables(cls, scope: str) -> List[str]:
        """Get list of valid variable names for a scope."""
        if scope == 'layer':
            return list(cls.LAYER_VARIABLES.keys())
        elif scope == 'constituent':
            return list(cls.CONSTITUENT_VARIABLES.keys())
        elif scope == 'kinetics':
            return list(cls.KINETICS_VARIABLES.keys())
        else:
            return []

    @classmethod
    def _get_layer_value(cls, layer, attr_name: str) -> float:
        """Get value from Layer object using explicit attribute access."""
        if attr_name == 'homeostatic_pressure':
            return layer.homeostatic_pressure
        elif attr_name == 'homeostatic_flow_rate':
            return layer.homeostatic_flow_rate
        elif attr_name == 'homeostatic_axial_stretch':
            return layer.homeostatic_axial_stretch
        elif attr_name == 'homeostatic_wss':
            return layer.homeostatic_wss
        else:
            raise ValueError(f"Unknown layer attribute: {attr_name}")

    @classmethod
    def _get_constituent_value(cls, constituent, attr_name: str) -> float:
        """Get value from Constituent object using explicit attribute access."""
        if attr_name == 'deposition_stretch':
            return constituent.deposition_stretch
        else:
            raise ValueError(f"Unknown constituent attribute: {attr_name}")


class PerturbationManager:
    """Manages multiple perturbations for different homeostatic variables.
    
    Responsibilities:
    - Store perturbations for each variable (pressure, flow, etc.)
    - Apply perturbations at given time
    - Track which perturbations are active
    
    Data Ownership:
    - Layer owns PerturbationManager
    - PerturbationManager owns individual Perturbation objects
    """
    
    def __init__(self, scope: str):
        """Initialize perturbation manager.
        
        Args:
            scope: Scope of perturbations ('layer', 'constituent', etc.)
        """
        self.scope = scope
        self.perturbations: Dict[str, Perturbation] = {}
    
    def add_perturbation(self, variable_name: str, perturbation: Perturbation):
        """Register a perturbation for a specific variable.
        
        Args:
            variable_name: Name of variable to perturb ('pressure', 'flow_rate', etc.)
            perturbation: Perturbation object
        """
        # Validate variable name for this scope (e.g., pressure, flow_rate, etc.)
        PerturbableVariableRegistry.validate(self.scope, variable_name)
        self.perturbations[variable_name] = perturbation
    
    def get_perturbed_value(self, variable_name: str, current_time: float, 
                           homeostatic_value: float) -> float:
        """Get perturbed value for a variable at current time.
        
        Args:
            variable_name: Name of variable ('pressure', 'flow_rate', etc.)
            current_time: Current simulation time (days)
            homeostatic_value: Baseline homeostatic value
        
        Returns:
            Perturbed value (or homeostatic value if no perturbation)
        """
        perturbation = self.perturbations.get(variable_name)
        
        if perturbation is None:
            return homeostatic_value
        
        return perturbation.apply(current_time, homeostatic_value)

    def compute_perturbed_values(self, class_object, current_time: float) -> Dict[str, float]:
        """Compute perturbed values.
        
        Args:
            class_object: Object to get values from (Layer, Constituent, etc.)
            current_time: Current simulation time (days)
        
        Returns:
            Dict of variable_name -> perturbed_value
        """
        perturbed_values = {}
        for variable_name, perturbation in self.perturbations.items():
            homeostatic_value = PerturbableVariableRegistry.get_value(
                self.scope, variable_name, class_object
            )
            perturbed_value = perturbation.apply(current_time, homeostatic_value)
            perturbed_values[variable_name] = perturbed_value
        
        return perturbed_values

    def is_any_active(self, current_time: float) -> bool:
        """Check if any perturbation is currently active."""
        for perturbation in self.perturbations.values():
            if perturbation.is_active(current_time):
                return True
        
        return False
    
    def get_active_perturbations(self, current_time: float) -> Dict[str, Perturbation]:
        """Get all currently active perturbations."""
        active = {}
        for var_name, perturbation in self.perturbations.items():
            if perturbation.is_active(current_time):
                active[var_name] = perturbation
        
        return active
    
    def __repr__(self) -> str:
        if len(self.perturbations) == 0:
            return f"PerturbationManager(scope='{self.scope}', no perturbations)"
        
        items = []
        for var_name, pert in self.perturbations.items():
            item_str = f"{var_name}: {pert}"
            items.append(item_str)
        
        items_joined = ', '.join(items)
        return f"PerturbationManager(scope='{self.scope}', {items_joined})"
    
    #TODO: add perturbation plotting functions