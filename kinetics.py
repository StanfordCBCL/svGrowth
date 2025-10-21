from abc import ABC, abstractmethod
import math

class DegradationRateFunction(ABC):
    """Abstract base class for constituent degradation rate functions (k_alpha)."""
    
    @abstractmethod
    def compute_k_alpha(self, current_timestep, deposition_timestep, state):
        """Compute k_alpha degradation rate.
        
        Args:
            current_timestep: Current simulation timestep
            deposition_timestep: When the material was deposited (tau)
            state: Dict containing state variables (stress, etc.)
        
        Returns:
            k_alpha value (degradation rate)
        """
        pass


class ConstantDegradationRate(DegradationRateFunction):
    """Constant degradation rate: k_alpha = k_h"""
    
    def __init__(self, k_h):
        self.k_h = k_h
    
    def compute_k_alpha(self, current_timestep, deposition_timestep, state):
        return self.k_h


class QuadraticDegradationRate(DegradationRateFunction):
    """Quadratic stimulus-modulated degradation: k_alpha = k_h * upsilon
    
    Where upsilon = 1 + sum(K_i * (delta_i)^2) for various stimuli.
    
    Supported stimuli:
    - sigma: Intramural/Cauchy stress
    - tauw: Wall shear stress
    - inflammation: Inflammatory markers (e.g., cytokines)
    
    Args:
        k_h: Baseline degradation rate
        stimulus_params: Dict of stimulus names to sensitivity parameters
                        e.g., {'sigma': 0.5, 'tauw': 0.1, 'inflammation': 0.3}
        use_relative_delta: If True, compute (x - x_h)/x_h; else x - x_h
    """
    
    def __init__(self, k_h, stimulus_params=None, use_relative_delta=True):
        self.k_h = k_h
        self.stimulus_params = stimulus_params or {}
        self.use_relative_delta = use_relative_delta
    
    def compute_k_alpha(self, current_timestep, deposition_timestep, state):
        upsilon = 1.0
        
        # Add contribution from each stimulus
        for stimulus_name, K_value in self.stimulus_params.items():
            current = state.get(stimulus_name, 0.0)
            homeostatic = state.get(f"{stimulus_name}_h", 0.0)
            
            # Compute delta
            if self.use_relative_delta:
                delta = (current - homeostatic) / homeostatic if homeostatic != 0 else 0.0
            else:
                delta = current - homeostatic
            
            # Add quadratic contribution
            upsilon += K_value * (delta ** 2)
        
        return self.k_h * upsilon


class DataDrivenDegradationRate(DegradationRateFunction):
    """K_alpha values loaded from data file."""
    
    def __init__(self, data_file_path):
        self.data = self._load_data(data_file_path)
    
    def _load_data(self, file_path):
        """Load k_alpha data from file.
        
        TODO: Implement data loading from .dat file
        Should return dict or array indexed by timestep
        """
        raise NotImplementedError("Data-driven degradation not yet implemented")
    
    def compute_k_alpha(self, current_timestep, deposition_timestep, state):
        return self.data.get(current_timestep, 0.0)


class SurvivalFunction(ABC):
    """Abstract base class for survival functions q(tau, t)."""
    
    @abstractmethod
    def compute_survival(self, deposition_timestep, current_timestep, degradation_function, state_history, dt):
        """Compute survival fraction from deposition to current time.
        
        Args:
            deposition_timestep: When material was deposited (tau)
            current_timestep: Current time (t)
            degradation_function: DegradationRateFunction instance
            state_history: Function or dict that returns state at any timestep
            dt: Time increment
        
        Returns:
            Survival fraction q in [0, 1]
        """
        pass


class ExponentialSurvival(SurvivalFunction):
    """Exponential survival: q(tau,t) = exp(-integral(k_alpha, tau, t))"""
    
    def compute_survival(self, deposition_timestep, current_timestep, degradation_function, state_history, dt):
        if current_timestep <= deposition_timestep:
            return 1.0
        
        # Integrate k_alpha from deposition to current time
        integral = 0.0
        for t in range(deposition_timestep + 1, current_timestep + 1):
            state = state_history(t)
            k_val = degradation_function.compute_k_alpha(t, deposition_timestep, state)
            integral += k_val * dt
        
        return math.exp(-integral)


class Kinetics:
    """Main kinetics class that combines degradation and survival functions."""
    
    def __init__(self, degradation_function, survival_function):
        self.degradation_function = degradation_function
        self.survival_function = survival_function
    
    def compute_survival(self, deposition_timestep, current_timestep, state_history, dt):
        """Compute survival from deposition to current time."""
        return self.survival_function.compute_survival(
            deposition_timestep, 
            current_timestep, 
            self.degradation_function, 
            state_history, 
            dt
        )
    
    @classmethod
    def from_parameters(cls, params):
        """Factory method to create Kinetics from parameter dictionary."""
        
        # Parse degradation function
        degradation_params = params.get('degradation', {})
        degradation_type = degradation_params.get('type', 'constant')
        
        if degradation_type == 'constant':
            k_h = degradation_params.get('k_h', 0.0)
            degradation_function = ConstantDegradationRate(k_h) 
            
        elif degradation_type == 'quadratic':  
            k_h = degradation_params.get('k_h', 0.0)
            stimulus_params = degradation_params.get('stimulus_params', {})  
            use_relative = degradation_params.get('use_relative_delta', True)
            degradation_function = QuadraticDegradationRate(k_h, stimulus_params, use_relative)  

        elif degradation_type == 'data_driven':
            data_file = degradation_params.get('data_file')
            degradation_function = DataDrivenDegradationRate(data_file)

        else:
            raise ValueError(f"Unknown degradation type: {degradation_type}")
        
        # Parse survival function
        survival_params = params.get('survival', {})
        survival_type = survival_params.get('type', 'exponential')
        
        if survival_type == 'exponential':
            survival_function = ExponentialSurvival()
        else:
            raise ValueError(f"Unknown survival function type: {survival_type}")
        
        return cls(degradation_function, survival_function)