import math
from abc import ABC, abstractmethod
from typing import Optional


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
    """Factory for creating perturbation objects from configuration."""
    
    @staticmethod
    def create(perturbation_type: str, perturbation_time: float, 
               perturbation_percentage: Optional[float] = None,
               duration: Optional[float] = None) -> Perturbation:
        """Create perturbation from parameters.
        
        Args:
            perturbation_type: Type of perturbation ("step", "linear", "latorre2018")
            perturbation_time: Time at which perturbation begins (days)
            perturbation_percentage: Magnitude as percentage (required for step/linear)
            duration: Duration for linear ramp (required for linear only)
            
        Returns:
            Perturbation object
            
        Raises:
            ValueError: If required parameters are missing or invalid type specified
        """
        perturbation_type = perturbation_type.lower()
        
        if perturbation_type == "step":
            if perturbation_percentage is None:
                raise ValueError("Perturbation percentage must be specified for step perturbation")
            return StepPerturbation(perturbation_time, perturbation_percentage)
        
        elif perturbation_type == "linear":
            if perturbation_percentage is None:
                raise ValueError("Perturbation percentage must be specified for linear perturbation")
            if duration is None:
                raise ValueError("Duration must be specified for linear perturbation")
            return LinearPerturbation(perturbation_time, perturbation_percentage, duration)
        
        elif perturbation_type == "latorre2018":
            return Latorre2018Perturbation(perturbation_time)
        
        else:
            raise ValueError(
                f"Invalid perturbation type: '{perturbation_type}'. "
                f"Choose 'step', 'linear', or 'latorre2018'."
            )


# Backward compatibility function
def apply_perturbation(current_time: float, perturbation_time: float, 
                      homeostatic_value: float, perturbation_type: str = "step",
                      perturbation_percentage: Optional[float] = None,
                      duration: Optional[float] = None) -> float:
    """Legacy function for backward compatibility.
    
    Creates a perturbation object and applies it. For new code, use the
    Perturbation classes directly.
    
    Args:
        current_time: Current simulation time (days)
        perturbation_time: Time at which perturbation begins (days)
        homeostatic_value: Baseline homeostatic value
        perturbation_type: Type of perturbation ("step", "linear", "latorre2018")
        perturbation_percentage: Magnitude as percentage
        duration: Duration for linear ramp (days)
        
    Returns:
        Perturbed value
    """
    perturbation = PerturbationFactory.create(
        perturbation_type, perturbation_time, 
        perturbation_percentage, duration
    )
    
    return perturbation.apply(current_time, homeostatic_value)