"""
Kinetics module for G&R simulations.

Handles:
1. Degradation rate functions (k_alpha)
2. Survival functions (q)
3. Production rate functions (mR_alpha)
4. Numerical integration of heredity integrals (q, mR_alpha, mR_alpha*q)
"""

from abc import ABC, abstractmethod
from typing import List
import math
from kinetics_interface import KineticsDataNotAvailableError
from integrators import SurvivalFunctionComputationFactory, IntegratorFactory


# =============================================================================
# MODULE-LEVEL HELPER FUNCTIONS
# =============================================================================

def _compute_delta(current, homeostatic):
    """Compute deviation from homeostatic value.
    
    Automatically handles both relative and absolute delta:
    - If homeostatic > 0: returns (current - homeostatic) / homeostatic (relative)
    - If homeostatic = 0: returns current (absolute)
    
    Args:
        current: Current stimulus value
        homeostatic: Homeostatic (reference) stimulus value
        
    Returns:
        Delta (dimensionless)
        
    Example:
        >>> _compute_delta(150.0, 100.0)
        0.5  # 50% increase
        >>> _compute_delta(75.0, 100.0)
        -0.25  # 25% decrease
    """
    if homeostatic == 0.0:
        return current  # Absolute delta for zero-homeostatic stimuli (e.g., inflammation)
    
    return (current - homeostatic) / homeostatic


# =============================================================================
# DEGRADATION RATE FUNCTIONS
# =============================================================================

class DegradationRateFunction(ABC):
    """Abstract base class for degradation rate functions k_alpha."""
    
    @abstractmethod
    def compute_k_alpha(self, context, current_timestep: int) -> float:
        """Compute degradation rate at current timestep.
        
        Args:
            context: KineticsContext providing data access
            current_timestep: Current simulation timestep
        
        Returns:
            k_alpha value (degradation rate, 1/day)
        """
        pass

    #TODO: validate inputs function

class ConstantDegradationRate(DegradationRateFunction):
    """Constant degradation rate: k_alpha = k_alpha_h"""
    
    def __init__(self, k_alpha_h):
        self.k_alpha_h = k_alpha_h
    
    def compute_k_alpha(self, context, current_timestep):
        return self.k_alpha_h


class QuadraticDegradationRate(DegradationRateFunction):
    """Quadratic stimulus-modulated degradation.
    
    Formula: k_alpha = k_alpha_h * (1 + Σ K_i * δ_i²)
    where δ_i = (stimulus_i - homeostatic_i) / homeostatic_i
    
    Supported stimuli (accessed via context):
    - intramural_stress: Layer-level total stress
    - constituent_stress: Constituent-level stress
    - wss: Wall shear stress
    - inflammation: Inflammatory markers
    """
    
    def __init__(self, k_alpha_h, gain_params=None):
        """Initialize quadratic degradation.
        
        Args:
            k_alpha_h: Homeostatic degradation rate (1/day)
            gain_params: Dict mapping stimulus names to gain coefficients K_i
                        Example: {'intramural_stress': 1.0, 'wss': 0.5}
        """
        self.k_alpha_h = k_alpha_h
        self.gain_params = gain_params or {}
    
    def compute_k_alpha(self, context, current_timestep):
        """Compute degradation rate: k_alpha = k_alpha_h * (1 + Σ K_i * δ_i²)"""
        upsilon = 1.0
        
        # Add contribution from each stimulus
        for stimulus_name, K_value in self.gain_params.items():
            try:
                current = context.get_stimulus(stimulus_name, current_timestep)
                homeostatic = context.get_stimulus_homeostatic(stimulus_name)
                
                delta = _compute_delta(current, homeostatic)
                upsilon += K_value * (delta ** 2)

            except KineticsDataNotAvailableError as e:
                print(f"[ERROR] Stimulus '{stimulus_name}' not available: {e}")
                pass
        
        return self.k_alpha_h * upsilon


# =============================================================================
# SURVIVAL FUNCTIONS
# =============================================================================

class SurvivalFunction(ABC):
    """Abstract base class for survival functions q(s, tau)."""
    
    @abstractmethod
    def compute_survival_function(
        self,
        k_alpha_history: list[float],
        tau_min: int,
        dt: float,
        current_timestep: int,
        integration_method: str,
        survival_strategy: str
    ) -> list[float]: 
        """Compute survival function q(s, τ) for all cohorts.
        
        The survival function q(s, τ) represents the fraction of material
        deposited at time τ that survives until time s.
        
        Args:
            k_alpha_history: List of degradation rates [k_alpha(0), ..., k_alpha(s)]
            tau_min: Minimum deposition time (birth time of constituent)
            dt: Time step size (days)
            current_timestep: Current time index s
            integration_method: Integration method ('simpson' or 'trapezoidal')
            survival_strategy: Strategy for survival computation ('backward', 'naive')
            
        Returns:
            List of survival values [q(s, tau_min), q(s, tau_min+1), ..., q(s, s-1), q(s, s)]
            where q(s, s) = 1.0 (just deposited material has 100% survival)
            
        Example:
            >>> # At timestep s=5, with tau_min=0
            >>> survival_values = compute_survival_function(
            ...     k_alpha_history, 0, 1.0, 5, 'simpson', 'backward'   
            ... )
            >>> # Returns: [q(5,0), q(5,1), q(5,2), q(5,3), q(5,4), q(5,5)]
            >>> #           [0.65,   0.72,   0.79,   0.86,   0.93,   1.00]
        """
        pass


class ExponentialSurvival(SurvivalFunction):
    """Exponential survival: q(s, τ) = exp(-∫[τ to s] k_α dt
    
    This implementation delegates to SurvivalFunctionComputationFactory to select
    the appropriate computational strategy (backward or naive).
    """
    
    def __init__(self, degradation_function):
        """Initialize exponential survival function.
        
        Args:
            degradation_function: DegradationRateFunction instance (for potential future use)
        """
        self.degradation_function = degradation_function
        
    def compute_survival_function(
        self,
        k_alpha_history: list[float],
        tau_min: int,
        dt: float,
        current_timestep: int,
        integration_method: str,
        survival_function_computation: str
    ) -> list[float]:
        """Compute survival function q(s, τ) for all cohorts.
        
        Args:
            k_alpha_history: Degradation rate history [k(0), ..., k(s)]
            tau_min: Earliest deposition time
            dt: Time step size (days)
            current_timestep: Current time s
            integration_method: Integration method ('simpson' or 'trapezoidal')
            survival_function_computation: Strategy to use ('backward', 'naive')
            
        Returns:
            Survival values [q(s, τ_min), ..., q(s, s)]
        """
        # Create strategy instance based on user's choice
        strategy = SurvivalFunctionComputationFactory.create(survival_function_computation)

        # Delegate to strategy
        return strategy.compute_survival(
            k_alpha_history, 
            tau_min, 
            dt, 
            current_timestep, 
            integration_method
        )


# =============================================================================
# PRODUCTION RATE FUNCTIONS
# =============================================================================

class ProductionRateFunction(ABC):
    """Abstract base class for production rate functions mR_alpha."""
    
    @abstractmethod
    def compute_production_rate(
        self,
        context,
        current_timestep: int,
        k_alpha: float,
        rhoR_alpha: float
    ) -> float:
        """Compute production rate at current timestep.
        
        Args:
            context: KineticsContext providing data access
            current_timestep: Current simulation timestep
            k_alpha: Current degradation rate (1/day)
            rhoR_alpha: Current referential mass density (kg/m³)
        
        Returns:
            mR_alpha value (production rate, kg/(m³·day))
        """
        pass


class LinearProductionRate(ProductionRateFunction):
    """Linear stimulus-modulated production.
    
    Formula: mR_alpha = rhoR_alpha * k_alpha * (1 + Σ K_i * δ_i)
    
    where:
    - rhoR_alpha: Current referential mass density
    - k_alpha: Current degradation rate
    - δ_i = (stimulus_i - homeostatic_i) / homeostatic_i
    
    At homeostasis (all δ_i = 0):
        mR_alpha = rhoR_alpha * k_alpha
    which maintains steady state: production = degradation.
    """
    
    def __init__(self, gain_params=None):
        """Initialize linear production.
        
        Args:
            gain_params: Dict mapping stimulus names to gain coefficients K_i
                        Example: {'intramural_stress': 2.0, 'wss': 1.0}
        """
        self.gain_params = gain_params or {}
    
    def compute_production_rate(self, context, current_timestep, k_alpha, rhoR_alpha):
        """Compute production rate: mR_alpha = rhoR_alpha * k_alpha * (1 + Σ K_i * δ_i)"""
        # Base production rate (maintains homeostasis when γ = 1)
        m_base = rhoR_alpha * k_alpha
        
        # Modulation factor
        gamma = 1.0
        
        # Add contribution from each stimulus
        for stimulus_name, K_value in self.gain_params.items():
            try:
                current = context.get_stimulus(stimulus_name, current_timestep)
                homeostatic = context.get_stimulus_homeostatic(stimulus_name)
                
                delta = _compute_delta(current, homeostatic)
                gamma += K_value * delta
                
            except KineticsDataNotAvailableError as e:
                print(f"[ERROR] Stimulus '{stimulus_name}' not available: {e}")
                pass
            
        return m_base * gamma

    #TODO: return the closed form solution of production rate
    # def __repr__(self):


# =============================================================================
# MAIN KINETICS CLASS
# =============================================================================

class Kinetics:
    """Main kinetics orchestrator combining degradation, survival, production, and integration."""
    
    def __init__(
        self,
        degradation_function: DegradationRateFunction,
        survival_function: SurvivalFunction,
        production_function: ProductionRateFunction
    ) -> None:
        self.degradation_function = degradation_function
        self.survival_function = survival_function
        self.production_function = production_function
    
    def compute_k_alpha(self, context, current_timestep):
        """Compute degradation rate at current timestep."""
        return self.degradation_function.compute_k_alpha(context, current_timestep)
    
    def compute_production_rate(self, context, current_timestep, k_alpha, rhoR_alpha):
        """Compute production rate at current timestep.
        
        Args:
            context: KineticsContext for stimulus access
            current_timestep: Current time s
            k_alpha: Degradation rate at current timestep
            rhoR_alpha: Referential mass density at current timestep
            
        Returns:
            Production rate mR_alpha (kg/(m³·day))
        """
        return self.production_function.compute_production_rate(
            context, current_timestep, k_alpha, rhoR_alpha
        )
    
    def compute_survival_function(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        current_timestep: int,
        integration_method: str,  
        survival_function_computation: str
    ) -> List[float]:
        """Compute survival function q(s, τ) for all cohorts.
        
        Returns:
            List [q(s, τ_min), ..., q(s, s)]
        """
        return self.survival_function.compute_survival_function(
            k_alpha_history, 
            tau_min, 
            dt, 
            current_timestep,
            integration_method,
            survival_function_computation
        )
    
    def compute_heredity_integral(
        self,
        mR_history: List[float],
        q_values: List[float],
        tau_min: int,
        dt: float,
        current_timestep: int,
        integration_method: str,
        rhoR_alpha_homeostatic: float
    ) -> float:
        """Compute heredity integral: ρ = ∫ mR(τ) × q(s,τ) dτ
        
        This includes:
        1. Pre-existing homeostatic material (τ <= 0, "founding cohort")
        2. Newly deposited material (τ > 0)
        
        Args:
            mR_history: Production rate history [mR(τ_min), ..., mR(s)]
            q_values: Survival function [q(s,τ_min), ..., q(s,s)]
            tau_min: Earliest deposition time
            dt: Time step size
            current_timestep: Current time s
            integration_method: Integration method ('simpson' or 'trapezoidal')
            rhoR_alpha_homeostatic: Homeostatic mass density
            
        Returns:
            Mass density rhoR_alpha
        """
        
        # STEP 1: Compute contribution from newly deposited material (τ > tau_min)
        # Form product mR × q for new material
        mq_values = [mR_history[tau_min + i] * q_values[i] 
                     for i in range(len(q_values))]
        
        n_points = len(mq_values)
        
        # TODO: add tau_min dependency
        integrator = IntegratorFactory.create(
            integration_method, dt, n_points - 1, 0
        )
        
        rhoR_deposited = integrator.integrate(mq_values)
        
        # STEP 2: Add contribution from pre-existing homeostatic material (τ <= 0)
        # This material was present at t=0 and survives with q(s, -∞)
        if tau_min == 0:
            q_initial = q_values[0]  
            rhoR_initial = rhoR_alpha_homeostatic * q_initial
        else:
            rhoR_initial = 0.0
        
        rhoR_alpha = rhoR_initial + rhoR_deposited

        return rhoR_alpha
    
    @classmethod
    def from_parameters(cls, params):
        """Factory method to create Kinetics from YAML parameters.
        
        Expected YAML structure:
        kinetics:
          production:
            stimulus_function_form: "linear"
            gain_params:
              intramural_stress: 2.0
              wss: -1.0
          degradation:
            survival_function_form: "exponential"
            deg_rate:
              type: "quadratic"
              k_alpha_h: 0.071
              gain_params:
                intramural_stress: 1.0
        """
        
        # Parse degradation
        degradation_params = params.get('degradation', {})
        
        # Extract deg_rate parameters
        deg_rate_params = degradation_params.get('deg_rate', {})
        deg_rate_type = deg_rate_params.get('type', 'constant')

        if deg_rate_type == 'constant':
            k_alpha_h = deg_rate_params.get('k_alpha_h', 0.0)
            deg_rate_function = ConstantDegradationRate(k_alpha_h)

        elif deg_rate_type == 'quadratic':
            k_alpha_h = deg_rate_params.get('k_alpha_h', 0.0)
            gain_params = deg_rate_params.get('gain_params', {})
            deg_rate_function = QuadraticDegradationRate(k_alpha_h, gain_params)
        else:
            raise ValueError(f"Unknown degradation rate type: {deg_rate_type}")
        
        # Parse survival
        survival_params = degradation_params
        survival_type = survival_params.get('survival_function_form', 'exponential')
        
        if survival_type == 'exponential':
            survival_function = ExponentialSurvival(deg_rate_function)
        else:
            raise ValueError(f"Unknown survival type: {survival_type}")
        
        # Parse production
        production_params = params.get('production', {})
        prod_type = production_params.get('stimulus_function_form', 'linear')
        
        if prod_type == 'linear':
            gain_params = production_params.get('gain_params', {})
            production_function = LinearProductionRate(gain_params)
        else:
            raise ValueError(f"Unknown production type: {prod_type}")
        
        return cls(deg_rate_function, survival_function, production_function)