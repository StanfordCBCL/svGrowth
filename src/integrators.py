from typing import List, Tuple
from abc import ABC, abstractmethod
import math

class Integrator(ABC):
    """
    Abstract base class for numerical integrators.

    Attributes:
        dt (float): Time step size between discrete points.
        start (int): Starting index for integration (inclusive).
        stop (int): Ending index for integration (inclusive).
    """

    def __init__(self, dt: float, start: int, stop: int) -> None:
        """
        Initialize the integrator.

        Args:
            dt (float): Time step size.
            start (int): Start index for integration (inclusive).
            stop (int): Stop index for integration (inclusive).
        """
        self.dt = dt
        self.start = start
        self.stop = stop

    @abstractmethod
    def integrate(self, f: List[float]) -> float:
        """
        Compute the numerical integral over the given discrete values.

        Args:
            f (List[float]): Discrete function values over which to integrate.

        Returns:
            float: Numerical approximation of the integral.
        """
        raise NotImplementedError


class TrapezoidIntegrator(Integrator):
    """
    Trapezoidal rule numerical integrator.

    Integrates discrete values using the trapezoidal rule
    from start to stop indices (counting backward).
    """

    def integrate(self, f: List[float]) -> float:
        """
        Perform trapezoidal integration.

        Args:
            f (List[float]): Discrete function values.

        Returns:
            float: Approximate integral value.
        
        Raises:
            ValueError: If fewer than two points are available.
        """
        n = (self.start - self.stop)/self.dt + 1 
        if n < 2:
            raise ValueError("At least two points are required for trapezoidal integration.")

        integral_value = 0.0
        # Integration backward in time
        for i in range(self.start, self.stop, -1):
            integral_value += (f[i] + f[i - 1]) * self.dt / 2

        return integral_value


class SimpsonIntegrator(Integrator):
    """
    Simpson's rule numerical integrator.

    Integrates discrete values using Simpson's 1/3 rule.
    If number of intervals is even, applies trapezoidal
    rule on the last interval.
    """

    def integrate(self, f: List[float]) -> float:
        """
        Perform Simpson's rule integration.

        Args:
            f (List[float]): Discrete function values.

        Returns:
            float: Approximate integral value.
        
        Raises:
            ValueError: If fewer than two points are available.
        """
        n = (self.start - self.stop)/self.dt + 1  # number of integration points

        if n < 2:
            raise ValueError("At least two points required for integration with Simpson's rule.")
        
        if n == 2:
            trapezoid_rule = TrapezoidIntegrator(self.dt, self.start, self.stop)
            return trapezoid_rule.integrate(f)

        simpson_stop = self.stop  # Stop index for the main Simpson's rule loop

        # If number of integration points is even, apply trapezoidal rule on last interval
        trapezoid_value = 0.0
        if n % 2 == 0:
            trapezoid_rule = TrapezoidIntegrator(self.dt, self.stop + 1, self.stop)
            trapezoid_value += trapezoid_rule.integrate(f)
            simpson_stop += 1  # Exclude last interval from Simpson’s rule

        # Add first and last terms
        simpson_value = f[self.start] + f[simpson_stop]

        # Add 4 * f(odd indices)
        for i in range(self.start - 1, simpson_stop, -2):
            simpson_value += 4 * f[i]

        # Add 2 * f(even indices)
        for i in range(self.start - 2, simpson_stop + 1, -2):
            simpson_value += 2 * f[i]

        simpson_value *= self.dt / 3

        return trapezoid_value + simpson_value


# =============================================================================
# INTEGRATOR FACTORY - Single place to manage all integrators
# =============================================================================

class IntegratorFactory:
    """Factory for creating integrators by name.
    
    Benefits:
    - Single place to register new integrators
    - No if-statements scattered throughout code
    - Easy to see all available methods
    - Extensible: just add to registry
    """
    
    # Registry of available integrators
    _registry = {
        'simpson': SimpsonIntegrator,
        'trapezoidal': TrapezoidIntegrator,
        # Add new methods here:
        # 'adaptive_simpson': AdaptiveSimpsonIntegrator,
        # 'gauss_quadrature': GaussQuadratureIntegrator,
    }
    
    @classmethod
    def create(cls, method: str, dt: float, start: int, stop: int) -> Integrator:
        """Create an integrator instance.
        
        Args:
            method: Integration method name ('simpson', 'trapezoidal', etc.)
            dt: Time step size
            start: Start index (inclusive)
            stop: Stop index (inclusive)
            
        Returns:
            Integrator instance
            
        Raises:
            ValueError: If method is not recognized
            
        Example:
            >>> integrator = IntegratorFactory.create('simpson', dt=0.1, start=10, stop=0)
            >>> result = integrator.integrate(data)
        """
        method_lower = method.lower()
        
        if method_lower not in cls._registry:
            available = ', '.join(cls._registry.keys())
            raise ValueError(
                f"Unknown integration method: '{method}'. "
                f"Available methods: {available}"
            )
        
        integrator_class = cls._registry[method_lower]
        return integrator_class(dt, start, stop)
    
    @classmethod
    def register(cls, name: str, integrator_class):
        """Register a new integrator type.
        
        Allows plugins/extensions to add custom integrators:
        
        Example:
            >>> class MyCustomIntegrator(Integrator):
            ...     def integrate(self, f):
            ...         # custom implementation
            ...         pass
            >>> 
            >>> IntegratorFactory.register('custom', MyCustomIntegrator)
        """
        cls._registry[name.lower()] = integrator_class
    
    @classmethod
    def available_methods(cls) -> List[str]:
        """Get list of available integration methods."""
        return list(cls._registry.keys())


# Convenience function for quick access
def create_integrator(method: str, dt: float, start: int, stop: int) -> Integrator:
    """Convenience function to create integrators.
    
    Example:
        >>> from integrators import create_integrator
        >>> integrator = create_integrator('simpson', 0.1, 10, 0)
    """
    return IntegratorFactory.create(method, dt, start, stop)


# =============================================================================
# SURVIVAL FUNCTION COMPUTATION STRATEGIES
# =============================================================================

class SurvivalFunctionComputation(ABC):
    """Abstract base class for survival function computation strategies.
    
    A survival strategy defines HOW to compute q(s, τ) values, independent
    of WHICH basic integrator (trapezoidal/Simpson's) is used.
    """
    
    @abstractmethod
    def compute_survival(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        current_time: int,
        integrator_method: str
    ) -> List[float]:
        """Compute survival function q(s, τ) for all cohorts.
        
        Args:
            k_alpha_history: Degradation rate history [k(0), ..., k(s)]
            tau_min: Earliest deposition time
            dt: Time step size
            current_time: Current time index s
            integrator_method: Which basic integrator to use ('simpson'/'trapezoidal')
            
        Returns:
            List of survival values [q(s, τ_min), ..., q(s, s)]
        """
        pass


class NaiveSurvivalFunctionComputation(SurvivalFunctionComputation):
    """Naive O(n²) survival computation.
    
    Computes each q(s, τ) independently by integrating k_α from τ to s.
    
    Simple and clear, but slow for large problems. Best for:
    - Debugging and validation
    - Small problems (< 50 cohorts)
    - Comparing against optimized methods
    
    Uses the specified integrator_method for each cohort's integral.
    """
    
    def compute_survival(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        current_time: int,
        integrator_method: str
    ) -> List[float]:
        """Compute survival using independent integration for each cohort.
        
        For each cohort τ:
            q(s, τ) = exp(-∫[τ to s] k_α(t) dt)
        
        Complexity: O(n²) where n = number of cohorts
        
        Example with 5 points, s=4:
            q(4,0) = exp(-∫[0→4] k_α)  ← integrate 4 intervals
            q(4,1) = exp(-∫[1→4] k_α)  ← integrate 3 intervals
            q(4,2) = exp(-∫[2→4] k_α)  ← integrate 2 intervals
            q(4,3) = exp(-∫[3→4] k_α)  ← integrate 1 interval
            q(4,4) = 1.0               ← no integration
            
        Total: 4+3+2+1 = 10 integrations for 5 cohorts
        """
        s = current_time
        q_values = []
        
        for tau in range(tau_min, s + 1):
            if tau == s:
                q_values.append(1.0)
            else:
                # Create integrator for this cohort's full lifetime [τ, s]
                integrator = IntegratorFactory.create(
                    integrator_method, dt, s, tau
                )
                integral = integrator.integrate(k_alpha_history)
                q_values.append(math.exp(-integral))
        
        return q_values


class BackwardSurvivalFunctionComputation(SurvivalFunctionComputation):
    """Optimized O(n) backward iteration survival computation.
    
    Computes survival values incrementally using the telescoping property:
        q(s, τ) = q(s, τ+1) × exp(-∫[τ to τ+1] k_α)
    
    This reduces complexity from O(n²) to O(n) by only integrating
    single intervals (or pairs for Simpson's) and using previous results.
    
    Best for:
    - Production simulations with many cohorts (> 50)
    - Long-term G&R simulations (years of daily time steps)
    
    Note on accuracy for Simpson's method:
    - Reference survival values (segment boundaries): 4th order accuracy
    - Intermediate survival values (midpoints): 2nd order accuracy
    - This is acceptable because intermediate values are only used
      as weights for the final heredity integral computation.
    """
    
    def compute_survival(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        current_time: int,
        integrator_method: str
    ) -> List[float]:
        """Compute survival using backward iteration.
        
        Complexity: O(n) where n = number of cohorts
        
        Example with 5 points, s=4, trapezoidal:
            q[4] = 1.0
            q[3] = q[4] × exp(-∫[3→4])  ← 1 integration
            q[2] = q[3] × exp(-∫[2→3])  ← 1 integration
            q[1] = q[2] × exp(-∫[1→2])  ← 1 integration
            q[0] = q[1] × exp(-∫[0→1])  ← 1 integration
            
        Total: 4 integrations for 5 cohorts
        """
        s = current_time
        n_cohorts = s - tau_min + 1
        
        q_values = [0.0] * n_cohorts
        q_values[-1] = 1.0  # q(s, s) = 1.0
        
        if n_cohorts == 1:
            return q_values
        
        # Choose algorithm based on integrator method
        # TODO: get rid of string check
        if integrator_method == 'simpson':
            return self._compute_backward_simpson(
                k_alpha_history, tau_min, dt, s, n_cohorts, q_values
            )
        else:  # trapezoidal
            return self._compute_backward_trapezoidal(
                k_alpha_history, tau_min, dt, s, n_cohorts, q_values
            )
    
    def _compute_backward_trapezoidal(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        s: int,
        n_cohorts: int,
        q_values: List[float]
    ) -> List[float]:
        """Backward iteration with trapezoidal rule.
        
        Sequential: q(s, τ) = q(s, τ+1) × exp(-∫[τ to τ+1] k_α)
        All q values have consistent 2nd order accuracy.
        """
        for i in range(n_cohorts - 2, -1, -1):
            tau = tau_min + i
            
            # Integrate single interval [τ, τ+1]
            #TODO: stop using 'trapezoidal' to create object.
            integrator = IntegratorFactory.create('trapezoidal', dt, tau + 1, tau)
            integral = integrator.integrate(k_alpha_history)
            
            # Incremental update
            q_values[i] = q_values[i + 1] * math.exp(-integral)
        
        return q_values
    
    def _compute_backward_simpson(
        self,
        k_alpha_history: List[float],
        tau_min: int,
        dt: float,
        s: int,
        n_cohorts: int,
        q_values: List[float]
    ) -> List[float]:
        """Backward iteration with Simpson's rule.
        
        Processes pairs of intervals for Simpson's rule efficiency.
        
        Accuracy notes:
        - Reference points (segment boundaries): 4th order (Simpson's rule)
        - Midpoints: 2nd order (trapezoidal)
        
        The midpoint accuracy is sufficient because these values are only
        used as intermediate weights when the constituent computes the
        heredity integral ∫ mR × q.
        
        Example with 5 points, s=4:
            Iter 1: Process [2,3,4]
                q[2] = exp(-∫Simpson[2→4])  ← 4th order
                q[3] = exp(-∫Trap[3→4])     ← 2nd order (midpoint)
            Iter 2: Process [0,1,2]
                q[0] = exp(-∫Simpson[0→2])  ← 4th order
                q[1] = exp(-∫Trap[1→2])     ← 2nd order (midpoint)
        """
        even_n = (n_cohorts % 2 == 0)

        # Track the reference point for incremental computation
        tau_ref = s
        q_ref = 1.0
        idx = n_cohorts - 1
        
        # Process pairs of intervals
        tau = s - 1
        while tau > tau_min:
            # Compute midpoint q(s, τ) using trapezoidal [τ, τ_ref]
            integrator_mid = IntegratorFactory.create('trapezoidal', dt, tau_ref, tau)
            integral_mid = integrator_mid.integrate(k_alpha_history)
            q_values[idx - 1] = math.exp(-integral_mid) * q_ref
            
            # Compute reference q(s, τ-1) using Simpson's [τ-1, τ_ref]
            integrator_old = IntegratorFactory.create('simpson', dt, tau_ref, tau - 1)
            integral_old = integrator_old.integrate(k_alpha_history)
            q_values[idx - 2] = math.exp(-integral_old) * q_ref
            
            # Shift reference
            tau_ref = tau - 1
            q_ref = q_values[idx - 2]
            idx -= 2
            tau -= 2
        
        # Handle remaining interval if even number of points
        if even_n:
            integrator = IntegratorFactory.create('trapezoidal', dt, tau_ref, tau_min)
            integral = integrator.integrate(k_alpha_history)
            q_values[0] = math.exp(-integral) * q_ref
        
        return q_values


# =============================================================================
# SURVIVAL FUNCTION COMPUTATION FACTORY
# =============================================================================

class SurvivalFunctionComputationFactory:
    """Factory for creating survival computation strategies.
    
    Allows users to choose between different algorithms for computing
    survival functions, independent of the underlying integration method.
    """
    
    _registry = {
        'naive': NaiveSurvivalFunctionComputation,
        'backward': BackwardSurvivalFunctionComputation,
    }
    
    @classmethod
    def create(cls, strategy: str) -> SurvivalFunctionComputation:
        """Create a survival strategy instance.
        
        Args:
            strategy: Strategy name ('naive' or 'backward')
            
        Returns:
            SurvivalFunctionComputation instance
            
        Raises:
            ValueError: If strategy is not recognized
            
        Example:
            >>> strategy = SurvivalFunctionComputationFactory.create('backward')
            >>> q_vals = strategy.compute_survival(k_alpha, 0, 0.1, 100, 'simpson')
        """
        strategy_lower = strategy.lower()
        
        if strategy_lower not in cls._registry:
            available = ', '.join(cls._registry.keys())
            raise ValueError(
                f"Unknown survival function computation: '{strategy}'. "
                f"Available strategies: {available}"
            )
        
        strategy_class = cls._registry[strategy_lower]
        return strategy_class()
    
    @classmethod
    def register(cls, name: str, strategy_class):
        """Register a new survival strategy.
        
        Example:
            >>> class AdaptiveSurvivalFunctionComputation(SurvivalFunctionComputation):
            ...     def compute_survival(self, ...):
            ...         # adaptive implementation
            ...         pass
            >>> 
            >>> SurvivalFunctionComputationFactory.register('adaptive', AdaptiveSurvivalFunctionComputation)
        """
        cls._registry[name.lower()] = strategy_class
    
    @classmethod
    def available_strategies(cls) -> List[str]:
        """Get list of available survival strategies."""
        return list(cls._registry.keys())