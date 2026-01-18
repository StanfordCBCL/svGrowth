"""
Generic solver module for G&R simulations.

Provides root-finding algorithms decoupled from specific problems.
Follows adapter pattern - Solver operates on SolverContext interface.
"""

from abc import ABC, abstractmethod
from typing import Tuple, Optional
from scipy.optimize import root_scalar
import numpy as np
from custom_logging import get_logger


# =============================================================================
# EXCEPTIONS
# =============================================================================

class SolverError(Exception):
    """Base exception for solver errors."""
    pass


class ConvergenceError(SolverError):
    """Raised when solver fails to converge."""
    pass


# =============================================================================
# SOLVER CONTEXT - Abstract interface
# =============================================================================

class SolverContext(ABC):
    """Abstract interface for solver objective functions.
    
    Decouples Solver from specific problem details (vessel geometry, stress, etc.)
    Similar to KineticsContext - provides adapter pattern.
    """
    
    @abstractmethod
    def evaluate_objective(self, trial_value: float) -> float:
        """Evaluate objective function at trial value.
        
        Args:
            trial_value: Trial value for unknown variable
            
        Returns:
            Residual (should be zero at solution)
        """
        pass
    
    @abstractmethod
    def get_search_bounds(self) -> Tuple[float, float]:
        """Get bounds for bracketing search.
        
        Returns:
            (lower_bound, upper_bound)
        """
        pass
    
    def get_initial_guess(self) -> float:
        """Get initial guess for iterative methods.
        
        Default: midpoint of search bounds.
        Override for better initial guess.
        
        Returns:
            Initial guess value
        """
        lower, upper = self.get_search_bounds()
        return (lower + upper) / 2.0
    
    def evaluate_jacobian(self, trial_value: float) -> Optional[float]:
        """Evaluate Jacobian (derivative) of objective function.
        
        Optional - only needed for Newton methods.
        If not implemented, finite differences will be used.
        
        Args:
            trial_value: Trial value for unknown variable
            
        Returns:
            Derivative df/dx, or None to use finite differences
        """
        return None  # Default: use finite differences

# =============================================================================
# GENERIC EQUILIBRIUM CONTEXT - For layer equilibrium problems
# =============================================================================

class LayerEquilibriumContext(SolverContext):
    """Generic equilibrium context for layer equilibrium problems.
    
    Works for any kinematics type (thin-wall, thick-wall, etc.).
    Delegates to Layer.compute_equilibrium_residual().
    """
    
    def __init__(self,
                 layer,
                 timestep: int,
                 dt: float,
                 integration_method: str,
                 survival_function_computation: str):
        """Initialize equilibrium context.
        
        Args:
            layer: Layer object to solve equilibrium for
            timestep: Current timestep
            dt: Time step size
            integration_method: Integration method for heredity integrals
            survival_function_computation: Survival computation method
        """
        self.layer = layer
        self.timestep = timestep
        self.dt = dt
        self.integration_method = integration_method
        self.survival_function_computation = survival_function_computation
        
        # Get reference value from kinematics (for bracketing)
        self.reference_value = layer.kinematics.get_reference_geometry_value(layer)
    
    def evaluate_objective(self, trial_geometry_value: float) -> float:
        """Evaluate equilibrium residual.
        
        Args:
            trial_value: Trial value for geometry (not used directly,
                        passed to layer for geometry update)
        
        Returns:
            Residual: σ_θθ(mixture) - σ_θθ(theoretical)
        """
        return self.layer.compute_equilibrium_residual(
            trial_geometry_value,
            self.timestep,
            self.dt,
            self.integration_method,
            self.survival_function_computation,
        )
    
    def get_search_bounds(self) -> Tuple[float, float]:
        """Get search bounds: ±5% bracket around reference value."""
        bracket = 0.6
        return (
            self.reference_value * (1.0 - bracket),
            self.reference_value * (1.0 + bracket)
        )
    
    def get_initial_guess(self) -> float:
        """Get initial guess from previous timestep."""
        if self.timestep > 0:
            # Get geometry variable name from kinematics
            geometry_var = self.layer.kinematics.get_equilibrium_variable()
            
            if geometry_var == 'mid_radius':
                return self.layer.get_mid_radius(self.timestep - 1)
            elif geometry_var == 'inner_radius':
                return self.layer.get_inner_radius(self.timestep - 1)
        
        return self.reference_value
    
# =============================================================================
# SOLVER - Generic root-finding algorithms
# =============================================================================

class Solver:
    """Generic root-finding solver.
    
    Stateless - operates on SolverContext (adapter pattern).
    Never stores problem-specific data.
    
    Similar to Kinetics/Mechanics classes - provides algorithms,
    context provides data access.
    """
    
    def __init__(self, 
                 method: str = 'brentq', 
                 tolerance: float = 1e-12, 
                 max_iterations: int = 100):
        """Initialize solver with method and tolerances.
        
        Args:
            method: Root-finding method
                   - 'brentq': Brent's method (default, robust)
                   - 'toms748': TOMS Algorithm 748 (more robust than Brent)
                   - 'bisect': Bisection method (slowest, most robust)
                   - 'newton_fd': Newton with finite difference Jacobian (future)
            tolerance: Convergence tolerance for residual
            max_iterations: Maximum iterations allowed
        """
        self.method = method
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.logger = get_logger(__name__)
        
        # Validate method
        valid_methods = ['brentq', 'toms748', 'bisect', 'newton_fd']
        if method not in valid_methods:
            raise ValueError(f"Unknown method '{method}'. Valid: {valid_methods}")
    
    def solve(self, context: SolverContext, verbose: bool = False) -> dict:
        """Solve for root using configured method.
        
        Args:
            context: Problem-specific adapter implementing objective function
            verbose: Print iteration details
            
        Returns:
            Dictionary with:
                - 'solution': Root value (float)
                - 'residual': Final residual value (float)
                - 'converged': Boolean convergence flag
                - 'iterations': Number of iterations (int)
                - 'message': Status message (str)
        """
        if self.method in ['brentq', 'toms748', 'bisect']:
            return self._solve_bracketing(context, verbose)
        elif self.method == 'newton_fd':
            return self._solve_newton_fd(context, verbose)
        else:
            raise ValueError(f"Method '{self.method}' not implemented")
    
    def _solve_bracketing(self, context: SolverContext, verbose: bool) -> dict:
        """Solve using bracketing method (Brent, TOMS748, bisection).
        
        These methods guarantee convergence if root is bracketed.
        """
        # Get search bounds
        lower, upper = context.get_search_bounds()
        
        self.logger.debug(
            f"Solver: {self.method}, bracket=[{lower:.6e}, {upper:.6e}]"
        )

        # Define objective function wrapper
        def objective(x):
            return context.evaluate_objective(x)
        
        try:
            # Check that root is bracketed
            # TODO: to verify: check bracketing inside root_scalar?
            # f_lower = objective(lower)
            # f_upper = objective(upper)
            
            # if verbose:
            #     print(f"    f(lower) = {f_lower:.6e}")
            #     print(f"    f(upper) = {f_upper:.6e}")
            
            # if f_lower * f_upper > 0:
            #     # Same sign - not bracketed
            #     raise ConvergenceError(
            #         f"Root not bracketed: f({lower:.6e})={f_lower:.6e}, "
            #         f"f({upper:.6e})={f_upper:.6e}"
            #     )
            
            # Solve using scipy
            result = root_scalar(
                objective,
                bracket=[lower, upper],
                method=self.method,
                rtol=self.tolerance, #TODO: expose to yaml file
                maxiter=self.max_iterations #TODO: expose to yaml file
            )
            
            self.logger.debug(
                f"Converged: x = {result.root:.6e}, "
                f"f(x) = {result.function_calls:.6e}, "
                f"iterations = {result.iterations}"
            )
            
            return {
                'solution': result.root,
                'residual': objective(result.root),
                'converged': result.converged,
                'iterations': result.iterations,
                'message': result.flag
            }
            
        except ValueError as e:
            # Bracketing failed or other scipy error
            return {
                'solution': None,
                'residual': None,
                'converged': False,
                'iterations': 0,
                'message': f'Solver error: {str(e)}'
            }
        except ConvergenceError as e:
            return {
                'solution': None,
                'residual': None,
                'converged': False,
                'iterations': 0,
                'message': str(e)
            }
    
    def _solve_newton_fd(self, context: SolverContext, verbose: bool) -> dict:
        """Solve using Newton method with finite difference Jacobian.
        
        Faster convergence than bracketing, but less robust.
        Requires good initial guess.
        """
        # Get initial guess
        x = context.get_initial_guess()
        epsilon = 1e-6  # Finite difference step
        
        if verbose:
            print(f"  Solver: Newton (finite differences), x0 = {x:.6e}")
        
        for iteration in range(self.max_iterations):
            # Evaluate function
            f_x = context.evaluate_objective(x)
            
            # Check convergence
            if abs(f_x) < self.tolerance:
                if verbose:
                    print(f"    Converged: x = {x:.6e}, f(x) = {f_x:.6e}, "
                          f"iterations = {iteration+1}")
                return {
                    'solution': x,
                    'residual': f_x,
                    'converged': True,
                    'iterations': iteration + 1,
                    'message': 'Converged'
                }
            
            # TODO: can implement analytical Jacobian here
            df_dx = context.evaluate_jacobian(x)
            
            if df_dx is None:
                # Use finite differences
                f_x_plus = context.evaluate_objective(x + epsilon)
                df_dx = (f_x_plus - f_x) / epsilon
            
            # Check for zero derivative
            if abs(df_dx) < 1e-12:
                return {
                    'solution': None,
                    'residual': f_x,
                    'converged': False,
                    'iterations': iteration + 1,
                    'message': 'Zero derivative encountered'
                }
            
            # Newton update
            x_new = x - f_x / df_dx
            
            if verbose and iteration % 10 == 0:
                print(f"    Iteration {iteration}: x = {x:.6e}, f(x) = {f_x:.6e}")
            
            x = x_new
        
        # Max iterations reached
        return {
            'solution': x,
            'residual': context.evaluate_objective(x),
            'converged': False,
            'iterations': self.max_iterations,
            'message': f'Max iterations ({self.max_iterations}) reached'
        }