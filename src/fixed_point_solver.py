"""
Fixed-point iteration solver for mass-geometry coupling in G&R simulations.

Solves the coupled system:
- Mass production depends on stress/geometry
- Geometry depends on mass (via incompressibility J = ρ_h/ρ)
- Stress depends on geometry

This creates a fixed-point problem requiring iterative refinement.
"""

from typing import Dict, Any, Optional
from custom_logging import get_logger
from configuration import Configuration


# =============================================================================
# EXCEPTIONS
# =============================================================================

class FixedPointError(Exception):
    """Base exception for fixed-point solver errors."""
    pass


class ConvergenceError(FixedPointError):
    """Raised when fixed-point iteration fails to converge."""
    pass


# =============================================================================
# CONVERGENCE CHECKER
# =============================================================================

class ConvergenceChecker:
    """Checks convergence criteria for fixed-point iteration.
    
    Responsibilities:
    - Compute convergence metrics (relative change in density)
    - Determine if iteration has converged
    - Provide diagnostic information
    """
    
    def __init__(self, tolerance: float = 1e-6):
        """Initialize convergence checker.
        
        Args:
            tolerance: Relative tolerance for convergence (dimensionless)
        """
        self.tolerance = tolerance
        self.logger = get_logger(__name__)
    
    def check_convergence(
        self,
        configuration: Configuration,
        timestep: int,
        rho_old: Dict[str, float]
    ) -> tuple[bool, float, Dict[str, float]]:
        """Check convergence for all layers.
        
        Args:
            configuration: Configuration object with layers
            timestep: Current timestep
            rho_old: Dict mapping layer name to old density
            
        Returns:
            Tuple of (all_converged, max_change, layer_changes)
            where layer_changes maps layer name to relative change
        """
        all_converged = True
        max_relative_change = 0.0
        layer_changes = {}
        
        for layer in configuration.layers:
            rho_new = layer.get_density(timestep)
            rho_prev = rho_old[layer.name]
            
            # Compute relative change
            if abs(rho_prev) > 1e-10:
                relative_change = abs(rho_new - rho_prev) / abs(rho_prev)
            else:
                relative_change = abs(rho_new - rho_prev)
            
            layer_changes[layer.name] = relative_change
            max_relative_change = max(max_relative_change, relative_change)
            
            # Check if this layer converged
            if relative_change > self.tolerance:
                all_converged = False
            
            self.logger.debug(
                f"  {layer.name}: ρ={rho_new:.2f} kg/m³, "
                f"Δρ/ρ={relative_change:.2e} "
                f"({'✓' if relative_change <= self.tolerance else '✗'})"
            )
        
        return all_converged, max_relative_change, layer_changes


# =============================================================================
# FIXED-POINT SOLVER
# =============================================================================

class FixedPointSolver:
    """Solves coupled mass-geometry system via fixed-point iteration.
    
    Algorithm:
        0. Initial solve (better starting point than guesses)
        1. Compute mass density (depends on stress/geometry from previous iteration)
        2. Solve geometric equilibrium (uses updated mass via incompressibility)
        3. Check convergence on mass density
        4. Repeat until converged or max iterations reached
    
    Responsibilities:
    - Orchestrate iteration loop
    - Perform initial solve
    - Check convergence
    - Provide diagnostics
    
    Delegates to:
    - Configuration: mass computation, equilibrium solving
    - ConvergenceChecker: convergence criteria
    """
    
    def __init__(
        self,
        tolerance: float = 1e-6,
        max_iterations: int = 50,
        perform_initial_solve: bool = True,
        verbose: bool = False
    ):
        """Initialize fixed-point solver.
        
        Args:
            tolerance: Convergence tolerance (relative change in density)
            max_iterations: Maximum iterations allowed
            perform_initial_solve: Whether to perform initial solve before iteration
            verbose: Print iteration details
        """
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.perform_initial_solve = perform_initial_solve
        self.verbose = verbose
        self.logger = get_logger(__name__)
        
        # Convergence checker
        self.convergence_checker = ConvergenceChecker(tolerance)
    
    def solve(
        self,
        configuration: Configuration,
        timestep: int,
        dt: float,
        integration_method: str,
        survival_function_computation: str,
        equilibrium_solver_params: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """Solve coupled mass-geometry system via fixed-point iteration.
        
        Args:
            configuration: Configuration object
            timestep: Current timestep
            dt: Time step size
            integration_method: Integration method for heredity integrals
            survival_function_computation: Survival computation strategy
            equilibrium_solver_params: Parameters for equilibrium solver
                - solver_method: Root-finding method (default: 'brentq')
                - tolerance: Equilibrium tolerance (default: 1e-5)
                - verbose: Print equilibrium details (default: False)
            
        Returns:
            Dictionary with:
                - 'converged': bool
                - 'iterations': int
                - 'final_residual': float (max Δρ/ρ)
                - 'rho_history': List[Dict[str, float]] (density per layer per iteration)
                - 'message': str
        """
        # Default equilibrium solver parameters
        eq_params = equilibrium_solver_params or {}
        solver_method = eq_params.get('solver_method', 'brentq')
        eq_tolerance = eq_params.get('tolerance', 1e-5)
        eq_verbose = eq_params.get('verbose', False)
        
        self.logger.info("Starting fixed-point iteration")
        
        # Initialize history
        rho_history = []
        max_relative_change = 0.0
        
        # =====================================================================
        # STEP 0: Initial Solve (Better Starting Point)
        # =====================================================================
        if self.perform_initial_solve:
            self.logger.debug("Performing initial solve")
            
            # Compute mass with guessed geometry/stress
            configuration.compute_all_rhoR(
                timestep,
                dt,
                integration_method,
                survival_function_computation
            )
            
            # Solve equilibrium with computed mass
            result = configuration.solve_equilibrium_geometry(
                timestep=timestep,
                dt=dt,
                integration_method=integration_method,
                survival_function_computation=survival_function_computation,
                solver_method=solver_method,
                tolerance=eq_tolerance,
                verbose=eq_verbose
            )
            
            if not result['all_converged']:
                return {
                    'converged': False,
                    'iterations': 0,
                    'final_residual': None,
                    'rho_history': rho_history,
                    'message': 'Initial equilibrium solve failed'
                }
            
            # Store initial densities
            initial_densities = {
                layer.name: layer.get_density(timestep)
                for layer in configuration.layers
            }
            rho_history.append(initial_densities)
            
            self.logger.debug("Initial solve complete")
        
        # =====================================================================
        # STEP 1-N: Fixed-Point Iteration Loop
        # =====================================================================
        iteration = 0
        converged = False
        
        while iteration < self.max_iterations:
            iteration += 1
            
            self.logger.debug(f"Fixed-point iteration {iteration}")
            
            # Store old densities BEFORE updating
            rho_old = {}
            for layer in configuration.layers:
                rho_old[layer.name] = layer.get_density(timestep)
            
            # STEP A: Compute mass densities (heredity integrals)
            configuration.compute_all_rhoR(
                timestep,
                dt,
                integration_method,
                survival_function_computation
            )
            
            # STEP B: Solve geometric equilibrium
            #TODO: consider calling solver here. Then solver should call configuration.solve_equilibrium_geometry() and pass in the parameters.
            result = configuration.solve_equilibrium_geometry(
                timestep=timestep,
                dt=dt,
                integration_method=integration_method,
                survival_function_computation=survival_function_computation,
                solver_method=solver_method,
                tolerance=eq_tolerance,
                verbose=eq_verbose
            )
            
            if not result['all_converged']:
                return {
                    'converged': False,
                    'iterations': iteration,
                    'final_residual': None,
                    'rho_history': rho_history,
                    'message': f"Equilibrium solver failed at iteration {iteration}"
                }
            
            # STEP C: Check convergence
            all_converged, max_change, layer_changes = \
                self.convergence_checker.check_convergence(
                    configuration, timestep, rho_old
                )
            
            # Store iteration history
            current_densities = {
                layer.name: layer.get_density(timestep)
                for layer in configuration.layers
            }
            rho_history.append(current_densities)
            max_relative_change = max_change
            
            # STEP D: Check stopping criteria
            if all_converged:
                converged = True
                self.logger.info(
                    f"✓ Fixed-point solver converged in {iteration} iterations "
                    f"(max Δρ/ρ = {max_change:.2e})"
                )
                break
        
        # =====================================================================
        # FINALIZE: Return Results
        # =====================================================================
        if not converged:
            message = (
                f"Fixed-point iteration did not converge in {self.max_iterations} iterations. "
                f"Max Δρ/ρ = {max_relative_change:.2e} > {self.tolerance:.2e}"
            )
            self.logger.warning(f"⚠️  {message}")
            
            return {
                'converged': False,
                'iterations': iteration,
                'final_residual': max_relative_change,
                'rho_history': rho_history,
                'message': message
            }
        
        return {
            'converged': True,
            'iterations': iteration,
            'final_residual': max_relative_change,
            'rho_history': rho_history,
            'message': f'Converged in {iteration} iterations'
        }