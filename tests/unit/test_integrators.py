"""Unit tests for numerical integrators.

This module tests all integrator classes using analytical test functions
where exact solutions are known. The framework is designed to be extensible:
adding a new integrator only requires creating a new test class.

Test Strategy:
- Analytical functions with known integrals
- Multiple integration bounds (including edge cases)
- Multiple time steps to verify convergence
- Parametrized tests for comprehensive coverage
"""
import pytest
import numpy as np
import math
from typing import Callable, Tuple, List, Dict, Any

from src.integrators import (
    TrapezoidIntegrator,
    SimpsonIntegrator,
    IntegratorFactory,
    Integrator,
)

# =============================================================================
# INTEGRATOR TEST REGISTRY
# =============================================================================

INTEGRATOR_TEST_REGISTRY = {
    'trapezoidal': {
        'class': TrapezoidIntegrator,
        'name': 'trapezoidal',
        'exact_up_to_degree': 1,  # Exact for polynomials of degree ≤ 1
        'min_points': 2,
        'convergence_order': 2,  # O(h²)
    },
    'simpson': {
        'class': SimpsonIntegrator,
        'name': 'simpson',
        'exact_up_to_degree': 3,  # Exact for polynomials of degree ≤ 3
        'min_points': 3,
        'convergence_order': lambda n_intervals=None: (
            4 if (n_intervals is None or n_intervals % 2 == 0)  # Even intervals (odd points): O(h⁴)
            else 2  # Odd intervals (even points): O(h²) - uses trapezoidal on last interval
        )
    },
    # Future integrators can be added here:
    # 'gauss_quadrature': { ... },
}


# =============================================================================
# FIXTURES: Analytical Test Functions
# =============================================================================

@pytest.fixture
def polynomial_functions():
    """Polynomial functions for exactness testing.
    
    Returns functions of increasing degree with their exact integrals.
    Used to test that integrators are exact up to their specified degree.
    """
    return [
        {
            'degree': 0,
            'name': 'constant',
            'f': lambda x: 5.0,
            'F': lambda a, b: 5.0 * (b - a),
            'description': 'f(x) = 5',
        },
        {
            'degree': 1,
            'name': 'linear',
            'f': lambda x: 3.0 * x + 2.0,
            'F': lambda a, b: 1.5 * (b**2 - a**2) + 2.0 * (b - a),
            'description': 'f(x) = 3x + 2',
        },
        {
            'degree': 2,
            'name': 'quadratic',
            'f': lambda x: x**2 - 2*x + 1,
            'F': lambda a, b: (b**3 - a**3)/3 - (b**2 - a**2) + (b - a),
            'description': 'f(x) = x² - 2x + 1',
        },
        {
            'degree': 3,
            'name': 'cubic',
            'f': lambda x: 2*x**3 - x**2 + 3*x - 1,
            'F': lambda a, b: (
                0.5 * (b**4 - a**4) 
                - (b**3 - a**3) / 3 
                + 1.5 * (b**2 - a**2) 
                - (b - a)
            ),
            'description': 'f(x) = 2x³ - x² + 3x - 1',
        },
        {
            'degree': 4,
            'name': 'quartic',
            'f': lambda x: x**4,
            'F': lambda a, b: (b**5 - a**5) / 5,
            'description': 'f(x) = x⁴',
        },
    ]


@pytest.fixture
def smooth_functions():
    """Smooth non-polynomial functions for convergence testing."""
    return [
        {
            'name': 'sine',
            'f': lambda x: math.sin(x),
            'F': lambda a, b: -math.cos(b) + math.cos(a),
            'description': 'f(x) = sin(x)',
        },
        {
            'name': 'cosine',
            'f': lambda x: math.cos(x),
            'F': lambda a, b: math.sin(b) - math.sin(a),
            'description': 'f(x) = cos(x)',
        },
        {
            'name': 'exponential',
            'f': lambda x: math.exp(x),
            'F': lambda a, b: math.exp(b) - math.exp(a),
            'description': 'f(x) = e^x',
        },
    ]


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def discretize_function(f: Callable, a: float, b: float, n_intervals: int) -> Tuple[List[float], float, int, int]:
    """Discretize a continuous function over [a, b] with n_intervals.
    
    Args:
        f: Function to discretize
        a: Lower bound
        b: Upper bound
        n_intervals: Number of intervals (n_points = n_intervals + 1)
    
    Returns:
        Tuple of (function_values, dt, start_index, stop_index)
        
    Note:
        - n_intervals determines the discretization
        - n_points = n_intervals + 1
        - Even n_intervals → odd n_points (good for Simpson's: O(h⁴))
        - Odd n_intervals → even n_points (Simpson's degrades to O(h²))
    """
    # Compute dt from number of intervals (no rounding errors!)
    dt = (b - a) / n_intervals
    
    # Create grid points
    n_points = n_intervals + 1
    x_values = [a + i * dt for i in range(n_points)]
    
    # Evaluate function at grid points
    f_values = [f(x) for x in x_values]
    
    # For backward integration: start = last index, stop = 0
    start_idx = len(f_values) - 1
    stop_idx = 0
    
    return f_values, dt, start_idx, stop_idx


def get_integrator_metadata(integrator_name: str) -> Dict[str, Any]:
    """Get metadata for an integrator from the registry.
    
    Args:
        integrator_name: Name of integrator (e.g., 'trapezoidal', 'simpson')
    
    Returns:
        Dictionary of integrator metadata
    
    Raises:
        KeyError: If integrator not in registry
    """
    if integrator_name not in INTEGRATOR_TEST_REGISTRY:
        available = ', '.join(INTEGRATOR_TEST_REGISTRY.keys())
        raise KeyError(
            f"Integrator '{integrator_name}' not in test registry. "
            f"Available: {available}"
        )
    return INTEGRATOR_TEST_REGISTRY[integrator_name]

def get_convergence_order(integrator_name: str, n_intervals: int = None) -> int:
    """Get convergence order for an integrator.
    
    Args:
        integrator_name: Name of integrator
        n_intervals: Number of intervals (for variable-order methods)
    
    Returns:
        Convergence order (e.g., 2 for O(h²), 4 for O(h⁴))
    """
    metadata = get_integrator_metadata(integrator_name)
    order = metadata['convergence_order']
    
    # Handle callable convergence order (e.g., Simpson's varies with n_intervals)
    if callable(order):
        return order(n_intervals)
    return order

def compute_expected_error(integrator_name: str, dt: float, 
                          interval_length: float, n_intervals: int = None) -> float:
    """Compute expected integration error.
    
    Args:
        integrator_name: Name of integrator
        dt: Time step size
        interval_length: Length of integration interval
        n_intervals: Number of intervals (for variable-order methods)
    
    Returns:
        Expected error magnitude: interval_length * dt^order
    
    Note:
        Error scales as O(h^p) where p is the convergence order.
        The interval_length factor accounts for error accumulation.
    """
    order = get_convergence_order(integrator_name, n_intervals)
    return interval_length * dt**order


def create_integrator(integrator_name: str, dt: float, start: int, stop: int) -> Integrator:
    """Create an integrator instance from the registry.
    
    Args:
        integrator_name: Name of integrator
        dt: Time step
        start: Start index
        stop: Stop index
    
    Returns:
        Integrator instance
    """
    metadata = get_integrator_metadata(integrator_name)
    integrator_class = metadata['class']
    return integrator_class(dt=dt, start=start, stop=stop)



# =============================================================================
# TEST CATEGORY 1: EXACTNESS TESTS
# =============================================================================

class TestExactness:
    """Test that integrators are exact for polynomials up to specified degree.
    
    Each integrator should integrate polynomials exactly up to a certain degree:
    - Trapezoidal: exact for degree ≤ 1
    - Simpson's: exact for degree ≤ 3
    """
    
    @pytest.mark.parametrize("integrator_name", INTEGRATOR_TEST_REGISTRY.keys())
    def test_exactness_for_polynomials(self, integrator_name, polynomial_functions):
        """Test that integrator is exact for polynomials up to its degree.
        
        For each integrator, test all polynomials up to and including
        its exactness degree. Should be exact within floating-point precision.
        """
        metadata = get_integrator_metadata(integrator_name)
        exact_degree = metadata['exact_up_to_degree']
        min_points = metadata['min_points']
        
        # Test all polynomials up to and including exact degree
        for func_dict in polynomial_functions:
            if func_dict['degree'] > exact_degree:
                continue  # Skip polynomials beyond exactness degree
            
            f = func_dict['f']
            F = func_dict['F']
            degree = func_dict['degree']
            
            # Use fine discretization
            a, b = 0.0, 1.0
            n_intervals = 100
            tolerance = 10e-12
            
            f_values, dt, start_idx, stop_idx = discretize_function(f, a, b, n_intervals)
            n_points = len(f_values)
            
            if n_points < min_points:
                pytest.skip(f"Only {n_points} points, need {min_points}")
            
            # Integrate
            integrator = create_integrator(integrator_name, dt, start_idx, stop_idx)
            result = integrator.integrate(f_values)
            
            # Compare to analytical
            expected = F(a, b)
            error = abs(result - expected)
            
            # Should be exact (within machine precision)
            assert error < tolerance, (
                f"{metadata['name']} should be exact for degree {degree} polynomial:\n"
                f"  Function: {func_dict['description']}\n"
                f"  Intervals: {n_intervals} (points: {n_points})\n"
                f"  Expected: {expected:.15f}\n"
                f"  Got:      {result:.15f}\n"
                f"  Error:    {error:.2e}\n"
                f"  Tolerance: {tolerance:.2e}\n"
                f"  (Exactness guaranteed up to degree {exact_degree})"
            )


# =============================================================================
# TEST CATEGORY 2: CONVERGENCE TESTS
# =============================================================================

class TestConvergence:
    """Test that integrators achieve expected convergence rates.
    
    Verifies that error decreases at the theoretical rate O(h^p) as dt → 0,
    where p is the convergence order.
    """
    @pytest.mark.parametrize("integrator_name", INTEGRATOR_TEST_REGISTRY.keys())
    def test_convergence_rate(self, integrator_name, smooth_functions):
        """Test convergence rate for smooth non-polynomial functions.
        
        For a given integrator:
        1. Compute errors at multiple dt values
        2. Verify error decreases monotonically
        3. Verify convergence rate matches theory: error(h/2) / error(h) ≈ 2^order
        
        Note: This test naturally produces odd number of points with the chosen
        dt values and interval [0, π].
        """
        metadata = get_integrator_metadata(integrator_name)
        min_points = metadata['min_points']
        
        # Use smooth function (exp) for convergence testing
        func_dict = smooth_functions[2]  # exp
        f = func_dict['f']
        F = func_dict['F']
        name = func_dict['name']
        
        # Integration interval
        a, b = 0.0, math.pi
        expected_integral = F(a, b)
        
        # Multiple refinement levels
        n_intervals_values = [32, 64, 128, 256]
        errors = []
        dt_values = []
        
        for n_intervals in n_intervals_values:
            f_values, dt, start_idx, stop_idx = discretize_function(f, a, b, n_intervals)
            n_points = len(f_values)
            
            if n_points < min_points:
                continue
            
            integrator = create_integrator(integrator_name, dt, start_idx, stop_idx)
            result = integrator.integrate(f_values)
            error = abs(result - expected_integral)
            
            errors.append(error)
            dt_values.append(dt)
        
        # Test 1: Error should decrease monotonically
        for i in range(len(errors) - 1):
            n_intervals = n_intervals_values[i]
            n_points = n_intervals + 1
            parity = "odd" if n_points % 2 == 1 else "even"
            
            assert errors[i] > errors[i + 1], (
                f"{metadata['name']} error did not decrease monotonically:\n"
                f"  n_intervals={n_intervals_values[i]}: error={errors[i]:.3e} "
                f"(n_points={n_points}, {parity}, dt={dt_values[i]:.6f})\n"
                f"  n_intervals={n_intervals_values[i+1]}: error={errors[i+1]:.3e} "
                f"(n_points={n_intervals_values[i+1]+1}, "
                f"{'odd' if (n_intervals_values[i+1]+1) % 2 == 1 else 'even'}, "
                f"dt={dt_values[i+1]:.6f})\n"
                f"  Function: {name}"
            )
        
        # Test 2: Convergence rate should match theory
        tolerance = 0.01
        for i in range(len(errors) - 1): 
            dt_ratio = dt_values[i] / dt_values[i + 1]
            error_ratio = errors[i] / errors[i + 1]
            observed_order = math.log(error_ratio) / math.log(dt_ratio)
            expected_order = get_convergence_order(integrator_name, n_intervals_values[i])
            
            assert abs(observed_order - expected_order) < tolerance, (
                f"{metadata['name']} convergence order does not match theory:\n"
                f"  n_intervals={n_intervals_values[i]} → {n_intervals_values[i+1]} "
                f"(n_points={n_points}, {parity}, dt={dt_values[i]:.6f} → {dt_values[i+1]:.6f})\n"
                f"  Expected order: O(h^{expected_order})\n"
                f"  Observed order: O(h^{observed_order:.2f})\n"
                f"  Error ratio: {error_ratio:.2f} (expected ≈ {dt_ratio**expected_order:.2f})\n"
                f"  Errors: {errors[i]:.3e} → {errors[i+1]:.3e}\n"
                f"  Deviation: {abs(observed_order - expected_order):.2f} (tolerance: {tolerance})\n"
                f"  Function: {name}\n"
                f"  Note: Convergence rates are asymptotic; coarser grids may deviate"
            )


# TODO: Currently Simpson's rule is not tested for odd number of intervals

# =============================================================================
# TEST CATEGORY 3: EDGE CASE TESTS
# =============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    @pytest.mark.parametrize("integrator_name", INTEGRATOR_TEST_REGISTRY.keys())
    def test_insufficient_points_raises_error(self, integrator_name):
        """Test that integrator raises error with too few points."""
        dt = 0.1
        f_values = [1.0]  
        
        integrator = create_integrator(integrator_name, dt, start=0, stop=0)
        
        with pytest.raises(ValueError, match="At least .* points"):
            integrator.integrate(f_values)
    
    @pytest.mark.parametrize("integrator_name", INTEGRATOR_TEST_REGISTRY.keys())
    def test_zero_length_interval(self, integrator_name):
        """Test handling of zero-length interval (start == stop)."""
        dt = 0.1
        f_values = [1.0, 2.0, 3.0]
        
        integrator = create_integrator(integrator_name, dt, start=1, stop=1)
        with pytest.raises(ValueError, match="At least .* points"):
            integrator.integrate(f_values)
    

# =============================================================================
# TEST CATEGORY 5: FACTORY TESTS
# =============================================================================

class TestIntegratorFactory:
    """Test the IntegratorFactory."""
    
    @pytest.mark.parametrize("integrator_name", INTEGRATOR_TEST_REGISTRY.keys())
    def test_factory_creates_correct_class(self, integrator_name):
        """Test that factory creates correct integrator class."""
        metadata = get_integrator_metadata(integrator_name)
        expected_class = metadata['class']
        
        integrator = IntegratorFactory.create(integrator_name, dt=0.1, start=10, stop=0)
        
        assert isinstance(integrator, expected_class), (
            f"Factory should create {expected_class.__name__}, "
            f"got {type(integrator).__name__}"
        )
        assert integrator.dt == 0.1
        assert integrator.start == 10
        assert integrator.stop == 0
    
    def test_factory_unknown_method_raises_error(self):
        """Test that unknown method raises helpful error."""
        with pytest.raises(ValueError, match="Unknown integration method"):
            IntegratorFactory.create('unknown_method', dt=0.1, start=10, stop=0)
    
    def test_factory_lists_available_methods(self):
        """Test that factory provides list of available methods."""
        methods = IntegratorFactory.available_methods()
        
        # Should include all registered methods
        for integrator_name in INTEGRATOR_TEST_REGISTRY.keys():
            assert integrator_name in methods, (
                f"Factory should list '{integrator_name}' as available method"
            )
















