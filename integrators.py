from typing import List, Optional
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
        n = self.start - self.stop + 1 
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
        n = self.start - self.stop + 1 # number of integration points

        if n < 2:
            raise ValueError("At least two points required for integration with Simpson's rule.")
        
        if n == 2:
            trapezoid_rule = TrapezoidIntegrator(self.dt, self.start, self.stop)
            return trapezoid_rule.integrate(f)

        integral_value = 0.0
        simpson_stop = self.stop  # Stop index for the main Simpson's rule loop

        # If number of integration points is even, apply trapezoidal rule on last interval
        if n % 2 == 0:
            trapezoid_rule = TrapezoidIntegrator(self.dt, self.stop + 1, self.stop)
            integral_value += trapezoid_rule.integrate(f)
            simpson_stop += 1  # Exclude last interval from Simpson’s rule

        # Add first and last terms
        integral_value += f[self.start] + f[simpson_stop]

        # Add 4 * f(odd indices)
        for i in range(self.start - 1, simpson_stop, -2):
            integral_value += 4 * f[i]

        # Add 2 * f(even indices)
        for i in range(self.start - 2, simpson_stop + 1, -2):
            integral_value += 2 * f[i]

        integral_value *= self.dt / 3
        return integral_value