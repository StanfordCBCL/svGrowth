"""
Custom exceptions for svGrowth G&R simulations.

Provides domain-specific exceptions for clear error handling.
"""

class GrowthRemodelingError(Exception):
    """Base exception for all G&R simulation errors."""
    pass


# ============================================================================
# Geometry Errors
# ============================================================================

class GeometryError(GrowthRemodelingError):
    """Raised for invalid or problematic geometry."""
    pass


class ThinWallViolationError(GeometryError):
    """Raised when thin-wall assumption (h/a < 0.2) is violated."""
    pass


# ============================================================================
# Solver Errors
# ============================================================================

class SolverError(GrowthRemodelingError):
    """Base exception for solver errors."""
    pass


class ConvergenceError(SolverError):
    """Raised when iterative solver fails to converge."""
    pass


class EquilibriumError(SolverError):
    """Raised when equilibrium cannot be satisfied."""
    pass


# ============================================================================
# Physics Errors
# ============================================================================

class PhysicsViolationError(GrowthRemodelingError):
    """Raised when physical constraints are violated."""
    pass


class IncompressibilityError(PhysicsViolationError):
    """Raised when incompressibility constraint (J ≈ 1) is violated."""
    pass


class StressError(PhysicsViolationError):
    """Raised when stress values are physically unreasonable."""
    pass


# ============================================================================
# Data Errors
# ============================================================================

class DataError(GrowthRemodelingError):
    """Base exception for data-related errors."""
    pass


class KineticsDataNotAvailableError(DataError):
    """Raised when kinetics requests unavailable data."""
    pass


class MechanicsDataNotAvailableError(DataError):
    """Raised when mechanics requests unavailable data."""
    pass


# ============================================================================
# Parameter Errors
# ============================================================================

class ParameterError(GrowthRemodelingError):
    """Raised for invalid input parameters."""
    pass


class ValidationError(ParameterError):
    """Raised when parameter validation fails."""
    pass