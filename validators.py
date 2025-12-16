"""
Validation utilities for svGrowth G&R simulations.

Provides physics-based validation with configurable severity levels.
"""

import logging
from enum import Enum
from typing import Optional
from exceptions import ValidationError as ValidationException

logger = logging.getLogger(__name__)


class ValidationSeverity(Enum):
    """Severity level for validation checks."""
    WARNING = "warning"  # Log warning, continue execution
    ERROR = "error"      # Raise exception, stop execution


class Validator:
    """Physics-based validation for G&R simulations.
    
    All methods are static for easy use without instantiation.
    Validation failures can either warn (log) or error (raise exception).
    """
    
    # ========================================================================
    # Basic Value Validation
    # ========================================================================
    
    @staticmethod
    def validate_positive(value: float, name: str, 
                         severity: ValidationSeverity = ValidationSeverity.ERROR) -> bool:
        """Validate that value is positive (> 0).
        
        Args:
            value: Value to validate
            name: Descriptive name for error messages
            severity: WARNING or ERROR
            
        Returns:
            True if valid, False if invalid (only when severity=WARNING)
            
        Raises:
            ValidationException: If invalid and severity=ERROR
        """
        if value <= 0:
            msg = f"{name} must be positive, got {value}"
            
            if severity == ValidationSeverity.ERROR:
                logger.error(msg)
                raise ValidationException(msg)
            else:
                logger.warning(msg)
                return False
        
        return True
    
    @staticmethod
    def validate_range(value: float, name: str,
                      min_val: Optional[float] = None,
                      max_val: Optional[float] = None,
                      severity: ValidationSeverity = ValidationSeverity.ERROR) -> bool:
        """Validate that value is within specified range.
        
        Args:
            value: Value to validate
            name: Descriptive name for error messages
            min_val: Minimum allowed value (inclusive)
            max_val: Maximum allowed value (inclusive)
            severity: WARNING or ERROR
            
        Returns:
            True if valid, False if invalid (only when severity=WARNING)
            
        Raises:
            ValidationException: If invalid and severity=ERROR
        """
        if min_val is not None and value < min_val:
            msg = f"{name} = {value} is below minimum {min_val}"
            
            if severity == ValidationSeverity.ERROR:
                logger.error(msg)
                raise ValidationException(msg)
            else:
                logger.warning(msg)
                return False
        
        if max_val is not None and value > max_val:
            msg = f"{name} = {value} exceeds maximum {max_val}"
            
            if severity == ValidationSeverity.ERROR:
                logger.error(msg)
                raise ValidationException(msg)
            else:
                logger.warning(msg)
                return False
        
        return True
    
    # ========================================================================
    # Geometry Validation
    # ========================================================================
    
    @staticmethod
    def validate_geometry(a: float, h: float, name: str = "") -> bool:
        """Validate geometric parameters.
        
        Checks:
        - Inner radius is positive
        - Thickness is positive
        - Thin-wall assumption (h/a < 0.2) holds
        
        Args:
            a: Inner radius (m)
            h: Thickness (m)
            name: Optional prefix for error messages
            
        Returns:
            True if valid
            
        Raises:
            ValidationException: If geometry is invalid
        """
        prefix = f"{name} " if name else ""
        
        # Positive values
        Validator.validate_positive(a, f"{prefix}inner_radius")
        Validator.validate_positive(h, f"{prefix}thickness")
        
        # Thin-wall assumption
        ratio = h / a
        if ratio > 0.2:
            logger.warning(
                f"{prefix}h/a ratio = {ratio:.3f} > 0.2. "
                f"Thin-wall assumption may not hold."
            )
        
        return True
    
    # ========================================================================
    # Physics Validation
    # ========================================================================
    
    @staticmethod
    def validate_incompressibility(J: float,
                                   warn_tolerance: float = 0.05,
                                   error_tolerance: float = 0.20) -> bool:
        """Validate incompressibility constraint: J ≈ 1.
        
        Args:
            J: Volume ratio (det(F) or ρ_h/ρ)
            warn_tolerance: Threshold for warning (deviation from 1.0)
            error_tolerance: Threshold for error (deviation from 1.0)
            
        Returns:
            True if valid
            
        Raises:
            ValidationException: If deviation exceeds error_tolerance
        """
        deviation = abs(J - 1.0)
        
        # Warning for moderate deviation
        if deviation > warn_tolerance:
            logger.warning(
                f"Incompressibility deviation: J = {J:.6f} "
                f"(deviation = {deviation*100:.2f}%, warn threshold = {warn_tolerance*100:.1f}%)"
            )
        
        # Error for severe deviation
        if deviation > error_tolerance:
            msg = (
                f"Severe incompressibility violation: J = {J:.6f} "
                f"(deviation = {deviation*100:.2f}%, error threshold = {error_tolerance*100:.1f}%)"
            )
            logger.error(msg)
            raise ValidationException(msg)
        
        return True
    
    @staticmethod
    def validate_stress(sigma: float, name: str,
                       warn_threshold: float = 500e3,
                       error_threshold: float = 1000e3) -> bool:
        """Validate stress magnitude is physically reasonable.
        
        Typical arterial stresses: 10-200 kPa
        Warning threshold: 500 kPa
        Error threshold: 1000 kPa (1 MPa)
        
        Args:
            sigma: Stress value (Pa)
            name: Stress component name (e.g., "σ_θ")
            warn_threshold: Threshold for warning (Pa)
            error_threshold: Threshold for error (Pa)
            
        Returns:
            True if valid
            
        Raises:
            ValidationException: If stress exceeds error_threshold
        """
        abs_sigma = abs(sigma)
        
        # Check for warning
        if abs_sigma > warn_threshold:
            logger.warning(
                f"{name} = {sigma/1000:.2f} kPa exceeds typical range "
                f"(warn threshold: {warn_threshold/1000:.0f} kPa)"
            )
        
        # Check for error
        if abs_sigma > error_threshold:
            msg = (
                f"{name} = {sigma/1000:.2f} kPa is unreasonably large "
                f"(error threshold: {error_threshold/1000:.0f} kPa)"
            )
            logger.error(msg)
            raise ValidationException(msg)
        
        return True
    
    @staticmethod
    def validate_strain_constraint(sigma_r: float, 
                                   tolerance: float = 1e3) -> bool:
        """Validate thin-wall kinematic constraint: σ_r ≈ 0.
        
        Args:
            sigma_r: Radial stress component (Pa)
            tolerance: Maximum allowed deviation (Pa)
            
        Returns:
            True if valid
        """
        if abs(sigma_r) > tolerance:
            logger.warning(
                f"Radial stress constraint violated: σ_r = {sigma_r:.2e} Pa "
                f"(should be ≈ 0, tolerance = {tolerance:.2e} Pa)"
            )
        
        return True