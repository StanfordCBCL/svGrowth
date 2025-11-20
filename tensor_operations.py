import numpy as np
from numba import njit
from typing import Tuple

# =============================================================================
# CORE TENSOR OPERATIONS (NumPy)
# =============================================================================

class TensorOperations:
    """Pure NumPy tensor operations.
    
    These are the building blocks used throughout the code.
    """

    @staticmethod
    def trace(A: np.ndarray) -> float:
        """Compute trace: tr(A) = sum of diagonal elements."""
        return float(np.trace(A))
        
    @staticmethod
    def right_cauchy_green(F: np.ndarray) -> np.ndarray:
        """C = F^T @ F"""
        return F.T @ F
    
    @staticmethod
    def left_cauchy_green(F: np.ndarray) -> np.ndarray:
        """B = F @ F^T"""
        return F @ F.T
    
    @staticmethod
    def green_strain(F: np.ndarray) -> np.ndarray:
        """E = 0.5 (C - I)"""
        C = F.T @ F
        I = np.eye(F.shape[0])
        return 0.5 * (C - I)
    
    @staticmethod
    def principal_stretches(F: np.ndarray) -> Tuple[float, float, float]:
        """Compute principal stretches from F.
        
        For diagonal F (as in thin-wall), just extract diagonal.
        For general F, use eigenvalue decomposition of C.
        """
        if np.allclose(F - np.diag(np.diagonal(F)), 0):
            # F is diagonal - just extract
            return F[0, 0], F[1, 1], F[2, 2]
        else:
            # F is not diagonal - eigenvalue decomposition
            C = F.T @ F
            eigenvals = np.linalg.eigvalsh(C)  # eigvalsh is faster for symmetric
            return tuple(np.sqrt(eigenvals))
    
    @staticmethod
    def invariants(C: np.ndarray) -> Tuple[float, float, float]:
        """Compute strain invariants I1, I2, I3.
        
        I1 = tr(C)
        I2 = 0.5 * (I1^2 - tr(C^2))
        I3 = det(C)
        """
        I1 = np.trace(C)
        I2 = 0.5 * (I1**2 - np.trace(C @ C))
        I3 = np.linalg.det(C)
        return I1, I2, I3
    
    @staticmethod
    def jacobian(F: np.ndarray) -> float:
        """J = det(F)"""
        return np.linalg.det(F)