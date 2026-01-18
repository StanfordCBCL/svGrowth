import numpy as np
from numba import njit
from typing import Tuple

# =============================================================================
# CORE TENSOR OPERATIONS (NumPy)
# =============================================================================

class TensorOperations:
    """Tensor operations for 3x3 matrices using NumPy.
    """

    @staticmethod
    def is_diagonal_3x3(A: np.ndarray) -> bool:
        """Check if 3x3 matrix is diagonal."""
        return A[0, 1] == 0 and A[0, 2] == 0 and A[1, 0] == 0 and \
               A[1, 2] == 0 and A[2, 0] == 0 and A[2, 1] == 0

    @staticmethod
    def inverse(A: np.ndarray) -> np.ndarray:
        """Compute inverse of 3x3 matrix."""
        if TensorOperations.is_diagonal_3x3(A):
            inv_A = np.zeros((3, 3))
            inv_A[0, 0] = 1.0 / A[0, 0]
            inv_A[1, 1] = 1.0 / A[1, 1]
            inv_A[2, 2] = 1.0 / A[2, 2]
            return inv_A
        else:
            return np.linalg.inv(A)
    
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
    def jacobian(F: np.ndarray) -> float:
        """J = det(F)"""
        return np.linalg.det(F)
    
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