"""
Mechanics module for G&R simulations.

Handles:
1. Stress computation for constituents (σ_hat_α)
2. Stress computation for layers (σ_total = Σ σ_α)
3. Numerical integration of heredity integrals (for constituent stress)

Follows same adapter pattern as Kinetics:
- Mechanics class is stateless (no data storage)
- Operates on MechanicsContext (adapter for data access)
- Layer/Constituent own the stress histories
"""

from abc import ABC, abstractmethod
from typing import List, Optional
import numpy as np
from mechanics_interface import MechanicsContext
from tensor_operations import TensorOperations

# =============================================================================
# MAIN MECHANICS CLASS - Stateless stress computations
# =============================================================================

class Mechanics:
    """Main mechanics orchestrator for stress computations.
    
    Stateless - operates on MechanicsContext, never stores data.
    Results are returned to caller (Constituent/Layer) for storage.
    
    Uses TensorOperations for all tensor algebra.
    """
    
    def __init__(self):
        """Initialize with tensor operations utility."""
        self.tensor_ops = TensorOperations()
       
    def compute_stress_trace(self, context: MechanicsContext, timestep: int) -> float:
        """Compute trace of stress tensor (for intramural stress).
        
        Delegates to TensorOperations.
        
        Args:
            stress: Stress tensor (3×3)
            
        Returns:
            tr(σ) = σ_rr + σ_θθ + σ_zz
        """
        stress_tensor = context.get_stress_tensor(timestep)
        return self.tensor_ops.trace(stress_tensor)