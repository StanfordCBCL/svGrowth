"""
Mechanics interface - Adapter layer for mechanics computations.

Provides concrete implementations of MechanicsContext for:
- Constituent-level stress access
- Layer-level stress access
"""

from abc import ABC, abstractmethod
import numpy as np


class MechanicsDataNotAvailableError(Exception):
    """Raised when mechanics requests data that is not available."""
    pass


class MechanicsContext(ABC):
    """Abstract interface for accessing mechanical data in stress computations.
    
    Decouples Mechanics from Layer/Constituent data structures.
    """
    
    @abstractmethod
    def get_stress_tensor(self, timestep: int) -> np.ndarray:
        """Get stress tensor at timestep (3×3 array)."""
        pass
    
    @abstractmethod
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        """Get deformation gradient F at timestep (3×3 array)."""
        pass
    
    @abstractmethod
    def get_mass_density(self, timestep: int) -> float:
        """Get mass density at timestep (kg/m³)."""
        pass


class ConstituentMechanicsContext(MechanicsContext):
    """Concrete implementation for constituent-level stress access."""
    
    def __init__(self, constituent):
        """Initialize with constituent reference."""
        self.constituent = constituent
        self.layer = constituent.layer
        
        if self.layer is None:
            raise ValueError(
                f"Constituent '{constituent.name}' has no layer reference. "
                "Add constituent to layer before creating mechanics context."
            )
    
    def get_stress_tensor(self, timestep: int) -> np.ndarray:
        """Get constituent stress tensor (3×3 array)."""
        if not hasattr(self.constituent, 'stress_history'):
            raise MechanicsDataNotAvailableError(
                f"Constituent '{self.constituent.name}' has no stress_history"
            )
        
        if not (0 <= timestep < len(self.constituent.stress_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range [0, {len(self.constituent.stress_history)-1}]"
            )
        
        return self.constituent.stress_history[timestep]
    
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        """Get F from layer (constituents deform with layer)."""
        if not self.layer or not hasattr(self.layer, 'F_history'):
            raise MechanicsDataNotAvailableError("Layer F_history not available")
        
        if not (0 <= timestep < len(self.layer.F_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range [0, {len(self.layer.F_history)-1}]"
            )
        
        return self.layer.F_history[timestep]
    
    def get_mass_density(self, timestep: int) -> float:
        """Get constituent mass density (kg/m³)."""
        if not hasattr(self.constituent, 'rhoR_alpha_history'):
            raise MechanicsDataNotAvailableError("Constituent mass density not available")
        
        if not (0 <= timestep < len(self.constituent.rhoR_alpha_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range for mass density"
            )
        
        return self.constituent.rhoR_alpha_history[timestep]


class LayerMechanicsContext(MechanicsContext):
    """Concrete implementation for layer-level stress access."""
    
    def __init__(self, layer):
        """Initialize with layer reference."""
        self.layer = layer
    
    def get_stress_tensor(self, timestep: int) -> np.ndarray:
        """Get layer stress tensor (3×3 array)."""
        if not hasattr(self.layer, 'stress_history'):
            raise MechanicsDataNotAvailableError("Layer has no stress_history")
        
        if not (0 <= timestep < len(self.layer.stress_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range [0, {len(self.layer.stress_history)-1}]"
            )
        
        return self.layer.stress_history[timestep]
    
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        """Get layer F (3×3 array)."""
        if not hasattr(self.layer, 'F_history'):
            raise MechanicsDataNotAvailableError("Layer has no F_history")
        
        if not (0 <= timestep < len(self.layer.F_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range [0, {len(self.layer.F_history)-1}]"
            )
        
        return self.layer.F_history[timestep]
    
    def get_mass_density(self, timestep: int) -> float:
        """Get layer mass density (kg/m³)."""
        if not hasattr(self.layer, 'rhoR_history'):
            raise MechanicsDataNotAvailableError("Layer has no rhoR_history")
        
        if not (0 <= timestep < len(self.layer.rhoR_history)):
            raise MechanicsDataNotAvailableError(
                f"Timestep {timestep} out of range [0, {len(self.layer.rhoR_history)-1}]"
            )
        
        return self.layer.rhoR_history[timestep]