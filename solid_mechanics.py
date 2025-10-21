import math
import numpy as np
from abc import ABC, abstractmethod
from enum import Enum

class IsotropyType(Enum):
    """Enumeration for material isotropy types."""
    ISOTROPIC = "isotropic"
    ANISOTROPIC = "anisotropic"

class ConstitutiveModel(ABC):
    """Abstract base class for constitutive models."""
    
    def __init__(self, isotropy_type, parameters):
        self.isotropy_type = IsotropyType(isotropy_type)
        self.parameters = parameters
        self.fiber_orientation = None
        
        # Validate parameters
        self._validate_parameters()
    
    @abstractmethod
    def _validate_parameters(self):
        """Validate model-specific parameters."""
        pass
    
    @abstractmethod
    def compute_stress(self, deformation_gradient, mass_density=1.0):
        """Compute Cauchy stress from deformation gradient."""
        pass
    
    @abstractmethod
    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        pass
    
    @classmethod
    def from_parameters(cls, model_data):
        """Factory method to create constitutive model from parameters."""
        model_type = model_data['type']
        isotropy = model_data['isotropy']
        parameters = model_data['parameters']
        
        # Create appropriate model
        if model_type == "neo_hookean":
            model = NeoHookeanModel(isotropy, parameters)
        elif model_type == "fung_exponential":
            model = FungExponentialModel(isotropy, parameters)
        elif model_type == "holzapfel_ogden":
            model = HolzapfelOgdenModel(isotropy, parameters)
        else:
            raise ValueError(f"Unknown constitutive model type: {model_type}")
        
        # Set fiber orientation if anisotropic
        if isotropy == "anisotropic" and 'fiber_orientation' in model_data:
            angle_degrees = model_data['fiber_orientation']
            model.set_fiber_orientation(angle_degrees)
        
        return model
    
    def set_fiber_orientation(self, angle_degrees):
        """Set fiber orientation for anisotropic materials."""
        if self.isotropy_type != IsotropyType.ANISOTROPIC:
            raise ValueError("Fiber orientation only applies to anisotropic materials")
        
        self.fiber_orientation = math.radians(angle_degrees)
        print(f"        Fiber orientation: {angle_degrees}° → {self.fiber_orientation:.3f} rad")
    
    def is_isotropic(self):
        """Check if material is isotropic."""
        return self.isotropy_type == IsotropyType.ISOTROPIC
    
    def is_anisotropic(self):
        """Check if material is anisotropic."""
        return self.isotropy_type == IsotropyType.ANISOTROPIC

# Registry of available models (for validation and documentation)
AVAILABLE_MODELS = {
    "neo_hookean": {
        "description": "Neo-Hookean model for rubber-like materials",
        "isotropy": ["isotropic"],
        "parameters": ["c"],
        "units": {"c": "Pa"}
    },
    "fung_exponential": {
        "description": "Fung exponential model for biological tissues",
        "isotropy": ["isotropic", "anisotropic"],
        "parameters": ["c1", "c2"],
        "units": {"c1": "Pa", "c2": "dimensionless"}
    },
    "holzapfel_ogden": {
        "description": "Holzapfel-Ogden model for arterial walls",
        "isotropy": ["anisotropic"],
        "parameters": ["c", "k1", "k2"],
        "units": {"c": "Pa", "k1": "Pa", "k2": "dimensionless"}
    }
}

class NeoHookeanModel(ConstitutiveModel):
    """Neo-Hookean constitutive model for isotropic materials."""
    
    def __init__(self, isotropy_type, parameters):
        if isotropy_type != "isotropic":
            raise ValueError("Neo-Hookean model is inherently isotropic")
        
        super().__init__(isotropy_type, parameters)
        self.c = parameters['c']  # Material constant in Pa
    
    def _validate_parameters(self):
        """Validate Neo-Hookean parameters."""
        if 'c' not in self.parameters:
            raise ValueError("Neo-Hookean model requires parameter 'c'")
        if self.parameters['c'] <= 0:
            raise ValueError("Neo-Hookean parameter 'c' must be positive")
    
    def compute_stress(self, deformation_gradient, mass_density=1.0):
        """Compute Cauchy stress for Neo-Hookean material."""
        F = deformation_gradient  # Deformation gradient (assumed scalar for 1D)
        
        # Left Cauchy-Green deformation tensor B = F * F^T
        B = F @ F.T
        
        # Neo-Hookean stress for incompressible material
        stress = self.c * B
       
        return stress
    
    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        F = deformation_gradient
        if isinstance(F, (int, float)):
            # W = c/2 * (F² + 1/F² - 2) for incompressible Neo-Hookean
            energy = 0.5 * self.c * (F*F + 1.0/(F*F) - 2.0)
        else:
            energy = 0.5 * self.c  # Placeholder
        
        return energy
    
    def __str__(self):
        return f"NeoHookean(c={self.c:.1f} Pa, isotropic)"
    
class FungExponentialModel(ConstitutiveModel):
    """Fung exponential constitutive model for anisotropic materials."""
    
    def __init__(self, isotropy_type, parameters):
        super().__init__(isotropy_type, parameters)
        self.c1 = parameters['c1']  # Material constant in Pa
        self.c2 = parameters['c2']  # Dimensionless exponential parameter
    
    def _validate_parameters(self):
        """Validate Fung exponential parameters."""
        required_params = ['c1', 'c2']
        for param in required_params:
            if param not in self.parameters:
                raise ValueError(f"Fung model requires parameter '{param}'")
        
        if self.parameters['c1'] <= 0:
            raise ValueError("Fung parameter 'c1' must be positive")
        if self.parameters['c2'] <= 0:
            raise ValueError("Fung parameter 'c2' must be positive")
    
    def compute_stress(self, deformation_gradient, mass_density=1.0):
        """Compute Cauchy stress for Fung exponential material."""
        F = deformation_gradient
        
        if isinstance(F, (int, float)):
            # Simplified Fung stress (would need proper fiber direction in 3D)
            strain = 0.5 * (F*F - 1.0)  # Green strain
            exponential_term = math.exp(self.c2 * strain * strain)
            stress = self.c1 * mass_density * strain * exponential_term
        else:
            # Placeholder for tensor implementation
            stress = self.c1 * mass_density * 1.2
        
        return stress
    
    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        F = deformation_gradient
        
        if isinstance(F, (int, float)):
            strain = 0.5 * (F*F - 1.0)
            energy = (self.c1 / (2.0 * self.c2)) * (math.exp(self.c2 * strain * strain) - 1.0)
        else:
            energy = self.c1 / (2.0 * self.c2)  # Placeholder
        
        return energy
    
    def __str__(self):
        isotropy_str = self.isotropy_type.value
        if self.fiber_orientation is not None:
            angle_deg = math.degrees(self.fiber_orientation)
            return f"FungExponential(c1={self.c1:.1f} Pa, c2={self.c2:.1f}, {isotropy_str}, θ={angle_deg:.1f}°)"
        else:
            return f"FungExponential(c1={self.c1:.1f} Pa, c2={self.c2:.1f}, {isotropy_str})"