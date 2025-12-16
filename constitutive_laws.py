import math
import numpy as np
from abc import ABC, abstractmethod
from enum import Enum
from tensor_operations import TensorOperations

class IsotropyType(Enum):
    """Enumeration for material isotropy types."""
    ISOTROPIC = "isotropic"
    ANISOTROPIC = "anisotropic"

class ConstitutiveModel(ABC):
    """Abstract base class for constitutive models."""
    
    # Subclasses must define this
    ISOTROPY_TYPE = None
    
    def __init__(self, parameters):
        self.parameters = parameters
        #TODO: move these to Fung explonential only? Not general for all laws.
        self.fiber_orientation = None
        self.fiber_angle_degrees = None
        
        # Set isotropy from class definition
        if self.ISOTROPY_TYPE is None:
            raise NotImplementedError(f"{self.__class__.__name__} must define ISOTROPY_TYPE")
        self.isotropy_type = self.ISOTROPY_TYPE
        
        # Validate parameters
        self._validate_parameters()
    
    @abstractmethod
    def _validate_parameters(self):
        """Validate model-specific parameters."""
        pass
    
    @abstractmethod
    def compute_PK2_stress(self, deformation_gradient):
        """Compute deviatoric 2nd Piola–Kirchhoff stress from deformation gradient."""
        pass
    
    @abstractmethod
    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        pass
    
    @classmethod
    def from_parameters(cls, model_data):
        """Factory method to create constitutive model from parameters."""
        model_type = model_data['type']
        parameters = model_data['parameters']

        kPa_to_Pa = 1000.0
        
        # Create appropriate model
        if model_type == "neo_hookean":
            if 'c' in parameters:
                parameters['c'] = float(parameters['c']) * kPa_to_Pa
            model = NeoHookeanModel(parameters)
        elif model_type == "fung_exponential":
            if 'c1' in parameters:
                parameters['c1'] = float(parameters['c1']) * kPa_to_Pa
            model = FungExponentialModel(parameters)
        elif model_type == "holzapfel_ogden":
            if 'c' in parameters:
                parameters['c'] = float(parameters['c']) * kPa_to_Pa
            if 'k1' in parameters:
                parameters['k1'] = float(parameters['k1']) * kPa_to_Pa
            model = HolzapfelOgdenModel(parameters)
        else:
            raise ValueError(f"Unknown constitutive model type: {model_type}")
        
        # Set fiber orientation if provided (anisotropic models only)
        if 'fiber_orientation' in model_data:
            angle_degrees = model_data['fiber_orientation']
            model.set_fiber_orientation(angle_degrees)
        
        return model
    
    def set_fiber_orientation(self, angle_degrees):
        """Set fiber orientation for anisotropic material.
        
        Overrides base class to store as numpy array (unit vector).
        
        Args:
            angle_degrees: Fiber angle measured from axial direction (degrees)
                          0° = axial, 90° = circumferential
        """
        if self.isotropy_type != IsotropyType.ANISOTROPIC:
            raise ValueError(
                f"{self.__class__.__name__} is {self.isotropy_type.value}. "
                "Fiber orientation only applies to anisotropic materials."
            )
        
        # TODO: to either generalize or add communication to kinematics for coordinate system and 
        # balance equations? For now, assume cylindrical coords (r, θ, z).
        self.fiber_angle_degrees = angle_degrees
        theta = math.radians(angle_degrees)
        
        # Fiber direction in cylindrical coords (r, θ, z)
        # Convention: angle measured from axial (z) direction
        a0 = np.array([
            0.0,              # radial component (none for thin-wall)
            np.sin(theta),    # circumferential component
            np.cos(theta)     # axial component
        ], dtype=float)
        
        # Store as unit vector
        self.fiber_orientation = a0 / np.linalg.norm(a0)
        
        print(f"        Fiber orientation: {angle_degrees}° → "
              f"a0 = [{a0[0]:.3f}, {a0[1]:.3f}, {a0[2]:.3f}]")
    
    def is_isotropic(self):
        """Check if material is isotropic."""
        return self.isotropy_type == IsotropyType.ISOTROPIC
    
    def is_anisotropic(self):
        """Check if material is anisotropic."""
        return self.isotropy_type == IsotropyType.ANISOTROPIC

    # ------------------------------------------------------------------
    # Static helper methods
    # ------------------------------------------------------------------
    
    @staticmethod
    def compute_fiber_stretch(F: np.ndarray, a0: np.ndarray) -> float:
        """Compute fiber stretch λ_f = sqrt(a0 · C · a0).
        
        Args:
            F: Deformation gradient (3×3)
            a0: Fiber direction unit vector (3,)
            
        Returns:
            Fiber stretch λ_f (dimensionless)
        """
        C = TensorOperations.right_cauchy_green(F)
        return math.sqrt(a0 @ C @ a0)
    
    @staticmethod
    def compute_fiber_strain(F: np.ndarray, a0: np.ndarray) -> float:
        """Compute fiber Green-Lagrange strain E_f = a0 · E · a0.
        
        Args:
            F: Deformation gradient (3×3)
            a0: Fiber direction unit vector (3,)
            
        Returns:
            Fiber strain E_f (dimensionless)
        """
        C = TensorOperations.right_cauchy_green(F)
        I = np.eye(3)
        E = 0.5 * (C - I) # TODO: tensor operations function for this?
        return float(a0 @ E @ a0)
    
# Registry of available models (for validation and documentation)
AVAILABLE_MODELS = {
    "neo_hookean": {
        "description": "Neo-Hookean model for rubber-like materials",
        "isotropy": "isotropic",
        "parameters": ["c"],
        "units": {"c": "Pa"}
        #TODO: add checks to perform
    },
    "fung_exponential": {
        "description": "Fung exponential model for biological tissues",
        "isotropy": "anisotropic",
        "parameters": ["c1", "c2"],
        "optional_parameters": ["allow_compression"],
        "required": ["fiber_orientation"],  # Must specify fiber direction
        "units": {"c1": "Pa", "c2": "dimensionless"},
        "defaults": {"allow_compression": False}
    },
    "holzapfel_ogden": {
        "description": "Holzapfel-Ogden model for arterial walls",
        "isotropy": "anisotropic",  
        "parameters": ["c", "k1", "k2"],
        "optional_parameters": ["allow_compression"],
        "required": ["fiber_orientation"],
        "units": {"c": "Pa", "k1": "Pa", "k2": "dimensionless"},
        "defaults": {"allow_compression": False}
    }
}


class NeoHookeanModel(ConstitutiveModel):
    """Neo-Hookean constitutive model for isotropic materials.
    
    This model is inherently isotropic - suitable for elastin and other
    rubber-like biological materials.
    """
    
    ISOTROPY_TYPE = IsotropyType.ISOTROPIC
    
    def __init__(self, parameters):
        super().__init__(parameters)
        self.c = parameters['c']  # Material constant in Pa
    
    def _validate_parameters(self):
        """Validate Neo-Hookean parameters."""
        if 'c' not in self.parameters:
            raise ValueError("Neo-Hookean model requires parameter 'c'")
        if self.parameters['c'] <= 0:
            raise ValueError("Neo-Hookean parameter 'c' must be positive")
    
    def compute_PK2_stress(self, deformation_gradient):
        """Compute deviatoric 2nd Piola–Kirchhoff stress for incompressible Neo-Hookean material with J=1."""
        F = deformation_gradient
        I = np.eye(3)

        # TODO: generalize if needed
        S_pk2 = self.c * I

        return S_pk2

    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        F = deformation_gradient
        energy = 1 # PLACEHOLDER. TODO: Implement properly
        
        return energy
    
    def __str__(self):
        return f"NeoHookean(c={self.c:.1f} Pa, isotropic)"


class FungExponentialModel(ConstitutiveModel):
    """Fung exponential constitutive model for anisotropic materials.
    
    This model is inherently anisotropic - suitable for collagen fibers,
    smooth muscle cells, and other fibrous biological tissues.
    """
    
    ISOTROPY_TYPE = IsotropyType.ANISOTROPIC
    
    def __init__(self, parameters):
        # TODO: Consider restructuring how anisotropic fibers are handled
        # Should it be a parameter in self?
        super().__init__(parameters)
        self.c1 = parameters['c1']  # Material constant in Pa
        self.c2 = parameters['c2']  # Dimensionless exponential parameter
        self.allow_compression = parameters.get('allow_compression', False) # TODO: think of how to handle defaults better
    
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
        
        if 'allow_compression' in self.parameters:
            if not isinstance(self.parameters['allow_compression'], bool):
                raise ValueError("Parameter 'allow_compression' must be boolean (True/False)")
    
    def compute_PK2_stress(self, deformation_gradient: np.ndarray) -> np.ndarray:
        """Compute deviatoric 2nd Piola-Kirchhoff stress for Fung model.
        
        Formula:
            S_pk2 = c1 · c2 · E_f · exp(c2 · E_f²) · (a0 ⊗ a0)
        
        where:
            E_f = a0 · E · a0  (scalar fiber strain)
            a0 = fiber direction in reference configuration
            E = 0.5 (C - I)   (Green-Lagrange strain tensor)
        
        Args:
            deformation_gradient: F (3×3)
            
        Returns:
            S (3×3) - 2nd Piola-Kirchhoff stress tensor (Pa)
        """
        # TODO: revisit variable names for clarity        
        F = deformation_gradient
        a0 = self.fiber_orientation
        
        # Compute fiber strain
        lambda_f = self.compute_fiber_stretch(F, a0)

        # Fiber is compressed - buckles and carries no load
        if not self.allow_compression and lambda_f < 1:
            lambda_f = 1

        # Fiber invariant
        Q1 = lambda_f**2 - 1

        # Exponential term
        Q2 = self.c2 * Q1**2

        # Stress coefficient
        coeff = self.c1 * Q1 * np.exp(Q2)
        
        # Structural tensor: a0 ⊗ a0 (dyadic product)
        A = np.outer(a0, a0)
        
        # 2nd PK stress: S = coeff · (a0 ⊗ a0)
        S_pk2 = coeff * A
        
        return S_pk2
    
    def compute_strain_energy(self, deformation_gradient: np.ndarray) -> float:
        """Compute strain energy density.
        
        Formula:
            Ψ = (c1 / 2) · [exp(c2 · E_f²) - 1]
        
        Args:
            deformation_gradient: F (3×3)
            
        Returns:
            Strain energy density (Pa or J/m³)
        """        
        F = deformation_gradient
        # a0 = self.fiber_orientation
        
        # # Compute fiber strain
        # lambda_f = self.compute_fiber_stretch(F, a0)

        # if not self.allow_compression and lambda_f < 1:
        #     # TODO: confirm correct behavior for strain energy, 0 as placeholder
        #     return 0
                
        # # Strain energy
        # energy = (self.c1 / 2.0) * (np.exp(self.c2 * E_f**2) - 1.0)
        
        energy = 1  # PLACEHOLDER. TODO: Implement properly
        
        return float(energy)
    
    def __str__(self) -> str:
        """String representation showing parameters in Pa."""
        model_str = f"FungExponential(c1={self.c1:.1f} Pa, c2={self.c2:.2f}, {self.ISOTROPY_TYPE.value}"
        
        if self.fiber_angle_degrees is not None:
            model_str += f", fiber={self.fiber_angle_degrees:.1f}°"

        if not self.allow_compression:
            model_str += ", no-compression"

        model_str += ")"
        return model_str


class HolzapfelOgdenModel(ConstitutiveModel):
    """Holzapfel-Ogden constitutive model for arterial walls.
    
    This model is inherently anisotropic - designed specifically for
    fiber-reinforced arterial tissue.
    """
    
    ISOTROPY_TYPE = IsotropyType.ANISOTROPIC
    
    def __init__(self, parameters):
        super().__init__(parameters)
        self.c = parameters['c']
        self.k1 = parameters['k1']
        self.k2 = parameters['k2']
        self.allow_compression = parameters.get('allow_compression', False) # TODO: think of how to handle defaults better
    
    def _validate_parameters(self):
        """Validate Holzapfel-Ogden parameters."""
        required_params = ['c', 'k1', 'k2']
        for param in required_params:
            if param not in self.parameters:
                raise ValueError(f"Holzapfel-Ogden model requires parameter '{param}'")
        
        if self.parameters['c'] <= 0:
            raise ValueError("Holzapfel-Ogden parameter 'c' must be positive")
        if self.parameters['k1'] <= 0:
            raise ValueError("Holzapfel-Ogden parameter 'k1' must be positive")
        if self.parameters['k2'] <= 0:
            raise ValueError("Holzapfel-Ogden parameter 'k2' must be positive")
    
    def compute_stress(self, deformation_gradient, mass_density=1.0):
        """Compute Cauchy stress for Holzapfel-Ogden material."""
        # Placeholder implementation
        return self.c * mass_density
    
    def compute_strain_energy(self, deformation_gradient):
        """Compute strain energy density."""
        # Placeholder implementation
        return self.c
    
    # def __str__(self):
    #     if self.fiber_orientation is not None:
    #         angle_deg = math.degrees(self.fiber_orientation)
    #         return f"HolzapfelOgden(c={self.c:.1f} Pa, k1={self.k1:.1f} Pa, k2={self.k2:.1f}, anisotropic, θ={angle_deg:.1f}°)"
    #     else:
    #         return f"HolzapfelOgden(c={self.c:.1f} Pa, k1={self.k1:.1f} Pa, k2={self.k2:.1f}, anisotropic, θ=not set)"