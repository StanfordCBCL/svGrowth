import math
import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple
import numpy as np
from enum import Enum
from tensor_operations import TensorOperations
from custom_logging import get_logger

# =============================================================================
# COORDINATE SYSTEMS
# =============================================================================

class CoordinateSystem(Enum):
    """Coordinate system types."""
    CARTESIAN = "cartesian"
    CYLINDRICAL = "cylindrical"
    SPHERICAL = "spherical"


# =============================================================================
# DEFORMATION KINEMATICS - Computes F and derived quantities
# =============================================================================

class DeformationKinematics(ABC):
    """Abstract base class for computing deformation gradient and strain measures.
    
    Different implementations:
    - ThinWallKinematics: Simplified thin-wall vessel
    - ThickWallKinematics: Full through-thickness integration
    - FEMKinematics: Finite element-based (future)
    
    All implementations use TensorOperations for consistent tensor algebra.
    """
    
    def __init__(self):
        """Initialize with tensor operations utility."""
        self.logger = get_logger(__name__)
        self.tensor_ops = TensorOperations()
    
    @abstractmethod
    def compute_deformation_gradient(
        self,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> np.ndarray:
        """Compute deformation gradient F.
        
        Args:
            current_geometry: Current configuration (a, h, lambda_z, etc.)
            reference_geometry: Reference configuration (A, H, Lambda_z, etc.)
            
        Returns:
            Deformation gradient tensor F (3x3 for 3D)
        """
        pass
    
    @abstractmethod
    def compute_stress_from_equilibrium(
        self,
        pressure: float,
        inner_radius: float,
        thickness: float
    ) -> float:
        """Compute stress required by equilibrium given geometry and pressure.
        
        This computes the PRIMARY stress components needed for solving equilibrium.
        For cylinder: circumferential stress σ_θ = P·a/h
        
        Args:
            pressure: Intramural (internal) pressure (Pa)
            inner_radius: Current inner radius (m)
            thickness: Current thickness (m)

        Returns:
            Primary stress components from equilibrium.
            For cylinder: {'sigma_theta': σ_θ}
        """
        pass
    
    @abstractmethod
    def compute_wss(
        self,
        flow_rate: float,
        inner_radius: float,
        blood_viscosity: float
    ) -> float:
        """Compute wall shear stress from flow rate.
        
        This is geometry-specific:
        - Thin-wall cylinder: Poiseuille flow τ_w = 4μQ/(πa³)
        
        Args:
            flow_rate: Volumetric flow rate (m³/s)
            inner_radius: Inner radius (m)
            blood_viscosity: Dynamic viscosity (Pa·s)
            
        Returns:
            Wall shear stress (Pa)
        """
        pass

    def compute_stretches(self, F: np.ndarray) -> Tuple[float, float, float]:
        """Extract principal stretches from F.
        
        Delegates to TensorOperations for consistent computation.
        
        Returns:
            (lambda_r, lambda_theta, lambda_z) or (lambda_1, lambda_2, lambda_3)
        """
        return self.tensor_ops.principal_stretches(F)
    
    def compute_right_cauchy_green(self, F: np.ndarray) -> np.ndarray:
        """Compute right Cauchy-Green tensor C = F^T F.
        
        Delegates to TensorOperations.
        """
        return self.tensor_ops.right_cauchy_green(F)
    
    def compute_left_cauchy_green(self, F: np.ndarray) -> np.ndarray:
        """Compute left Cauchy-Green tensor B = F F^T.
        
        Delegates to TensorOperations.
        """
        return self.tensor_ops.left_cauchy_green(F)
    
    def compute_green_strain(self, F: np.ndarray) -> np.ndarray:
        """Compute Green-Lagrange strain E = 0.5 (C - I).
        
        Delegates to TensorOperations.
        """
        return self.tensor_ops.green_strain(F)
    
    def compute_invariants(self, F: np.ndarray) -> Tuple[float, float, float]:
        """Compute strain invariants I1, I2, I3 from F.
        
        Delegates to TensorOperations (computes via C = F^T F).
        
        Returns:
            (I1, I2, I3)
        """
        C = self.compute_right_cauchy_green(F)
        return self.tensor_ops.invariants(C)
    
    def compute_jacobian(self, F: np.ndarray) -> float:
        """Compute Jacobian (volume ratio) J = det(F).
        
        Delegates to TensorOperations.
        """
        return self.tensor_ops.jacobian(F)

    @abstractmethod
    def compute_volume_ratio(
        self,
        F_history: List[np.ndarray],
        rhoR_history: List[float],
        rhoR_h: float,
        timestep: int
    ) -> float:
        """Compute volume ratio J at specified timestep.
        
        Kinematics-specific logic determines whether to use:
        - F (fundamental): J = det(F)
        - Density (constraint): J = ρ_h / ρ
        
        Args:
            F_history: Deformation gradient history
            rhoR_history: Mass density history
            rhoR_h: Homeostatic density
            timestep: Target timestep
            
        Returns:
            Volume ratio J
        """
        pass
        
    def verify_incompressibility(
        self,
        F: np.ndarray,
        tolerance: float = 1e-6
    ) -> bool:
        """Verify incompressibility constraint det(F) = 1.
        
        Args:
            F: Deformation gradient
            tolerance: Tolerance for constraint violation
            
        Returns:
            True if |det(F) - 1| < tolerance
        """
        J = self.compute_jacobian(F)
        return abs(J - 1.0) < tolerance
    
    @abstractmethod
    def compute_lagrange_multiplier(self, sigma_material: np.ndarray) -> np.ndarray:
        """Compute Lagrange multiplier to enforce kinematic constraint.
        
        The Lagrange multiplier is a reaction stress that enforces the
        constraint specific to this kinematics (e.g., σ_r = 0 for thin-wall).
        
        Args:
            sigma_material: Unconstrained material stress tensor (3×3)
            
        Returns:
            Lagrange multiplier tensor (3×3) to be subtracted from σ_material
            
        Note: Layer is responsible for applying the subtraction.
        """
        pass

    @abstractmethod
    def get_component_name(self, index: int) -> str:
        """Get component name for tensor index.
        
        Maps tensor indices to coordinate system component names.
        For example, cylindrical: 0→'r', 1→'theta', 2→'z'
        
        Args:
            index: Tensor component index (0, 1, or 2)
            
        Returns:
            Component name string
            
        Raises:
            ValueError: If index out of range
        """
    pass

    @abstractmethod
    def get_equilibrium_variable(self) -> str:
        #TODO: refactor this to a class?
        """Get name of geometry variable to solve for in equilibrium.
        
        Different kinematics solve for different variables:
        - Thin-wall: 'mid_radius'
        - Thick-wall: 'inner_radius'
        
        Returns:
            Variable name (string)
        """
        pass

    @abstractmethod
    def get_reference_geometry_value(self, layer) -> float:
        """Get reference value for equilibrium variable (for bracketing).
        
        Typically the homeostatic value of the variable.
        
        Args:
            layer: Layer object to get value from
            
        Returns:
            Reference value (float, in meters)
        """
        pass

# =============================================================================
# THIN-WALL KINEMATICS (Current Implementation)
# =============================================================================

class ThinWallKinematics(DeformationKinematics):
    """Kinematics for thin-walled cylindrical vessels.
    
    Assumptions:
    - h << a (thickness << radius)
    - Uniform stress/strain through thickness
    - Incompressibility (det(F) = 1)
    - Axisymmetric deformation
    
    This is the simplest implementation for current G&R simulations.
    
    Uses TensorOperations for all tensor algebra, ensuring consistency
    with potential future thick-wall or FEM implementations.
    """
    
    def __init__(self, coord_system: CoordinateSystem = CoordinateSystem.CYLINDRICAL):
        """Initialize thin-wall kinematics.
        
        Args:
            coord_system: Coordinate system for interpretation
        """
        super().__init__()
        self.coord_system = coord_system
    
    def compute_deformation_gradient(
        self,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> np.ndarray:
        """Compute F for thin-wall vessel at mid-wall.
        
        Args:
            current_geometry: {
                'inner_radius': a,      # Current inner radius (m)
                'thickness': h,         # Current thickness (m)
                'axial_stretch': lambda_z  # Axial stretch (-)
            }
            reference_geometry: {
                'inner_radius': A,      # Reference inner radius (m)
                'thickness': H,         # Reference thickness (m)
                'axial_stretch': Lambda_z  # Reference axial stretch (-)
            }
            
        Returns:
            F as 3x3 diagonal matrix [λ_r, λ_θ, λ_z] in cylindrical coordinates
        """
        # Extract current geometry
        a = current_geometry['inner_radius']
        h = current_geometry['thickness']
        lambda_z = current_geometry['axial_stretch']
        
        # Extract reference geometry
        A = reference_geometry['inner_radius']
        H = reference_geometry['thickness']
        
        # Mid-wall radii (thin-wall approximation: evaluate at midpoint)
        r_mid = a + 0.5 * h
        R_mid = A + 0.5 * H
        
        # Circumferential stretch (at mid-wall)
        lambda_theta = r_mid / R_mid
        
        # Radial stretch (from incompressibility: λ_r · λ_θ · λ_z = 1)
        lambda_r = 1.0 / (lambda_theta * lambda_z)
        
        # F in principal directions (cylindrical coordinates)
        # F = diag(λ_r, λ_θ, λ_z)
        F = np.diag([lambda_r, lambda_theta, lambda_z])
        
        return F

    def compute_lagrange_multiplier(self, sigma: np.ndarray) -> np.ndarray:
        """Compute Lagrange multiplier for thin-wall constraint: σ_r = 0.
        
        For thin-wall, the constraint is that radial stress vanishes.
        The Lagrange multiplier is a uniform hydrostatic stress equal
        to the radial component:
        
        λ = σ_r · I  (diagonal matrix with σ_r on diagonal)
        
        After subtraction: σ_constrained = σ - λ
        - σ_r becomes 0 (constraint satisfied)
        - σ_θ, σ_z shift by -σ_r (maintains deviatoric stress)
        
        Args:
            sigma: Material stress before constraint (3×3)
            
        Returns:
            Lagrange multiplier tensor (3×3 diagonal)
        """
        # Extract radial stress (index 0 in cylindrical coords)
        # TODO: refactor this tight coupling to sigma[0,0]
        sigma_r = sigma[0, 0]
        
        # Lagrange multiplier is uniform hydrostatic stress
        lagrange = np.diag([sigma_r, sigma_r, sigma_r])
        
        return lagrange
        
    def compute_circumferential_stretch(
        self,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> float:
        """Compute circumferential stretch λ_θ.
        
        λ_θ = r_mid / R_mid
        
        Args:
            current_geometry: Current configuration
            reference_geometry: Reference configuration
            
        Returns:
            Circumferential stretch λ_θ (-)
        """
        a = current_geometry['inner_radius']
        h = current_geometry['thickness']
        A = reference_geometry['inner_radius']
        H = reference_geometry['thickness']
        
        r_mid = a + 0.5 * h
        R_mid = A + 0.5 * H
        
        return r_mid / R_mid
    
    def compute_radial_stretch(
        self,
        J: float,
        lambda_theta: float,
        lambda_z: float
    ) -> float:
        """Compute radial stretch from volume ratio.
        
        λ_r = J / (λ_θ · λ_z)
        
        Args:
            J: Volume ratio (-)
            lambda_theta: Circumferential stretch
            lambda_z: Axial stretch
            
        Returns:
            Radial stretch λ_r (-)
        """
        #TODO: if J not provided, then calculate it here.
        return J / (lambda_theta * lambda_z)
    
    def compute_thickness_from_incompressibility(
        self,
        J: float,
        reference_thickness: float,
        F: Optional[np.ndarray] = None,
        lambda_theta: Optional[float] = None,
        lambda_z: Optional[float] = None
    ) -> float:
        """Compute current thickness from volume ratio.
        
        From λ_r = h/H and λ_r = J/(λ_θ · λ_z):
        h = (J / (λ_θ · λ_z)) * H
        """
        # Extract stretches from F if provided (thin-wall: F is diagonal)
        # TODO: use component recognition to get lambda_theta, lambda_z
        if F is not None:
            # For thin-wall, F is diagonal: [λ_r, λ_θ, λ_z]
            lambda_theta = F[1, 1]
            lambda_z = F[2, 2]
        elif lambda_theta is not None and lambda_z is not None:
            # Use provided stretches
            pass
        else:
            raise ValueError(
                "Must provide either F or both lambda_theta and lambda_z"
            )
        
        h = (J / (lambda_theta * lambda_z)) * reference_thickness
        return h
    
    def compute_all_kinematic_quantities(
        self,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> Dict[str, any]:
        """Compute all kinematic quantities in one call.
        
        This is useful for debugging and output.
        
        Args:
            current_geometry: Current configuration
            reference_geometry: Reference configuration
            
        Returns:
            Dictionary with:
                - F: Deformation gradient (3x3)
                - C: Right Cauchy-Green tensor (3x3)
                - E: Green-Lagrange strain (3x3)
                - lambda_r, lambda_theta, lambda_z: Principal stretches
                - I1, I2, I3: Strain invariants
                - J: Jacobian (should be 1.0)
        """
        # Compute F
        F = self.compute_deformation_gradient(current_geometry, reference_geometry)
        
        # Extract stretches using TensorOperations
        lambda_r, lambda_theta, lambda_z = self.compute_stretches(F)
        
        # Compute derived quantities using TensorOperations
        C = self.compute_right_cauchy_green(F)
        E = self.compute_green_strain(F)
        I1, I2, I3 = self.compute_invariants(F)
        J = self.compute_jacobian(F)
        
        return {
            'F': F,
            'C': C,
            'E': E,
            'lambda_r': lambda_r,
            'lambda_theta': lambda_theta,
            'lambda_z': lambda_z,
            'I1': I1,
            'I2': I2,
            'I3': I3,
            'J': J
        }
    
    def compute_stress_from_equilibrium(
        self,
        pressure: float,
        inner_radius: float,
        thickness: float
    ) -> float:
        """Compute circumferential stress from equilibrium.
        
        This is the PRIMARY equilibrium equation used during solving:
        σ_θ = P · a / h
        
        Workflow:
        1. Simulation uses this to get target σ_θ
        2. Iteration updates geometry to match this stress
        3. After convergence, compute_axial_force() can be called for completeness
        
        Args:
            pressure: Internal pressure (Pa)
            geometry: Must contain 'inner_radius' (a) and 'thickness' (h)
            
        Returns:
            {'sigma_theta': σ_θ, 'sigma_r': σ_r}
            Note: Does NOT return σ_z 
        """
        a = inner_radius
        h = thickness

        # Check thin-wall validity (h/a < 0.1)
        # TODO: move this at the end of the simulation step
        # if h/a > 0.1:
        #     self.logger.warning(
        #         f"Thin-wall assumption may not be valid: h/a = {h/a:.3f} > 0.1. "
        #         f"[WARNING] Results may be inaccurate."
        #     )
        
        # Circumferential equilibrium
        sigma_theta = pressure * a / h
        
        return sigma_theta
    
    def compute_axial_force(
        self,
        pressure: float,
        geometry: Dict[str, float],
        sigma_z: float
    ) -> float:
        """Compute required axial force (post-solve).
        
        After solving for geometry and computing material stress σ_z,
        determine what axial force f_z maintains equilibrium.
        
        Axial equilibrium:
        σ_z =  f_z / A_wall 
        
        Therefore:
        f_z = σ_z · A_wall
        
        where A_wall = π · h · (2a + h)
        
        Interpretation (in the case of external force not included here):
        - f_z > 0: Requires external tension (e.g., tethering)
        - f_z = 0: Closed-end vessel (σ_z ≈ P·a/(2h))
        - f_z < 0: Requires external compression

        Args:
            pressure: Internal pressure (Pa)
            geometry: Must contain 'inner_radius' (a) and 'thickness' (h)
            sigma_z: Axial stress from material constitutive law (Pa)
            
        Returns:
            Required axial force f_z (N)
        """
        a = geometry['inner_radius']
        h = geometry['thickness']
        
        # Wall cross-sectional area
        A_wall = math.pi * h * (2.0 * a + h)
        
        # Required axial force from equilibrium
        f_z = sigma_z * A_wall

        return f_z
    
    def compute_wss(
        self,
        flow_rate: float,
        inner_radius: float,
        blood_viscosity: float
    ) -> float:
        """Compute wall shear stress for cylindrical pipe (Poiseuille flow).
        
        For steady, laminar, fully-developed flow in a circular pipe:
        
        WSS = τ_w = (4 * μ * Q) / (π * a³)
        
        Derivation:
        - Poiseuille flow: u(r) = (ΔP/(4μL)) * (a² - r²)
        - Wall velocity gradient: du/dr|_wall = -(ΔP/(2μL)) * a
        - WSS: τ_w = μ * du/dr|_wall
        - Flow rate: Q = (π * ΔP * a⁴) / (8μL)
        - Eliminate ΔP: τ_w = (4μQ) / (πa³)
        
        Args:
            flow_rate: Volumetric flow rate Q (m³/s)
            inner_radius: Inner radius a (m)
            blood_viscosity: Dynamic viscosity μ (Pa·s)
            
        Returns:
            Wall shear stress τ_w (Pa)
        """
        # Check for valid inputs
        if flow_rate < 0:
            raise ValueError(f"Flow rate must be non-negative, got {flow_rate}")
        if inner_radius <= 0:
            raise ValueError(f"Inner radius must be positive, got {inner_radius}")
        if blood_viscosity <= 0:
            raise ValueError(f"Blood viscosity must be positive, got {blood_viscosity}")
        
        # Poiseuille WSS formula
        # TODO: unlcear why blood viscosity not used in Latorre 2018. Rethink units.
        # wss = (4.0 * blood_viscosity * flow_rate) / (math.pi * inner_radius**3)
        wss = flow_rate / inner_radius**3

        return wss

    def compute_volume_ratio(
        self,
        F_history: List[np.ndarray],
        rhoR_history: List[float],
        rhoR_h: float,
        timestep: int
    ) -> float:
        """Compute volume ratio for thin-wall kinematics.
        
        Strategy:
        1. Simplified calculation from incompressibility: J = ρ_h / ρ 
        2. Fallback to F-based: J = det(F) (if density not available)
        
        Args:
            F_history: Deformation gradient history
            rhoR_history: Mass density history
            rhoR_h: Homeostatic density
            timestep: Target timestep
            
        Returns:
            Volume ratio J
        """

        # Try density first (if available)
        if 0 <= timestep < len(rhoR_history):
            rhoR = rhoR_history[timestep]
            if not np.isnan(rhoR):
                J = rhoR / rhoR_h
                return  J

        # Fallback to F
        if 0 <= timestep < len(F_history):
            F = F_history[timestep]
            J = self.compute_jacobian(F)
            return J
            
        # Neither available
        raise ValueError(
            f"Cannot compute J at timestep {timestep}: "
            f"neither F nor density available"
        )

    def get_component_name(self, index: int) -> str:
        """Get component name for given tensor index.
        
        For cylindrical coordinates:
            0 → 'r' (radial)
            1 → 'theta' (circumferential)
            2 → 'z' (axial)
        
        Args:
            index: Tensor component index (0, 1, or 2)
            
        Returns:
            Component name string
        """
        component_names = ['r', 'theta', 'z']
        if not (0 <= index < 3):
            raise ValueError(f"Index {index} out of range [0, 2]")
        return component_names[index]
                
    def compute_inner_radius(
        self,
        mid_radius: float,
        thickness: float
    ) -> float:
        """Compute inner radius from mid-radius and thickness.
        
        For thin-wall cylinder:
        r_mid = a + h/2
        Therefore: a = r_mid - h/2
        
        Args:
            mid_radius: Mid-wall radius (m)
            thickness: Wall thickness (m)
            
        Returns:
            Inner radius a (m)
        """
        #TODO: abstract to functions like update_geometry() that will be common across kinematics. Currently inner radius too specific?
        return mid_radius - 0.5 * thickness

    def get_equilibrium_variable(self) -> str:
        """Thin-wall solves for mid_radius."""
        return 'mid_radius'

    def get_reference_geometry_value(self, layer) -> float:
        #TODO: refactor
        """Homeostatic mid_radius for bracketing."""
        return layer.get_mid_radius(0)

# =============================================================================
# THICK-WALL KINEMATICS (Future Implementation)
# =============================================================================

class ThickWallKinematics(DeformationKinematics):
    """Kinematics for thick-walled cylindrical vessels.
    
    Accounts for:
    - Variation of stress/strain through thickness
    - Non-uniform radial deformation r(R)
    - Exact incompressibility constraint
    
    This is needed for:
    - Accurate residual stress computation
    - Vessels with h/a > 0.1
    - Multi-layer vessels with significant thickness
    
    Uses TensorOperations for consistency with thin-wall implementation.
    """
    
    def __init__(self, n_integration_points: int = 5):
        """Initialize thick-wall kinematics.
        
        Args:
            n_integration_points: Number of points through thickness for integration
        """
        super().__init__()
        self.n_points = n_integration_points
        self._setup_integration_points()
    
    def _setup_integration_points(self) -> None:
        """Setup Gauss quadrature points for through-thickness integration."""
        # Gauss-Legendre quadrature points and weights on [-1, 1]
        # Map to [0, 1] for through-thickness coordinate
        
        if self.n_points == 2:
            # 2-point Gauss
            xi = np.array([-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)])
            w = np.array([1.0, 1.0])
        elif self.n_points == 3:
            # 3-point Gauss
            xi = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
            w = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])
        elif self.n_points == 5:
            # 5-point Gauss
            xi = np.array([
                -np.sqrt((5.0 + 2.0*np.sqrt(10.0/7.0))/9.0),
                -np.sqrt((5.0 - 2.0*np.sqrt(10.0/7.0))/9.0),
                0.0,
                np.sqrt((5.0 - 2.0*np.sqrt(10.0/7.0))/9.0),
                np.sqrt((5.0 + 2.0*np.sqrt(10.0/7.0))/9.0)
            ])
            w = np.array([
                (322.0 - 13.0*np.sqrt(70.0))/900.0,
                (322.0 + 13.0*np.sqrt(70.0))/900.0,
                128.0/225.0,
                (322.0 + 13.0*np.sqrt(70.0))/900.0,
                (322.0 - 13.0*np.sqrt(70.0))/900.0
            ])
        else:
            raise ValueError(f"Unsupported number of integration points: {self.n_points}")
        
        # Map from [-1, 1] to [0, 1]
        self.integration_points = 0.5 * (xi + 1.0)
        self.integration_weights = 0.5 * w
    
    def compute_deformation_gradient(
        self,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> np.ndarray:
        """Compute F accounting for through-thickness variation.
        
        For thick-wall, F varies through thickness. This method returns
        F at a representative point (e.g., mid-wall).
        
        For full thick-wall analysis, use compute_deformation_gradient_at_point().
        """
        raise NotImplementedError(
            "Thick-wall kinematics requires through-thickness integration. "
            "Use compute_deformation_gradient_at_point() or ThinWallKinematics "
            "for current simulations."
        )
    
    def compute_deformation_gradient_at_point(
        self,
        R: float,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> np.ndarray:
        """Compute F at specific radial position R in reference configuration.
        
        Solves for r(R) from incompressibility:
        r² - a² = (R² - A²) · λ_z
        
        Then:
        λ_r = dr/dR
        λ_θ = r/R
        λ_z = given
        
        Args:
            R: Radial position in reference configuration (m)
            current_geometry: Current configuration
            reference_geometry: Reference configuration
            
        Returns:
            F at position R
        """
        a = current_geometry['inner_radius']
        A = reference_geometry['inner_radius']
        lambda_z = current_geometry['axial_stretch']
        
        # Solve for r(R) from incompressibility
        r_squared = a**2 + (R**2 - A**2) * lambda_z
        r = np.sqrt(r_squared)
        
        # Stretches
        lambda_theta = r / R
        # dr/dR from chain rule on r² - a² = (R² - A²) · λ_z
        # 2r · dr/dR = 2R · λ_z
        lambda_r = (R * lambda_z) / r
        
        # F in cylindrical coordinates
        F = np.diag([lambda_r, lambda_theta, lambda_z])
        
        return F
    
    def integrate_through_thickness(
        self,
        integrand_func,
        current_geometry: Dict[str, float],
        reference_geometry: Dict[str, float]
    ) -> float:
        """Integrate a function through wall thickness using Gauss quadrature.
        
        ∫[A to B] f(R) dR ≈ Σ w_i · f(R_i)
        
        Args:
            integrand_func: Function f(R, F) that takes radial position and F
            current_geometry: Current configuration
            reference_geometry: Reference configuration
            
        Returns:
            Integrated value
        """
        A = reference_geometry['inner_radius']
        H = reference_geometry['thickness']
        B = A + H  # Outer radius
        
        integral = 0.0
        
        for i in range(self.n_points):
            # Map integration point to [A, B]
            xi = self.integration_points[i]
            R = A + xi * H
            
            # Compute F at this point
            F = self.compute_deformation_gradient_at_point(
                R, current_geometry, reference_geometry
            )
            
            # Evaluate integrand
            value = integrand_func(R, F)
            
            # Add to integral with weight
            integral += self.integration_weights[i] * value
        
        # Scale by thickness
        integral *= H
        
        return integral
    
    def compute_stress_from_equilibrium(
        self,
        pressure: float,
        geometry: Dict[str, float]
    ) -> Dict[str, float]:
        """Thick-wall equilibrium requires radial integration.
        
        Radial equilibrium: dσ_r/dr + (σ_r - σ_θ)/r = 0
        
        This cannot be solved in closed form - requires iteration
        or integration through thickness.
        """
        raise NotImplementedError(
            "Thick-wall equilibrium requires iterative solution. "
            "Use ThinWallKinematics for current simulations."
        )
    
    def compute_axial_force(
        self,
        pressure: float,
        geometry: Dict[str, float],
        sigma_z: float
    ) -> float:
        """Thick-wall axial force computation."""
        raise NotImplementedError(
            "Thick-wall axial force requires through-thickness integration."
        )
    

# =============================================================================
# SIMPLER FACTORY (FOR YAML INTEGRATION)
# =============================================================================

class KinematicsFactory:
    """Factory for creating kinematics instances."""
    
    _registry = {
        'thin_wall_cylinder': ThinWallKinematics,
        # Add more as we implement them:
        # 'thin_wall_sphere': ThinWallSphereKinematics,
        # 'thick_wall_cylinder': ThickWallKinematics,
    }
    
    @classmethod
    def create(cls, geometry_type: str):
        """Create kinematics from geometry type.
        
        Args:
            geometry_type: Type string from YAML (e.g., 'thin_wall_cylinder')
            
        Returns:
            DeformationKinematics instance
            
        Raises:
            ValueError: If geometry_type not recognized
        """
        if geometry_type not in cls._registry:
            available = ', '.join(cls._registry.keys())
            raise ValueError(
                f"Unknown geometry type: '{geometry_type}'. "
                f"Available: {available}"
            )
        
        kinematics_class = cls._registry[geometry_type]
        return kinematics_class()
    
    @classmethod
    def available_types(cls) -> list[str]:
        """Get list of available geometry types."""
        return list(cls._registry.keys())