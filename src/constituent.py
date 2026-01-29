import math
import warnings
import numpy as np
from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod
from kinetics import Kinetics  
from kinetics_interface import ConstituentKineticsContext
from constitutive_laws import ConstitutiveModel
from mechanics import Mechanics
from mechanics_interface import ConstituentMechanicsContext

# Forward reference for Layer (avoids circular import)
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from layer import Layer

# =========================================================================
# ABSTRACT BASE CLASS
# =========================================================================

class Constituent(ABC):
    """Abstract base class for all constituents."""
    
    def __init__(self, name):
        self.name = name
        self.layer: Optional['Layer'] = None  # Reference to parent layer (set by Layer.add_constituent)
        self.kinetics: Optional[Kinetics] = None
        self.constitutive_model: Optional[ConstitutiveModel] = None # TODO: consider making this mandatory for SingleConstituent. Never initialized for MultiFiberFamilyConstituent.
        
        self.homeostatic_referential_density: Optional[float] = None
        self.deposition_stretch_scalar: Optional[float] = None
        self.deposition_stretch: Optional[np.ndarray] = None

        # Active stress properties (optional - only for contractile constituents)
        self.is_active = False
        self.active_properties = None
        
        # Active stress histories (only allocated if is_active=True)
        self.active_radius_history = []  # [a_act(0), a_act(1), ...]
        self.active_survival_history = []  # [[q_act(0,0)], [q_act(1,0), q_act(1,1)], ...]
        self.sigma_hat_act_history = []  # [σ^_act(0), σ^_act(1), ...] (tensors)
        self.active_stress_history = []  # [σ_act(0), σ_act(1), ...] (tensors)

        # Referential mass density evolution (rhoR_alpha) - timestep-indexed history
        self.rhoR_alpha_history = []  # List of mass densities over time
        self.survival_history = [] # List of survival fractions over time (for each cohort)
        self.k_alpha_history = []  # List of degradation rates of survival function over time
        self.mR_alpha_history = []  # List of production rates over time
        self.stimulus_function_history = []  # List of stimulus function values over time #TODO: only stored for plot currently, streamline with Plotting class

        self.F_alpha_history: List[List[np.ndarray]] = []
        self.sigma_hat_history: List[List[np.ndarray]] = []  # List of constituent-based partial stress over time TODO: perhaps make temporary array?
        self.stress_history: List[np.ndarray] = []

        self.tau_min = 0

    # =========================================================================
    # INITIALIZATION
    # =========================================================================

    @staticmethod 
    def from_parameters(name: str, properties: Dict[str, Any], layer: 'Layer' = None) -> 'Constituent':
        """Factory method to create appropriate constituent type."""
        constituent_type = properties.get('constituent_type', 'single')
        
        if constituent_type == 'multi_fiber_family':
            return MultiFiberFamilyConstituent.from_parameters(name, properties, layer)
        else:
            return SingleConstituent.from_parameters(name, properties, layer)
    
    # =========================================================================
    # PROPERTY ACCESSORS
    # =========================================================================

    def get_rhoR_alpha(self, timestep):
        """Get referential mass density for specified timestep."""
        if timestep < 0 or timestep >= len(self.rhoR_alpha_history):
            raise IndexError(f"Timestep {timestep} out of range [0, {len(self.rhoR_alpha_history)-1}]")
        
        return self.rhoR_alpha_history[timestep]
    
    def set_rhoR_alpha(self, timestep, new_density):
        """Update mass density for specified timestep."""
        if timestep < 0 or timestep >= len(self.rhoR_alpha_history):
            raise IndexError(f"Timestep {timestep} out of range [0, {len(self.rhoR_alpha_history)-1}]")
        
        self.rhoR_alpha_history[timestep] = new_density

    # Helper method to check if timestep exists
    # TODO: generalize for any history
    def _timestep_exists(self, timestep):
        """Check if a specific timestep exists in history."""
        return 0 <= timestep < len(self.rhoR_alpha_history)

    # =========================================================================
    # UTILITIES: HISTORY MANAGEMENT
    # =========================================================================

    def _update_history(self, history: List, timestep: int, value):
        """Update history at timestep (append if new, update if exists).
        
        This method handles the fixed-point iteration case where the same
        timestep may be computed multiple times during refinement.
        
        TODO: Future: When pre-allocated arrays are used, this will simply be:
            history[timestep] = value
        
        Args:
            history: History list to update (e.g., rhoR_alpha_history)
            timestep: Timestep index
            value: Value to store
        """
        if timestep < len(history):
            # Timestep exists - update it (fixed-point refinement)
            history[timestep] = value
        else:
            # New timestep - append
            # Special case: if history is empty, we can append any timestep
            # (handles cases where t=0 wasn't stored, e.g., elastin F_alpha_history)
            if len(history) == 0:
                history.append(value)
            # Otherwise, verify we're appending to the next sequential timestep
            elif timestep != len(history):
                raise ValueError(
                    f"Cannot append to timestep {timestep}: "
                    f"expected {len(history)} (history length)"
                )
            else:
                history.append(value)

    # =========================================================================
    # ABSTRACT METHODS (must be implemented by subclasses)
    # =========================================================================

    #TODO: Consider moving this logic to the numerics/simulation layer, and use set_rhoR_alpha function for the prediction.   
    @abstractmethod
    def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Compute referential mass density of constituent for target timestep."""
        pass

    @abstractmethod
    def compute_rhoR_alpha(self, target_timestep, dt, integration_method, survival_strategy):
        """Compute referential mass density of constituent for target timestep."""
        pass

    @abstractmethod
    def compute_sigma_alpha(self, target_timestep: int, dt: float, 
                       integration_method: str,
                       survival_function_computation: str) -> np.ndarray:
        """Compute constituent stress via heredity integral.
    
        σ_α(s) = ∫₀ˢ [mR(τ)·q(s,τ)/ρ_h] · σ̂_α(s,τ) dτ
    
        where:
            σ̂_α(s,τ) - stress from cohort deposited at τ, evaluated at s
            mR(τ) - production rate at deposition time τ
            q(s,τ) - survival function (fraction surviving from τ to s)
            ρ_h - homeostatic mass density
    
        Returns:
            3×3 stress tensor in Pa
        """
        pass

# =========================================================================
# SINGLE CONSTITUENT IMPLEMENTATION
# =========================================================================

class SingleConstituent(Constituent):
    """Single constituent (elastin, muscle, etc.)."""
    
    def __init__(self, name):
        super().__init__(name)
        self.params = {}
        self.mechanics = Mechanics()

    # =========================================================================
    # INITIALIZATION
    # =========================================================================
         
    @classmethod
    def from_parameters(cls, name: str, properties: Dict[str, Any], layer: 'Layer' = None) -> 'SingleConstituent':
        """Create and fully initialize single constituent from parameters."""
        print(f"    Initializing single constituent: {name}")
        
        constituent = cls(name)

        # Set layer reference immediately if provided
        if layer is not None:
            constituent.layer = layer
        
        # Initialize homeostatic mass density (t=0)
        constituent.homeostatic_referential_density = properties['mass_fraction'] * 1050.0
        constituent.rhoR_alpha_history.append(constituent.homeostatic_referential_density)
        print(f"        Homeostatic mass density (rhoR_alpha at t=0): {constituent.homeostatic_referential_density:.2f} kg/m³")

        # Initialize constitutive model
        if 'constitutive_model' in properties:
            constituent.constitutive_model = ConstitutiveModel.from_parameters(properties['constitutive_model'])
            print(f"        Constitutive model: {constituent.constitutive_model}")
        
        # Initialize active properties if present
        if 'active_properties' in properties:
            constituent._initialize_active_properties(properties['active_properties'])

        # Initialize deposition stretch and deformation gradient
        deposition_stretch_input = properties['deposition_stretch']
        constituent._initialize_deposition_stretch(deposition_stretch_input)

        # Initialize kinetics
        if 'kinetics' in properties:
            constituent.kinetics = Kinetics.from_parameters(properties['kinetics'])
            constituent._initialize_homeostatic_kinetics()
        else:
            print(f"        No kinetics (elastin-like constituent)")
        
        # TODO: Set tau_min based on constituent type

        # Store all parameters
        constituent.params = properties.copy()
        
        return constituent

    def _initialize_deposition_stretch(self, deposition_stretch_input) -> None:
        """Initialize deposition stretch tensor G_α from YAML input.
        
        Supports three input formats:
        1. Scalar (isotropic): deposition_stretch: 1.40 #TODO: correct
        → G_α = 1.40 · I
        
        2. List (anisotropic diagonal): deposition_stretch: [1.2, 1.5, 1.1]
        → G_α = diag(1.2, 1.5, 1.1)
        
        3. Dict (full tensor - future): deposition_stretch: {type: "tensor", values: [[...]]}
        → G_α = full 3×3 tensor
        
        Args:
            deposition_stretch_input: Value from YAML (scalar, list, or dict)
        
        Raises:
            ValueError: If input format is invalid
        """        
        # Only accept scalar input
        if not isinstance(deposition_stretch_input, (int, float)):
            raise ValueError(
                f"deposition_stretch must be a scalar, got {type(deposition_stretch_input)}"
            )
        
        self.deposition_stretch_scalar = float(deposition_stretch_input)

        has_fiber = (
            self.constitutive_model is not None and 
            hasattr(self.constitutive_model, 'fiber_angle_degrees') and
            self.constitutive_model.fiber_angle_degrees is not None
        )
    
        if has_fiber:
            # MODE 2: FIBER-ONLY DEPOSITION
            fiber_angle = self.constitutive_model.fiber_angle_degrees
            
            # Determine fiber direction component based on angle
            # TODO: Create mapping class in deformation kinematics to be used in the code
            if abs(fiber_angle - 0.0) < 1e-6:
                # 0° = Axial (z-direction, index 2)
                fiber_component = 2
                coord_name = 'z (axial)'
                
            elif abs(fiber_angle - 90.0) < 1e-6:
                # 90° = Circumferential (θ-direction, index 1)
                fiber_component = 1
                coord_name = 'θ (circumferential)'
                
            elif abs(fiber_angle - 45.0) < 1e-6:
                # 45° = Helical (not yet implemented)
                raise NotImplementedError(
                    f"Helical fiber deposition (45°) not yet implemented. "
                    f"Only axial (0°) and circumferential (90°) are supported."
                )
                
            else:
                # General angle - not yet implemented
                raise NotImplementedError(
                    f"Fiber deposition at angle {fiber_angle:.1f}° not yet implemented. "
                    f"Only axial (0°) and circumferential (90°) are supported. "
                    f"For helical (45°), implementation is needed."
                )
            
            # Create tensor with ones except fiber direction
            self.deposition_stretch = np.eye(3)
            self.deposition_stretch[fiber_component, fiber_component] = self.deposition_stretch_scalar
            
            print(f"        Deposition stretch: λ_dep = {self.deposition_stretch_scalar:.2f} (fiber-only)")
            print(f"        Fiber angle: {fiber_angle:.1f}° → {coord_name}")
            print(f"        G_α = diag({self.deposition_stretch[0,0]:.2f}, "
                f"{self.deposition_stretch[1,1]:.2f}, {self.deposition_stretch[2,2]:.2f})")
            
        else:
            # MODE 1: ISOTROPIC 3D DEPOSITION #TODO: rename
            # Incompressibility: λ_r · λ_θ · λ_z = 1
            # If λ_θ = λ_z = λ, then λ_r = 1/λ²
            lambda_r = 1.0 / (self.deposition_stretch_scalar ** 2)
            self.deposition_stretch = np.diag([lambda_r, self.deposition_stretch_scalar, self.deposition_stretch_scalar])
            
            print(f"        Deposition stretch: λ_dep = {self.deposition_stretch_scalar:.2f} (isotropic 3D)")
            print(f"        G_α = diag({lambda_r:.4f}, {self.deposition_stretch_scalar:.2f}, {self.deposition_stretch_scalar:.2f})")

    def _initialize_homeostatic_kinetics(self) -> None:
        """Initialize kinetics histories at t=0 (homeostatic state).
        
        Assumes rhoR_alpha_history[0] already initialized.
        
        At homeostasis:
        - Production = Degradation (steady state)
        - mR_h = rhoR_h * k_alpha_h
        - q(0,0) = 1.0 (just deposited material survives 100%)
        """
        print(f"        Homeostatic kinetics:")
               
        # Get homeostatic degradation rate from degradation function
        k_alpha_h = self.kinetics.degradation_function.k_alpha_h
        
        # At homeostasis: production = degradation
        # mR_h = rhoR_h * k_alpha_h (steady state condition)
        mR_alpha_h = self.homeostatic_referential_density * k_alpha_h
        
        # Initialize kinetics histories with t=0 values
        self.k_alpha_history.append(k_alpha_h)
        self.mR_alpha_history.append(mR_alpha_h)
        self.stimulus_function_history.append(1.0) #TODO: compute from kinetics instead of hardcode
        
        # Initialize survival history
        # At t=0, we have one cohort deposited at τ=0 with q(0,0) = 1.0
        self.survival_history.append([1.0])
        
        # Print summary
        print(f"          k_alpha(0) = {k_alpha_h:.6f} 1/day")
        print(f"          mR_alpha(0) = {mR_alpha_h:.6f} kg/(m³·day)")
        print(f"          q(0,0) = 1.0")

    def _initialize_active_properties(self, active_props: Dict[str, float]) -> None:
        """Initialize active properties for contractile constituents (e.g., SMCs).
        
        Args:
            active_props: Dictionary with:
                - T_act_h: Maximum active stress (kPa)
                - k_act: Activation rate constant (1/day)
                - lambda_0: Minimum stretch for activation (-)
                - lambda_m: Optimal stretch for activation (-)
                - CB: Basal calcium concentration (-)
                - CS: Shear stress sensitivity (optional, default=0)
        """
        print("        Active properties:")
        
        self.is_active = True
        
        # Unit conversions
        # TODO: centralize unit conversions in a utility module
        kPa_to_Pa = 1.0e3
        
        # Store active properties
        self.active_properties = {
            'T_act': active_props['T_act_h'] * kPa_to_Pa,  # Convert to Pa
            'k_act': active_props['k_act'],                # 1/day
            'lambda_0': active_props['lambda_0'],          # Dimensionless
            'lambda_m': active_props['lambda_m'],          # Dimensionless
            'CB': active_props['CB'],                      # Dimensionless
            'CS': active_props['CS']                       # Dimensionless
        }
        
        # Initialize active histories at t=0
        # At homeostasis: a_act(0) = a_h (homeostatic inner radius)
        a_h = self.layer.get_homeostatic_inner_radius()
        self.active_radius_history.append(a_h)
        
        # Active survival: q_act(0,0) = 1.0
        self.active_survival_history.append([1.0])
        
        # Active stress at homeostasis (compute using homeostatic conditions)
        sigma_hat_act_h = self._compute_active_stress_at_homeostasis()
        self.sigma_hat_act_history.append(sigma_hat_act_h)
        
        # Print summary
        print(f"          T_act = {active_props['T_act_h']:.1f} kPa")
        print(f"          k_act = {self.active_properties['k_act']:.6f} 1/day")
        print(f"          λ_0 = {self.active_properties['lambda_0']:.2f}, λ_m = {self.active_properties['lambda_m']:.2f}")
        print(f"          CB = {self.active_properties['CB']:.4f}")
        print(f"          a_act(0) = {a_h*1000:.4f} mm")
        print(f"          σ_hat_act(0) = {sigma_hat_act_h/1000:.2f} kPa")

    # =========================================================================
    # COMPUTATION: MASS DENSITY
    # =========================================================================

    def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess referential mass density for target timestep."""
        #TODO: Consider moving this logic to the numerics/simulation layer, and use set_rhoR_alpha function for the prediction.   

        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        
        if target_timestep == 0:
            raise ValueError("Cannot guess timestep 0: should be initialized with homeostatic value")
        
        if guess_method == "from_previous_timestep":
            #TODO: When pre-allocating with 0 values, this check needs to be adjusted.
            if self._timestep_exists(target_timestep):
                warnings.warn(f"Target timestep {target_timestep} already exists. Returning existing value.")
                return self.get_rhoR_alpha(target_timestep)
            
            previous_timestep = target_timestep - 1
            
            if not self._timestep_exists(previous_timestep):
                raise ValueError(f"Cannot guess timestep {target_timestep}: previous timestep {previous_timestep} does not exist")
            
            #TODO: Change this part for pre-allocated arrays, no longer append but set value.
            previous_density = self.rhoR_alpha_history[previous_timestep]
            self.rhoR_alpha_history.append(previous_density)
        else:
            raise ValueError(f"Unknown guess method: {guess_method}")
        
        return self.get_rhoR_alpha(target_timestep)

    def compute_rhoR_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute referential mass density of constituent at target timestep.

        Args:
            target_timestep: Target time index
            dt: Time step size (from simulation)
            integration_method: Integration method ('simpson' or 'trapezoidal')
            survival_function_computation: Survival strategy ('naive', 'backward')
        """
        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        
        if target_timestep == 0:
            raise ValueError("Cannot compute timestep 0: should be initialized with homeostatic value")

        # If no degradation/production - maintain previous referential mass density
        if self.kinetics is None:
            previous_rhoR_alpha = self.rhoR_alpha_history[target_timestep - 1]
            self._update_history(self.rhoR_alpha_history, target_timestep, previous_rhoR_alpha)
            return self.get_rhoR_alpha(target_timestep)
        
        # Pass necessary data from constituent to kinetics class
        context = ConstituentKineticsContext(self)

        # STEP 1: Compute k_alpha at target timestep (Constituent asks Kinetics for computation)
        k_alpha = self.kinetics.compute_k_alpha(context, target_timestep)
        self._update_history(self.k_alpha_history, target_timestep, k_alpha)

        # STEP 2: Compute production rate
        # TODO: improve fail-safe if we don't have a guess for rhoR_alpha.
        # Consider how this interacts with guess_rhoR_alpha (which is needed for this step).
        rhoR_alpha = self.get_rhoR_alpha(target_timestep) # comes from guess
        mR_alpha, stimulus_function = self.kinetics.compute_production_rate(context, target_timestep, k_alpha, rhoR_alpha)
        self._update_history(self.mR_alpha_history, target_timestep, mR_alpha)
        self._update_history(self.stimulus_function_history, target_timestep, stimulus_function)

        # STEP 3: Compute survival function q(s, tau) for all cohorts
        survival_values = self.kinetics.compute_survival_function(
            k_alpha_history=self.k_alpha_history,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method,
            survival_function_computation=survival_function_computation
        )
        self._update_history(self.survival_history, target_timestep, survival_values)

        # Step 4: Compute heredity integral (just integrate mR × q!)
        rhoR_alpha = self.kinetics.compute_heredity_integral(
            self.mR_alpha_history,
            survival_values,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method,
            rhoR_alpha_homeostatic=self.homeostatic_referential_density
        )
        
        self.rhoR_alpha_history[target_timestep] = rhoR_alpha

        return self.get_rhoR_alpha(target_timestep)

    # =========================================================================
    # COMPUTATION: STRESS
    # =========================================================================

    def get_stress(self):
        """Get intramural stress contribution."""
        return self.stress_contribution

    def compute_sigma_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute referential mass density of constituent at target timestep.

        Args:
            target_timestep: Target time index
            dt: Time step size (from simulation)
            integration_method: Integration method ('simpson' or 'trapezoidal')
        """
        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        
        if target_timestep == 0:
            raise ValueError("Cannot compute timestep 0: should be initialized with homeostatic value")
            
        # Pass necessary data from constituent to kinetics class
        context = ConstituentMechanicsContext(self)

        if self.kinetics is None:
            #print(f"        {self.name}: No kinetics - single cohort (τ=0)")

            # STEP 1: Compute F_α for single cohort (τ=0)
            F_alpha = self.mechanics.compute_F_alpha_for_cohort(
                context, target_timestep, deposition_timestep=0
            )
            self._update_history(self.F_alpha_history, target_timestep, [F_alpha])

            # STEP 2: Compute S_hat_α (PK2 stress) from constitutive model
            S_hat_alpha = self.mechanics.compute_S_hat_alpha_for_cohort(
                context, F_alpha
            )

            # STEP 3: Compute σ̂_α (Cauchy stress)
            J = self.layer.get_volume_ratio(target_timestep)
            sigma_hat_alpha = self.mechanics.compute_sigma_hat_alpha_for_cohort(
                J, F_alpha, S_hat_alpha
            )
            self._update_history(self.sigma_hat_history, target_timestep, [sigma_hat_alpha])
            
            # STEP 4: No heredity integral - σ_α = (ρR_α/ρR) σ̂_α
            # TODO: abstract with mass fraction?
            rhoR_alpha = self.get_rhoR_alpha(target_timestep)
            rhoR_h = self.layer.get_homeostatic_density()
            sigma_alpha = (rhoR_alpha / rhoR_h) * sigma_hat_alpha
            self._update_history(self.stress_history, target_timestep, sigma_alpha)
            
            return sigma_alpha

        #print(f"        {self.name}: With kinetics - multiple cohorts")
        # STEP 1: Compute F_α(s,τ) for all cohorts τ ∈ [τ_min, s]
        F_alpha_cohorts = self.mechanics.compute_F_alpha_for_all_cohorts(
            context, target_timestep
        )
        self._update_history(self.F_alpha_history, target_timestep, F_alpha_cohorts)
        
        # STEP 2: Compute S_hat_α(s,τ) (PK2 stress) for all cohorts
        S_hat_alpha_cohorts = self.mechanics.compute_S_hat_alpha_for_all_cohorts(
            context, F_alpha_cohorts
        )
        
        # STEP 3: Compute σ̂_α(s,τ) (Cauchy stress) for all cohorts
        J = self.layer.get_volume_ratio(target_timestep)
        sigma_hat_alpha_cohorts = self.mechanics.compute_sigma_hat_alpha_for_all_cohorts(
            J, F_alpha_cohorts, S_hat_alpha_cohorts
        )
        self._update_history(self.sigma_hat_history, target_timestep, sigma_hat_alpha_cohorts)
        
        # STEP 4: Integrate heredity integral: σ_α = ∫ (mR·q/ρ_h) · σ̂_α dτ
        # TODO: consider passing whole survival history to integrate_constituent_stress
        # to improve readability? 
        # TODO: break into helper function
        survival_values = self.survival_history[target_timestep]
        rhoR_homeostatic = self.layer.get_homeostatic_density()
        rhoR_alpha_homeostatic = self.get_rhoR_alpha(0)
        sigma_alpha = self.mechanics.integrate_constituent_stress(
            sigma_hat_alpha_cohorts,
            self.mR_alpha_history,
            survival_values,
            rhoR_homeostatic,
            self.tau_min,
            target_timestep,
            dt,
            integration_method,
            rhoR_alpha_homeostatic
        )

        # STEP 2: Compute active stress (if constituent is contractile)
        if self.is_active:
            #print(f"        {self.name}: Computing active stress")
            
            # Compute active radius a_act(s)
            a_act = self.compute_active_radius(
                target_timestep,
                dt,
                integration_method,
                survival_function_computation 
            )
            
            # Compute σ̂_act (circumferential component only)
            sigma_hat_act = self.compute_active_stress_component(target_timestep)
            self._update_history(self.sigma_hat_act_history, target_timestep, sigma_hat_act)
            
            # Compute σ_act = (ρR_α / (J·ρR_h)) · σ̂_act 
            sigma_alpha_active = np.zeros((3, 3))
            J = self.layer.get_volume_ratio(target_timestep)
            rhoR_alpha = self.get_rhoR_alpha(target_timestep)
            rhoR_h = self.layer.get_homeostatic_density()
            active_stress_scalar = (rhoR_alpha / (J * rhoR_h)) * sigma_hat_act
            sigma_alpha_active[1, 1] = active_stress_scalar  # Circumferential
            
            # TODO: Add to debug
            # print(f"          a_act = {a_act*1000:.4f} mm")
            # print(f"          σ̂_act = {sigma_hat_act/1000:.2f} kPa")
            # print(f"          σ_act = {active_stress_scalar/1000:.2f} kPa")
            
            # STEP 3: Sum passive + active
            sigma_alpha = sigma_alpha + sigma_alpha_active
            
        # Store total stress
        self._update_history(self.stress_history, target_timestep, sigma_alpha)
            
        return sigma_alpha

    # =========================================================================
    # COMPUTATION: ACTIVE STRESS
    # =========================================================================

    def _compute_active_stress_at_homeostasis(self) -> float:
        """Compute homeostatic active stress (at t=0).
        
        At homeostasis:
        - a = a_h, a_act = a_h → λ_act = 1.0
        - τ_w = τ_w_h → C = CB
        
        Returns:
            σ_act_h (Pa) - circumferential active stress
        """
        # At homeostasis
        lambda_act = 1.0  # a_h / a_act_h = 1.0
        C = self.active_properties['CB']  # No shear deviation at homeostasis
        
        # Parabolic activation function
        lm = self.active_properties['lambda_m']
        l0 = self.active_properties['lambda_0']
        parab_act = 1.0 - ((lm - lambda_act) / (lm - l0))**2
        
        # Active stress
        T_act = self.active_properties['T_act']
        sigma_hat_act_h = T_act * (1.0 - np.exp(-C**2)) * lambda_act * parab_act
        
        return sigma_hat_act_h

    def compute_active_radius(
        self,
        target_timestep: int,
        dt: float,
        integration_method: str,
        survival_function_computation: str,
    ) -> float:
        """Compute active radius via heredity integral.
        
        a_act(s) = ∫₀ˢ k_act · a(τ) · q_act(s,τ) dτ
        
        where:
            a(τ) - layer inner radius at deposition time τ
            q_act(s,τ) = exp(-k_act·(s-τ)) - active material survival
            k_act - activation rate constant
        
        Args:
            target_timestep: Current time s
            dt: Time step size (days)
            integration_method: Integration method ('simpson' or 'trapezoidal')
            survival_function_computation: Survival strategy ('naive', 'backward')
            
        Returns:
            Active radius a_act (m)
        """
        # TODO: Refactor
        if not self.is_active:
            raise RuntimeError(f"Constituent '{self.name}' is not active")
        
        if target_timestep == 0:
            # Already initialized in _initialize_active_properties
            return self.active_radius_history[0]
        
        # STEP 1: Compute active survival function q_act(s,τ) = exp(-k_act·(s-τ))
        # This is an EXPONENTIAL survival with CONSTANT degradation rate
        k_act = self.active_properties['k_act']
        
        # Build constant k_act history [k_act, k_act, ..., k_act]
        k_act_history = [k_act] * (target_timestep + 1)
        
        # Compute survival using exponential survival (reuse kinetics!)
        from kinetics import ExponentialSurvival, ConstantDegradationRate
        
        deg_function = ConstantDegradationRate(k_act)
        survival_function = ExponentialSurvival(deg_function)
        
        q_act_values = survival_function.compute_survival_function(
            k_alpha_history=k_act_history,
            tau_min=0,  # Active material starts at t=0
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method,
            survival_function_computation=survival_function_computation
        )
        
        self._update_history(self.active_survival_history, target_timestep, q_act_values)
        
        # STEP 2: Build integrand: k_act · a(τ) · q_act(s,τ)
        integrand = []
        for tau in range(0, target_timestep + 1):
            a_tau = self.layer.get_inner_radius(tau)
            integrand.append(k_act * a_tau * q_act_values[tau])
        
        # STEP 3: Integrate
        from integrators import IntegratorFactory
        
        n_points = len(integrand)
        integrator = IntegratorFactory.create(
            integration_method,
            dt,
            n_points - 1,  # Last index
            0              # First index
        )
        
        a_act_deposited = integrator.integrate(integrand)

        # STEP 4: Add initial active material contribution
        # TODO: clean up step 4
        a_act_homeostatic = self.active_radius_history[0]  # a_act(0) from homeostasis
        q_act_initial = q_act_values[0]  # q_act(s, 0) - survival of initial active material
        a_act_initial = a_act_homeostatic * q_act_initial

        a_act = a_act_initial + a_act_deposited
        
        self._update_history(self.active_radius_history, target_timestep, a_act)
        
        return a_act

    def compute_active_stress_component(
        self,
        target_timestep: int
    ) -> float:
        """Compute circumferential active stress σ̂_act.
        
        Formula (from Latorre 2018):
            σ̂_act = T_act · [1 - exp(-C²)] · λ_act · parab_act
        
        where:
            C = CB - CS·(τ_w/τ_w_h - 1) - vasomotor control
            λ_act = a/a_act - active stretch
            parab_act = 1 - [(λ_m - λ_act)/(λ_m - λ_0)]² - activation curve
        
        Args:
            target_timestep: Current time s
            
        Returns:
            σ̂_act (Pa) - circumferential active stress
        """
        if not self.is_active:
            return 0.0
        
        # Get current geometry
        a = self.layer.get_inner_radius(target_timestep)
        a_act = self.active_radius_history[target_timestep]
        
        # Get wall shear stress (for vasomotor control)
        tau_w = self.layer.wss_history[target_timestep]
        tau_w_h = self.layer.homeostatic_wss
        
        # Active properties
        T_act = self.active_properties['T_act']
        CB = self.active_properties['CB']
        CS = self.active_properties['CS']
        lambda_m = self.active_properties['lambda_m']
        lambda_0 = self.active_properties['lambda_0']
        
        # STEP 1: Vasomotor control
        C = CB - CS * (tau_w / tau_w_h - 1.0)
        
        # STEP 2: Active stretch
        lambda_act = a / a_act if a_act > 0 else 0.0
        
        # STEP 3: Parabolic activation function
        # parab_act = 1 - [(λ_m - λ_act)/(λ_m - λ_0)]²
        denominator = (lambda_m - lambda_0) if abs(lambda_m - lambda_0) > 1e-10 else 1.0
        parab_act = 1.0 - ((lambda_m - lambda_act) / denominator)**2
        
        # Clip parab_act to [0, 1] (physical constraint)
        parab_act = max(0.0, min(1.0, parab_act))
        
        # STEP 4: Active stress (hat quantity - per cohort)
        sigma_hat_act = T_act * (1.0 - np.exp(-C**2)) * lambda_act * parab_act
        
        return sigma_hat_act

    # =========================================================================
    # OUTPUT: DATA EXPORT
    # =========================================================================

    def to_dict(self, timestep: int, time: float = None) -> Dict[str, Any]:
        """Convert constituent state to dictionary for serialization.
        
        Args:
            timestep: Timestep index
            time: Simulation time (optional)
            
        Returns:
            Dictionary containing constituent state
        """
        data = {
            'time': time if time is not None else float(timestep),
            'rho_alpha': self.get_rhoR_alpha(timestep)
        }
        
        # Add kinetics data if present
        if hasattr(self, 'k_alpha_history') and len(self.k_alpha_history) > timestep:
            data['k_alpha'] = self.k_alpha_history[timestep]
        else:
            data['k_alpha'] = None
            
        if hasattr(self, 'mR_alpha_history') and len(self.mR_alpha_history) > timestep:
            data['mR_alpha'] = self.mR_alpha_history[timestep]
        else:
            data['mR_alpha'] = None
        
        if hasattr(self, 'stimulus_function_history') and len(self.stimulus_function_history) > timestep:
            data['stimulus_function'] = self.stimulus_function_history[timestep]
        else:
            data['stimulus_function'] = None

        # Add stress if available
        # TODO: streamline with a loop
        if len(self.stress_history) > timestep:
            data['stress_r'] = self.stress_history[timestep][0,0]
            data['stress_theta'] = self.stress_history[timestep][1,1]
            data['stress_z'] = self.stress_history[timestep][2,2]
        else:
            data['stress_r'] = None
            data['stress_theta'] = None
            data['stress_z'] = None

        if len(self.sigma_hat_history) > timestep:
            data['sigma_hat_r'] = self.sigma_hat_history[timestep][0][0,0]
            data['sigma_hat_theta'] = self.sigma_hat_history[timestep][0][1,1]
            data['sigma_hat_z'] = self.sigma_hat_history[timestep][0][2,2]
        else:
            data['sigma_hat_r'] = None
            data['sigma_hat_theta'] = None
            data['sigma_hat_z'] = None
               
        if len(self.sigma_hat_act_history) > timestep and self.is_active:
            data['a_act'] = self.active_radius_history[timestep]
            data['sigma_hat_act'] = self.sigma_hat_act_history[timestep]
        else:
            data['a_act'] = None
            data['sigma_hat_act'] = None

        return data

# =========================================================================
# MULTI-FIBER FAMILY CONSTITUENT IMPLEMENTATION
# =========================================================================

class MultiFiberFamilyConstituent(Constituent):
    #TODO: To clean up after code works with single constituent.
    """Multi-fiber family constituent (e.g., collagen with multiple orientations)."""
    
    def __init__(self, name):
        super().__init__(name)
        self.fiber_families = []
        self.shared_properties = {}
        self.total_mass_fraction = 0.0

    # =========================================================================
    # INITIALIZATION
    # =========================================================================
        
    @classmethod
    def from_parameters(cls, name, properties):
        """Create and fully initialize multi-fiber family constituent from parameters."""
        print(f"    Initializing multi-fiber family constituent: {name}")
        
        constituent = cls(name)
        
        # Store shared properties
        constituent.shared_properties = properties['shared_properties'].copy()
        constituent.total_mass_fraction = constituent.shared_properties['total_mass_fraction']
        
        # Validate mass fraction ratios sum to 1.0
        constituent._validate_mass_fraction_ratios(properties['fiber_families'])
        
        # Create individual fiber families
        # (Each fiber family initializes its own rhoR and kinetics via SingleConstituent)
        for family_name, family_props in properties['fiber_families'].items():
            fiber_family = constituent._create_fiber_family(family_name, family_props)
            constituent.fiber_families.append(fiber_family)
        
        print(f"        Created {len(constituent.fiber_families)} fiber families")
        
        # Initialize total mass density history (sum of all families at t=0)
        # Assumes each fiber family has rhoR_alpha_history[0] already initialized
        total_homeostatic_density = sum(
            family.rhoR_alpha_history[0] for family in constituent.fiber_families
        )
        constituent.rhoR_alpha_history.append(total_homeostatic_density)
        print(f"        Total rhoR(0) = {total_homeostatic_density:.2f} kg/m³")
        
        # Initialize total kinetics histories (if fiber families have kinetics)
        constituent._initialize_total_kinetics_histories()
        
        return constituent
    
    def _validate_mass_fraction_ratios(self, fiber_families):
        """Validate that mass fraction ratios sum to 1.0."""
        total_ratio = sum(family['mass_fraction_ratio'] for family in fiber_families.values())
        
        if abs(total_ratio - 1.0) > 1e-6:
            raise ValueError(f"Fiber family mass fraction ratios must sum to 1.0, got {total_ratio}")
        
        print(f"        Mass fraction validation: {total_ratio:.6f} ≈ 1.0 ✓")
    
    def _create_fiber_family(self, family_name, family_props):
        """Create individual fiber family."""
        # Calculate absolute mass fraction for this family
        family_mass_fraction = (self.total_mass_fraction * 
                              family_props['mass_fraction_ratio'])
        
        # Build properties for this fiber family (shared + individual)
        family_properties = self.shared_properties.copy()
        family_properties['mass_fraction'] = family_mass_fraction
        family_properties['constitutive_model'] = family_props['constitutive_model']
        
        # Create as single constituent
        full_name = f"{self.name}_{family_name}"
        fiber_family = SingleConstituent.from_parameters(full_name, family_properties)
        
        print(f"          {family_name}: mass fraction = {family_mass_fraction:.4f}")
        
        return fiber_family

    def _initialize_total_kinetics_histories(self) -> None:
        """Initialize total kinetics histories by summing across fiber families.
        
        Assumes:
        - Each fiber family has already initialized its own kinetics
        - rhoR_alpha_history[0] exists for this multi-fiber constituent
        
        The multi-fiber constituent's total kinetics values are the sum
        of individual fiber family kinetics at t=0.
        """
        # Check if any fiber family has kinetics
        if not self.fiber_families or self.fiber_families[0].kinetics is None:
            print(f"        No kinetics for multi-fiber constituent")
            return
        
        print(f"        Total homeostatic kinetics:")
        
        # Sum kinetics values from all fiber families at t=0
        total_k_alpha_h = sum(
            family.k_alpha_history[0] for family in self.fiber_families
        )
        total_mR_alpha_h = sum(
            family.mR_alpha_history[0] for family in self.fiber_families
        )
        
        # Initialize histories
        self.k_alpha_history.append(total_k_alpha_h)
        self.mR_alpha_history.append(total_mR_alpha_h)
        
        # Survival history for multi-fiber is aggregate (same as fiber families)
        self.survival_history.append([1.0])
        
        print(f"          k_alpha(0) = {total_k_alpha_h:.6f} 1/day")
        print(f"          mR_alpha(0) = {total_mR_alpha_h:.6f} kg/(m³·day)")
        print(f"          q(0,0) = 1.0")

    # =========================================================================
    # COMPUTATION: MASS DENSITY
    # =========================================================================
    
    def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess referential mass density for target timestep using specified method."""
        
        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        
        if target_timestep == 0:
            raise ValueError("Cannot guess timestep 0: should be initialized with homeostatic value")
        
        if guess_method == "from_previous_timestep":
            # Check if target timestep already exists
            if self._timestep_exists(target_timestep):
                warnings.warn(f"Target timestep {target_timestep} already exists. Returning existing value.")
                return self.get_rhoR_alpha(target_timestep)
            
            previous_timestep = target_timestep - 1
            
            # Check if previous timestep exists for this multi-fiber constituent
            if not self._timestep_exists(previous_timestep):
                raise ValueError(f"Cannot guess timestep {target_timestep}: previous timestep {previous_timestep} does not exist")
            
            # Guess for each fiber family first
            for family in self.fiber_families:
                family.guess_rhoR_alpha(target_timestep, guess_method)
            
            # Now update our total density (append, not set)
            total_density = sum(family.get_rhoR_alpha(target_timestep) for family in self.fiber_families)
            self.rhoR_alpha_history.append(total_density)
            
        else:
            raise ValueError(f"Unknown guess method: {guess_method}")
        
        return self.get_rhoR_alpha(target_timestep)

    def compute_rhoR_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute rhoR_alpha for all fiber families and update total."""
        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        if target_timestep == 0:
            raise ValueError("Cannot compute timestep 0: should be initialized with homeostatic value")

        # Compute for each fiber family
        for family in self.fiber_families:
            family.compute_rhoR_alpha(
                target_timestep, 
                dt, 
                integration_method,
                survival_function_computation
            )

        # Sum total and append
        total_density = sum(family.get_rhoR_alpha(target_timestep) for family in self.fiber_families)
        self.rhoR_alpha_history.append(total_density)

        return self.get_rhoR_alpha(target_timestep)

    # =========================================================================
    # COMPUTATION: STRESS
    # =========================================================================

    def get_stress(self):
        """Get total intramural stress from all fiber families."""
        return sum(family.get_stress() for family in self.fiber_families)
        
    def compute_stress(self, current_timestep):
        """Compute total stress from all fiber families."""
        total_stress = sum(family.compute_stress(current_timestep) for family in self.fiber_families)
        return total_stress

    # =========================================================================
    # PROPERTY ACCESSORS
    # =========================================================================
      
    def get_fiber_families(self):
        """Return list of individual fiber families."""
        return self.fiber_families