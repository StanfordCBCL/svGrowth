import math
import warnings
from typing import List, Dict, Any, Optional
from solid_mechanics import ConstitutiveModel 
from abc import ABC, abstractmethod
from kinetics import Kinetics  
from kinetics_interface import ConstituentKineticsContext
from solid_mechanics import ConstitutiveModel
from mechanics_interface import ConstituentMechanicsContext

# Forward reference for Layer (avoids circular import)
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from layer import Layer

class Constituent(ABC):
    """Abstract base class for all constituents."""
    
    def __init__(self, name):
        self.name = name
        self.layer: Optional['Layer'] = None  # Reference to parent layer (set by Layer.add_constituent)
        self.kinetics: Optional[Kinetics] = None
        self.constitutive_model: Optional[ConstitutiveModel] = None # TODO: consider making this mandatory for SingleConstituent. Never initialized for MultiFiberFamilyConstituent.
        
        # Homeostatic properties
        self.homeostatic_referential_density: Optional[float] = None

        # Referential mass density evolution (rhoR_alpha) - timestep-indexed history
        self.q_history = []  # List of survival function values over time (for each cohort)
        self.rhoR_alpha_history = []  # List of mass densities over time
        self.sigma_hat_history = []  # List of constituent-based partial stress over time TODO: perhaps make temporary array?
        self.stress_history = []  # List of stress over time
        self.wss_history = []  # List of strain energy density over time
        self.survival_history = [] # List of survival fractions over time (for each cohort)
        self.k_alpha_history = []  # List of degradation rates of survival function over time
        self.mR_alpha_history = []  # List of production rates over time

    @classmethod
    def from_parameters(cls, name, properties):
        """Factory method to create appropriate constituent type."""
        constituent_type = properties.get('constituent_type', 'single')
        
        if constituent_type == 'multi_fiber_family':
            return MultiFiberFamilyConstituent.from_parameters(name, properties)
        else:
            return SingleConstituent.from_parameters(name, properties)
    
    # UTILITY/CONVENIENCE METHODS - Same implementation for all constituents
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

    # ABSTRACT METHODS - Must be implemented by subclasses
    @abstractmethod
    def get_stress(self):
        """Get intramural stress contribution."""
        pass

    @abstractmethod
    def compute_stress(self, current_timestep):
        """Compute total stress contribution."""
        pass

    #TODO: Consider moving this logic to the numerics/simulation layer, and use set_rhoR_alpha function for the prediction.   
    @abstractmethod
    def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Compute referential mass density of constituent for target timestep."""
        pass

    @abstractmethod
    def compute_rhoR_alpha(self, target_timestep, dt, integration_method, survival_strategy):
        """Compute referential mass density of constituent for target timestep."""
        pass

class SingleConstituent(Constituent):
    """Single constituent (elastin, muscle, etc.)."""
    
    def __init__(self, name):
        super().__init__(name)
        self.params = {}
        self.tau_min = 0
    
    @classmethod
    def from_parameters(cls, name, properties):
        """Create and fully initialize single constituent from parameters."""
        print(f"    Initializing single constituent: {name}")
        
        constituent = cls(name)
        
        # Initialize homeostatic mass density (t=0)
        constituent.homeostatic_referential_density = properties['mass_fraction'] * 1050.0
        constituent.rhoR_alpha_history.append(constituent.homeostatic_referential_density)
        print(f"        Homeostatic mass density (rhoR_alpha at t=0): {constituent.homeostatic_referential_density:.2f} kg/m³")

        # Initialize kinetics
        if 'kinetics' in properties:
            constituent.kinetics = Kinetics.from_parameters(properties['kinetics'])
            constituent._initialize_homeostatic_kinetics()
        else:
            print(f"        No kinetics (elastin-like constituent)")
        
        # TODO: Set tau_min based on constituent type
        
        # Initialize constitutive model
        if 'constitutive_model' in properties:
            constituent.constitutive_model = ConstitutiveModel.from_parameters(properties['constitutive_model'])
            print(f"        Constitutive model: {constituent.constitutive_model}")
        
        # Initialize active properties if present
        if 'active_properties' in properties:
            constituent._initialize_active_properties(properties['active_properties'])
        
        # Store all parameters
        constituent.params = properties.copy()
        
        return constituent
    
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
        
        # Initialize survival history
        # At t=0, we have one cohort deposited at τ=0 with q(0,0) = 1.0
        self.survival_history.append([1.0])
        
        # Print summary
        print(f"          k_alpha(0) = {k_alpha_h:.6f} 1/day")
        print(f"          mR_alpha(0) = {mR_alpha_h:.6f} kg/(m³·day)")
        print(f"          q(0,0) = 1.0")

    def _initialize_active_properties(self, active_props):
        """Initialize active properties."""
        print("        Active properties initialized")
        self.is_active = True
        self.T_act_h = active_props['T_act_h'] * 1000.0  # kPa -> Pa
        # ... other active properties
    
    def get_mass_density(self, current_timestep):
        """Get current mass density (convenience method)."""
        return self.get_rhoR_alpha(current_timestep)

    def get_stress(self):
        """Get intramural stress contribution."""
        return self.stress_contribution
        
    def compute_stress(self, current_timestep):
        """Compute stress using constitutive model."""
        if self.constitutive_model:
            deformation_gradient = 1.1  # Mock value
            current_density = self.get_rhoR_alpha(current_timestep)
            self.stress_contribution = self.constitutive_model.compute_stress(
                deformation_gradient, current_density
            )
        else:
            self.stress_contribution = 0.0
        
        return self.stress_contribution
    
    def update_kinetics(self, dt):
        """Update mass evolution."""
        self._update_production_rate()
        dmass_dt = self.production_rate - self.degradation_rate * self.mass_density
        self.mass_density += dmass_dt * dt
        self.history['mass_density'].append(self.mass_density)
    
    def _update_production_rate(self):
        """Update production rate based on stimuli."""
        stimulus = 1.0  # Mock
        base_production = self.params.get('stress_production_gain', 0.0)
        self.production_rate = base_production * stimulus

    #TODO: Consider moving this logic to the numerics/simulation layer, and use set_rhoR_alpha function for the prediction.   
    def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess referential mass density for target timestep."""
          
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
            self.rhoR_alpha_history.append(previous_rhoR_alpha)
            return self.get_rhoR_alpha(target_timestep)
        
        # Pass necessary data from constituent to kinetics class
        context = ConstituentKineticsContext(self)

        # STEP 1: Compute k_alpha at target timestep (Constituent asks Kinetics for computation)
        k_alpha = self.kinetics.compute_k_alpha(context, target_timestep)
        self.k_alpha_history.append(k_alpha)

        # STEP 2: Compute production rate
        # TODO: improve fail-safe if we don't have a guess for rhoR_alpha.
        # Consider how this interacts with guess_rhoR_alpha (which is needed for this step).
        rhoR_alpha = self.get_rhoR_alpha(target_timestep) # comes from guess
        mR_alpha = self.kinetics.compute_production_rate(context, target_timestep, k_alpha, rhoR_alpha)
        self.mR_alpha_history.append(mR_alpha)

        # STEP 3: Compute survival function q(s, tau) for all cohorts
        survival_values = self.kinetics.compute_survival_function(
            k_alpha_history=self.k_alpha_history,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method,
            survival_function_computation=survival_function_computation
        )
        self.q_history.append(survival_values)

        # Step 4: Compute heredity integral (just integrate mR × q!)
        rhoR_alpha = self.kinetics.compute_heredity_integral(
            self.mR_alpha_history,
            survival_values,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method
        )
        
        self.rhoR_alpha_history[target_timestep] = rhoR_alpha

        return self.get_rhoR_alpha(target_timestep)

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

        # STEP 0: Compute survival function q(s, tau) for all cohorts
        # TODO: This assums we already have k_alpha and mR_alpha. Figure out where to store q.
        # Generally, q should be known as we first compute_rhoR_alpha before compute_sigma_alpha.
        # This is computed on-the-fly from k_alpha_history, not stored
        survival_values = self.kinetics.compute_survival_function(
            k_alpha_history=self.k_alpha_history,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method,
            survival_function_computation=survival_function_computation
        )

         # Step 1: Compute sigma hat_alpha(s, τ) for all cohorts
        sigma_hat_alpha = self.mechanics.compute_sigma_hat_alpha(
            self.mR_alpha_history,
            survival_values,
            sigma_hat_alpha = self.sigma_hat_history,
            tau_min=self.tau_min,
            dt=dt,  
            current_timestep=target_timestep,
            integration_method=integration_method
        )
        self.sigma_hat_alpha_history.append(sigma_hat_alpha)

        # Step 2: Compute heredity integral (just integrate (mR/J) × q x sigma_hat_alpha!)
        sigma_alpha = self.mechanics.compute_heredity_integral(
            self.mR_alpha_history / J,
            survival_values,
            sigma_hat_alpha = self.sigma_hat_history,
            tau_min=self.tau_min,
            dt=dt,
            current_timestep=target_timestep,
            integration_method=integration_method
        )
        self.sigma_alpha_history.append(sigma_alpha)

        return self.get_sigma_alpha(target_timestep)


class MultiFiberFamilyConstituent(Constituent):
    #TODO: To clean up after code works with single constituent.
    """Multi-fiber family constituent (e.g., collagen with multiple orientations)."""
    
    def __init__(self, name):
        super().__init__(name)
        self.fiber_families = []
        self.shared_properties = {}
        self.total_mass_fraction = 0.0
    
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
    
    def get_mass_density(self, current_timestep):
        """Get total mass density from all fiber families."""
        return sum(family.get_rhoR_alpha(current_timestep) for family in self.fiber_families)

    def get_stress(self):
        """Get total intramural stress from all fiber families."""
        return sum(family.get_stress() for family in self.fiber_families)
        
    def compute_stress(self, current_timestep):
        """Compute total stress from all fiber families."""
        total_stress = sum(family.compute_stress(current_timestep) for family in self.fiber_families)
        return total_stress
    
    def add_timestep_with_guess(self, current_timestep):
        """Add next timestep for all fiber families and update total."""
        # Add timestep for each fiber family
        for family in self.fiber_families:
            family.add_timestep_with_guess(current_timestep)
        
        # Update total mass density history
        total_density = sum(family.get_rhoR_alpha(current_timestep + 1) for family in self.fiber_families)
        self.rhoR_alpha_history.append(total_density)
    
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
        
    def get_fiber_families(self):
        """Return list of individual fiber families."""
        return self.fiber_families

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