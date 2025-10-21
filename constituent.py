import math
import warnings
from typing import List, Dict, Any, Optional
from solid_mechanics import ConstitutiveModel 
from abc import ABC, abstractmethod
from kinetics import Kinetics  

class Constituent(ABC):
    """Abstract base class for all constituents."""
    
    def __init__(self, name):
        self.name = name
        self.layer = None  # Reference to parent layer (set by Layer.add_constituent)
        self.constitutive_model = None
        self.kinetics = None 
        
        # Referential mass density evolution (rhoR_alpha) - timestep-indexed history
        self.rhoR_alpha_history = []  # List of mass densities over time
        self.sigma_hat_history = []  # List of constituent-based partial stress over time TODO: perhaps make temporary array?
        self.stress_history = []  # List of stress over time
        self.wss_history = []  # List of strain energy density over time
        self.survival_history = [] # List of survival fractions over time (for each cohort)
        self.k_alpha_history = []  # List of degradation rates over time

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
    
    def _compute_and_store_k_alpha(self, target_timestep):
        """Compute k_alpha and store in history."""
        # Constituent provides STATE to Kinetics
        state = self._get_current_state(target_timestep)
        
        # Kinetics computes using its degradation function
        k_alpha = self.kinetics.compute_k_alpha(
            current_timestep=target_timestep,
            deposition_timestep=target_timestep,  # For newly deposited material
            state=state
        )
        
        # Constituent stores the result
        self.k_alpha_history.append(k_alpha)
        return k_alpha

    def _update_survival_history(self, target_timestep):
        """Update survival for all cohorts."""
        dt = 1.0  # TODO: Get from simulation
        
        # Update survival for each existing cohort (tau)
        for tau in range(len(self.survival_history)):
            # Constituent provides state history accessor
            state_history_func = self._get_state_history_function()
            
            # Kinetics computes survival using its survival function
            q = self.kinetics.compute_survival(
                deposition_timestep=tau,
                current_timestep=target_timestep,
                state_history=state_history_func,
                dt=dt
            )
            
            # Constituent stores the result
            self.survival_history[tau].append(q)
        
        # Initialize new cohort deposited at target_timestep
        self.survival_history.append([1.0])

    def _compute_and_store_production_rate(self, target_timestep, k_alpha):
        """Compute production rate and store in history."""
        # Constituent provides STATE to Kinetics
        state = self._get_current_state(target_timestep)
        
        # Kinetics computes using its production function
        mR_alpha = self.kinetics.compute_production_rate(
            current_timestep=target_timestep,
            k_alpha=k_alpha,  # Pass pre-computed k_alpha
            state=state
        )
        
        # Constituent stores the result
        self.mR_alpha_history.append(mR_alpha)
        return mR_alpha

     def _integrate_mass_evolution(self, target_timestep, mR_alpha):
        """Integrate mass evolution (heredity integral + production)."""
        # This is CONSTITUENT's responsibility - it knows its own history
        
        # Heredity integral: sum of survived mass from all cohorts
        heredity_integral = 0.0
        for tau in range(len(self.survival_history)):
            if tau <= target_timestep:
                rho_tau = self.rhoR_alpha_history[tau]
                survival_index = target_timestep - tau
                q = self.survival_history[tau][survival_index]
                heredity_integral += rho_tau * q
        
        # Add production (simplified - you may need to integrate this too)
        # new_density = heredity_integral + mR_alpha * dt
        new_density = heredity_integral  # For now, just heredity
        
        return new_density
       
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
    
    def compute_rhoR_alpha(self, target_timestep):
        """Compute referential mass density of constituent for target timestep."""

        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")

        if target_timestep == 0:
            raise ValueError("Cannot compute timestep 0: should be initialized with homeostatic value")

        # If no degradation/production - maintain previous referential mass density
        if self.kinetics is None:
            previous_rhoR_alpha = self.rhoR_alpha_history[target_timestep - 1]
            self.rhoR_alpha_history.append(previous_rhoR_alpha)
            return self.get_rhoR_alpha(target_timestep)
        
        # STEP 1: Compute k_alpha once (Constituent asks Kinetics for computation)
        k_alpha = self._compute_and_store_k_alpha(target_timestep)
        
        # STEP 2: Update survival history (Constituent asks Kinetics for survival)
        q = self._update_survival_history(target_timestep)
        
        # STEP 3: Compute production rate (Constituent asks Kinetics for production)
        mR_alpha = self._compute_and_store_production_rate(target_timestep, k_alpha)
        
        # STEP 4: Integrate (Constituent does the integration itself)
        rhoR_alpha = self._integrate_mass_evolution(target_timestep, mR_alpha)
        
        # STEP 5: Store result
        self.rhoR_alpha_history.append(rhoR_alpha)

        return self.get_rhoR_alpha(target_timestep)
    
    # ABSTRACT METHODS - Must be implemented by subclasses
    @abstractmethod
    def get_stress(self):
        """Get intramural stress contribution."""
        pass

    @abstractmethod
    def compute_stress(self, current_timestep):
        """Compute total stress contribution."""
        pass
    

class SingleConstituent(Constituent):
    """Single constituent (elastin, muscle, etc.)."""
    
    def __init__(self, name):
        super().__init__(name)
        self.params = {}
        self.constitutive_model = None
        self.degradation_rate = 0.0
        self.production_rate = 0.0
        self.stress_contribution = 0.0
    
    @classmethod
    def from_parameters(cls, name, properties):
        """Create and fully initialize single constituent from parameters."""
        print(f"    Initializing single constituent: {name}")
        
        constituent = cls(name)
        
        # Initialize constitutive model
        if 'constitutive_model' in properties:
            constituent.constitutive_model = ConstitutiveModel.from_parameters(properties['constitutive_model'])
            print(f"        Constitutive model: {constituent.constitutive_model}")
        
        # Initialize kinetics
        if 'kinetics' in properties:
            constituent.kinetics = Kinetics.from_parameters(properties['kinetics'])
        
        # Initialize homeostatic mass density (t=0)
        homeostatic_density = properties['mass_fraction'] * 1050.0
        constituent.rhoR_alpha_history.append(homeostatic_density)
        
        # Initialize active properties if present
        if 'active_properties' in properties:
            constituent._initialize_active_properties(properties['active_properties'])
        
        # Store all parameters
        constituent.params = properties.copy()
        
        return constituent
       
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


class MultiFiberFamilyConstituent(Constituent):
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
        for family_name, family_props in properties['fiber_families'].items():
            fiber_family = constituent._create_fiber_family(family_name, family_props)
            constituent.fiber_families.append(fiber_family)
        
        print(f"        Created {len(constituent.fiber_families)} fiber families")
        
        # Initialize total mass density history (sum of all families at t=0)
        total_homeostatic_density = sum(family.get_rhoR_alpha(0) for family in constituent.fiber_families)
        constituent.rhoR_alpha_history.append(total_homeostatic_density)
        
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

    def compute_rhoR_alpha(self, target_timestep):
        """Compute rhoR_alpha for all fiber families and update total."""
        if target_timestep < 0:
            raise ValueError(f"Target timestep {target_timestep} must be >= 0")
        if target_timestep == 0:
            raise ValueError("Cannot compute timestep 0: should be initialized with homeostatic value")

        # Compute for each fiber family
        for family in self.fiber_families:
            family.compute_rhoR_alpha(target_timestep)

        # Sum total and append
        total_density = sum(family.get_rhoR_alpha(target_timestep) for family in self.fiber_families)
        self.rhoR_alpha_history.append(total_density)

        return self.get_rhoR_alpha(target_timestep)
        
    def get_fiber_families(self):
        """Return list of individual fiber families."""
        return self.fiber_families