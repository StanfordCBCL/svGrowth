from turtle import title
import numpy as np
from custom_logging import get_logger, CustomLogger
from typing import List, Dict, Any, Optional
from constituent import Constituent
from deformation_kinematics import DeformationKinematics, KinematicsFactory
from mechanics import Mechanics
from mechanics_interface import LayerMechanicsContext
from solver import Solver, LayerEquilibriumContext
from perturbations import PerturbationManager, PerturbationFactory

class Layer:
    """Represents a tissue layer as a constrained mixture with structurally significant constituents."""

    def __init__(self, name):
        self.name = name
        self.constituents: List[Constituent] = []
        self.kinematics: Optional[DeformationKinematics] = None
        self.mechanics = Mechanics()
        self.perturbation_manager = PerturbationManager(scope='layer')
        self.logger: CustomLogger = get_logger(__name__)
        
        # Homeostatic geometry (= reference configuration in G&R)
        self.homeostatic_inner_radius: Optional[float] = None
        self.homeostatic_thickness: Optional[float] = None
        self.homeostatic_axial_stretch: Optional[float] = None
        self.homeostatic_density: Optional[float] = None
        self.homeostatic_pressure: Optional[float] = None
        self.homeostatic_flow_rate: Optional[float] = None
        self.homeostatic_wss: Optional[float] = None 
        self.blood_viscosity: Optional[float] = None

        # Derived homeostatic values
        self.homeostatic_mid_radius: Optional[float] = None
        self.homeostatic_outer_radius: Optional[float] = None

        # Homeostatic stress and deformation
        self.homeostatic_F: Optional[np.ndarray] = None  
        self.homeostatic_stress: Optional[np.ndarray] = None

        # Time histories
        self.rhoR_history: List[float] = []
        self.inner_radius_history: List[float] = []
        self.mid_radius_history: List[float] = []
        self.thickness_history: List[float] = []
        self.axial_stretch_history: List[float] = []
        self.J_history: List[float] = []  
        
        self.F_history: List[np.ndarray] = [] 
        self.stress_history: List[np.ndarray] = []
        
        self.pressure_history: List[float] = []
        self.flow_rate_history: List[float] = []
        self.wss_history: List[float] = []



    # =========================================================================
    # INITIALIZATION
    # =========================================================================
        
    @classmethod
    def from_parameters(cls, layer_data: Dict[str, Any]) -> 'Layer':
        """Create and fully initialize Layer from parameter dictionary.
        
        This is the main factory method that orchestrates layer creation:
        1. Create empty layer
        2. Initialize constituents
        3. Set kinematics from geometry type
        4. Initialize homeostatic geometry and stress


        Args:
            layer_data: Layer section from YAML
            params: Full parameter dictionary
            
        Returns:
            Fully initialized Layer instance
        """      
        # Step 1: Create empty layer
        layer = cls(name=layer_data['layer_name'])

        # STEP 2: Set kinematics from geometry type
        geometry = layer_data['geometry']
        geometry_type = geometry['type']
        layer.set_kinematics(geometry_type)

        # STEP 3: Initialize homeostatic geometry from parameters
        homeostatic_loading = layer_data['loading_variables']
        layer._initialize_homeostatic_geometry(
            a_h=geometry['a_h'],
            h_h=geometry['h_h'],
            lambda_z_h=geometry['lambda_z_h'],
            rhoR_h=layer_data['rhoR_h'],
            pressure_h=homeostatic_loading['P_h'],
            flow_rate_h=homeostatic_loading['Q_h'],
            blood_viscosity=homeostatic_loading['blood_viscosity']
        )

        # Step 4: Initialize constituents for this layer
        constituents = layer_data['constituents']        
        for name, properties in constituents.items():
            constituent = Constituent.from_parameters(name, properties, layer)
            layer.add_constituent(constituent)

        # Step 5: Initialize perturbations if any
        if 'perturbations' in layer_data:
            layer._initialize_perturbations(layer_data['perturbations'])
        
        return layer
    
    def set_kinematics(self, geometry_type: str) -> None:
        """Set kinematics for this layer."""
        # TODO: validation should happen at the inputfile leve, not when setting kinematics.
        try:
            self.kinematics = KinematicsFactory.create(geometry_type)
        except ValueError as e:
            available_types = KinematicsFactory.available_types()
            raise ValueError(
                f"Invalid geometry type: '{geometry_type}'\n"
                f"Available types: {', '.join(available_types)}\n"
            ) from e

    def _initialize_homeostatic_geometry(
        self,
        a_h: float,
        h_h: float,
        lambda_z_h: float,
        rhoR_h: float,
        pressure_h: float,
        flow_rate_h: float,
        blood_viscosity: float,
    ) -> None:
        """Initialize homeostatic (reference) geometry and start time histories.
        
        This method:
        1. Converts units
        2. Sets homeostatic values
        3. Computes derived values (mid_radius, outer_radius)
        4. Initializes time history arrays with t=0 values
        
        Args:
            a_h: Homeostatic inner radius (mm)
            h_h: Homeostatic thickness (mm)
            lambda_z_h: Homeostatic axial stretch (dimensionless)
            rhoR_h: Homeostatic mass density (kg/m³)
        """
        self._ensure_kinematics_set()
    
        # Unit conversions
        mm_to_m = 1.0e-3
        kPa_to_Pa = 1.0e3
        Pa_s_to_Pa_day = 86400.0  # Pa·s → Pa·day (86400 s/day)
        
        # Set homeostatic values
        self.homeostatic_inner_radius = a_h * mm_to_m
        self.homeostatic_thickness = h_h * mm_to_m
        self.homeostatic_axial_stretch = lambda_z_h
        self.homeostatic_density = rhoR_h
        self.homeostatic_pressure = pressure_h * kPa_to_Pa
        self.homeostatic_flow_rate = flow_rate_h # Store in m³/day (as input)
        self.blood_viscosity = blood_viscosity  # Pa·s

        # Compute derived homeostatic values
        self.homeostatic_mid_radius = (
            self.homeostatic_inner_radius + 0.5 * self.homeostatic_thickness
        )
        self.homeostatic_outer_radius = (
            self.homeostatic_inner_radius + self.homeostatic_thickness
        )
        
        # Homeostatic F = I (reference configuration)
        self.homeostatic_F = np.eye(3)

        # Initialize stress as zero tensor (will be computed from constituents later)
        # TODO: Compute from constituent stress.
        self.homeostatic_stress = np.zeros((3, 3))
        
        # Initialize WSS
        self.homeostatic_wss = self.kinematics.compute_wss(
            flow_rate_h,
            self.homeostatic_inner_radius,
            self.blood_viscosity*Pa_s_to_Pa_day
        )

        # Initialize time history arrays with t=0 (homeostatic) values
        self.rhoR_history = [self.homeostatic_density]
        self.inner_radius_history = [self.homeostatic_inner_radius]
        self.mid_radius_history = [self.homeostatic_mid_radius]
        self.thickness_history = [self.homeostatic_thickness]
        self.axial_stretch_history = [self.homeostatic_axial_stretch]
        self.J_history = [1.0]  # J = 1 at homeostasis (rhoR_0/rhoR_0)
        self.F_history = [self.homeostatic_F.copy()]
        self.stress_history = [self.homeostatic_stress.copy()]
        self.pressure_history = [self.homeostatic_pressure]
        self.flow_rate_history = [self.homeostatic_flow_rate]
        self.wss_history = [self.homeostatic_wss]

    def _initialize_perturbations(self, perturbations_config: Dict[str, Any]):
        """Initialize perturbations from configuration.
        
        Args:
            perturbations_config: Dict mapping variable names to perturbation configs
                Example: {'pressure': {'type': 'step', ...}, 'flow_rate': {...}}
        """        
        for variable_name, pert_config in perturbations_config.items():
            perturbation = PerturbationFactory.create(pert_config)
            self.perturbation_manager.add_perturbation(variable_name, perturbation)
       
    def _ensure_kinematics_set(self) -> None:
        """Ensure kinematics has been set before use."""
        if self.kinematics is None:
            raise RuntimeError(
                f"Kinematics not set for layer '{self.name}'. "
                f"Call set_kinematics() or use from_parameters() factory."
            )

    def _ensure_homeostatic_initialized(self) -> None:
        """Ensure homeostatic values have been initialized."""
        if self.homeostatic_inner_radius is None:
            raise RuntimeError(
                f"Homeostatic geometry not initialized for layer '{self.name}'. "
                f"Call initialize_homeostatic_geometry() first."
            )

    # ========================================================================
    # CONSTITUENT MANAGEMENT
    # ========================================================================

    def add_constituent(self, constituent):
        """Add constituent to layer and set back-reference."""
        self.constituents.append(constituent)

    def get_all_fiber_families(self) -> List[Constituent]:
        """Get all individual fiber families from all constituents."""
        all_families = []
        
        for constituent in self.constituents:
            if hasattr(constituent, 'get_fiber_families'):
                # Multi-fiber family constituent
                all_families.extend(constituent.get_fiber_families())
            else:
                # Single constituent - treat as single family
                all_families.append(constituent)
        
        return all_families
    
    def get_constituent_summary(self) -> Dict[str, Dict[str, Any]]:
        """Get summary of constituents and fiber families."""
        summary = {}
        for constituent in self.constituents:
            if hasattr(constituent, 'fiber_families'):
                # Multi-fiber family
                summary[constituent.name] = {
                    'type': 'multi_fiber_family',
                    'total_mass_fraction': constituent.total_mass_fraction,
                    'fiber_families': len(constituent.fiber_families)
                }
            else:
                # Single constituent
                summary[constituent.name] = {
                    'type': 'single',
                    'mass_fraction': constituent.params.get('mass_fraction', 0.0)
                }
        return summary

    # =========================================================================
    # PROPERTY ACCESSORS: HOMEOSTATIC STATE
    # =========================================================================
    
    def get_homeostatic_inner_radius(self) -> float:
        """Get homeostatic inner radius."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_inner_radius
    
    def get_homeostatic_thickness(self) -> float:
        """Get homeostatic thickness."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_thickness

    def get_homeostatic_mid_radius(self) -> float:
        """Get homeostatic mid-radius."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_mid_radius

    def get_homeostatic_axial_stretch(self) -> float:
        """Get homeostatic axial stretch."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_axial_stretch
    
    def get_homeostatic_density(self) -> float:
        """Get homeostatic density."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_density

    def get_homeostatic_deformation_gradient(self) -> np.ndarray:
        """Get homeostatic deformation gradient (identity matrix)."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_F.copy()

    def get_homeostatic_geometry(self) -> Dict[str, float]:
        """Return homeostatic (reference) layer geometry."""
        self._ensure_homeostatic_initialized()
        
        return {
            'inner_radius': self.homeostatic_inner_radius,
            'mid_radius': self.homeostatic_mid_radius,
            'outer_radius': self.homeostatic_outer_radius,
            'thickness': self.homeostatic_thickness,
            'axial_stretch': self.homeostatic_axial_stretch,
        }

    # ========================================================================
    # PROPERTY ACCESSORS: CURRENT STATE
    # ========================================================================
    
    def get_inner_radius(self, timestep: int) -> float:
        """Get inner radius at specified timestep."""
        return self.inner_radius_history[timestep]
    
    def get_thickness(self, timestep: int) -> float:
        """Get thickness at specified timestep."""
        return self.thickness_history[timestep]

    def get_mid_radius(self, timestep: int) -> float:
        """Get mid-wall radius at specified timestep."""
        return self.get_inner_radius(timestep) + 0.5 * self.get_thickness(timestep)
    
    def get_outer_radius(self, timestep: int) -> float:
        """Get outer radius at specified timestep."""
        return self.get_inner_radius(timestep) + self.get_thickness(timestep)
        
    def get_axial_stretch(self, timestep: int) -> float:
        """Get axial stretch at specified timestep."""
        return self.axial_stretch_history[timestep]
    
    def get_density(self, timestep: int) -> float:
        """Get mass density at specified timestep."""
        return self.rhoR_history[timestep]
    
    def get_stress(self, timestep: int) -> np.ndarray:
        """Get stress tensor at specified timestep."""
        return self.stress_history[timestep].copy()
    
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        """Get deformation gradient F at specified timestep."""
        return self.F_history[timestep].copy()

    def get_volume_ratio(self, timestep: int) -> float:
        """Get volume ratio J at specified timestep from history."""
        return self.J_history[timestep]
    
    def get_pressure(self, timestep: int) -> float:
        """Get pressure at specified timestep."""
        return self.pressure_history[timestep]
    
    def get_flow_rate(self, timestep: int) -> float:
        """Get flow rate at specified timestep."""
        return self.flow_rate_history[timestep]
    
    def get_wss(self, timestep: int) -> float:
        """Get wall shear stress at specified timestep."""
        return self.wss_history[timestep]

    def get_geometry(self, timestep: int) -> Dict[str, float]:
        """Return layer geometry at specified timestep."""
        return {
            'inner_radius': self.get_inner_radius(timestep),
            'mid_radius': self.get_mid_radius(timestep),
            'outer_radius': self.get_outer_radius(timestep),
            'thickness': self.get_thickness(timestep),
            'axial_stretch': self.get_axial_stretch(timestep)
        }

    # ========================================================================
    # INITIAL GUESSES
    # ========================================================================

    def guess_all_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess mass densities for all constituents at target timestep using specified method."""
        #print(f"  Guessing rhoR_alpha for layer '{self.name}' at timestep {target_timestep} using '{guess_method}' method")
        
        rhoR = 0.0
        for constituent in self.constituents:
            rhoR_alpha = constituent.guess_rhoR_alpha(target_timestep, guess_method)
            rhoR += rhoR_alpha
        
        self._update_history(self.rhoR_history, target_timestep, rhoR)
        #print(f"    Total guessed ρ = {rhoR:.2f} kg/m³")

        return rhoR

    def guess_loading_variables(self, timestep: int) -> None:
        """Guess loading variables for next timestep."""
        self._update_history(self.pressure_history, timestep, self.get_pressure(timestep-1))
        self._update_history(self.flow_rate_history, timestep, self.get_flow_rate(timestep-1))

    def guess_geometry(self, timestep: int, 
                  guess_method: str = "from_previous_timestep") -> None:
        """Guess geometry for next timestep AND compute initial deformation gradient.
    
        Guessing strategy:
        1. Guess mid_radius and axial_stretch from previous timestep
        2. Copy thickness as initial guess (temporary - will be refined)
        3. Derive inner_radius = r_mid - h/2
        4. Compute initial F from guessed geometry
    
        The thickness guess will be refined by refine_geometry_from_density()
        using the incompressibility constraint: J = rhoR_h / rhoR.
    
        Args:
            timestep: Target timestep to guess geometry for
            guess_method: How to make the guess ("from_previous_timestep")
        """
        self._ensure_kinematics_set()
    
        if timestep == 0:
            raise ValueError("Cannot guess geometry at timestep 0 (use homeostatic)")
    
        #print(f"    Guessing geometry for layer '{self.name}' at timestep {timestep}")
        #print(f"      Method: {guess_method}")
    
        # STEP 1: Make initial guess based on method
        if guess_method == "from_previous_timestep":
            previous_timestep = timestep - 1
        
            # Guess independent variables
            guessed_mid_radius = self.mid_radius_history[previous_timestep]
            guessed_axial_stretch = self.axial_stretch_history[previous_timestep]
        
            # Copy thickness (temporary - will be refined by incompressibility)
            guessed_thickness = self.thickness_history[previous_timestep]
        
            # Derive inner_radius from mid_radius
            guessed_inner_radius = guessed_mid_radius - 0.5 * guessed_thickness

            for constituent in self.constituents:
                if constituent.is_active:
                    guessed_active_radius = constituent.active_radius_history[previous_timestep]
                    constituent.active_radius_history.append(guessed_active_radius)
        
        else:
            raise ValueError(
                f"Unsupported guess_method: '{guess_method}'. "
                f"Supported: 'from_previous_timestep'"
            )
    
        # STEP 2: Store guessed geometry in histories
        self.inner_radius_history.append(guessed_inner_radius)
        self.mid_radius_history.append(guessed_mid_radius)
        self.thickness_history.append(guessed_thickness)
        self.axial_stretch_history.append(guessed_axial_stretch)
    
        #print(f"      Guessed: r_mid={guessed_mid_radius*1000:.4f} mm, "
        #      f"λ_z={guessed_axial_stretch:.4f}")
        #print(f"      Temporary: h={guessed_thickness*1000:.4f} mm "
        #      f"(will be refined by incompressibility)")
        #print(f"      Derived: a={guessed_inner_radius*1000:.4f} mm")
    
        # STEP 3: Compute initial deformation gradient F from guessed geometry
        current_geometry = {
            'inner_radius': guessed_inner_radius,
            'mid_radius': guessed_mid_radius,
            'thickness': guessed_thickness,
            'axial_stretch': guessed_axial_stretch
        }
        
        reference_geometry = self.get_homeostatic_geometry()
        
        F_initial = self.kinematics.compute_deformation_gradient(
            current_geometry, 
            reference_geometry
        )
        self.F_history.append(F_initial)
    
        #print(f"      Initial F computed (det(F) = {np.linalg.det(F_initial):.6f})")

        # STEP 4: Compute J from F (fundamental definition)
        #TODO: refactor?
        J_initial = self.kinematics.compute_volume_ratio(
            self.F_history,
            self.rhoR_history,
            self.homeostatic_density,
            timestep
        )
        self.J_history.append(J_initial)

        #print(f"      Initial F computed (det(F) = J = {J_initial:.6f})")
        
    def guess_stress_and_wss(self, target_timestep: int) -> None:
        """Guess stress and wall shear stress at target timestep from previous timestep.

        Used to provide initial stress estimate for kinetics computation.
        Will be refined after mass density is computed.
        
        Args:
            target_timestep: Timestep to guess stress for stress and WSS
        """
        if target_timestep == 0:
            raise ValueError("Cannot guess stress at timestep 0 (use homeostatic)")
        
        # Copy stress from previous timestep as initial guess
        previous_stress = self.stress_history[target_timestep-1].copy()
        self._update_history(self.stress_history, target_timestep, previous_stress)
        
        previous_wss = self.wss_history[target_timestep-1]
        # TODO: potential update to WSS from perturbations
        self._update_history(self.wss_history, target_timestep, previous_wss)

        #print(f"    Guessed stress and WSS for layer '{self.name}' at timestep {target_timestep}")

    # =========================================================================
    # COMPUTATION: MASS DENSITY
    # =========================================================================

    def compute_all_rhoR_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute mass densities for all constituents at target timestep."""
        #print(f"  Computing rhoR_alpha for layer '{self.name}' at timestep {target_timestep}")

        rhoR = 0.0
        for constituent in self.constituents:
            rhoR_alpha = constituent.compute_rhoR_alpha(
                target_timestep, 
                dt, 
                integration_method, 
                survival_function_computation
            )
            rhoR += rhoR_alpha      

        # Update total layer density (replace guessed value)
        self.rhoR_history[target_timestep] = rhoR
        #print(f"    Total computed ρ = {rhoR:.2f} kg/m³")

        return rhoR

    # =========================================================================
    # COMPUTATION: STRESS
    # =========================================================================

    def compute_all_stress(self, target_timestep: int, dt: float,
                      integration_method: str,
                      survival_function_computation: str) -> np.ndarray:  
        """Compute total layer stress from constituent heredity integrals.

        Algorithm:
        1. For each constituent α:
           a. Compute σ̂_α(s,τ) for all cohorts τ (constitutive model)
           b. Integrate: σ_α = ∫ (mR·q/ρ_h) · σ̂_α dτ
           c. Include active stress if constituent has active properties
        2. Sum constituent stresses: σ_constituents = Σ σ_α (tensor addition)
        3. Compute Lagrange multiplier to enforce kinematic constraint
        4. Apply constraint: σ = σ_constituents - λ·I

        Note: Each σ_α includes both passive (structural) and active (contractile)
        contributions from that constituent's constitutive model.

        Args:
            timestep: Target timestep to compute stress for
            dt: Time step size (days)
            integration_method: Integration method ('simpson' or 'trapezoidal')
            survival_function_computation: Survival strategy ('backward', 'naive')
        """
        self._ensure_kinematics_set()
        
        #print(f"  Computing stress for layer '{self.name}' at timestep {target_timestep}")
        
        # STEP 1: Sum stress from all constituents (passive + active)
        # Each constituent's stress includes:
        # - Passive: from constitutive law (all constituents)
        # - Active: from muscle contraction (example: smooth muscle)
        sigma_constituents = np.zeros((3, 3))
        
        for constituent in self.constituents:
            sigma_alpha = constituent.compute_sigma_alpha(
                target_timestep,
                dt,
                integration_method,
                survival_function_computation
            )
            # σ_constituents += σ_α (includes both passive and active)
            sigma_constituents += sigma_alpha

        # STEP 2: Compute Lagrange multiplier to enforce kinematic constraint
        # (e.g., thin-wall: σ_r = 0, thick-wall: equilibrium)
        lagrange_multiplier = self.kinematics.compute_lagrange_multiplier(
            sigma_constituents
        )

        # STEP 3: Apply constraint by subtracting Lagrange multiplier
        sigma = sigma_constituents - lagrange_multiplier

        # STEP 4: Store constrained stress (replace guessed stress)
        self.stress_history[target_timestep] = sigma

        # STEP 5: Print summary using coordinate-agnostic component names
        # TODO: consider removing this print from here or streamlining
        # print(f"    Total stress (diagonal components):")
        # for i in range(3):
        #     component_name = self.kinematics.get_component_name(i)
        #     stress_kPa = sigma[i, i] / 1000
        #     print(f"      σ_{component_name[0]} = {stress_kPa:.2f} kPa")

        # Verify constraint satisfaction (for thin-wall: σ_r should be ≈0)
        if abs(sigma[0, 0]) > 1e-6:
            print(f"    WARNING: Radial stress not zero! σ_r = {sigma[0,0]:.2e} Pa")
            
        return sigma

    def compute_homeostatic_stress_direct(self) -> None:
        """Compute homeostatic stress directly from constituent properties.

        At homeostasis, bypass heredity integral and compute directly:
        σ_α = (ρ_α/ρ_h) · σ̂_α(G_α)

        This avoids integration issues at t=0.
        """
        sigma_total = np.zeros((3, 3))

        for constituent in self.constituents:        
            # Get constituent properties
            rho_alpha = constituent.homeostatic_referential_density
            rho_h = self.homeostatic_density
            G_alpha = constituent.deposition_stretch
        
            # Compute F_α(0,0) = G_α
            F_alpha = G_alpha
        
            # Compute C_α
            C_alpha = F_alpha.T @ F_alpha
        
            # Get constitutive model
            model = constituent.constitutive_model
        
            # ✅ NEW: Compute S_hat_α from constitutive law
            S_hat_alpha = model.compute_PK2_stress(F_alpha)
        
            # ✅ NEW: Compute σ̂_α (push-forward)
            # TODO: compute from the ratio of J from layer. Currently hardcoded to 1 for testing.
            # J_alpha = np.linalg.det(F_alpha)
            J_alpha = 1
            
            sigma_hat_alpha = (1.0 / J_alpha) * (F_alpha @ S_hat_alpha @ F_alpha.T)

            # Store sigma_hat_alpha at homeostasis for use in heredity integrals
            constituent.sigma_hat_history.append([sigma_hat_alpha])

            # Weight by mass fraction (PASSIVE STRESS)
            sigma_alpha_passive = (rho_alpha / rho_h) * sigma_hat_alpha
 
            # Start with passive stress
            sigma_alpha = sigma_alpha_passive.copy()
        
            # Add active stress if present
            if constituent.is_active:
                sigma_act = constituent._compute_active_stress_at_homeostasis()
                
                # Weight by mass fraction
                sigma_act_weighted = (rho_alpha / rho_h) * sigma_act
                
                # Add to circumferential component
                sigma_alpha[1, 1] += sigma_act_weighted

            # ✅ NEW: Store constituent stress at timestep 0
            constituent.stress_history.append(sigma_alpha)

            # Accumulate
            sigma_total += sigma_alpha

        # Apply kinematic constraint
        lagrange = self.kinematics.compute_lagrange_multiplier(sigma_total)
        
        sigma_constrained = sigma_total - lagrange

        # Store
        self.homeostatic_stress = sigma_constrained
        self.stress_history[0] = sigma_constrained.copy()

    def get_stress_trace(self, timestep: int) -> float:
        """Get trace of stress tensor (intramural stress) at timestep."""
        context = LayerMechanicsContext(self)
        return self.mechanics.compute_stress_trace(context, timestep)

    # =========================================================================
    # COMPUTATION: VOLUME RATIO
    # =========================================================================

    def compute_volume_ratio(self, target_timestep: int) -> float:
        """Compute volume ratio J at specified timestep.
        
        Delegates to kinematics, which decides whether to use F or density
        based on data availability and kinematics-specific logic.
        
        Args:
            target_timestep: Timestep at which to compute J
            
        Returns:
            Volume ratio J (dimensionless)
        """
        #TODO: consider removing this method and calling kinematics directly?
        self._ensure_kinematics_set()
        
        return self.kinematics.compute_volume_ratio(
            F_history=self.F_history,
            rhoR_history=self.rhoR_history,
            rhoR_h=self.homeostatic_density,
            timestep=target_timestep
        )

    # =========================================================================
    # PERTURBATIONS
    # =========================================================================

    def apply_perturbations(self, timestep, current_time):
        """Apply perturbations to the layer."""
            # Early exit if no active perturbations
        if not self.perturbation_manager.is_any_active(current_time):
            return

        perturbed_values = self.perturbation_manager.compute_perturbed_values(self, current_time)
        #TODO: refactor how we update histories. Think if we should include this info in perturbations.py.
        self._update_histories_with_perturbed_values(timestep, current_time, perturbed_values)
 
    def _update_histories_with_perturbed_values(self, timestep: int, current_time: float,
                                                perturbed_values: Dict[str, float]):
        """Update history arrays with perturbed values.
        
        Layer knows its structure and how to update its histories.
        Explicit mapping: variable_name -> history_array.
        
        Args:
            timestep: Timestep index
            current_time: Physical time (for logging)
            perturbed_values: Dict of variable_name -> perturbed_value
        """
        # Pressure
        if 'pressure' in perturbed_values:
            while len(self.pressure_history) <= timestep:
                self.pressure_history.append(self.homeostatic_pressure)
            self.pressure_history[timestep] = perturbed_values['pressure']
        
        # Flow rate
        if 'flow_rate' in perturbed_values:
            while len(self.flow_rate_history) <= timestep:
                self.flow_rate_history.append(self.homeostatic_flow_rate)
            self.flow_rate_history[timestep] = perturbed_values['flow_rate']
        
        # Axial stretch
        if 'axial_stretch' in perturbed_values:
            while len(self.axial_stretch_history) <= timestep:
                self.axial_stretch_history.append(self.homeostatic_axial_stretch)
            self.axial_stretch_history[timestep] = perturbed_values['axial_stretch']
        
        # Wall shear stress
        if 'wss' in perturbed_values:
            while len(self.wss_history) <= timestep:
                self.wss_history.append(self.homeostatic_wss)
            self.wss_history[timestep] = perturbed_values['wss']

    # =========================================================================
    # GEOMETRY: EQUILIBRIUM SOLVER INTERFACE
    # =========================================================================

    def compute_equilibrium_residual(
        self,
        trial_value: float,
        timestep: int,
        dt: float,
        integration_method: str,
        survival_function_computation: str,
    ) -> float:
        """Compute equilibrium residual for trial geometry value.
        
        Executes complete equilibrium workflow:
        1. Update geometry from incompressibility (uses existing method)
        2. Update WSS (depends on inner radius)
        3. Compute mixture stress (via constituent heredity integrals)
        4. Compute theoretical stress (from equilibrium condition)
        5. Return residual: σ_θθ(mixture) - σ_θθ(theoretical)
        
        This method is called by solver for each trial value.
        
        Args:
            trial_value: Trial value for geometry variable
            timestep: Current timestep
            dt: Time step size
            integration_method: Integration method for heredity integrals
            survival_function_computation: Survival computation method
            
        Returns:
            Residual (float) - should be zero at equilibrium
        """
            
        # STEP 1: Set trial value for the geometry variable
        # TODO: Consider moving out of this function to avoid redundant updates
        self.update_geometry_from_trial(
            trial_value, 
            timestep, 
            trial_geometry_variable=self.kinematics.get_equilibrium_variable()
        )
        
        # STEP 2: Update WSS (depends on inner radius)
        self.update_wss(timestep)
        
        # STEP 3: Compute mixture stress from constituents
        sigma_mixture =self.compute_all_stress(
            timestep, dt, integration_method, survival_function_computation
        )
        sigma_theta_mixture = sigma_mixture[1, 1]  # Circumferential component

        # STEP 4: Compute theoretical stress from equilibrium
        P = self.get_pressure(timestep)  
        a = self.get_inner_radius(timestep)
        h = self.get_thickness(timestep)
        
        sigma_theoretical = self.kinematics.compute_stress_from_equilibrium(P, a, h)
        
        # print(f"    [RESIDUAL DEBUG] σ_theoretical = {sigma_theoretical:.2f} Pa, σ_mixture = {sigma_theta_mixture:.2f} Pa")
        # print(f"    [RESIDUAL DEBUG] Residual = {sigma_theta_mixture - sigma_theoretical:.2f} Pa")
        
        return sigma_theta_mixture - sigma_theoretical     

    def solve_equilibrium_geometry(
        self,
        timestep: int,
        dt: float,
        integration_method: str,
        survival_function_computation: str,
        solver_method: str = 'brentq',
        tolerance: float = 1e-12,
        verbose: bool = False
    ) -> dict:
        """Solve for equilibrium geometry where mixture stress = theoretical stress.
        
        Creates generic context and delegates to solver.
        
        Args:
            timestep: Target timestep
            dt: Time step size
            integration_method: Integration method for heredity integrals
            survival_function_computation: Survival computation method
            solver_method: Root-finding method ('brentq', 'toms748', etc.)
            tolerance: Convergence tolerance
            verbose: Print solver details
            
        Returns:
            Dictionary with:
                - 'solution': Equilibrium geometry value
                - 'converged': Boolean convergence flag
                - 'iterations': Number of iterations
                - 'residual': Final residual
                - 'message': Status message
        """

        
        self._ensure_kinematics_set()
        
        if timestep == 0:
            raise ValueError(
                "Cannot solve equilibrium at timestep 0 (homeostatic state)"
            )
        
        # Create context (defined in solver.py)
        context = LayerEquilibriumContext(
            layer=self,
            timestep=timestep,
            dt=dt,
            integration_method=integration_method,
            survival_function_computation=survival_function_computation
        )
        
        # Solve
        solver = Solver(method=solver_method, tolerance=tolerance)
        result = solver.solve(context, verbose=verbose)
        
        if not result['converged']:
            print(f"  ⚠️  Layer '{self.name}': Equilibrium solver failed - "
                f"{result['message']}")
        
        return result   

    # =========================================================================
    # GEOMETRY: UPDATES
    # =========================================================================

    def update_geometry_from_volume_ratio(self, timestep: int, guess_variable: str = "mid_radius") -> None:
        """Update geometry based on incompressibility constraint using volume ratio J.
        
        Uses J = ρ_h / ρ(t) to refine geometry guess after mass density computation.
        
        Args:
            timestep: Target timestep to update geometry for
            guess_variable: Which variable to hold constant during refinement
                        - "mid_radius": Hold mid-radius constant, update inner_radius and thickness
                        - "inner_radius": [FUTURE] Hold inner radius constant
        """
        self._ensure_kinematics_set()
        
        #TODO: Refactor with factory pattern, consider moving to kinematics?
        #TODO: Consider exposing guess_variable value as input?
        if timestep == 0:
            raise ValueError("Cannot update geometry at timestep 0 (homeostatic state)")
        
        if guess_variable == "mid_radius":
            J = self.get_volume_ratio(timestep)  
            F = self.get_deformation_gradient(timestep)
            homeostatic_thickness = self.get_homeostatic_thickness()
            mid_radius = self.get_mid_radius(timestep)

            thickness = self.kinematics.compute_thickness_from_incompressibility(
                J=J,
                reference_thickness=homeostatic_thickness,
                F=F
            )
            inner_radius = self.kinematics.compute_inner_radius(mid_radius, thickness)
                
            # Update geometry histories
            self.inner_radius_history[timestep] = inner_radius
            self.thickness_history[timestep] = thickness

            # TODO: add conversion factor class
            inner_radius_old = self.get_inner_radius(timestep)
            thickness_old = self.get_thickness(timestep)
            print(f"      h: {thickness_old*1000:.4f} → {thickness*1000:.4f} mm")
            print(f"      a: {inner_radius_old*1000:.4f} → {inner_radius*1000:.4f} mm")

        elif guess_variable == "inner_radius":
            # TODO: [FUTURE] Implement: hold inner radius constant, update thickness
            raise NotImplementedError(
                f"Geometry update with guess_variable='{guess_variable}' not yet implemented. "
                f"Currently only 'mid_radius' is supported."
            )
        else:
            raise ValueError(
                f"Unknown guess_variable: '{guess_variable}'. "
                f"Supported: 'mid_radius', 'inner_radius', 'outer_radius'"
            )
        
    def update_geometry_from_trial(
        self,
        trial_value: float, # this is trial mid_radius value # TODO: rename to trial_geometry_value
        timestep: int,
        trial_geometry_variable: str = 'mid_radius'
    ) -> None:
        """Step 1: Update geometry based on trial value.
        
        Delegates to kinematics for geometry computation.
        """
        if trial_geometry_variable == 'mid_radius':
            # TODO: couple to deformation kinematics to generalize to other geometries eg thick wall or 3D
            # Update stretches from trial mid_radius
            homeostatic_mid_radius = self.get_homeostatic_mid_radius()
            lambda_theta_trial = trial_value / homeostatic_mid_radius
            lambda_z = self.get_axial_stretch(timestep)
            J = self.get_volume_ratio(timestep)
            lambda_r_trial = self.kinematics.compute_radial_stretch(
                J=J,
                lambda_theta=lambda_theta_trial,
                lambda_z=lambda_z
            )

            # Update thickness and inner radius from trial mid_radius
            h_h = self.get_homeostatic_thickness()
            thickness = self.kinematics.compute_thickness_from_incompressibility(
                J=J, reference_thickness=h_h, lambda_theta=lambda_theta_trial, lambda_z=lambda_z
            )
            inner_radius = self.kinematics.compute_inner_radius(
                mid_radius=trial_value,
                thickness=thickness
            )

            # Update geometry histories
            self.mid_radius_history[timestep] = trial_value
            self.inner_radius_history[timestep] = inner_radius
            self.thickness_history[timestep] = thickness

            # Update deformation gradient from trial mid_radius
            F_trial = np.diag([lambda_r_trial, lambda_theta_trial, lambda_z])
            self.F_history[timestep] = F_trial
            #self.update_deformation_gradient(timestep)

        else:
            raise NotImplementedError(
                f"Geometry variable '{trial_geometry_variable}' not supported"
            )            

    def update_wss(self, timestep: int) -> None:
        """Update wall shear stress at given timestep based on current geometry.
        
        Recomputes WSS using current inner radius and flow rate.
        Uses Poiseuille flow formula: WSS = (4·μ·Q)/(π·a³)
        
        Args:
            timestep: Target timestep to update WSS for
        """
        self._ensure_kinematics_set()
        
        if timestep == 0:
            raise ValueError("Cannot update WSS at timestep 0 (homeostatic state already set)")
        
        inner_radius = self.get_inner_radius(timestep)
        flow_rate = self.get_flow_rate(timestep)
        
        wss = self.kinematics.compute_wss(
            flow_rate=flow_rate,
            inner_radius=inner_radius,
            blood_viscosity=self.blood_viscosity
        )
        
        self.wss_history[timestep] = wss
        #print(f"    Updated WSS for layer '{self.name}': {wss:.6f} Pa (a={inner_radius*1000:.6f} mm)")

    def update_deformation_gradient(self, timestep: int) -> None:
        """Update deformation gradient at given timestep based on current geometry.
        
        Recomputes deformation gradient using current inner radius and flow rate.
        """
        self._ensure_kinematics_set()
        
        if timestep == 0:
            raise ValueError("Cannot update deformation gradient at timestep 0 (homeostatic state already set)")

        current_geometry={'inner_radius': self.get_inner_radius(timestep), 
                        'thickness': self.get_thickness(timestep), 
                        'axial_stretch': self.get_axial_stretch(timestep)}

        reference_geometry={'inner_radius': self.get_homeostatic_inner_radius(), 
                            'thickness': self.get_homeostatic_thickness(), 
                            'axial_stretch': self.get_homeostatic_axial_stretch()}

        F = self.kinematics.compute_deformation_gradient(
            current_geometry=current_geometry,
            reference_geometry=reference_geometry
        )
        self.F_history[timestep] = F

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
            history: History list to update (e.g., rhoR_history)
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

    # ========================================================================
    # OUTPUT: DATA EXPORT
    # ========================================================================

    def to_dict(self, timestep: int, time: float) -> Dict[str, Any]:
        """Convert layer state to dictionary for serialization.
        
        Args:
            timestep: Timestep index
            time: Simulation time (days)
            
        Returns:
            Dictionary containing layer state
        """
        F = self.F_history[timestep]
        sigma = self.stress_history[timestep]
        
        return {
            'time': time,
            'a': self.get_inner_radius(timestep),
            'h': self.get_thickness(timestep),
            'rho': self.get_density(timestep),
            'lambda_z': self.get_axial_stretch(timestep),
            'J': self.get_volume_ratio(timestep),
            'F_radial': F[0,0],
            'F_circ': F[1,1],
            'F_axial': F[2,2],
            'stress': self.get_stress_trace(timestep),
            'stress_radial': sigma[0,0],
            'stress_circ': sigma[1,1],
            'stress_axial': sigma[2,2],
            'pressure': self.get_pressure(timestep),
            'flow_rate': self.get_flow_rate(timestep),
            'wss': self.get_wss(timestep)
        }

    # =========================================================================
    # OUTPUT: PRINT SUMMARY
    # =========================================================================

    def print_state(self, timestep: int = None, logger: Optional[CustomLogger] = None):
        """Print layer state.
        
        Args:
            timestep: Timestep to print. If None, prints homeostatic state.
            logger: Logger to use for printing.
                If None, uses own logger (standalone mode).
                If provided, uses that logger (nested mode).
            
        Examples:
            layer.print_state()                    # Homeostatic, standalone
            layer.print_state(timestep=10)         # At timestep 10, standalone
            layer.print_state(timestep=10, logger=config.logger)  # Nested in config
        """
        # Use provided logger (for nesting) or own logger (for standalone)
        if logger is None:
            logger = self.logger
    
        if timestep is None:
            box_title = f"Layer: {self.name}"
        else:
            box_title = f"Layer: {self.name} (timestep={timestep})"
    
        with logger.box(box_title):
            self._print_content(timestep)
    
    def _print_content(self, timestep: int):
        """Print layer content (shared by standalone and nested modes)."""
        self.logger.box_line(f"Kinematics: {self.kinematics.__class__.__name__}")
        self.logger.box_line(f"Constituents: {len(self.constituents)}")

        self._print_mass(timestep)       
        self._print_geometry(timestep)
        self._print_loading(timestep)
        self._print_stress(timestep)
        
        for constituent in self.constituents:
            self.logger.box_line()
            constituent.print_state(timestep=timestep, logger=self.logger)
            

    def _print_geometry(self, timestep: int):
        """Print geometry section using kinematics metadata.
        
        TODO: Refactor kinematics to use registry pattern for variable metadata.
        For now, this is a simplified implementation that will be replaced when
        deformation_kinematics.py is restructured to define component symbols,
        variable names, units, etc. in a systematic registry.
        """
        # Determine title
        if timestep is None:
            box_title = "Homeostatic Geometry"
            geometry = self.get_homeostatic_geometry()
        else:
            box_title = f"Geometry (t={timestep})"
            geometry = self.get_geometry(timestep)
        
        self.logger.box_section(box_title)
        
        # TODO: Replace with kinematics.get_geometry_display_info() once registry is implemented
        # For now, hardcode thin-wall variables
        
        # Inner radius
        a_mm = geometry['inner_radius'] * 1e3
        h_mm = geometry['thickness'] * 1e3
        lambda_z = geometry['axial_stretch']

        if timestep is None:
            self.logger.box_item_aligned("Inner radius (a_h):", f"{a_mm:.4f} mm", label_width=25)
            self.logger.box_item_aligned("Thickness (h_h):", f"{h_mm:.4f} mm", label_width=25)
            self.logger.box_item_aligned("Axial stretch (λ_z_h):", f"{lambda_z:.4f}", label_width=25)
        else:
            self.logger.box_item_aligned("Inner radius (a):", f"{a_mm:.4f} mm", label_width=25)
            self.logger.box_item_aligned("Thickness (h):", f"{h_mm:.4f} mm", label_width=25)
            self.logger.box_item_aligned("Axial stretch (λ_z):", f"{lambda_z:.4f}", label_width=25)
            
    def _print_loading(self, timestep: int):
        """Print loading variables."""
        # TODO: refactor
        if timestep is None:
            P_kPa = self.homeostatic_pressure * 1e-3
            Q = self.homeostatic_flow_rate
            wss_kPa = self.homeostatic_wss * 1e-3
            box_title = "Homeostatic Loading"
            suffix = "_h"
        else:
            P_kPa = self.get_pressure(timestep) * 1e-3
            Q = self.get_flow_rate(timestep)
            wss_kPa = self.get_wss(timestep) * 1e-3
            box_title = f"Loading (timestep={timestep})"
            suffix = ""
        
        self.logger.box_section(box_title)
        self.logger.box_item_aligned(f"Pressure (P{suffix}):", f"{P_kPa:.2f} kPa", label_width=25)
        self.logger.box_item_aligned(f"Flow rate (Q{suffix}):", f"{Q:.2f} m³/day", label_width=25)
        # TODO: check wss computation and units (currently computing shear rate for Lattore2018)
        #self.logger.box_item_aligned(f"WSS (τ_w{suffix}):", f"{wss_kPa:.2f} kPa", label_width=25)

    def _print_mass(self, timestep: int):
        """Print mass density."""
        if timestep is None:
            rho = self.homeostatic_density
            box_title = "Homeostatic Mass Density"
            suffix = "_h"
        else:
            rho = self.get_density(timestep)
            box_title = f"Mass Density (timestep={timestep})"
            suffix = ""
        
        self.logger.box_section(box_title)
        self.logger.box_item_aligned(f"Referential mass density (ρ{suffix}):", f"{rho:.2f} kg/m³", label_width=25)
        
    def _print_stress(self, timestep: int):
        """Print stress tensor."""
        #TODO: Refactor to use kinematics.get_component_symbol(i) once registry is implemented
        if timestep is None:
            sigma = self.homeostatic_stress
            box_title = "Homeostatic Intramural Stress (diagonal)"
        else:
            sigma = self.get_stress(timestep)
            box_title = f"Intramural Stress (t={timestep}, diagonal)"
        
        self.logger.box_section(box_title)
        for i in range(3):
            component_symbol = ["rr", "θθ", "zz"][i]
            stress_kPa = sigma[i, i] / 1000
            self.logger.box_item_aligned(f"σ_{component_symbol}:", f"{stress_kPa:.2f} kPa", label_width=25)

    # def compute_homeostatic_stress_direct(self) -> None:
    #     """Compute homeostatic stress directly from constituent properties.
    
    #     At homeostasis, bypass heredity integral and compute directly:
    #     σ_α = (ρ_α/ρ_h) · σ̂_α(G_α)
    
    #     This avoids integration issues at t=0.
    #     """
    #     print(f"  Computing homeostatic stress for layer '{self.name}' (direct method)")
    
    #     sigma_total = np.zeros((3, 3))
    
    #     for constituent in self.constituents:
    #         print(f"    {constituent.name}:")
        
    #         # Get constituent properties
    #         rho_alpha = constituent.homeostatic_referential_density
    #         rho_h = self.homeostatic_density
    #         G_alpha = constituent.deposition_stretch
        
    #         # Compute F_α(0,0) = G_α
    #         F_alpha = G_alpha
        
    #         # Get constitutive model
    #         model = constituent.constitutive_model
        
    #         # Compute S_hat_α from constitutive law
    #         S_hat_alpha = model.compute_PK2_stress(F_alpha)
        
    #         # Compute σ̂_α (push-forward)
    #         J_alpha = np.linalg.det(F_alpha)
    #         sigma_hat_alpha = (1.0 / J_alpha) * (F_alpha @ S_hat_alpha @ F_alpha.T)
        
    #         # Weight by mass fraction
    #         sigma_alpha = (rho_alpha / rho_h) * sigma_hat_alpha
        
    #         # Add active stress if present
    #         if constituent.is_active:
    #             sigma_act = constituent._compute_active_stress_at_homeostasis()
    #             sigma_alpha[1, 1] += (rho_alpha / rho_h) * sigma_act
    #             print(f"      σ_act(0) = {sigma_act/1000:.2f} kPa")
        
    #         # Accumulate
    #         sigma_total += sigma_alpha
        
    #         print(f"      σ_α(0) diagonal = [{sigma_alpha[0,0]/1000:.2f}, "
    #               f"{sigma_alpha[1,1]/1000:.2f}, {sigma_alpha[2,2]/1000:.2f}] kPa")
    
    #     # Apply kinematic constraint
    #     lagrange = self.kinematics.compute_lagrange_multiplier(sigma_total)
    #     sigma_constrained = sigma_total - lagrange
    
    #     # Store
    #     self.homeostatic_stress = sigma_constrained
    #     self.stress_history[0] = sigma_constrained.copy()
    
    #     #TODO: add get_component_name and get_component_symbol methods to kinematics
    #     # print(f"    Homeostatic stress (constrained):")
    #     # for i in range(3):
    #     #     component_name = self.kinematics.get_component_name(i)
    #     #     stress_kPa = sigma_constrained[i, i] / 1000
    #     #     print(f"      σ_{component_name[0]}(0) = {stress_kPa:.2f} kPa")