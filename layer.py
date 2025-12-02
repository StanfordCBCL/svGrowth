from typing import List, Dict, Any, Optional
from constituent import Constituent
from deformation_kinematics import DeformationKinematics, KinematicsFactory
from mechanics import Mechanics
from mechanics_interface import LayerMechanicsContext
import numpy as np

class Layer:
    def __init__(self, name):
        self.name = name
        self.constituents: List[Constituent] = []
        self.kinematics: Optional[DeformationKinematics] = None
        self.mechanics = Mechanics()
        
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
        
        self.wss_history: List[float] = []
        self.flow_rate_history: List[float] = []
    
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
        print(f"\nInitializing layer: {layer_data['layer_name']}")
        
        # Step 1: Create empty layer (just name)
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
        print(f"  Constituents ({len(constituents)} total):")
        
        for name, properties in constituents.items():
            constituent = Constituent.from_parameters(name, properties, layer)
            layer.add_constituent(constituent)
        
        return layer
    
    def set_kinematics(self, geometry_type: str) -> None:
        """Set kinematics for this layer."""
        try:
            self.kinematics = KinematicsFactory.create(geometry_type)
            print(f"    Kinematics set: {type(self.kinematics).__name__}")
        except ValueError as e:
            available_types = KinematicsFactory.available_types()
            raise ValueError(
                f"Invalid geometry type: '{geometry_type}'\n"
                f"Available types: {', '.join(available_types)}\n"
            ) from e
    
    def _ensure_kinematics_set(self) -> None:
        """Ensure kinematics has been set before use."""
        if self.kinematics is None:
            raise RuntimeError(
                f"Kinematics not set for layer '{self.name}'. "
                f"Call set_kinematics() or use from_parameters() factory."
            )

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
        
        print("  Homeostatic geometry:")
        
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
        self.wss_history = [self.homeostatic_wss]
        self.flow_rate_history = [self.homeostatic_flow_rate]

        # Print summary
        print(f"    Inner radius (a_h): {a_h} mm → {self.homeostatic_inner_radius:.6f} m")
        print(f"    Thickness (h_h): {h_h} mm → {self.homeostatic_thickness:.6f} m")
        print(f"    Axial stretch (λ_z_h): {self.homeostatic_axial_stretch}")
        print(f"    Reference density (ρ_h): {self.homeostatic_density} kg/m³")
        print(f"    Loading:")
        print(f"      P_h = {pressure_h} kPa → {self.homeostatic_pressure:.2f} Pa")
        print(f"      Q_h = {flow_rate_h} m³/day")
        print(f"      μ = {self.blood_viscosity*1000:.2f} mPa·s")        
        print(f"    Wall shear stress (computed):")
        print(f"      WSS_h = {self.homeostatic_wss:.4f} Pa")

    # ========================================================================
    # PROPERTY ACCESSORS (Access by timestep)
    # ========================================================================
    
    def get_inner_radius(self, timestep: int) -> float:
        """Get inner radius at specified timestep."""
        return self.inner_radius_history[timestep]
    
    def get_thickness(self, timestep: int) -> float:
        """Get thickness at specified timestep."""
        return self.thickness_history[timestep]
    
    def get_axial_stretch(self, timestep: int) -> float:
        """Get axial stretch at specified timestep."""
        return self.axial_stretch_history[timestep]
    
    def get_density(self, timestep: int) -> float:
        """Get mass density at specified timestep."""
        return self.rhoR_history[timestep]
    
    def get_stress(self, timestep: int) -> np.ndarray:
        """Get stress tensor at timestep (3×3 array)."""
        return self.stress_history[timestep].copy()
    
    def get_mid_radius(self, timestep: int) -> float:
        """Get mid-wall radius at specified timestep."""
        return self.get_inner_radius(timestep) + 0.5 * self.get_thickness(timestep)
    
    def get_outer_radius(self, timestep: int) -> float:
        """Get outer radius at specified timestep."""
        return self.get_inner_radius(timestep) + self.get_thickness(timestep)
    
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        """Get deformation gradient F at specified timestep."""
        return self.F_history[timestep].copy()

    # ========================================================================
    # HOMEOSTATIC VALUE ACCESSORS
    # ========================================================================
    
    def get_homeostatic_inner_radius(self) -> float:
        """Get homeostatic inner radius."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_inner_radius
    
    def get_homeostatic_thickness(self) -> float:
        """Get homeostatic thickness."""
        self._ensure_homeostatic_initialized()
        return self.homeostatic_thickness
    
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
        return self.homeostatic_F.copy()
    
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

    def guess_all_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess mass densities for all constituents at target timestep using specified method."""
        print(f"  Guessing rhoR_alpha for layer '{self.name}' at timestep {target_timestep} using '{guess_method}' method")
        
        rhoR = 0.0
        for constituent in self.constituents:
            rhoR_alpha = constituent.guess_rhoR_alpha(target_timestep, guess_method)
            rhoR += rhoR_alpha

        return rhoR

    def compute_all_rhoR_alpha(self, target_timestep, dt, integration_method, survival_function_computation):
        """Compute mass densities for all constituents at target timestep."""
        print(f"  Computing rhoR_alpha for layer '{self.name}' at timestep {target_timestep}")

        rhoR = 0.0
        for constituent in self.constituents:
            rhoR_alpha = constituent.compute_rhoR_alpha(
                target_timestep, 
                dt, 
                integration_method, 
                survival_function_computation
            )
            rhoR += rhoR_alpha      

        return rhoR

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
    
        print(f"    Guessing geometry for layer '{self.name}' at timestep {timestep}")
        print(f"      Method: {guess_method}")
    
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
    
        print(f"      Guessed: r_mid={guessed_mid_radius*1000:.4f} mm, "
              f"λ_z={guessed_axial_stretch:.4f}")
        print(f"      Temporary: h={guessed_thickness*1000:.4f} mm "
              f"(will be refined by incompressibility)")
        print(f"      Derived: a={guessed_inner_radius*1000:.4f} mm")
    
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
    
        print(f"      Initial F computed (det(F) = {np.linalg.det(F_initial):.6f})")

    def add_constituent(self, constituent):
        """Add constituent to layer and set back-reference."""
        self.constituents.append(constituent)

    # ========================================================================
    # GEOMETRY QUERIES
    # ========================================================================
    
    def get_geometry(self, timestep: int) -> Dict[str, float]:
        """Return layer geometry at specified timestep.
        
        Args:
            timestep: Timestep at which to get geometry
            
        Returns:
            Dictionary with geometry values
        """
        return {
            'inner_radius': self.get_inner_radius(timestep),
            'thickness': self.get_thickness(timestep),
            'mid_radius': self.get_mid_radius(timestep),
            'outer_radius': self.get_outer_radius(timestep),
            'axial_stretch': self.get_axial_stretch(timestep),
            'density': self.get_density(timestep),
            'stress': self.get_stress(timestep)
        }
    
    def get_homeostatic_geometry(self) -> Dict[str, float]:
        """Return homeostatic (reference) layer geometry."""
        self._ensure_homeostatic_initialized()
        
        return {
            'inner_radius': self.homeostatic_inner_radius,
            'thickness': self.homeostatic_thickness,
            'mid_radius': self.homeostatic_mid_radius,
            'outer_radius': self.homeostatic_outer_radius,
            'axial_stretch': self.homeostatic_axial_stretch,
            'density': self.homeostatic_density
        }
    
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

        
    # ========================================================================
    # STRESS ESTIMATION
    # ========================================================================

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
        self.stress_history.append(previous_stress)
        
        previous_wss = self.wss_history[target_timestep-1]
        self.wss_history.append(previous_wss)

        print(f"    Guessed stress and WSS for layer '{self.name}' at timestep {target_timestep}")

    # ========================================================================
    # HELPER METHODS
    # ========================================================================

    def get_stress_trace(self, timestep: int) -> float:
        """Get trace of stress tensor (intramural stress) at timestep."""
        context = LayerMechanicsContext(self)
        return self.mechanics.compute_stress_trace(context, timestep)
    
    def compute_volume_ratio(self, timestep: int) -> float:
        """Compute volume ratio using kinematics assumption."""
        self._ensure_kinematics_set()
        
        rhoR_h = self.homeostatic_density
        rhoR = self.rhoR_history[timestep]
        
        # Delegate to kinematics (may have different formulas)
        return self.kinematics.compute_volume_ratio(rhoR, rhoR_h)

    def compute_all_stress(self, target_timestep: int, dt: float,
                      integration_method: str,
                      survival_function_computation: str) -> None:
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
        
        print(f"  Computing stress for layer '{self.name}' at timestep {timestep}")
        
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
        # TODO: consider removing this print from here
        print(f"    Total stress (diagonal components):")
        for i in range(3):
            component_name = self.kinematics.get_component_name(i)
            stress_kPa = sigma[i, i] / 1000
            print(f"      σ_{component_name[0]} = {stress_kPa:.2f} kPa")

        # Verify constraint satisfaction (for thin-wall: σ_r should be ≈0)
        if abs(sigma[0, 0]) > 1e-6:
            print(f"    WARNING: Radial stress not zero! σ_r = {sigma[0,0]:.2e} Pa")

    def _compute_active_stress(self, timestep: int) -> float:
        """Compute active stress contribution from smooth muscle.
    
        σ_act = (ρ_α/ρ_h) · T_act · (1 - exp(-C²)) · λ_act · parab(λ_act)
    
        where:
        - C = CB - CS · (τ_w/τ_w_h - 1)  (vasomotor control)
        - λ_act = a / a_act  (active stretch)
        - parab = 1 - [(λ_m - λ_act)/(λ_m - λ_0)]²
        """
        # TODO: Implement active stress logic
        # This needs:
        # - Active constituent identification
        # - a_act history computation (from heredity integral)
        # - Vasomotor control parameters (CB, CS)
        # - Active contraction parameters (T_act, λ_m, λ_0)
        
        return 0.0  # Placeholder

    def compute_homeostatic_stress_direct(self) -> None:
        """Compute homeostatic stress directly from constituent properties.

        At homeostasis, bypass heredity integral and compute directly:
        σ_α = (ρ_α/ρ_h) · σ̂_α(G_α)

        This avoids integration issues at t=0.
        """
        print(f"  Computing homeostatic stress for layer '{self.name}' (direct method)")

        sigma_total = np.zeros((3, 3))

        for constituent in self.constituents:
            print(f"    {constituent.name}:")
        
            # Get constituent properties
            rho_alpha = constituent.homeostatic_referential_density
            rho_h = self.homeostatic_density
            G_alpha = constituent.deposition_stretch
        
            print(f"      ρ_α = {rho_alpha:.2f} kg/m³, ρ_h = {rho_h:.2f} kg/m³")
            print(f"      Mass fraction = {rho_alpha/rho_h:.4f}")
            print(f"      G_α diagonal = [{G_alpha[0,0]:.3f}, {G_alpha[1,1]:.3f}, {G_alpha[2,2]:.3f}]")
        
            # Compute F_α(0,0) = G_α
            F_alpha = G_alpha
            print(f"      F_α diagonal = [{F_alpha[0,0]:.3f}, {F_alpha[1,1]:.3f}, {F_alpha[2,2]:.3f}]")
        
            # Compute C_α
            C_alpha = F_alpha.T @ F_alpha
            print(f"      C_α diagonal = [{C_alpha[0,0]:.3f}, {C_alpha[1,1]:.3f}, {C_alpha[2,2]:.3f}]")
        
            # Get constitutive model
            model = constituent.constitutive_model
        
            # ✅ NEW: Compute S_hat_α from constitutive law
            S_hat_alpha = model.compute_PK2_stress(F_alpha)
            print(f"      S_hat_α (PK2 stress) diagonal = [{S_hat_alpha[0,0]:.2f}, {S_hat_alpha[1,1]:.2f}, {S_hat_alpha[2,2]:.2f}] Pa")
        
            # ✅ NEW: Compute σ̂_α (push-forward)
            # TODO: compute from the ratio of J from layer. Currently hardcoded to 1 for testing.
            # J_alpha = np.linalg.det(F_alpha)
            J_alpha = 1
            print(f"      J_α = {J_alpha:.4f}")
            
            sigma_hat_alpha = (1.0 / J_alpha) * (F_alpha @ S_hat_alpha @ F_alpha.T)
            print(f"      σ̂_α (Cauchy stress) diagonal = [{sigma_hat_alpha[0,0]:.2f}, {sigma_hat_alpha[1,1]:.2f}, {sigma_hat_alpha[2,2]:.2f}] Pa")
        
            # Weight by mass fraction (PASSIVE STRESS)
            sigma_alpha_passive = (rho_alpha / rho_h) * sigma_hat_alpha
            print(f"      σ_α (passive, weighted) diagonal = [{sigma_alpha_passive[0,0]/1000:.2f}, {sigma_alpha_passive[1,1]/1000:.2f}, {sigma_alpha_passive[2,2]/1000:.2f}] kPa")
            
            # Start with passive stress
            sigma_alpha = sigma_alpha_passive.copy()
        
            # Add active stress if present
            if constituent.is_active:
                sigma_act = constituent._compute_active_stress_at_homeostasis()
                
                # Weight by mass fraction
                sigma_act_weighted = (rho_alpha / rho_h) * sigma_act
                
                print(f"      σ_act(0) = {sigma_act/1000:.2f} kPa (unweighted)")
                print(f"      σ_act(0, weighted) = {sigma_act_weighted/1000:.2f} kPa")
                
                # Add to circumferential component
                sigma_alpha[1, 1] += sigma_act_weighted
        
            # Accumulate
            sigma_total += sigma_alpha
        
            print(f"      σ_α(0) TOTAL diagonal = [{sigma_alpha[0,0]/1000:.2f}, "
                f"{sigma_alpha[1,1]/1000:.2f}, {sigma_alpha[2,2]/1000:.2f}] kPa")

        # ✅ NEW: Print total before constraint
        print(f"\n    Total stress BEFORE constraint:")
        print(f"      σ_total diagonal = [{sigma_total[0,0]/1000:.2f}, "
            f"{sigma_total[1,1]/1000:.2f}, {sigma_total[2,2]/1000:.2f}] kPa")

        # Apply kinematic constraint
        lagrange = self.kinematics.compute_lagrange_multiplier(sigma_total)
        
        # ✅ NEW: Print Lagrange multiplier
        print(f"    Lagrange multiplier:")
        print(f"      λ diagonal = [{lagrange[0,0]/1000:.2f}, "
            f"{lagrange[1,1]/1000:.2f}, {lagrange[2,2]/1000:.2f}] kPa")
        
        sigma_constrained = sigma_total - lagrange

        # ✅ NEW: Print total after constraint
        print(f"    Total stress AFTER constraint:")
        print(f"      σ_constrained diagonal = [{sigma_constrained[0,0]/1000:.2f}, "
            f"{sigma_constrained[1,1]/1000:.2f}, {sigma_constrained[2,2]/1000:.2f}] kPa")

        # Store
        self.homeostatic_stress = sigma_constrained
        self.stress_history[0] = sigma_constrained.copy()

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