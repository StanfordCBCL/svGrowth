================================================================================
svGrowth Function List by Class
================================================================================

================================================================================
1. Configuration
================================================================================

Core Functions:
- from_parameters() - Create and initialize from YAML parameters
- add_layer() - Add layer to configuration
- compute_all_rhoR_alpha() - Compute mass densities for all constituents
- compute_all_stress() - Compute stresses for all layers
- guess_all_rhoR_alpha() - Guess mass densities for next timestep
- guess_geometry() - Guess geometry for all layers
- guess_stress_and_wss() - Guess stress and WSS for all layers

Helper Functions:
- _validate_parameters() - Validate parameter structure
- _enforce_layer_interactions() - Enforce constraints between layers (future)

================================================================================
2. Layer
================================================================================

Core Functions:
- from_parameters() - Create and initialize from parameters
- set_kinematics() - Set kinematics type (thin-wall/thick-wall)
- add_constituent() - Add constituent to layer
- compute_homeostatic_stress_direct() - Compute homeostatic stress from constituent properties
- compute_all_rhoR_alpha() - Compute mass densities for all constituents
- compute_all_stress() - Compute total layer stress from constituent heredity integrals
- guess_all_rhoR_alpha() - Guess mass densities for constituents
- guess_geometry() - Guess geometry for next timestep
- guess_stress_and_wss() - Guess stress and WSS

Property Accessors:
- get_inner_radius() - Get inner radius at timestep
- get_thickness() - Get thickness at timestep
- get_axial_stretch() - Get axial stretch at timestep
- get_density() - Get mass density at timestep
- get_stress() - Get stress tensor at timestep
- get_mid_radius() - Get mid-wall radius at timestep
- get_outer_radius() - Get outer radius at timestep
- get_deformation_gradient() - Get F at timestep
- get_volume_ratio() - Get J at timestep
- get_homeostatic_inner_radius() - Get homeostatic inner radius
- get_homeostatic_thickness() - Get homeostatic thickness
- get_homeostatic_axial_stretch() - Get homeostatic axial stretch
- get_homeostatic_density() - Get homeostatic density
- get_homeostatic_deformation_gradient() - Get homeostatic F
- get_geometry() - Get complete geometry at timestep
- get_homeostatic_geometry() - Get complete homeostatic geometry

Helper Functions:
- _initialize_homeostatic_geometry() - Initialize homeostatic geometry and histories
- _ensure_kinematics_set() - Ensure kinematics initialized
- _ensure_homeostatic_initialized() - Ensure homeostatic values initialized
- get_stress_trace() - Get trace of stress tensor (intramural stress)
- compute_volume_ratio() - Compute volume ratio J
- get_all_fiber_families() - Get all fiber families (unused)
- get_constituent_summary() - Get constituent summary (unused)

================================================================================
3. Constituent (Abstract Base)
================================================================================

Core Functions:
- from_parameters() - Factory to create constituent type
- compute_rhoR_alpha() - Compute mass density via heredity integral
- compute_sigma_alpha() - Compute constituent stress via heredity integral
- guess_rhoR_alpha() - Guess mass density for next timestep
- get_stress() - Get stress contribution (abstract)
- compute_stress() - Compute stress using constitutive model (abstract)

Helper Functions:
- get_rhoR_alpha() - Get mass density at timestep
- set_rhoR_alpha() - Set mass density at timestep (unused)
- _timestep_exists() - Check if timestep exists in history

================================================================================
4. SingleConstituent
================================================================================

Core Functions:
- from_parameters() - Create and initialize single constituent
- compute_rhoR_alpha() - Compute mass density (with or without kinetics)
- compute_sigma_alpha() - Compute stress (passive + active if applicable)
- compute_active_radius() - Compute active radius via heredity integral
- compute_active_stress_component() - Compute circumferential active stress
- guess_rhoR_alpha() - Guess mass density for next timestep
- get_stress() - Get stress contribution
- compute_stress() - Compute stress using constitutive model

Helper Functions:
- _initialize_deposition_stretch() - Initialize G_alpha tensor
- _initialize_homeostatic_kinetics() - Initialize kinetics histories at t=0
- _initialize_active_properties() - Initialize active stress properties
- _compute_active_stress_at_homeostasis() - Compute sigma_act at t=0

================================================================================
5. MultiFiberFamilyConstituent
================================================================================

Core Functions:
- from_parameters() - Create and initialize multi-fiber family
- compute_rhoR_alpha() - Compute total mass density from all families
- guess_rhoR_alpha() - Guess total mass density
- get_stress() - Get total stress from all families
- compute_stress() - Compute total stress from all families
- get_fiber_families() - Return list of fiber families (unused)

Helper Functions:
- _validate_mass_fraction_ratios() - Validate fiber family ratios sum to 1
- _create_fiber_family() - Create individual fiber family
- _initialize_total_kinetics_histories() - Initialize total kinetics from families

================================================================================
6. Mechanics
================================================================================

Core Functions:
- compute_F_alpha_for_cohort() - Compute F_alpha(s,tau) for single cohort
- compute_F_alpha_for_all_cohorts() - Compute F_alpha for all cohorts
- compute_S_hat_alpha_for_cohort() - Compute PK2 stress for single cohort
- compute_S_hat_alpha_for_all_cohorts() - Compute PK2 stress for all cohorts
- compute_sigma_hat_alpha_for_cohort() - Compute Cauchy stress for single cohort
- compute_sigma_hat_alpha_for_all_cohorts() - Compute Cauchy stress for all cohorts
- integrate_constituent_stress() - Integrate stress heredity integral

Helper Functions:
- compute_stress_trace() - Compute tr(sigma) for intramural stress

================================================================================
7. Kinetics
================================================================================

Core Functions:
- from_parameters() - Factory to create kinetics from YAML
- compute_k_alpha() - Compute degradation rate
- compute_production_rate() - Compute production rate
- compute_survival_function() - Compute survival q(s,tau) for all cohorts
- compute_heredity_integral() - Compute rho = integral(mR*q dτ)

Degradation Rate Function Classes:
- ConstantDegradationRate.compute_k_alpha() - Constant k_alpha
- QuadraticDegradationRate.compute_k_alpha() - Quadratic stimulus-modulated k_alpha

Survival Function Classes:
- ExponentialSurvival.compute_survival_function() - Exponential survival

Production Rate Function Classes:
- LinearProductionRate.compute_production_rate() - Linear stimulus-modulated mR

Helper Functions:
- _compute_delta() - Module-level helper for stimulus deviation

================================================================================
8. DeformationKinematics (Abstract Base)
================================================================================

Core Functions:
- compute_deformation_gradient() - Compute F (abstract)
- compute_lagrange_multiplier() - Compute constraint multiplier (abstract)
- compute_volume_ratio() - Compute J (abstract)
- get_component_name() - Get coordinate component name (abstract)

Helper Functions:
- compute_stretches() - Extract principal stretches from F
- compute_right_cauchy_green() - Compute C = F^T*F
- compute_left_cauchy_green() - Compute B = F*F^T (unused)
- compute_green_strain() - Compute E = 0.5(C-I) (unused)
- compute_invariants() - Compute I1, I2, I3 (unused)
- compute_jacobian() - Compute J = det(F)
- verify_incompressibility() - Check if J approximately equals 1 (unused)

================================================================================
9. ThinWallKinematics
================================================================================

Core Functions:
- compute_deformation_gradient() - Compute F for thin-wall vessel
- compute_lagrange_multiplier() - Compute constraint lambda for sigma_r = 0
- compute_stress_from_equilibrium() - Compute sigma_theta from equilibrium
- compute_wss() - Compute wall shear stress (Poiseuille flow)
- compute_volume_ratio() - Compute J from density or F
- get_component_name() - Get cylindrical component name (r/theta/z)

Helper Functions:
- compute_circumferential_stretch() - Compute lambda_theta
- compute_radial_stretch() - Compute lambda_r from incompressibility
- compute_thickness_from_incompressibility() - Compute h from lambda_theta, lambda_z
- compute_all_kinematic_quantities() - Compute all quantities for debugging (unused)
- compute_axial_force() - Compute required axial force (unused)

================================================================================
10. ThickWallKinematics (Future Implementation)
================================================================================

Entire class unused - future thick-wall vessel implementation

================================================================================
11. ConstitutiveModel (Abstract Base)
================================================================================

Core Functions:
- from_parameters() - Factory to create model from YAML
- compute_PK2_stress() - Compute S (abstract)
- compute_strain_energy() - Compute psi (abstract)
- is_isotropic() - Check if isotropic
- is_anisotropic() - Check if anisotropic

Helper Functions:
- set_fiber_orientation() - Set fiber angle for anisotropic models
- project_stress_to_fiber() - Project stress onto fiber direction (static)
- rotate_tensor() - Rotate tensor by angle (static)

================================================================================
12. NeoHookeanModel
================================================================================

Core Functions:
- compute_PK2_stress() - Neo-Hookean S
- compute_strain_energy() - Neo-Hookean psi

Helper Functions:
- _validate_parameters() - Validate c parameter

================================================================================
13. FungExponentialModel
================================================================================

Core Functions:
- compute_PK2_stress() - Fung exponential S
- compute_strain_energy() - Fung exponential psi

Helper Functions:
- _validate_parameters() - Validate c1, c2 parameters

================================================================================
14. HolzapfelOgdenModel (Unused)
================================================================================

Entire class unused - future arterial wall model

================================================================================
15. TensorOperations
================================================================================

Core Functions:
- inverse() - Compute matrix inverse
- trace() - Compute trace
- right_cauchy_green() - Compute C = F^T*F
- jacobian() - Compute J = det(F)

Helper Functions (Mostly Unused):
- left_cauchy_green() - Compute B = F*F^T (unused)
- green_strain() - Compute E (unused)
- principal_stretches() - Extract stretches (unused directly)
- invariants() - Compute I1, I2, I3 (unused)

================================================================================
16. Integrators
================================================================================

Core Functions:
- IntegratorFactory.create() - Create integrator by name
- TrapezoidIntegrator.integrate() - Trapezoidal rule integration
- SimpsonIntegrator.integrate() - Simpson's rule integration

Helper Functions:
- IntegratorFactory.register() - Register new integrator (unused)
- IntegratorFactory.available_methods() - List methods (unused)
- create_integrator() - Convenience function

================================================================================
17. Survival Function Computation
================================================================================

Core Functions:
- SurvivalFunctionComputationFactory.create() - Create strategy
- BackwardSurvivalFunctionComputation.compute_survival() - Optimized backward computation
- NaiveSurvivalFunctionComputation.compute_survival() - Naive forward computation (unused)

Helper Functions:
- SurvivalFunctionComputationFactory.register() - Register strategy (unused)
- SurvivalFunctionComputationFactory.available_strategies() - List strategies (unused)

================================================================================
18. Simulation
================================================================================

Core Functions:
- run() - Run complete G&R simulation
- _initialize_timestep() - Advance one timestep

Helper Functions (Stubs/Unused):
- _setup_output() - Setup output files (stub)
- _write_step_data() - Write timestep data (stub)
- _cleanup_output() - Cleanup output (stub)
- save_final_configuration() - Save final state (unused)

================================================================================
19. IOHandler
================================================================================

Core Functions:
- load_parameters() - Load YAML parameter file
- save_parameters() - Save parameters to YAML
- setup_output_file() - Setup output file handle
- write_simulation_header() - Write header line
- write_simulation_step() - Write timestep data

Helper Functions:
- _check_file_extension() - Validate YAML extension

================================================================================
20. Mechanics/Kinetics Interface Classes
================================================================================

ConstituentMechanicsContext:
- get_stress_tensor() - Get constituent stress
- get_deformation_gradient() - Get layer F
- get_mass_density() - Get constituent density
- get_deposition_stretch() - Get G_alpha
- get_constitutive_model() - Get constitutive model
- get_homeostatic_density() - Get homeostatic density
- get_tau_min() - Get earliest deposition time
- get_production_rate() - Get mR at timestep
- get_survival_values() - Get q values

LayerMechanicsContext:
- get_stress_tensor() - Get layer stress
- get_deformation_gradient() - Get layer F
- get_mass_density() - Get layer density

ConstituentKineticsContext:
- get_stimulus() - Get stimulus value at timestep
- get_stimulus_homeostatic() - Get homeostatic stimulus
- get_rhoR_alpha() - Get constituent density
- get_survival() - Get survival value
- Helper methods for specific stimuli (stress, WSS, inflammation)

================================================================================
21. KinematicsFactory
================================================================================

- create() - Create kinematics by geometry type
- register() - Register new kinematics type (unused)
- available_types() - List available types (unused)

================================================================================
22. Perturbations Module (Entire module unused)
================================================================================

All perturbation classes and factory unused

================================================================================
END OF FUNCTION LIST
================================================================================