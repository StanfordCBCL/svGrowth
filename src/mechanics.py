"""
Mechanics module for G&R simulations.

Handles:
1. Deformation gradient computation for constituents (F_α)
2. Stress computation for constituents (σ̂_α, σ_α)
3. Numerical integration of heredity integrals

Follows adapter pattern like Kinetics:
- Mechanics class is stateless (no data storage)
- Operates on MechanicsContext (adapter for data access)
- Constituent owns the stress/deformation histories
"""

from abc import ABC, abstractmethod
from typing import List, Optional, Tuple
import numpy as np
from mechanics_interface import MechanicsContext, ConstituentMechanicsContext
from tensor_operations import TensorOperations
from integrators import IntegratorFactory

# =============================================================================
# MAIN MECHANICS CLASS - Stateless stress computations
# =============================================================================

class Mechanics:
    """Main mechanics orchestrator for deformation and stress computations.
    
    Stateless - operates on MechanicsContext, never stores data.
    Results are returned to caller (Constituent) for storage.
    
    Uses TensorOperations for all tensor algebra.
    """
    
    def __init__(self):
        """Initialize with tensor operations utility."""
        self.tensor_ops = TensorOperations()
    
    # =========================================================================
    # F_α COMPUTATION 
    # =========================================================================
    
    def compute_F_alpha_for_cohort(
        self,
        context: ConstituentMechanicsContext,
        current_timestep: int,
        deposition_timestep: int
    ) -> np.ndarray:
        """Compute F_α(s,τ) for single cohort deposited at τ.
        
        This is the BUILDING BLOCK function that both kinetics and no-kinetics constituents use.
        
        Formula:
            F_α(s,τ) = F(s) · F(τ)⁻¹ · G_α
        
        Special case (no kinetics, τ=0):
            F_α(s,0) = F(s) · I⁻¹ · G_α⁻¹ = F(s) · G_α
        
        Args:
            context: Constituent mechanics context (provides layer F, G_α)
            current_timestep: Current time s
            deposition_timestep: Deposition time τ
            
        Returns:
            F_α(s,τ) as 3×3 matrix
        """
        # Get current layer deformation
        F_current = context.get_deformation_gradient(current_timestep)
        
        # Get deposition stretch (inverse)
        G_alpha = context.get_deposition_stretch()
        
        # OPTIMIZATION: If τ=0, F(τ=0) = I, so skip inversion
        if deposition_timestep == 0:
            # F_α = F(s) · I⁻¹ · G_α⁻¹ = F(s) · G_α
            F_alpha = F_current @ G_alpha
        else:
            # General case: F_α = F(s) · F(τ)⁻¹ · G_α
            F_tau = context.get_deformation_gradient(deposition_timestep)
            F_tau_inv = TensorOperations.inverse(F_tau)
            F_alpha = F_current @ F_tau_inv @ G_alpha
        
        return F_alpha
    
    def compute_F_alpha_for_all_cohorts(
        self,
        context: ConstituentMechanicsContext,
        current_timestep: int
    ) -> List[np.ndarray]:
        """Compute F_α for all active cohorts.
        
        Reuses compute_F_alpha_for_cohort() for each cohort.
        
        Args:
            context: Constituent mechanics context
            current_timestep: Current time s
            
        Returns:
            List [F_α(s,τ_min), ..., F_α(s,s)] (length: s - τ_min + 1)
        """
        tau_min = context.get_tau_min()
        
        F_alpha_cohorts = []
        
        for tau in range(tau_min, current_timestep + 1):
            F_alpha_tau = self.compute_F_alpha_for_cohort(
                context, current_timestep, tau
            )
            F_alpha_cohorts.append(F_alpha_tau)
        
        return F_alpha_cohorts

    # =========================================================================
    # S_hat_α COMPUTATION (delegates to constitutive model)
    # =========================================================================
    
    def compute_S_hat_alpha_for_cohort(
        self,
        context: ConstituentMechanicsContext,
        F_alpha: np.ndarray
    ) -> np.ndarray:
        """Compute S_hat_α (PK2 stress) for single cohort given F_α.
        
        Delegates to constitutive model for material-specific computation.
        
        Formula:
            S_hat_α = ∂ψ_α/∂F_α  (from constitutive model)
        
        Args:
            context: Constituent mechanics context (provides constitutive model)
            F_alpha: Deformation gradient for this cohort (3×3)
            
        Returns:
            S_hat_α as 3×3 PK2 stress tensor (Pa)
        """
        # Get constitutive model from context
        constitutive_model = context.get_constitutive_model()
        
        # Delegate to constitutive model
        S_hat_alpha = constitutive_model.compute_PK2_stress(F_alpha)
        
        return S_hat_alpha
    
    def compute_S_hat_alpha_for_all_cohorts(
        self,
        context: ConstituentMechanicsContext,
        F_alpha_cohorts: List[np.ndarray]
    ) -> List[np.ndarray]:
        """Compute S_hat_α (PK2 stress) for all cohorts.
        
        Reuses compute_S_hat_alpha_for_cohort() for each F_α.
        
        Args:
            context: Constituent mechanics context
            F_alpha_cohorts: List of F_α matrices
            
        Returns:
            List [S_hat_α(s,τ_min), ..., S_hat_α(s,s)]
        """
        # TODO: Consider changing inputs to timestep instead?

        S_hat_alpha_cohorts = []
        
        for F_alpha_tau in F_alpha_cohorts:
            S_hat_alpha_tau = self.compute_S_hat_alpha_for_cohort(
                context, F_alpha_tau
            )
            S_hat_alpha_cohorts.append(S_hat_alpha_tau)
        
        return S_hat_alpha_cohorts
        
    # =========================================================================
    # σ̂_α COMPUTATION 
    # =========================================================================
    
    def compute_sigma_hat_alpha_for_cohort(
        self,
        J: float,
        F_alpha: np.ndarray,
        S_hat_alpha: np.ndarray
    ) -> np.ndarray:
        """Compute σ̂_α (partial constituent stress) for single cohort given F_α and S_hat_α.
        
        This is the PUSH-FORWARD operation from PK2 to partial constituent stress.
        
        Formula:
            σ̂_α = (1/J_α) F_α · S_hat_α · F_α^T 
        
        where J_α = det(F_α)
        
        Args:
            J: Layer volume ratio at current time J(s) = ρ_h/ρ(s)
            F_alpha: Deformation gradient for this cohort (3×3)
            S_hat_alpha: PK2 stress from constitutive model (3×3)
            
        Returns:
            σ̂_α(s,τ) as 3×3 Cauchy stress tensor (Pa)
        """
        # Compute sigma_hat
        # TODO: call constituent or kinematics function for volume ratio J_alpha calc
        # TODO: simplify calculation: J_alpha = J(s)/J(τ), from layer kinematics
        # TODO: to understand why J_alpha = np.linalg.det(F_alpha). Add J tau back in.
        sigma_hat_alpha = (1.0 / J) * (F_alpha @ S_hat_alpha @ F_alpha.T)
        
        return sigma_hat_alpha
    
    def compute_sigma_hat_alpha_for_all_cohorts(
        self,
        J: float,
        F_alpha_cohorts: List[np.ndarray],
        S_hat_alpha_cohorts: List[np.ndarray]
    ) -> List[np.ndarray]:
        """Compute σ̂_α (partial constituent stress) for all cohorts.
        
        Reuses compute_sigma_hat_alpha_for_cohort() for each (F_α, S_hat_α) pair.
        
        Args:
            J: Layer volume ratio at current time J(s) = ρ_h/ρ(s)
            F_alpha_cohorts: List of F_α matrices
            S_hat_alpha_cohorts: List of S_hat_α matrices
            
        Returns:
            List [σ̂_α(s,τ_min), ..., σ̂_α(s,s)]
        """
        if len(F_alpha_cohorts) != len(S_hat_alpha_cohorts):
            raise ValueError(
                f"Mismatch: {len(F_alpha_cohorts)} F_α but {len(S_hat_alpha_cohorts)} S_hat_α"
            )
        
        sigma_hat_alpha_cohorts = []
        for F_alpha_tau, S_hat_alpha_tau in zip(F_alpha_cohorts, S_hat_alpha_cohorts):
            sigma_hat_alpha_tau = self.compute_sigma_hat_alpha_for_cohort(
                J, F_alpha_tau, S_hat_alpha_tau
            )
            sigma_hat_alpha_cohorts.append(sigma_hat_alpha_tau)
        
        return sigma_hat_alpha_cohorts
    
    # =========================================================================
    # HEREDITY INTEGRAL (stress integration)
    # =========================================================================
    
    def integrate_constituent_stress(
        self,
        sigma_hat_alpha_cohorts: List[np.ndarray],
        mR_alpha_history: List[float],
        survival_values: List[float],
        rhoR_homeostatic: float,
        tau_min: int,
        current_timestep: int,
        dt: float,
        integration_method: str,
        rhoR_alpha_homeostatic: float
    ) -> np.ndarray:
        """Integrate stress heredity integral: σ_α = ∫ (mR·q/ρ_h) · σ̂_α dτ.
        
        This includes:
        1. Initial material stress (τ <= 0, pre-existing at homeostasis)
        2. Deposited material stress (τ > tau_min, tracked cohorts)

        This is a HELPER function called by Constituent after it has:
        1. Computed F_α for all cohorts
        2. Computed S_hat_α for all cohorts
        3. Computed σ̂_α for all cohorts
        
        Args:
            sigma_hat_alpha_cohorts: List [σ̂_α(τ_min), ..., σ̂_α(s)]
            mR_alpha_history: Production rate history [mR(τ_min), ..., mR(s)]
            survival_values: Survival function [q(s,τ_min), ..., q(s,s)]
            rhoR_homeostatic: Homeostatic mass density (kg/m³)
            tau_min: Earliest deposition time
            current_timestep: Current time s
            dt: Time step size
            integration_method: Integration method ('simpson'/'trapezoidal')
            rhoR_alpha_homeostatic: Homeostatic mass density 
            
        Returns:
            Total constituent stress σ_α (3×3)
        """
        # TODO: rename variables, remove len dependency?       

        # STEP 1: Integrate deposited material stress (τ >= tau_min)
        sigma_deposited = np.zeros((3, 3))
        # Integrate each stress component [i,j]
        for i in range(3):
            for j in range(3):
                # Extract σ̂_ij(τ) for all cohorts and weight by (mR·q/ρ_h)
                weighted_stress_ij = []
                
                # TODO: Currently tau_min not working with integration schemes
                # See relative vs absolute time indexing discussion
                for tau in range(tau_min, current_timestep+1):
                    
                    # Weight: (mR(τ) · q(s,τ)) / ρ_h
                    #TODO: to add J(tau) here
                    weight = (mR_alpha_history[tau] * survival_values[tau]) / rhoR_homeostatic
                    
                    # Weighted stress component
                    weighted_stress_ij.append(
                        weight * sigma_hat_alpha_cohorts[tau][i, j]
                    )
                
                # Integrate this component
                integrator = IntegratorFactory.create(
                    integration_method, 
                    dt, 
                    current_timestep,  # start
                    tau_min  # stop
                )
                sigma_deposited[i, j] = integrator.integrate(weighted_stress_ij)

        # STEP 2: Add initial material stress contribution (τ <= 0)
        if tau_min == 0:            
            q_initial = survival_values[0]  
            weight_initial = (rhoR_alpha_homeostatic / rhoR_homeostatic) * q_initial
            #TODO: potential source of errror?
            sigma_initial = weight_initial * sigma_hat_alpha_cohorts[0]

            # TODO: convert to DEBUG: Print initial material contribution
            # print(f"        [InitialMaterial] q(s,0) = {q_initial:.6f}")
            # print(f"        [InitialMaterial] rhoR_alpha_h = {rhoR_alpha_homeostatic:.2f} kg/m³")
            # print(f"        [InitialMaterial] weight = {weight_initial:.6f}")
            # print(f"        [InitialMaterial] sigma_initial[1,1] = {sigma_initial[1,1]/1000:.2f} kPa")
            # print(f"        [InitialMaterial] sigma_deposited[1,1] = {sigma_deposited[1,1]/1000:.2f} kPa")
        else:
            sigma_initial = np.zeros((3, 3))

        sigma_alpha = sigma_initial + sigma_deposited

        return sigma_alpha
    
    # =========================================================================
    # UTILITY METHODS
    # =========================================================================
    
    def compute_stress_trace(self, context: MechanicsContext, timestep: int) -> float:
        """Compute trace of stress tensor (for intramural stress).
        
        Delegates to TensorOperations.
        
        Args:
            context: Mechanics context
            timestep: Timestep index
            
        Returns:
            tr(σ) = σ_rr + σ_θθ + σ_zz
        """
        stress_tensor = context.get_stress_tensor(timestep)
        return self.tensor_ops.trace(stress_tensor)