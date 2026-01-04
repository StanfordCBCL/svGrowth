"""
Interface layer between data sources (Constituent, Layer, Configuration) and Kinetics.

This module defines:
1. Abstract interface (contract) for what Kinetics needs
2. Concrete adapters for different data sources
3. Factory for creating appropriate adapters
"""

from abc import ABC, abstractmethod
from typing import Optional


# =============================================================================
# EXCEPTIONS
# =============================================================================

class KineticsDataNotAvailableError(Exception):
    """Raised when kinetics requests data that is not available."""
    pass

# =============================================================================
# ABSTRACT INTERFACE - The contract that Kinetics depends on
# =============================================================================

class KineticsContext(ABC):
    """Abstract interface for accessing constituent/layer data in kinetics computations.
    
    This adapter pattern decouples kinetics computations from the specific
    data structures of Constituent and Layer classes.
    
    Implementations should:
    1. Map stimulus names (e.g., 'intramural_stress') to actual data sources
    2. Handle missing data gracefully (raise KineticsDataNotAvailableError)
    3. Provide both current and homeostatic values for stimuli
    """
    
    # =================================================================
    # ABSTRACT METHODS - Must be implemented by subclasses
    # =================================================================
    
    @abstractmethod
    def get_stimulus(self, stimulus_name: str, timestep: int) -> float:
        """Get current value of a stimulus at specified timestep.
        
        Args:
            stimulus_name: Name of stimulus (e.g., 'intramural_stress', 'wss')
            timestep: Timestep at which to retrieve value
            
        Returns:
            Current stimulus value
            
        Raises:
            ValueError: If stimulus_name is not recognized
            KineticsDataNotAvailableError: If data is not available
            
        Example:
            >>> context.get_stimulus('intramural_stress', 10)
            150000.0  # Pa
        """
        pass
    
    @abstractmethod
    def get_stimulus_homeostatic(self, stimulus_name: str) -> float:
        """Get homeostatic (reference) value of a stimulus.
        
        Args:
            stimulus_name: Name of stimulus (e.g., 'intramural_stress', 'wss')
            
        Returns:
            Homeostatic stimulus value
            
        Raises:
            ValueError: If stimulus_name is not recognized
            KineticsDataNotAvailableError: If data is not available
            
        Example:
            >>> context.get_stimulus_homeostatic('intramural_stress')
            100000.0  # Pa
        """
        pass
    
    # =================================================================
    # CONSTITUENT-SPECIFIC DATA ACCESS
    # =================================================================
    
    def get_rhoR_alpha(self, timestep: int) -> float:
        """Get constituent mass density at specified timestep.
        
        Args:
            timestep: Timestep at which to retrieve density
            
        Returns:
            Referential mass density (kg/m³)
            
        Raises:
            KineticsDataNotAvailableError: If data is not available
        """
        if not hasattr(self, 'constituent'):
            raise KineticsDataNotAvailableError(
                "Context does not have constituent reference"
            )
        
        return self.constituent.get_rhoR_alpha(timestep)
    
    def get_survival(self, deposition_timestep: int, current_timestep: int) -> float:
        """Get survival fraction for cohort deposited at deposition_timestep.
        
        Args:
            deposition_timestep: When cohort was deposited
            current_timestep: Current timestep
            
        Returns:
            Survival fraction q(τ, s) ∈ [0, 1]
            
        Raises:
            KineticsDataNotAvailableError: If data is not available
        """
        if not hasattr(self, 'constituent'):
            raise KineticsDataNotAvailableError(
                "Context does not have constituent reference"
            )
        
        cohort_age = current_timestep - deposition_timestep
        
        if deposition_timestep < 0 or deposition_timestep >= len(self.constituent.survival_history):
            raise KineticsDataNotAvailableError(
                f"No survival data for cohort deposited at t={deposition_timestep}"
            )
        
        if cohort_age < 0 or cohort_age >= len(self.constituent.survival_history[deposition_timestep]):
            raise KineticsDataNotAvailableError(
                f"No survival data for cohort age {cohort_age}"
            )
        
        return self.constituent.survival_history[deposition_timestep][cohort_age]


# =============================================================================
# CONSTITUENT ADAPTER - Provides constituent + layer data
# =============================================================================

class ConstituentKineticsContext(KineticsContext):
    """Concrete implementation: Maps stimulus names to Constituent/Layer data sources."""
    
    def __init__(self, constituent):
        """Initialize context with constituent reference.
        
        Args:
            constituent: Constituent instance (provides access to Layer via constituent.layer)
        """
        self.constituent = constituent
        self.layer = constituent.layer
    
    # =================================================================
    # IMPLEMENT ABSTRACT METHODS
    # =================================================================
    
    def get_stimulus(self, stimulus_name: str, timestep: int) -> float:
        """Map stimulus name to appropriate data source.
        
        Supported stimuli:
        - 'intramural_stress': Layer stress history
        - 'constituent_stress': Constituent stress history
        - 'wss': Wall shear stress from layer
        - 'inflammation': Inflammation (layer or constituent)
        """
        if stimulus_name == 'intramural_stress':
            return self._get_layer_stress(timestep)
        
        elif stimulus_name == 'constituent_stress':
            return self._get_constituent_stress(timestep)
        
        elif stimulus_name == 'wss':
            return self._get_wall_shear_stress(timestep)
        
        elif stimulus_name == 'inflammation':
            return self._get_inflammation(timestep)
        
        else:
            raise ValueError(
                f"Unknown stimulus: '{stimulus_name}'. "
                f"Supported stimuli: intramural_stress, constituent_stress, wss, inflammation"
            )
    
    def get_stimulus_homeostatic(self, stimulus_name: str) -> float:
        """Get homeostatic value for stimulus."""
        if stimulus_name == 'intramural_stress':
            return self._get_layer_stress_homeostatic()
        
        elif stimulus_name == 'constituent_stress':
            # Use layer homeostatic as proxy (constituent doesn't track its own)
            return self._get_layer_stress_homeostatic()
        
        elif stimulus_name == 'wss':
            return self._get_wall_shear_stress_homeostatic()
        
        elif stimulus_name == 'inflammation':
            return 0.0  # Homeostatic inflammation is zero
        
        else:
            raise ValueError(
                f"Unknown stimulus: '{stimulus_name}'. "
                f"Supported stimuli: intramural_stress, constituent_stress, wss, inflammation"
            )
    
    # =================================================================
    # PRIVATE HELPER METHODS - Actual data access
    # =================================================================
    
    def _get_layer_stress(self, timestep: int) -> float:
        """Get intramural stress from layer (with fallback to previous timestep)."""
        if not self.layer:
            raise KineticsDataNotAvailableError("Layer not available")
        
        try:
            # Try to get stress at requested timestep
            return self.layer.get_stress_trace(timestep)
        except (IndexError, ValueError, AttributeError):
            # Fallback to previous timestep
            if timestep > 0:
                try:
                    print(f"  [INFO]: Layer stress not available at t={timestep}, using t={timestep-1}")
                    return self.layer.get_stress_trace(timestep - 1)
                except (IndexError, ValueError, AttributeError):
                    pass
        
        raise KineticsDataNotAvailableError(
            f"Could not get layer stress at timestep {timestep} or {timestep-1}"
        )

    def _get_layer_stress_homeostatic(self) -> float:
        """Get homeostatic layer stress."""
        if not self.layer:
            raise KineticsDataNotAvailableError("Layer not available")
        try:
            return self.layer.get_stress_trace(0)  # Assume t=0 is homeostatic
            # TODO: Clean this up with separate function.
        except (IndexError, ValueError, AttributeError) as e:
            raise KineticsDataNotAvailableError(f"Could not get layer stress: {e}")        

    def _get_constituent_stress(self, timestep: int) -> float:
        """Get stress from this specific constituent (with fallback)."""
        if not hasattr(self.constituent, 'stress_history'):
            raise KineticsDataNotAvailableError("Constituent stress not available")
        
        try:
            # Try to get stress at requested timestep
            return self.constituent.stress_history[timestep]
        except (IndexError, ValueError, AttributeError):
            # Fallback to previous timestep
            if timestep > 0:
                try:
                    print(f"  Warning: Constituent stress not available at t={timestep}, using t={timestep-1}")
                    return self.constituent.stress_history[timestep - 1]
                except (IndexError, ValueError, AttributeError):
                    pass
        
        raise KineticsDataNotAvailableError(
            f"Could not get constituent stress at timestep {timestep} or {timestep-1}"
        )

    def _get_wall_shear_stress(self, timestep: int) -> float:
        """Get wall shear stress from layer (with fallback)."""
        if not self.layer or not hasattr(self.layer, 'wss_history'):
            raise KineticsDataNotAvailableError("Wall shear stress not available")
        
        try:
            # Try to get WSS at requested timestep
            return self.layer.wss_history[timestep]
        except (IndexError, ValueError, AttributeError):
            # Fallback to previous timestep
            if timestep > 0:
                try:
                    print(f"  [INFO]: WSS not available at t={timestep}, using t={timestep-1}")
                    return self.layer.wss_history[timestep - 1]
                except (IndexError, ValueError, AttributeError):
                    pass
        
        raise KineticsDataNotAvailableError(
            f"Could not get WSS at timestep {timestep} or {timestep-1}"
        )

    def _get_wall_shear_stress_homeostatic(self) -> float:
        """Get homeostatic wall shear stress."""
        if not self.layer or not hasattr(self.layer, 'homeostatic_wss'):
            raise KineticsDataNotAvailableError(
                f"Homeostatic wall shear stress not available"
            )
        return self.layer.homeostatic_wss
    
    def _get_inflammation(self, timestep: int) -> float:
        """Get inflammation (with fallback to previous timestep)."""
        # Try layer first (systemic inflammation)
        if self.layer and hasattr(self.layer, 'inflammation_history'):
            try:
                return self.layer.inflammation_history[timestep]
            except (IndexError, ValueError, AttributeError):
                if timestep > 0:
                    try:
                        print(f"  [INFO]: Layer inflammation not available at t={timestep}, using t={timestep-1}")
                        return self.layer.inflammation_history[timestep - 1]
                    except (IndexError, ValueError, AttributeError):
                        pass
    
        # Try constituent (local inflammation)
        if hasattr(self.constituent, 'inflammation_history'):
            try:
                return self.constituent.inflammation_history[timestep]
            except (IndexError, ValueError, AttributeError):
                if timestep > 0:
                    try:
                        print(f"  [INFO]: Constituent inflammation not available at t={timestep}, using t={timestep-1}")
                        return self.constituent.inflammation_history[timestep - 1]
                    except (IndexError, ValueError, AttributeError):
                        pass
                    
        raise KineticsDataNotAvailableError(
            f"Inflammation not available at timestep {timestep} or {timestep-1}"
        )