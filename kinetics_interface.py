"""
Interface layer between data sources (Constituent, Layer, Configuration) and Kinetics.

This module defines:
1. Abstract interface (contract) for what Kinetics needs
2. Concrete adapters for different data sources
3. Factory for creating appropriate adapters
"""

from abc import ABC, abstractmethod
from typing import Optional, Callable


# =============================================================================
# ABSTRACT INTERFACE - The contract that Kinetics depends on
# =============================================================================

class KineticsContext(ABC):
    """Abstract interface for providing state variables to kinetics computations.
    
    This decouples Kinetics from knowing about specific data structures.
    Any class can implement this interface to provide data to Kinetics.
    """
    
    @abstractmethod
    def get_cauchy_stress(self) -> float:
        """Get current Cauchy stress (total stress)."""
        pass
    
    @abstractmethod
    def get_cauchy_stress_homeostatic(self) -> float:
        """Get homeostatic Cauchy stress."""
        pass
    
    @abstractmethod
    def get_wall_shear_stress(self) -> float:
        """Get current wall shear stress."""
        pass
    
    @abstractmethod
    def get_wall_shear_stress_homeostatic(self) -> float:
        """Get homeostatic wall shear stress."""
        pass
    
    @abstractmethod
    def get_constituent_stress(self) -> float:
        """Get constituent's stress contribution (if applicable)."""
        pass
    
    @abstractmethod
    def get_inflammation(self) -> float:
        """Get inflammation level."""
        pass
    
    @abstractmethod
    def get_rhoR_alpha(self) -> float:
        """Get referential mass density (if applicable)."""
        pass
    
    @abstractmethod
    def get_timestep(self) -> int:
        """Get current timestep."""
        pass


# =============================================================================
# CONSTITUENT ADAPTER - Provides constituent + layer data
# =============================================================================

class ConstituentKineticsContext(KineticsContext):
    """Adapter: Provides constituent and layer data to kinetics.
    
    Use this when computing kinetics for a specific constituent.
    """
    
    def __init__(self, constituent, timestep: int):
        """Initialize context for constituent at specific timestep.
        
        Args:
            constituent: Constituent instance
            timestep: Which timestep to access
        """
        self.constituent = constituent
        self.layer = constituent.layer
        self.timestep = timestep
    
    def get_cauchy_stress(self) -> float:
        """Get total Cauchy stress from layer."""
        if self.layer and hasattr(self.layer, 'stress_history'):
            if self.timestep < len(self.layer.stress_history):
                return self.layer.stress_history[self.timestep]
        return 0.0
    
    def get_cauchy_stress_homeostatic(self) -> float:
        """Get homeostatic Cauchy stress from layer."""
        if self.layer and hasattr(self.layer, 'stress_homeostatic'):
            return self.layer.stress_homeostatic
        return 0.0
    
    def get_wall_shear_stress(self) -> float:
        """Get wall shear stress from layer."""
        if self.layer and hasattr(self.layer, 'wss_history'):
            if self.timestep < len(self.layer.wss_history):
                return self.layer.wss_history[self.timestep]
        return 0.0
    
    def get_wall_shear_stress_homeostatic(self) -> float:
        """Get homeostatic WSS from layer."""
        if self.layer and hasattr(self.layer, 'wss_homeostatic'):
            return self.layer.wss_homeostatic
        return 0.0
    
    def get_constituent_stress(self) -> float:
        """Get constituent's own stress contribution."""
        if hasattr(self.constituent, 'stress_history'):
            if self.timestep < len(self.constituent.stress_history):
                return self.constituent.stress_history[self.timestep]
        return 0.0
    
    def get_inflammation(self) -> float:
        """Get inflammation - try layer first, then constituent."""
        # Try layer (systemic)
        if self.layer and hasattr(self.layer, 'inflammation_history'):
            if self.timestep < len(self.layer.inflammation_history):
                return self.layer.inflammation_history[self.timestep]
        
        # Try constituent (local)
        if hasattr(self.constituent, 'inflammation_history'):
            if self.timestep < len(self.constituent.inflammation_history):
                return self.constituent.inflammation_history[self.timestep]
        
        return 0.0
    
    def get_rhoR_alpha(self) -> float:
        """Get constituent's referential mass density."""
        if self.timestep < len(self.constituent.rhoR_alpha_history):
            return self.constituent.rhoR_alpha_history[self.timestep]
        return 0.0
    
    def get_timestep(self) -> int:
        """Get current timestep."""
        return self.timestep


# =============================================================================
# LAYER ADAPTER - Provides layer-level data
# =============================================================================

class LayerKineticsContext(KineticsContext):
    """Adapter: Provides layer-level data to kinetics.
    
    Use this when computing layer-level kinetics (e.g., ECM remodeling).
    """
    
    def __init__(self, layer, timestep: int):
        """Initialize context for layer at specific timestep.
        
        Args:
            layer: Layer instance
            timestep: Which timestep to access
        """
        self.layer = layer
        self.timestep = timestep
    
    def get_cauchy_stress(self) -> float:
        """Get total layer stress."""
        if hasattr(self.layer, 'stress_history'):
            if self.timestep < len(self.layer.stress_history):
                return self.layer.stress_history[self.timestep]
        return 0.0
    
    def get_cauchy_stress_homeostatic(self) -> float:
        """Get homeostatic stress."""
        if hasattr(self.layer, 'stress_homeostatic'):
            return self.layer.stress_homeostatic
        return 0.0
    
    def get_wall_shear_stress(self) -> float:
        """Get wall shear stress."""
        if hasattr(self.layer, 'wss_history'):
            if self.timestep < len(self.layer.wss_history):
                return self.layer.wss_history[self.timestep]
        return 0.0
    
    def get_wall_shear_stress_homeostatic(self) -> float:
        """Get homeostatic WSS."""
        if hasattr(self.layer, 'wss_homeostatic'):
            return self.layer.wss_homeostatic
        return 0.0
    
    def get_constituent_stress(self) -> float:
        """Layer doesn't have constituent-specific stress."""
        return 0.0
    
    def get_inflammation(self) -> float:
        """Get layer-level inflammation."""
        if hasattr(self.layer, 'inflammation_history'):
            if self.timestep < len(self.layer.inflammation_history):
                return self.layer.inflammation_history[self.timestep]
        return 0.0
    
    def get_rhoR_alpha(self) -> float:
        """Get total layer mass density (sum of all constituents)."""
        if hasattr(self.layer, 'total_mass_history'):
            if self.timestep < len(self.layer.total_mass_history):
                return self.layer.total_mass_history[self.timestep]
        return 0.0
    
    def get_timestep(self) -> int:
        """Get current timestep."""
        return self.timestep


# =============================================================================
# CONFIGURATION ADAPTER - Provides configuration-level data
# =============================================================================

class ConfigurationKineticsContext(KineticsContext):
    """Adapter: Provides configuration/simulation-level data to kinetics.
    
    Use this for global kinetics parameters that don't depend on 
    specific constituents or layers (e.g., system-wide inflammation).
    """
    
    def __init__(self, configuration, timestep: int):
        """Initialize context from configuration.
        
        Args:
            configuration: Configuration instance
            timestep: Which timestep to access
        """
        self.configuration = configuration
        self.timestep = timestep
    
    def get_cauchy_stress(self) -> float:
        """Get average stress across all layers."""
        if hasattr(self.configuration, 'global_stress_history'):
            if self.timestep < len(self.configuration.global_stress_history):
                return self.configuration.global_stress_history[self.timestep]
        return 0.0
    
    def get_cauchy_stress_homeostatic(self) -> float:
        """Get global homeostatic stress."""
        if hasattr(self.configuration, 'global_stress_homeostatic'):
            return self.configuration.global_stress_homeostatic
        return 0.0
    
    def get_wall_shear_stress(self) -> float:
        """Get global WSS."""
        if hasattr(self.configuration, 'global_wss_history'):
            if self.timestep < len(self.configuration.global_wss_history):
                return self.configuration.global_wss_history[self.timestep]
        return 0.0
    
    def get_wall_shear_stress_homeostatic(self) -> float:
        """Get global homeostatic WSS."""
        if hasattr(self.configuration, 'global_wss_homeostatic'):
            return self.configuration.global_wss_homeostatic
        return 0.0
    
    def get_constituent_stress(self) -> float:
        """Not applicable at configuration level."""
        return 0.0
    
    def get_inflammation(self) -> float:
        """Get systemic inflammation level."""
        if hasattr(self.configuration, 'systemic_inflammation_history'):
            if self.timestep < len(self.configuration.systemic_inflammation_history):
                return self.configuration.systemic_inflammation_history[self.timestep]
        return 0.0
    
    def get_rhoR_alpha(self) -> float:
        """Not applicable at configuration level."""
        return 0.0
    
    def get_timestep(self) -> int:
        """Get current timestep."""
        return self.timestep


# =============================================================================
# CONTEXT FACTORY - Creates appropriate context based on source
# =============================================================================

class KineticsContextFactory:
    """Factory for creating kinetics contexts from different data sources."""
    
    def __init__(self, data_source):
        """Initialize factory with data source.
        
        Args:
            data_source: Can be Constituent, Layer, or Configuration
        """
        self.data_source = data_source
        self._context_type = self._determine_context_type()
    
    def _determine_context_type(self):
        """Determine which context type to use based on data source."""
        from constituent import Constituent
        from layer import Layer
        from configuration import Configuration
        
        if isinstance(self.data_source, Constituent):
            return ConstituentKineticsContext
        elif isinstance(self.data_source, Layer):
            return LayerKineticsContext
        elif isinstance(self.data_source, Configuration):
            return ConfigurationKineticsContext
        else:
            raise ValueError(f"Unknown data source type: {type(self.data_source)}")
    
    def create_context(self, timestep: int) -> KineticsContext:
        """Create context for specified timestep.
        
        Args:
            timestep: Timestep to access
        
        Returns:
            Appropriate KineticsContext instance
        """
        return self._context_type(self.data_source, timestep)
    
    def __call__(self, timestep: int) -> KineticsContext:
        """Allow factory to be called as function: factory(timestep)"""
        return self.create_context(timestep)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_context_for_constituent(constituent, timestep: int) -> ConstituentKineticsContext:
    """Convenience function to create constituent context."""
    return ConstituentKineticsContext(constituent, timestep)


def create_context_for_layer(layer, timestep: int) -> LayerKineticsContext:
    """Convenience function to create layer context."""
    return LayerKineticsContext(layer, timestep)


def create_context_for_configuration(configuration, timestep: int) -> ConfigurationKineticsContext:
    """Convenience function to create configuration context."""
    return ConfigurationKineticsContext(configuration, timestep)