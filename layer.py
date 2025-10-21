from typing import List, Dict, Any, Optional
from constituent import Constituent
import constituent

class Layer:
    def __init__(self, name, layer_type, params):
        self.name = name
        self.layer_type = layer_type
        self.params = params
        self.constituents: List[Constituent] = []

        # Initialize geometry from parameters
        self._initialize_geometry_from_parameters()
        
        # Mechanical state (computed from constituents)
        self.cauchy_stress = 0.0
        self.total_mass_density = 0.0
        self.deformation_gradient = 1.0
    
    @classmethod
    def from_parameters(cls, layer_data, params):
        """Create and initialize Layer from parameter dictionary."""
        print(f"\nInitializing layer: {layer_data['layer_name']}")
        
        # Create layer with homeostatic geometry
        layer = cls._initialize_layer_at_homeostasis(layer_data)
        
        # Initialize constituents for this layer
        constituents = layer_data['constituents']
        print(f"  Constituents ({len(constituents)} total):")
        
        for name, properties in constituents.items():
            # Create constituent (single or multi-fiber family)
            constituent = Constituent.from_parameters(name, properties)
            layer.add_constituent(constituent)
        
        return layer
    
    @classmethod
    def _initialize_layer_at_homeostasis(cls, layer_data):
        """Initialize layer homeostatic geometry (a_h, h_h, lambda_h, rho_h)."""
        print("  Homeostatic geometry:")
        
        # Unit conversion
        mm_to_m = 1.0e-3
        
        # Homeostatic geometry parameters with unit conversion
        homeostatic_params = {
            'inner_radius': layer_data['a_h'] * mm_to_m,      # a_h: mm → m
            'thickness': layer_data['h_h'] * mm_to_m,         # h_h: mm → m  
            'axial_stretch': layer_data['lambda_z_h'],        # lambda_z_h (dimensionless)
            'reference_density': layer_data['rhoR_h']         # rho_h in kg/m³
        }
        
        print(f"    Inner radius (a_h): {layer_data['a_h']} mm → {homeostatic_params['inner_radius']:.6f} m")
        print(f"    Thickness (h_h): {layer_data['h_h']} mm → {homeostatic_params['thickness']:.6f} m")
        print(f"    Axial stretch (λ_z_h): {homeostatic_params['axial_stretch']}")
        print(f"    Reference density (ρ_h): {homeostatic_params['reference_density']} kg/m³")
        
        # Create layer with homeostatic geometry
        layer = cls(layer_data['layer_name'], "biological", homeostatic_params)
        
        return layer
    
    def _initialize_geometry_from_parameters(self):
        """Initialize geometry parameters from params."""
        self.inner_radius = self.params.get('inner_radius', 0.0)
        self.thickness = self.params.get('thickness', 0.0)
        self.axial_stretch = self.params.get('axial_stretch', 1.0)
        self.reference_density = self.params.get('reference_density', 1050.0)
        self.mid_radius = self.inner_radius + 0.5 * self.thickness

    def update_constituent_kinetics(self, dt):
        """Update kinetics for all constituents in this layer."""
        print(f"Updating constituent kinetics for layer {self.name}")
        for constituent in self.constituents:
            constituent.update_kinetics(dt)
    
    def update_mechanical_state(self, dt):
        """Update layer mechanical state from constituent properties."""
        print(f"Computing mechanical state for layer {self.name}")
        
        # 1. Compute total stress from all constituents
        self.cauchy_stress = self._compute_total_stress()
        
        # 2. Update geometry based on mechanical equilibrium
        self._update_geometry_from_mechanics()
        
        # 3. Update mass density from constituent masses
        self._update_mass_density()
    
    def _compute_total_stress(self):
        """Compute total Cauchy stress from all constituents."""
        total_stress = 0.0
        for constituent in self.constituents:
            constituent_stress = constituent.compute_stress()
            total_stress += constituent_stress
        return total_stress
    
    def _update_geometry_from_mechanics(self):
        """Update layer geometry based on mechanical equilibrium."""
        # Example: Update thickness based on stress state
        # This would involve solving equilibrium equations
        # For now, mock implementation
        stress_factor = 1.0 + 0.001 * self.cauchy_stress  # Mock deformation
        self.thickness *= stress_factor
        self.mid_radius = self.inner_radius + 0.5 * self.thickness
    
    def _update_mass_density(self):
        """Update total mass density from constituent contributions."""
        self.total_mass_density = sum(
            constituent.get_mass_density() for constituent in self.constituents
        )
    
    def enforce_contact_with(self, inner_layer):
        """Enforce geometric contact constraint with inner layer."""
        # Ensure this layer's inner radius matches inner layer's outer radius
        inner_layer_outer_radius = inner_layer.inner_radius + inner_layer.thickness
        self.inner_radius = inner_layer_outer_radius
        self.mid_radius = self.inner_radius + 0.5 * self.thickness
    
    def get_total_intramural_stress(self):
        """Sum intramural stress from all constituents."""
        return sum(constituent.get_stress() for constituent in self.constituents)
    
    def get_cauchy_stress(self, fluid_pressure=0.0):
        """Convert to Cauchy stress by adding fluid pressure."""
        intramural_stress = self.get_total_intramural_stress()
        # For isotropic fluid pressure: σ_cauchy = σ_intramural - p*I
        return intramural_stress - fluid_pressure
    
    def get_total_mass(self):
        """Get total mass of this layer."""
        volume = 3.14159 * (
            (self.inner_radius + self.thickness)**2 - self.inner_radius**2
        ) * 1.0  # Assuming unit length
        return self.total_mass_density * volume

    def add_constituent(self, constituent):
        """Add constituent to layer and set back-reference."""
        constituent.layer = self  # Set back-reference to parent layer
        self.constituents.append(constituent)

    def get_geometry(self):
        """Return current layer geometry."""
        return {
            'inner_radius': self.inner_radius,
            'thickness': self.thickness,
            'mid_radius': self.mid_radius,
            'axial_stretch': self.axial_stretch,
            'reference_density': self.reference_density
        }
    
    def get_all_fiber_families(self):
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
    
    def get_constituent_summary(self):
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
    
    def guess_all_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
        """Guess mass densities for all constituents at target timestep using specified method."""
        print(f"  Guessing rhoR_alpha for layer '{self.name}' at timestep {target_timestep} using '{guess_method}' method")
        
        for constituent in self.constituents:
            constituent.guess_rhoR_alpha(target_timestep, guess_method)
    
    def compute_all_rhoR_alpha(self, target_timestep):
        """Compute mass densities for all constituents at target timestep."""
        print(f"  Computing rhoR_alpha for layer '{self.name}' at timestep {target_timestep}")
        
        for constituent in self.constituents:
            constituent.compute_rhoR_alpha(target_timestep)