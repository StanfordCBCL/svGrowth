# StressState Design Pattern - Reference Document

## Overview

This document describes the **StressState abstraction pattern** for coordinate-system-agnostic stress handling in pyGrowth.

## Core Principle

**Decouple stress storage from coordinate system knowledge using the Adapter Pattern**

```
┌─────────────────────────────────────────────────────────┐
│ StressState (Wrapper)                                   │
│  - Stores: np.ndarray([σ_1, σ_2, σ_3])                 │
│  - References: Kinematics (knows coordinate system)     │
│  - Provides: get_component(name) → asks kinematics      │
└─────────────────────────────────────────────────────────┘
                          ↓
              ┌───────────┴────────────┐
              ↓                        ↓
    ┌─────────────────┐      ┌─────────────────┐
    │ Layer           │      │ Constituent     │
    │  stress_history │      │  sigma_hat_hist │
    │  List[Stress    │      │  sigma_alpha    │
    │       State]    │      │  List[Stress    │
    └─────────────────┘      │       State]    │
                             └─────────────────┘
```

## Key Design Decisions

### 1. **Stress Storage: Vector-Only (For Now)**
- **Internal storage:** `np.ndarray(shape=(3,))` - principal stresses [σ_1, σ_2, σ_3]
- **Why vectors:** Simpler implementation, sufficient for thin-wall
- **Future:** Easy to extend to full 3×3 tensors without API changes

### 2. **Coordinate System: Owned by Kinematics**
- **NOT global:** Each Layer has kinematics instance
- **Set from config:** Via `geometry.type` in YAML
- **Access pattern:** StressState asks kinematics for component mappings

### 3. **No Hardcoded Indices**
```python
# ❌ BAD - Tight coupling
sigma_theta = stress_vector[1]  # What is index 1?

# ✅ GOOD - Coordinate-agnostic
sigma_theta = stress_state.get_component('theta')  # Asks kinematics
```

### 4. **Kinematics as Adapter**
Kinematics provides component name → index mapping:
```python
class DeformationKinematics(ABC):
    @abstractmethod
    def get_component_names(self) -> List[str]:
        """['r', 'theta', 'z'] for cylindrical"""
        pass
    
    @abstractmethod
    def get_component_index(self, component: str) -> int:
        """'theta' → 1 for cylindrical"""
        pass
    
    @abstractmethod
    def get_intramural_component_name(self) -> str:
        """'theta' for cylinder"""
        pass
```

## Architecture Analogy

**This follows the SAME pattern as KineticsContext:**

```
┌────────────────────────────────────────────────────────┐
│ KineticsContext Pattern (EXISTING)                    │
├────────────────────────────────────────────────────────┤
│ Kinetics (algorithm) ← KineticsContext (adapter) ←    │
│                        Constituent/Layer (data owner)  │
└────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────┐
│ StressState Pattern (PROPOSED)                        │
├────────────────────────────────────────────────────────┤
│ User code ← StressState (adapter) ← Kinematics        │
│                                      (coord system)    │
└────────────────────────────────────────────────────────┘
```

Both use **adapter pattern** to decouple algorithm/access from implementation details!

## Implementation Checklist

### Phase 1: Core StressState Class
- [ ] Create `stress_state.py` with vector storage
- [ ] Implement `get_component()` - asks kinematics for index
- [ ] Implement `get_intramural_stress()` - asks kinematics for component name
- [ ] Implement arithmetic ops: `__add__`, `__mul__`, `copy()`
- [ ] Add `get_all_components()` → Dict[str, float]

### Phase 2: Kinematics Interface
- [ ] Add `get_component_names()` to DeformationKinematics ABC
- [ ] Add `get_component_index(component: str)` to ABC
- [ ] Add `get_intramural_component_name()` to ABC
- [ ] Implement in ThinWallKinematics
- [ ] Update `compute_stress_from_equilibrium()` to return StressState

### Phase 3: Layer Integration
- [ ] Change `homeostatic_stress: Optional[float]` → `Optional[StressState]`
- [ ] Change `stress_history: List[float]` → `List[StressState]`
- [ ] Update `_initialize_homeostatic_geometry()` to use StressState
- [ ] Update `get_intramural_stress()` to extract from StressState
- [ ] Add `get_stress_component()`, `get_stress_vector()` accessors

### Phase 4: Constituent Integration
- [ ] Change `sigma_hat_history: List[??]` → `List[StressState]`
- [ ] Change `sigma_alpha_history: List[??]` → `List[StressState]`
- [ ] Update initialization to use StressState
- [ ] Update stress computations (Mechanics class - future)

### Phase 5: Interface Updates
- [ ] Update `kinetics_interface.py` accessors
  - `_get_layer_stress()` → extract from StressState
  - `_get_layer_stress_homeostatic()` → extract from StressState
- [ ] No API changes needed (still returns float)

## Benefits Checklist

After implementation, verify:
- [ ] **No hardcoded indices:** Search codebase for `stress[0]`, `stress[1]`, etc. → should be zero!
- [ ] **Coordinate agnostic:** Can switch cylindrical ↔ spherical without changing Layer/Constituent
- [ ] **Self-documenting:** `get_component('theta')` is clearer than `stress[1]`
- [ ] **Type safe:** `stress_history: List[StressState]` not `List[Union[float, np.ndarray, dict]]`
- [ ] **Easy tensor migration:** Only need to change StressState internals, not user code

## Code Examples to Remember

### Creating StressState
```python
# Kinematics computes and wraps
stress_vector = np.array([sigma_r, sigma_theta, sigma_z])
stress_state = StressState(stress_vector, kinematics=self)
```

### Accessing Components
```python
# By name (coordinate-agnostic)
sigma_theta = stress.get_component('theta')

# Primary intramural stress
sigma_intramural = stress.get_intramural_stress()

# All components
components = stress.get_all_components()
# → {'r': 0.0, 'theta': 130000.0, 'z': 65000.0}
```

### Arithmetic
```python
# Addition (e.g., constituent stresses)
sigma_total = sigma_elastin + sigma_collagen

# Scaling (e.g., survival weighting)
sigma_weighted = q * sigma_hat

# Accumulation
sigma_total = StressState.zeros(kinematics)
for constituent in constituents:
    sigma_total = sigma_total + constituent.stress_history[t]
```

### Storage in History
```python
# Initialize
self.stress_history = [homeostatic_stress.copy()]

# Append
new_stress = kinematics.compute_stress(...)
self.stress_history.append(new_stress)
```

## Migration Path to Tensors (Future)

When thick-wall or FEM is needed:

1. **Update StressState internal storage:**
```python
# Change one line:
self._data = data  # Can now be shape (3,) OR (3,3)
```

2. **Update component access:**
```python
def get_component(self, name: str) -> float:
    idx = self._kinematics.get_component_index(name)
    if self._data.ndim == 1:  # Vector
        return self._data[idx]
    else:  # Tensor
        return self._data[idx, idx]  # Diagonal
```

3. **That's it!** All user code unchanged.

## Testing Strategy

### Unit Tests for StressState
```python
def test_component_access_cylindrical():
    kin = ThinWallKinematics()
    stress = StressState(np.array([0, 130e3, 65e3]), kin)
    assert stress.get_component('theta') == 130e3
    assert stress.get_component('r') == 0
    assert stress.get_component('z') == 65e3

def test_intramural_stress():
    kin = ThinWallKinematics()
    stress = StressState(np.array([0, 130e3, 65e3]), kin)
    # Cylinder: intramural is 'theta'
    assert stress.get_intramural_stress() == 130e3

def test_stress_addition():
    kin = ThinWallKinematics()
    s1 = StressState(np.array([1, 2, 3]), kin)
    s2 = StressState(np.array([4, 5, 6]), kin)
    s_sum = s1 + s2
    assert np.allclose(s_sum.as_vector(), [5, 7, 9])
```

### Integration Test
```python
def test_layer_stress_initialization():
    layer = Layer.from_parameters(yaml_data)
    
    # Should be StressState
    assert isinstance(layer.homeostatic_stress, StressState)
    
    # Should have correct value
    sigma_theta = layer.homeostatic_stress.get_component('theta')
    expected = P * a / h  # Thin-wall
    assert abs(sigma_theta - expected) / expected < 0.01
    
    # History should contain StressState
    assert isinstance(layer.stress_history[0], StressState)
```

## Questions to Ask When Implementing

1. **Am I hardcoding an index?** → Use `get_component()` instead
2. **Am I checking coordinate system?** → Let kinematics handle it
3. **Am I storing raw np.ndarray?** → Wrap in StressState
4. **Am I extracting scalar from stress?** → Use accessor methods
5. **Am I doing stress arithmetic?** → Use `__add__`, `__mul__` operators

## Files to Modify

1. **New file:** `stress_state.py` (~150 lines)
2. **Update:** `deformation_kinematics.py` (add component methods)
3. **Update:** `layer.py` (stress storage → StressState)
4. **Update:** `constituent.py` (stress storage → StressState)
5. **Update:** `kinetics_interface.py` (extract from StressState)
6. **Future:** `mechanics.py` (new file, similar to `kinetics.py`)

## Key Insight to Remember

> **The abstraction (StressState) provides benefits independent of storage format (vector vs tensor).**
> 
> Benefits come from:
> - Decoupling from coordinate system (via kinematics adapter)
> - Type safety (uniform interface)
> - Self-documenting code (named components)
> - Easy future migration (change internals, not API)

This is the **same philosophy as KineticsContext** - adapter pattern to decouple!

---

## Related Architecture Documents
- See `ARCHITECTURE.md` for adapter pattern philosophy
- See `kinetics_interface.py` for similar adapter implementation
- Compare with `KineticsContext` design - same pattern!