# Mechanics Interface Design Discussion: Factory vs Adapter Pattern

**Date:** November 24, 2025  
**Context:** Discussion about whether to replace the current Adapter pattern in `mechanics_interface.py` with a Factory pattern  
**Participants:** Development team discussion for software engineering review

---

## Table of Contents
1. [Background](#background)
2. [Current Implementation (Adapter Pattern)](#current-implementation-adapter-pattern)
3. [Proposed Alternative (Factory Pattern)](#proposed-alternative-factory-pattern)
4. [Comparative Analysis](#comparative-analysis)
5. [Recommendation](#recommendation)
6. [Alternative Approaches Considered](#alternative-approaches-considered)
7. [Appendix: Code Examples](#appendix-code-examples)

---

## Background

The `mechanics_interface.py` module provides an abstraction layer between:
- **Data sources**: `Constituent` and `Layer` classes (store stress, deformation, density data)
- **Computation**: `Mechanics` class (performs mechanical calculations)

**Current architecture uses the Adapter Pattern** to decouple `Mechanics` from knowing the internal data structures of `Constituent` and `Layer`.

### Question Under Discussion
> "Is it possible to use a Factory object instead to initialize mechanics to bypass the need for `mechanics_interface.py`? This factory object should contain all the mapping between the variables of the class that calls it (e.g., constituent or layer class) and the function of mechanics. Would this approach work?"

---

## Current Implementation (Adapter Pattern)

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    Current Architecture                      │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ┌─────────────┐                                           │
│  │ Constituent │                                           │
│  │   or Layer  │                                           │
│  └──────┬──────┘                                           │
│         │ creates                                          │
│         ▼                                                  │
│  ┌──────────────────────┐                                 │
│  │ MechanicsContext     │ (Adapter)                        │
│  │ - ConstituentContext │                                  │
│  │ - LayerContext       │                                  │
│  └──────┬───────────────┘                                 │
│         │ passed to                                        │
│         ▼                                                  │
│  ┌──────────────────────┐                                 │
│  │     Mechanics        │ (Stateless Computation)          │
│  │ - compute_stress()   │                                  │
│  └──────────────────────┘                                 │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### Code Flow

```python
# In Constituent class
context = ConstituentMechanicsContext(self)  # Create adapter
stress_trace = self.mechanics.compute_stress_trace(context, timestep)

# In Mechanics class
def compute_stress_trace(self, context: MechanicsContext, timestep: int) -> float:
    """Mechanics doesn't know about Constituent or Layer structure."""
    stress_tensor = context.get_stress_tensor(timestep)  # Abstract interface
    return self.tensor_ops.trace(stress_tensor)
```

### Key Characteristics

| Aspect | Implementation |
|--------|----------------|
| **Pattern** | Adapter Pattern |
| **Coupling** | Loose - Mechanics depends on abstract `MechanicsContext` |
| **Testability** | High - Easy to mock contexts |
| **Extensibility** | High - Add new contexts without changing Mechanics |
| **Complexity** | Low - Simple interface with clear responsibilities |

---

## Proposed Alternative (Factory Pattern)

### Conceptual Design

```python
class MechanicsFactory:
    """Factory that creates Mechanics instances pre-configured with data accessors."""
    
    @classmethod
    def for_constituent(cls, constituent):
        """Create Mechanics with constituent-specific data accessors."""
        mechanics = Mechanics()
        
        # Map constituent data to mechanics functions
        mechanics._get_stress = lambda t: constituent.stress_history[t]
        mechanics._get_F = lambda t: constituent.layer.F_history[t]
        mechanics._get_rhoR = lambda t: constituent.rhoR_alpha_history[t]
        
        return mechanics
    
    @classmethod
    def for_layer(cls, layer):
        """Create Mechanics with layer-specific data accessors."""
        mechanics = Mechanics()
        
        # Map layer data to mechanics functions
        mechanics._get_stress = lambda t: layer.stress_history[t]
        mechanics._get_F = lambda t: layer.F_history[t]
        mechanics._get_rhoR = lambda t: layer.rhoR_history[t]
        
        return mechanics
```

### Usage Pattern

```python
# In Constituent.__init__()
self.mechanics = MechanicsFactory.for_constituent(self)

# Later, in computation
stress_trace = self.mechanics.compute_stress_trace(timestep)  # No context needed
```

### Intended Benefits

1. **No interface file needed**: Eliminate `mechanics_interface.py`
2. **Simpler usage**: No need to create context objects
3. **Direct access**: Mechanics has direct accessors to data

---

## Comparative Analysis

### 1. Coupling Analysis

#### Current (Adapter) ✅

```python
# Mechanics is DECOUPLED from data structures
class Mechanics:
    def compute_stress_trace(self, context: MechanicsContext, timestep: int):
        # Knows NOTHING about Constituent or Layer internals
        stress = context.get_stress_tensor(timestep)  # Abstract interface
        return self.trace(stress)
```

**Dependency Graph:**
```
Mechanics → MechanicsContext (abstract)
             ↑
             └── ConstituentMechanicsContext (concrete, knows Constituent)
```

#### Factory Approach ❌

```python
# Mechanics is COUPLED to data structures via factory
class MechanicsFactory:
    def for_constituent(cls, constituent):
        # Factory MUST know about constituent.stress_history, layer.F_history, etc.
        mechanics._get_stress = lambda t: constituent.stress_history[t]
        mechanics._get_F = lambda t: constituent.layer.F_history[t]
```

**Dependency Graph:**
```
MechanicsFactory → Constituent (knows internal structure)
MechanicsFactory → Layer (knows internal structure)
Mechanics → MechanicsFactory (implicit coupling)
```

**Verdict:** Adapter maintains loose coupling; Factory introduces tight coupling.

---

### 2. Testability Analysis

#### Current (Adapter) ✅

```python
# Unit test for Mechanics (no real data structures needed)
def test_compute_stress_trace():
    # Easy to mock
    mock_context = Mock(spec=MechanicsContext)
    mock_context.get_stress_tensor.return_value = np.array([
        [100, 0, 0],
        [0, 50, 0],
        [0, 0, 25]
    ])
    
    mechanics = Mechanics()
    result = mechanics.compute_stress_trace(mock_context, 0)
    
    assert result == 175.0  # tr(σ) = 100 + 50 + 25
```

**Benefits:**
- No need to create `Constituent` or `Layer` objects
- Test focuses purely on mechanics logic
- Fast execution (no complex setup)

#### Factory Approach ❌

```python
# Unit test for Mechanics (requires complex mocking)
def test_compute_stress_trace():
    # Must mock entire data structures
    mock_constituent = Mock()
    mock_constituent.stress_history = [np.array([[100, 0, 0], [0, 50, 0], [0, 0, 25]])]
    mock_constituent.layer = Mock()
    mock_constituent.layer.F_history = [np.eye(3)]
    mock_constituent.rhoR_alpha_history = [1050.0]
    
    mechanics = MechanicsFactory.for_constituent(mock_constituent)
    result = mechanics.compute_stress_trace(0)
    
    assert result == 175.0
```

**Issues:**
- Must create complex mock hierarchies (`constituent.layer.F_history`)
- Test is coupled to data structure details
- Harder to maintain when data structures change

**Verdict:** Adapter enables cleaner, simpler unit tests.

---

### 3. Single Responsibility Principle (SRP)

#### Current (Adapter) ✅

**Separation of Concerns:**
- **`Mechanics`**: Pure computation (stateless algorithms)
- **`MechanicsContext`**: Data access (adapts data structures)
- **`Constituent/Layer`**: Data storage (manages history)

```python
# Each class has ONE responsibility

Mechanics.compute_stress_trace()
└── Responsibility: Compute trace of stress tensor

MechanicsContext.get_stress_tensor()
└── Responsibility: Retrieve stress data from storage

Constituent.stress_history
└── Responsibility: Store stress data over time
```

#### Factory Approach ❌

**Responsibilities Blur:**
- **`MechanicsFactory`**: Knows data structures + creates mechanics + maps accessors
- **`Mechanics`**: Computation + holds data accessor lambdas

```python
# Factory has MULTIPLE responsibilities

MechanicsFactory.for_constituent()
├── Responsibility: Know Constituent data structure
├── Responsibility: Create Mechanics instance
└── Responsibility: Map data accessors
```

**Verdict:** Adapter maintains clear SRP; Factory violates it.

---

### 4. Extensibility Analysis

#### Scenario: Add New Data Source (e.g., `ExternalSupport`)

**Current (Adapter) ✅**

```python
# Add new adapter - NO changes to Mechanics!
class ExternalSupportMechanicsContext(MechanicsContext):
    def __init__(self, external_support):
        self.support = external_support
    
    def get_stress_tensor(self, timestep):
        return self.support.reaction_stress[timestep]
    
    def get_deformation_gradient(self, timestep):
        return self.support.prescribed_F[timestep]
    
    def get_mass_density(self, timestep):
        return 0.0  # External support has no mass

# Usage (Mechanics unchanged!)
context = ExternalSupportMechanicsContext(my_support)
stress = mechanics.compute_stress_trace(context, timestep)
```

**Factory Approach ❌**

```python
# Must modify factory AND potentially Mechanics
class MechanicsFactory:
    @classmethod
    def for_external_support(cls, support):
        mechanics = Mechanics()
        # New mapping logic
        mechanics._get_stress = lambda t: support.reaction_stress[t]
        mechanics._get_F = lambda t: support.prescribed_F[t]
        mechanics._get_rhoR = lambda t: 0.0
        return mechanics

# If Mechanics interface changes, all factory methods need updates
```

**Verdict:** Adapter follows Open/Closed Principle; Factory requires modification.

---

### 5. Mapping Elimination Analysis

#### Does Factory Eliminate Mapping? ❌

**Current (Adapter):**
```python
# Mapping in ConstituentMechanicsContext
def get_stress_tensor(self, timestep):
    return self.constituent.stress_history[timestep]  # Mapping here

def get_deformation_gradient(self, timestep):
    return self.constituent.layer.F_history[timestep]  # Mapping here
```

**Factory:**
```python
# Mapping in MechanicsFactory
mechanics._get_stress = lambda t: constituent.stress_history[t]  # Mapping here
mechanics._get_F = lambda t: constituent.layer.F_history[t]       # Mapping here
```

**Conclusion:** Factory doesn't eliminate mapping - it just moves it from the adapter to the factory. The mapping must exist somewhere.

---

### 6. Consistency with Existing Architecture

The codebase **already uses Adapter Pattern** successfully in `kinetics_interface.py`:

```python
# Kinetics interface (parallel to mechanics interface)
class KineticsContext(ABC):
    @abstractmethod
    def get_stimulus(self, stimulus_name, timestep):
        pass

class ConstituentKineticsContext(KineticsContext):
    def __init__(self, constituent):
        self.constituent = constituent
    
    def get_stimulus(self, stimulus_name, timestep):
        if stimulus_name == 'intramural_stress':
            return self.layer.get_stress_trace(timestep)
        # ...

# Usage in Constituent
context = ConstituentKineticsContext(self)
k_alpha = self.kinetics.compute_k_alpha(context, timestep)
```

**Maintaining consistency:**
- ✅ Adapter: Keeps parallel structure with kinetics interface
- ❌ Factory: Creates architectural inconsistency

---

## Recommendation

### Keep the Current Adapter Pattern ✅

**Reasons:**

1. **Loose Coupling**: Mechanics remains decoupled from data structures
2. **High Testability**: Easy to unit test with mocked contexts
3. **Clear Responsibilities**: Each class has one well-defined job
4. **Extensible**: Add new contexts without modifying Mechanics
5. **Consistent**: Matches existing kinetics interface architecture
6. **Proven Design**: Adapter pattern is ideal for this exact use case

### When Would Factory Be Appropriate?

Factory pattern makes sense when:

| Scenario | Fits svGrowth? |
|----------|----------------|
| Complex object creation with multiple steps | ❌ No - contexts are simple (`ConstituentMechanicsContext(self)`) |
| Need to hide implementation details | ❌ No - contexts are intentionally explicit |
| Multiple configuration variants | ❌ No - only 2 contexts (constituent, layer) |
| Creating families of related objects | ❌ No - contexts are independent |

**Verdict:** Factory doesn't provide value for this use case.

---

## Alternative Approaches Considered

### Option 1: Factory + Context (Hybrid)

```python
class MechanicsContextFactory:
    """Factory to create contexts (but keep adapter pattern)."""
    
    @classmethod
    def for_constituent(cls, constituent):
        return ConstituentMechanicsContext(constituent)
    
    @classmethod
    def for_layer(cls, layer):
        return LayerMechanicsContext(layer)
```

**Analysis:**
- ❌ Doesn't save code: `MechanicsContextFactory.for_constituent(self)` vs `ConstituentMechanicsContext(self)`
- ❌ Adds unnecessary indirection
- ❌ Makes code less explicit

**Verdict:** Not beneficial.

---

### Option 2: Bind Context at Initialization

```python
class Constituent:
    def __init__(self, name):
        self.mechanics = Mechanics()
        self._mechanics_context = None  # Lazy initialization
    
    def _ensure_mechanics_context(self):
        """Create mechanics context if not already created."""
        if self._mechanics_context is None:
            if self.layer is None:
                raise RuntimeError(
                    "Cannot create mechanics context: constituent not added to layer"
                )
            self._mechanics_context = ConstituentMechanicsContext(self)
    
    def compute_stress_trace(self, timestep: int) -> float:
        """Compute stress trace using mechanics."""
        self._ensure_mechanics_context()
        return self.mechanics.compute_stress_trace(self._mechanics_context, timestep)
```

**Analysis:**
- ✅ Hides context creation from user
- ❌ Adds complexity (lazy initialization, error checking)
- ❌ Makes data flow less explicit
- ❌ Harder to understand control flow

**Verdict:** Trades simplicity for marginal convenience.

---

### Option 3: Helper Methods (Minimal Change)

```python
class Constituent:
    def compute_stress_trace_via_mechanics(self, timestep: int) -> float:
        """Helper: compute stress trace using mechanics adapter."""
        context = ConstituentMechanicsContext(self)
        return self.mechanics.compute_stress_trace(context, timestep)
```

**Analysis:**
- ✅ Reduces boilerplate for common operations
- ✅ Keeps adapter pattern intact
- ✅ Explicit and clear

**Verdict:** Best option if boilerplate is a concern.

---

## Appendix: Code Examples

### Current Implementation (Complete)

```python
# mechanics_interface.py
class MechanicsContext(ABC):
    """Abstract interface for mechanical data access."""
    
    @abstractmethod
    def get_stress_tensor(self, timestep: int) -> np.ndarray:
        pass
    
    @abstractmethod
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        pass
    
    @abstractmethod
    def get_mass_density(self, timestep: int) -> float:
        pass


class ConstituentMechanicsContext(MechanicsContext):
    """Adapter for constituent-level stress access."""
    
    def __init__(self, constituent):
        self.constituent = constituent
        self.layer = constituent.layer
    
    def get_stress_tensor(self, timestep: int) -> np.ndarray:
        return self.constituent.stress_history[timestep]
    
    def get_deformation_gradient(self, timestep: int) -> np.ndarray:
        return self.layer.F_history[timestep]
    
    def get_mass_density(self, timestep: int) -> float:
        return self.constituent.rhoR_alpha_history[timestep]


# constituent.py
class Constituent:
    def update_stress(self, timestep: int):
        """Update constituent stress using mechanics."""
        context = ConstituentMechanicsContext(self)
        stress_trace = self.mechanics.compute_stress_trace(context, timestep)
        # Use stress_trace...


# mechanics.py
class Mechanics:
    def compute_stress_trace(self, context: MechanicsContext, timestep: int) -> float:
        """Compute tr(σ) - works with ANY context."""
        stress_tensor = context.get_stress_tensor(timestep)
        return np.trace(stress_tensor)
```

---

### Factory Alternative (For Comparison)

```python
# mechanics_factory.py (would replace mechanics_interface.py)
class MechanicsFactory:
    @classmethod
    def for_constituent(cls, constituent):
        mechanics = Mechanics()
        mechanics._get_stress = lambda t: constituent.stress_history[t]
        mechanics._get_F = lambda t: constituent.layer.F_history[t]
        mechanics._get_rhoR = lambda t: constituent.rhoR_alpha_history[t]
        return mechanics


# constituent.py
class Constituent:
    def __init__(self, name):
        self.mechanics = MechanicsFactory.for_constituent(self)
    
    def update_stress(self, timestep: int):
        stress_trace = self.mechanics.compute_stress_trace(timestep)
        # Use stress_trace...


# mechanics.py (modified to use internal accessors)
class Mechanics:
    def compute_stress_trace(self, timestep: int) -> float:
        """Compute tr(σ) using internal accessors."""
        stress_tensor = self._get_stress(timestep)  # Uses lambda set by factory
        return np.trace(stress_tensor)
```

---

## Summary Table

| Criterion | Adapter (Current) | Factory (Proposed) |
|-----------|-------------------|-------------------|
| **Coupling** | ✅ Loose | ❌ Tight |
| **Testability** | ✅ Easy | ❌ Complex |
| **SRP** | ✅ Clear | ❌ Blurred |
| **Extensibility** | ✅ Open/Closed | ❌ Requires modification |
| **Mapping Elimination** | ❌ No | ❌ No (just moved) |
| **Consistency** | ✅ Matches kinetics | ❌ Diverges |
| **Complexity** | ✅ Simple | ❌ More complex |
| **Boilerplate** | ⚠️ One extra line | ✅ Slightly less |

**Overall:** Adapter pattern is superior for this use case.

---

## Questions for Software Engineering Discussion

1. **Performance**: Are there performance concerns with creating context objects? (Answer: No - context creation is O(1) and occurs infrequently)

2. **Boilerplate**: Is `context = ConstituentMechanicsContext(self)` too verbose? (Answer: It's explicit and clear - not a problem)

3. **Consistency**: Should we maintain parallel structure with `kinetics_interface.py`? (Answer: Yes - consistency aids understanding)

4. **Future**: What if we add more data sources (e.g., external support, fluid-structure interaction)? (Answer: Adapter handles this elegantly)

5. **Testing**: How important is unit test simplicity? (Answer: Critical for confidence in mechanics computations)

---

## Conclusion

**Recommendation: Keep the current Adapter Pattern implementation.**

The adapter pattern is the correct design choice for the mechanics interface because:
- It maintains loose coupling between computation and data structures
- It enables simple, focused unit tests
- It follows the Single Responsibility Principle
- It's consistent with the existing kinetics interface
- It's extensible without modification (Open/Closed Principle)

The Factory pattern would introduce tight coupling, complicate testing, and provide no meaningful benefits for this use case.

---

**Document Version:** 1.0  
**Last Updated:** November 24, 2025  
**Status:** Ready for team review