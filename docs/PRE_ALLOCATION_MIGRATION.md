# Pre-Allocation Migration Guide

## Overview

This document describes the migration from **dynamic append-based** history arrays to **pre-allocated static** arrays. This change will:

1. ✅ Eliminate fallback logic in kinetics/mechanics interfaces
2. ✅ Improve performance (no list resizing)
3. ✅ Prepare code for C++ translation
4. ✅ Make simulation state more explicit
5. ✅ Simplify debugging (can see full timeline with zeros for uncomputed values)

---

## Problem Statement

### Current Issue
When kinetics needs stress at timestep `t`, but stress hasn't been computed yet:

```python
# Kinetics runs BEFORE mechanics in simulation loop
context.get_stimulus('intramural_stress', t=5)
# ↓
# But stress_history only has [0,1,2,3,4] - timestep 5 doesn't exist!
# ↓ 
# Must add fallback logic: try t, then t-1, then error
```

This requires messy fallback logic scattered throughout interfaces:
```python
try:
    return self.layer.stress_history[timestep]
except IndexError:
    if timestep > 0:
        try:
            return self.layer.stress_history[timestep - 1]  # Fallback!
        except IndexError:
            pass
    raise Error(...)
```

### Root Cause
**Temporal coupling** between data existence and computation order. Data doesn't exist until computed, but computation order depends on data availability.

---

## Solution: Pre-Allocated Arrays

### Core Idea
**Allocate all history arrays at simulation start** with size `n_timesteps + 1`:

```python
# OLD: Append-based (dynamic)
self.stress_history = []
self.stress_history.append(homeostatic_stress)  # t=0
# Later...
self.stress_history.append(computed_stress)     # t=1

# NEW: Pre-allocated (static)
self.stress_history = np.zeros((n_timesteps + 1, 3, 3))  # All timesteps at once
self.stress_history[0] = homeostatic_stress     # t=0
# Later...
self.stress_history[5] = computed_stress        # t=5 (can write any timestep)
```

### Benefits

| Aspect | Append-Based | Pre-Allocated |
|--------|-------------|---------------|
| **Fallback logic** | Required (messy) | Not needed (clean) |
| **Performance** | O(n) resizing overhead | O(1) constant time |
| **Debugging** | Can't see future | Can see zeros = not computed |
| **C++ translation** | Difficult (dynamic) | Natural (static arrays) |
| **Code clarity** | `history[-1]` ambiguous | `history[t]` explicit |

---

## Architecture Changes

### 1. Flow of `n_steps` Through System

```
Simulation (knows n_steps)
    ↓ allocate_histories(n_steps + 1)
Configuration
    ↓ allocate_histories(n_steps + 1)
Layer
    ↓ allocate_histories(n_steps + 1)
Constituent
```

Each level allocates its own arrays and propagates call down.

### 2. Initialization Pattern

```python
# At simulation start (before any computation)
def allocate_histories(self, n_timesteps):
    """Allocate all history arrays upfront."""
    # Allocate empty arrays
    self.stress_history = np.zeros((n_timesteps, 3, 3))
    self.rhoR_history = np.zeros(n_timesteps)
    # ... all other histories ...
    
    # Initialize t=0 with homeostatic values
    self._set_homeostatic_values_at_t0()

def _set_homeostatic_values_at_t0(self):
    """Set known t=0 values."""
    self.stress_history[0] = self.homeostatic_stress.copy()
    self.rhoR_history[0] = self.homeostatic_density
    # ... all other t=0 values ...
```

### 3. Computation Pattern

```python
# OLD: Append
def compute_rhoR_alpha(self, target_timestep, ...):
    # ... compute rhoR_alpha ...
    self.rhoR_alpha_history.append(rhoR_alpha)  # ❌ Appends to end

# NEW: Index assignment
def compute_rhoR_alpha(self, target_timestep, ...):
    # ... compute rhoR_alpha ...
    self.rhoR_alpha_history[target_timestep] = rhoR_alpha  # ✅ Writes to index
```

### 4. Guess Pattern

```python
# OLD: Complex logic with existence checks
def guess_rhoR_alpha(self, target_timestep):
    if self._timestep_exists(target_timestep):  # Check if exists
        return self.get_rhoR_alpha(target_timestep)
    if not self._timestep_exists(target_timestep - 1):
        raise ValueError(...)
    previous = self.rhoR_alpha_history[-1]     # Get last
    self.rhoR_alpha_history.append(previous)   # Append
    return self.rhoR_alpha_history[-1]

# NEW: Simple copy
def guess_rhoR_alpha(self, target_timestep):
    """Copy previous value to current index."""
    self.rhoR_alpha_history[target_timestep] = \
        self.rhoR_alpha_history[target_timestep - 1]
    return self.rhoR_alpha_history[target_timestep]
```

---

## Implementation Plan

### Phase 1: Add Allocation Infrastructure (1-2 hours)

#### 1.1 Update `Simulation.__init__()`

```python
class Simulation:
    def __init__(self, configuration, ..., n_steps, ...):
        self.configuration = configuration
        self.n_steps = n_steps
        
        # NEW: Allocate all histories at start
        self.configuration.allocate_histories(n_steps + 1)
```

#### 1.2 Add `Configuration.allocate_histories()`

```python
class Configuration:
    def allocate_histories(self, n_timesteps):
        """Allocate history arrays for all layers and constituents."""
        print(f"Allocating histories for {n_timesteps} timesteps...")
        
        for layer in self.layers:
            layer.allocate_histories(n_timesteps)
```

#### 1.3 Add `Layer.allocate_histories()`

```python
class Layer:
    def allocate_histories(self, n_timesteps):
        """Pre-allocate history arrays."""
        self.n_timesteps = n_timesteps
        
        # Allocate layer histories (scalars)
        self.rhoR_history = np.zeros(n_timesteps)
        self.inner_radius_history = np.zeros(n_timesteps)
        self.thickness_history = np.zeros(n_timesteps)
        self.axial_stretch_history = np.zeros(n_timesteps)
        self.wss_history = np.zeros(n_timesteps)
        self.flow_rate_history = np.zeros(n_timesteps)
        
        # Allocate layer histories (tensors)
        self.stress_history = np.zeros((n_timesteps, 3, 3))
        self.F_history = np.zeros((n_timesteps, 3, 3))
        
        # Initialize t=0 with homeostatic values
        self._set_homeostatic_values_at_t0()
        
        # Allocate for constituents
        for constituent in self.constituents:
            constituent.allocate_histories(n_timesteps)
    
    def _set_homeostatic_values_at_t0(self):
        """Set t=0 values from homeostatic state."""
        self.rhoR_history[0] = self.homeostatic_density
        self.inner_radius_history[0] = self.homeostatic_inner_radius
        self.thickness_history[0] = self.homeostatic_thickness
        self.axial_stretch_history[0] = self.homeostatic_axial_stretch
        self.stress_history[0] = self.homeostatic_stress.copy()
        self.F_history[0] = self.homeostatic_F.copy()
        self.wss_history[0] = self.homeostatic_wss
        self.flow_rate_history[0] = self.homeostatic_flow_rate
```

#### 1.4 Add `Constituent.allocate_histories()`

```python
class SingleConstituent(Constituent):
    def allocate_histories(self, n_timesteps):
        """Pre-allocate history arrays."""
        self.n_timesteps = n_timesteps
        
        # Scalar histories
        self.rhoR_alpha_history = np.zeros(n_timesteps)
        self.k_alpha_history = np.zeros(n_timesteps)
        self.mR_alpha_history = np.zeros(n_timesteps)
        
        # Tensor histories
        self.stress_history = np.zeros((n_timesteps, 3, 3))
        
        # Survival history (special case - ragged array)
        self.survival_history = [None] * n_timesteps  # Will be lists
        
        # Initialize t=0 with homeostatic
        homeostatic_density = self.params['mass_fraction'] * 1050.0
        self.rhoR_alpha_history[0] = homeostatic_density
```

---

### Phase 2: Replace Append with Indexing (30 min)

#### 2.1 Find and Replace Pattern

**Search for:** `.append(`  
**Context:** Inside history updates

**Examples:**

```python
# Layer.py - Multiple replacements
# OLD:
self.rhoR_history.append(rhoR)
self.inner_radius_history.append(a)
self.stress_history.append(sigma)

# NEW:
self.rhoR_history[timestep] = rhoR
self.inner_radius_history[timestep] = a
self.stress_history[timestep] = sigma
```

```python
# Constituent.py - Multiple replacements
# OLD:
self.k_alpha_history.append(k_alpha)
self.mR_alpha_history.append(mR_alpha)
self.rhoR_alpha_history.append(rhoR_alpha)

# NEW:
self.k_alpha_history[target_timestep] = k_alpha
self.mR_alpha_history[target_timestep] = mR_alpha
self.rhoR_alpha_history[target_timestep] = rhoR_alpha
```

#### 2.2 Special Case: Survival History (Ragged)

Survival history is **ragged** (each cohort has different length):

```python
# OLD:
self.survival_history.append(survival_values)  # Append list

# NEW:
self.survival_history[target_timestep] = survival_values  # Assign to index
```

Pre-allocated as list of `None`:
```python
self.survival_history = [None] * n_timesteps
```

---

### Phase 3: Update Guess Methods (30 min)

#### 3.1 Simplify Guess Logic

```python
# Layer.py
def guess_rhoR_alpha(self, target_timestep):
    """Guess mass density at target timestep."""
    if target_timestep == 0:
        raise ValueError("Cannot guess timestep 0")
    
    # Simple copy (array already allocated)
    self.rhoR_history[target_timestep] = self.rhoR_history[target_timestep - 1]
    return self.rhoR_history[target_timestep]

# Same pattern for:
# - guess_inner_radius()
# - guess_thickness()
# - guess_stress()
```

```python
# Constituent.py
def guess_rhoR_alpha(self, target_timestep, guess_method="from_previous_timestep"):
    """Guess mass density at target timestep."""
    if target_timestep == 0:
        raise ValueError("Cannot guess timestep 0")
    
    if guess_method == "from_previous_timestep":
        # Simple copy (array already allocated)
        self.rhoR_alpha_history[target_timestep] = \
            self.rhoR_alpha_history[target_timestep - 1]
    else:
        raise ValueError(f"Unknown guess method: {guess_method}")
    
    return self.rhoR_alpha_history[target_timestep]
```

#### 3.2 Remove Helper Methods

Delete these methods (no longer needed):

```python
# ❌ DELETE - No longer needed
def _timestep_exists(self, timestep):
    """Check if a specific timestep exists in history."""
    return 0 <= timestep < len(self.rhoR_alpha_history)
```

---

### Phase 4: Update Validation (1 hour)

#### 4.1 Update Accessor Bounds Checks

```python
# Layer.py
def get_inner_radius(self, timestep: int) -> float:
    """Get inner radius at timestep."""
    # NEW: Check against n_timesteps (not len)
    if not (0 <= timestep < self.n_timesteps):
        raise IndexError(
            f"Timestep {timestep} out of range [0, {self.n_timesteps-1}]"
        )
    
    return self.inner_radius_history[timestep]

# Apply same pattern to:
# - get_thickness()
# - get_density()
# - get_stress()
# - get_axial_stretch()
# etc.
```

#### 4.2 Optional: Add Computed Check

If you want to verify a timestep has been computed (not just allocated):

```python
class Layer:
    def allocate_histories(self, n_timesteps):
        # ... existing code ...
        
        # Optional: Track which timesteps are computed
        self._computed_mask = np.zeros(n_timesteps, dtype=bool)
        self._computed_mask[0] = True  # t=0 is initialized
    
    def get_inner_radius(self, timestep: int) -> float:
        """Get inner radius at timestep."""
        if not (0 <= timestep < self.n_timesteps):
            raise IndexError(f"Timestep {timestep} out of range")
        
        # Optional: Check if computed
        if not self._computed_mask[timestep]:
            raise ValueError(f"Timestep {timestep} not yet computed")
        
        return self.inner_radius_history[timestep]
    
    def compute_geometry(self, timestep):
        # ... compute ...
        self.inner_radius_history[timestep] = a
        self._computed_mask[timestep] = True  # Mark as computed
```

**Note:** This is optional. For most use cases, trusting sequential computation order is sufficient.

---

### Phase 5: Remove Fallback Logic (30 min)

#### 5.1 Clean Up Kinetics Interface

```python
# kinetics_interface.py

# OLD: Complex fallback logic
def _get_layer_stress(self, timestep: int) -> float:
    if not self.layer:
        raise KineticsDataNotAvailableError("Layer not available")
    
    try:
        return self.layer.get_stress_trace(timestep)
    except (IndexError, ValueError, AttributeError):
        if timestep > 0:
            try:
                print(f"  [INFO]: Using t={timestep-1}")
                return self.layer.get_stress_trace(timestep - 1)
            except:
                pass
        raise KineticsDataNotAvailableError(...)

# NEW: Clean direct access
def _get_layer_stress(self, timestep: int) -> float:
    """Get intramural stress from layer (no fallback needed)."""
    if not self.layer:
        raise KineticsDataNotAvailableError("Layer not available")
    
    # Direct access - stress always exists after allocation
    return self.layer.get_stress_trace(timestep)
```

Apply same cleanup to:
- `_get_constituent_stress()`
- `_get_wall_shear_stress()`
- `_get_inflammation()`

#### 5.2 Clean Up Mechanics Interface

Same pattern:

```python
# mechanics_interface.py

# OLD: Bounds check with len()
if not (0 <= timestep < len(self.layer.stress_history)):
    raise MechanicsDataNotAvailableError(...)

# NEW: Bounds check with n_timesteps
if not (0 <= timestep < self.layer.n_timesteps):
    raise MechanicsDataNotAvailableError(
        f"Timestep {timestep} out of range [0, {self.layer.n_timesteps-1}]"
    )
```

---

## Testing Strategy

### Unit Tests

```python
def test_pre_allocation():
    """Test that histories are pre-allocated correctly."""
    config = Configuration.from_parameters(params)
    layer = config.layers[0]
    
    # Should allocate n_timesteps + 1
    n_timesteps = 101  # 100 steps + t=0
    layer.allocate_histories(n_timesteps)
    
    # Check allocation
    assert len(layer.stress_history) == n_timesteps
    assert len(layer.rhoR_history) == n_timesteps
    assert layer.stress_history.shape == (n_timesteps, 3, 3)
    
    # Check t=0 initialized
    assert layer.rhoR_history[0] == layer.homeostatic_density
    assert np.allclose(layer.stress_history[0], layer.homeostatic_stress)
    
    # Check rest is zeros
    assert layer.rhoR_history[1] == 0.0  # Not yet computed
    assert np.allclose(layer.stress_history[1], np.zeros((3,3)))

def test_index_assignment():
    """Test that index assignment works."""
    layer = Layer("test")
    layer.allocate_histories(10)
    
    # Assign to specific timestep
    layer.rhoR_history[5] = 1050.0
    assert layer.rhoR_history[5] == 1050.0
    
    # Other timesteps unaffected
    assert layer.rhoR_history[4] == 0.0  # Still zero
    assert layer.rhoR_history[6] == 0.0

def test_guess_simplification():
    """Test simplified guess method."""
    const = SingleConstituent("test")
    const.allocate_histories(10)
    const.rhoR_alpha_history[0] = 100.0  # Set t=0
    
    # Guess t=1 from t=0
    const.guess_rhoR_alpha(1)
    assert const.rhoR_alpha_history[1] == 100.0
    
    # Guess t=5 from t=4 (should fail - t=4 not set)
    const.rhoR_alpha_history[4] = 150.0
    const.guess_rhoR_alpha(5)
    assert const.rhoR_alpha_history[5] == 150.0
```

### Integration Test

```python
def test_full_simulation_with_pre_allocation():
    """Test that simulation works with pre-allocated arrays."""
    params = load_parameters("latorre2018_updated.yaml")
    config = Configuration.from_parameters(params)
    
    sim = Simulation(
        configuration=config,
        n_steps=10,  # Small test
        dt=1.0
    )
    
    # Should allocate 11 timesteps (0-10)
    layer = config.layers[0]
    assert len(layer.stress_history) == 11
    
    # Run simulation
    sim.run()
    
    # All timesteps should be filled
    for t in range(11):
        assert layer.rhoR_history[t] > 0  # Non-zero (computed)
        assert not np.allclose(layer.stress_history[t], np.zeros((3,3)))
```

---

## Migration Checklist

### Before Starting
- [ ] Create backup branch: `git checkout -b backup-before-pre-allocation`
- [ ] Run existing tests to establish baseline
- [ ] Document any custom modifications

### Phase 1: Allocation Infrastructure
- [ ] Add `Simulation.allocate_histories()` call in `__init__`
- [ ] Add `Configuration.allocate_histories()`
- [ ] Add `Layer.allocate_histories()` + `_set_homeostatic_values_at_t0()`
- [ ] Add `SingleConstituent.allocate_histories()`
- [ ] Add `MultiFiberFamilyConstituent.allocate_histories()`
- [ ] Test: Histories allocated with correct size

### Phase 2: Replace Append
- [ ] Find all `.append(` in `layer.py` → replace with `[timestep] =`
- [ ] Find all `.append(` in `constituent.py` → replace with `[target_timestep] =`
- [ ] Special case: `survival_history[target_timestep] = values`
- [ ] Test: Values written to correct indices

### Phase 3: Update Guess Methods
- [ ] Simplify `Layer.guess_*()` methods
- [ ] Simplify `Constituent.guess_rhoR_alpha()`
- [ ] Remove `_timestep_exists()` helper
- [ ] Test: Guess still works correctly

### Phase 4: Update Validation
- [ ] Update `Layer` accessor bounds checks
- [ ] Update `Constituent` accessor bounds checks
- [ ] Update `MechanicsContext` bounds checks
- [ ] Test: Proper errors for invalid timesteps

### Phase 5: Remove Fallbacks
- [ ] Clean up `kinetics_interface.py` (remove try/except fallbacks)
- [ ] Clean up `mechanics_interface.py` (update bounds checks)
- [ ] Test: No fallback logic triggered

### Final Verification
- [ ] Run full test suite
- [ ] Run short simulation (10 steps)
- [ ] Run longer simulation (100 steps)
- [ ] Check performance improvement
- [ ] Verify output matches baseline

---

## Expected Results

### Before (Append-Based)
```python
# Simulation loop creates data on-demand
for t in range(n_steps):
    guess()     # Creates new entries
    compute()   # Fills new entries
    # stress_history grows: [t=0] → [t=0, t=1] → [t=0, t=1, t=2] ...
```

**Issues:**
- ❌ Fallback logic needed everywhere
- ❌ List resizing overhead
- ❌ Can't see full timeline
- ❌ Complex bounds checking

### After (Pre-Allocated)
```python
# One-time allocation at start
allocate_histories(n_steps + 1)
# stress_history = [t=0=homeostatic, t=1=0, t=2=0, ..., t=100=0]

# Simulation loop writes to indices
for t in range(n_steps):
    guess()     # Copies to index t
    compute()   # Writes to index t
    # stress_history[t] updated in place
```

**Benefits:**
- ✅ No fallback logic needed
- ✅ O(1) constant time access
- ✅ Can see full timeline (zeros = not computed)
- ✅ Simple bounds checking (`0 <= t < n_timesteps`)

---

## Performance Comparison

### Append-Based (Current)
```python
# Worst case: O(n) list resizing every ~8 appends
for t in range(1000):
    self.history.append(value)  # Periodic O(n) resize
# Total: O(n²) in worst case
```

### Pre-Allocated (New)
```python
# Allocation: O(n) once
self.history = np.zeros(1000)

# Writes: O(1) each
for t in range(1000):
    self.history[t] = value  # Always O(1)
# Total: O(n)
```

**Expected speedup:** 2-5× for large simulations (> 365 days)

---

## C++ Translation Ready

### Python (After Migration)
```python
# Pre-allocated arrays
self.stress_history = np.zeros((n_timesteps, 3, 3))
self.stress_history[t] = sigma
```

### C++ (Direct Translation)
```cpp
// Static arrays
double stress_history[N_TIMESTEPS][3][3];
stress_history[t] = sigma;  // Same syntax!
```

**No complex translation needed** - direct array indexing maps 1:1 to C++.

---

## Rollback Plan

If issues arise:

1. **Revert to backup branch:**
   ```bash
   git checkout backup-before-pre-allocation
   ```

2. **Or fix forward:**
   - Check allocation size: `n_steps + 1` (include t=0)
   - Check index bounds: `0 <= t < n_timesteps`
   - Check t=0 initialization in `_set_homeostatic_values_at_t0()`
   - Check survival_history (special case - list of lists)

---

## Future Enhancements

After basic migration works:

### 1. Add Computed Mask (Optional)
```python
self._computed_mask = np.zeros(n_timesteps, dtype=bool)
```
Tracks which timesteps are actually computed vs just allocated.

### 2. Memory Optimization
```python
# For very long simulations, use memory-mapped files
self.stress_history = np.memmap('stress.dat', dtype=float, 
                                 shape=(n_timesteps, 3, 3))
```

### 3. Checkpoint/Restart
```python
def save_checkpoint(self, timestep):
    """Save arrays up to timestep."""
    np.savez(f'checkpoint_{timestep}.npz',
             stress=self.stress_history[:timestep+1],
             rhoR=self.rhoR_history[:timestep+1])

def load_checkpoint(self, filename):
    """Resume from checkpoint."""
    data = np.load(filename)
    self.stress_history[:len(data['stress'])] = data['stress']
    # ...
```

---

## Questions & Answers

**Q: What if I don't know `n_steps` at initialization?**  
A: This approach requires knowing simulation length upfront. If dynamic length is needed, consider chunked allocation (allocate in blocks of 100).

**Q: What about memory for very long simulations?**  
A: For n=10,000 steps × 10 constituents × 5 histories:  
`10,000 × 10 × 5 × 8 bytes = 4 MB` (negligible)  
Even 1 million steps = 400 MB (still manageable).

**Q: Can I still use append for debugging?**  
A: Yes, but only during development. Production code should use indexing.

**Q: What if a constituent is added mid-simulation?**  
A: New constituents added after `allocate_histories()` must call it themselves with remaining timesteps.

---

## Summary

### Key Changes
1. **Add allocation flow:** `Simulation → Configuration → Layer → Constituent`
2. **Replace `.append()` → `[timestep] =`** throughout
3. **Simplify guess methods** (no existence checks needed)
4. **Remove fallback logic** (data always available)
5. **Update bounds checks** (use `n_timesteps` not `len()`)

### Impact
- **Code:** -200 lines (delete fallback logic, helper methods)
- **Performance:** 2-5× faster for large simulations
- **Clarity:** Explicit timestep indexing vs implicit append
- **Maintenance:** Simpler, fewer edge cases

### Timeline
- **Phase 1-2:** 2 hours (infrastructure + find/replace)
- **Phase 3-4:** 1 hour (simplify + validation)
- **Phase 5:** 30 min (cleanup)
- **Testing:** 1-2 hours
- **Total:** ~5 hours for complete migration

---

**Ready to proceed?** Start with Phase 1, test incrementally, and migrate one phase at a time! 🚀