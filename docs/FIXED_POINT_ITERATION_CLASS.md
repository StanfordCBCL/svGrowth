# Fixed-Point Iteration Solver for Mass-Geometry Coupling

## Overview

The `FixedPointSolver` addresses a fundamental coupling in constrained mixture growth and remodeling (G&R) models: **mass production depends on stress, which depends on geometry, which depends on mass**. This creates a nonlinear coupled system that requires iterative refinement to solve.

## The Coupled Problem

### Physical Dependencies

```
Mass Production → Depends on stress/strain stimuli
       ↓
Total Mass Density (ρ)
       ↓
Incompressibility (J = ρ_h/ρ) → Determines geometry
       ↓
Geometry (a, h)
       ↓
Mixture Stress (via equilibrium) → Feeds back to mass production
       ↑_______________________________________________|
```

### Mathematical Formulation

At each time step `s`, we need to solve the coupled system:

**Mass Evolution:**
```
ρ_α(s) = ∫₀ˢ mR_α(τ, σ(s), a(s)) · q(s,τ) dτ
```
where `mR_α` (mass production) depends on current stress `σ(s)` and geometry `a(s)`.

**Geometric Equilibrium:**
```
σ_θθ(mixture) = P · a(s) / h(s)
```
where geometry is constrained by incompressibility:
```
J = ρ_h / ρ(s) = det(F)
```

### The Challenge

This is a **fixed-point problem** because:
1. To compute mass, we need stress (which requires geometry)
2. To compute geometry, we need mass (through incompressibility)
3. Neither can be solved independently!

---

## Solution Strategy: Fixed-Point Iteration

### Algorithm

```
Input: Initial guess for ρ(s), a(s)
Output: Converged ρ(s), a(s) satisfying both equations

1. Initialize: ρ⁽⁰⁾ = ρ(s-1), a⁽⁰⁾ = a(s-1)

2. For k = 1, 2, ..., max_iter:
   
   a) Compute mass production with current geometry/stress:
      ρ⁽ᵏ⁾ = ∫₀ˢ mR(τ, σ⁽ᵏ⁻¹⁾, a⁽ᵏ⁻¹⁾) · q(s,τ) dτ
   
   b) Solve geometric equilibrium with updated mass:
      Find a⁽ᵏ⁾, h⁽ᵏ⁾ such that:
        - σ_θθ(mixture) = P · a / h
        - J = ρ_h / ρ⁽ᵏ⁾
   
   c) Check convergence:
      If |ρ⁽ᵏ⁾ - ρ⁽ᵏ⁻¹⁾| / ρ⁽ᵏ⁻¹⁾ < tol:
         CONVERGED → Return ρ⁽ᵏ⁾, a⁽ᵏ⁾
   
   d) Update for next iteration:
      ρ⁽ᵏ⁻¹⁾ ← ρ⁽ᵏ⁾, a⁽ᵏ⁻¹⁾ ← a⁽ᵏ⁾

3. If max_iter reached without convergence:
   Return failure with diagnostics
```

### Convergence Criteria

**Primary:** Relative change in total mass density
```python
convergence_metric = |ρ_new - ρ_old| / ρ_old
converged = convergence_metric < tolerance
```

**Typical values:**
- `tolerance = 1e-6` (very strict)
- `tolerance = 1e-4` (reasonable for most cases)
- `max_iterations = 50` (sufficient for well-posed problems)

---

## Implementation Design

### Class Architecture

```python
class FixedPointSolver:
    """
    Solves coupled mass-geometry system via fixed-point iteration.
    
    Responsibilities:
    - Orchestrate iteration between mass computation and geometry solve
    - Check convergence on total mass density
    - Handle iteration failures gracefully
    - Provide diagnostics (iteration history, convergence metrics)
    
    NOT responsible for:
    - Computing mass production (delegates to Configuration)
    - Solving equilibrium (delegates to equilibrium solver)
    - Storing simulation state (delegates to Configuration)
    """
```

### Key Methods

#### `__init__(tolerance, max_iterations)`
Initialize solver with convergence parameters.

#### `solve(configuration, timestep, dt, ...)`
Main solve method:
- **Input:** Configuration object, timestep, numerical parameters
- **Output:** Dictionary with convergence results
- **Side effects:** Updates configuration state iteratively

#### Return Value Structure
```python
{
    'converged': bool,              # True if converged within tolerance
    'iterations': int,              # Number of iterations performed
    'final_residual': float,        # Final |Δρ|/ρ value
    'rho_history': [ρ₁, ρ₂, ...],  # Mass density at each iteration
    'message': str                  # Status/error message (if any)
}
```

---

## Integration with Existing Code

### Current Workflow (Missing Iteration)

```python
# simulation.py - CURRENT (INCOMPLETE)
def run(self):
    for step in range(self.n_steps):
        next_timestep = self.current_timestep + 1
        
        # Guess
        self.configuration.guess_all_rhoR_alpha(next_timestep)
        self.configuration.guess_geometry(next_timestep)
        
        # Compute once (NO ITERATION!)
        self.configuration.compute_all_rhoR(next_timestep, ...)
        self.configuration.solve_equilibrium_geometry(next_timestep, ...)
        
        # Problem: Mass and geometry may not be consistent!
        self.current_timestep = next_timestep
```

### With FixedPointSolver

```python
# simulation.py - WITH FIXED-POINT ITERATION
from fixed_point_solver import FixedPointSolver

class Simulation:
    def __init__(self, ...):
        self.fixed_point_solver = FixedPointSolver(
            tolerance=self.tolerance,
            max_iterations=self.max_iterations
        )
    
    def run(self):
        for step in range(self.n_steps):
            next_timestep = self.current_timestep + 1
            
            # Initialize guesses
            self.configuration.guess_all_rhoR_alpha(next_timestep)
            self.configuration.guess_geometry(next_timestep)
            self.configuration.guess_stress_and_wss(next_timestep)
            
            # Solve coupled system iteratively
            result = self.fixed_point_solver.solve(
                configuration=self.configuration,
                timestep=next_timestep,
                dt=self.dt,
                integration_method=self.integration_method,
                survival_function_computation=self.survival_function_computation,
                equilibrium_solver_params={
                    'solver_method': 'toms748',
                    'tolerance': 1e-3,
                    'verbose': self.verbose
                },
                verbose=self.verbose
            )
            
            if not result['converged']:
                raise RuntimeError(f"Fixed-point iteration failed: {result['message']}")
            
            self.current_timestep = next_timestep
```

---

## Required Helper Methods

### In `configuration.py`

```python
def get_total_density(self, timestep: int) -> float:
    """Get total referential mass density at timestep.
    
    Sums density across all layers (for multi-layer models) or
    returns single layer density.
    
    Args:
        timestep: Timestep index
        
    Returns:
        Total referential mass density (kg/m³)
    """
    total_rho = 0.0
    for layer in self.layers:
        total_rho += layer.get_density(timestep)
    return total_rho
```

---

## Typical Convergence Behavior

### Well-Posed Problem (Near Homeostasis)
```
Iteration 1: ρ = 1050.0417 kg/m³, Δρ/ρ = 3.97e-05
✓ Converged in 1 iteration
```

### Moderate Perturbation
```
Iteration 1: ρ = 1052.34 kg/m³, Δρ/ρ = 2.23e-03
Iteration 2: ρ = 1052.41 kg/m³, Δρ/ρ = 6.65e-05
✓ Converged in 2 iterations
```

### Large Perturbation
```
Iteration 1: ρ = 1065.12 kg/m³, Δρ/ρ = 1.43e-02
Iteration 2: ρ = 1066.89 kg/m³, Δρ/ρ = 1.66e-03
Iteration 3: ρ = 1067.01 kg/m³, Δρ/ρ = 1.12e-04
Iteration 4: ρ = 1067.02 kg/m³, Δρ/ρ = 9.37e-06
✓ Converged in 4 iterations
```

### Divergence Warning Signs
```
Iteration 1: ρ = 1050.00 kg/m³, Δρ/ρ = 5.00e-04
Iteration 2: ρ = 1050.05 kg/m³, Δρ/ρ = 4.76e-04
Iteration 3: ρ = 1050.10 kg/m³, Δρ/ρ = 4.76e-04
...
Iteration 50: ρ = 1052.38 kg/m³, Δρ/ρ = 4.55e-04
⚠️ Max iterations reached (oscillating or slow convergence)
```

---

## Comparison with Other Solvers

| Solver | Purpose | Input | Output |
|--------|---------|-------|--------|
| **FixedPointSolver** | Mass-geometry coupling | Configuration state | Converged state |
| **EquilibriumSolver** | Geometric equilibrium | Trial geometry | Equilibrium geometry |
| **Root-finding (scipy)** | Find zeros of f(x) | Function, bracket | Root x where f(x)=0 |

**Nesting structure:**
```
FixedPointSolver
  └─> calls Configuration.solve_equilibrium_geometry()
       └─> calls EquilibriumSolver
            └─> calls scipy.root_scalar
```

---

## Troubleshooting

### Problem: Slow Convergence (many iterations)

**Possible causes:**
- Strong coupling between mass production and stress
- Large time step `dt`
- Stiff production function (e.g., high stress sensitivity)

**Solutions:**
- Reduce time step
- Use under-relaxation: `ρ_new = ω·ρ_computed + (1-ω)·ρ_old` with `ω < 1`
- Increase tolerance slightly

### Problem: Oscillation (Δρ not decreasing)

**Possible causes:**
- Numerical instability
- Poor initial guess
- Incompatible equilibrium/mass constraints

**Solutions:**
- Check equilibrium solver is converging properly
- Verify mass production formula
- Try different initial guess strategy

### Problem: Max Iterations Reached

**Possible causes:**
- Tolerance too strict
- Problem genuinely doesn't converge (check physics)
- Numerical issues in sub-solvers

**Solutions:**
- Increase max_iterations
- Relax tolerance
- Add diagnostics to check residual trend

---

## Future Enhancements

### 1. Adaptive Relaxation
```python
# Automatically adjust relaxation parameter based on convergence
if oscillating:
    omega = max(0.5, omega * 0.9)  # Dampen updates
else:
    omega = min(1.0, omega * 1.1)  # Speed up
```

### 2. Multi-Level Convergence
```python
# Check convergence on multiple quantities
converged = (
    mass_converged and
    geometry_converged and
    stress_converged
)
```

### 3. Anderson Acceleration
Advanced technique to accelerate fixed-point iteration by combining multiple previous iterates.

### 4. Line Search
If updates are too large, use line search to find optimal step size.

---

## References

### Mathematical Background
- **Fixed-point iteration**: Numerical Analysis textbooks (e.g., Burden & Faires)
- **Constrained mixture theory**: Humphrey & Rajagopal (2002), Baek et al. (2006)

### Similar Implementations
- Latorre et al. (2018): "A mechanobiologically equilibrated constrained mixture model"
- Matlab/Fortran implementations in `other_implementations/` folder

---

## Implementation Checklist

When migrating from Option 1 to Option 2:

- [ ] Create `fixed_point_solver.py` with `FixedPointSolver` class
- [ ] Add `get_total_density()` method to `Configuration`
- [ ] Update `Simulation.__init__()` to instantiate solver
- [ ] Replace explicit iteration loop in `Simulation.run()` with solver call
- [ ] Add unit tests for `FixedPointSolver`
- [ ] Verify convergence behavior matches Option 1
- [ ] Update simulation YAML to include fixed-point parameters

---

## Notes

- **Thread safety:** Current design assumes single-threaded execution
- **Memory:** Stores full iteration history (minimal overhead for typical cases)
- **Performance:** Convergence typically within 1-5 iterations for well-posed problems
- **Robustness:** Gracefully handles non-convergence with detailed diagnostics

---

*Last updated: 2025-12-09*
*Document version: 1.0*