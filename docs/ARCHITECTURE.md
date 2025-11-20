# pyGrowth Architecture

## Data Ownership Principle

### Core Rule: Separation of Data Ownership and Data Access

**Data Storage (Ownership):**
- `Layer` owns layer-level histories:
  - `inner_radius_history[]`
  - `thickness_history[]`
  - `lambda_z_history[]`
  - `pressure_history[]`
  - `intramural_stress_history[]`

- `Constituent` owns constituent-level histories:
  - `rhoR_alpha_history[]`
  - `k_alpha_history[]`
  - `mR_alpha_history[]`
  - `q_history[]`
  - `sigma_hat_history[]`
  - `sigma_alpha_history[]`

**Data Access (via Context Interfaces):**
- Context objects provide access without transferring ownership
- `LayerMechanicalContext` - wraps Layer + Kinematics
- `ConstituentKineticsContext` - wraps Constituent
- `ConstitutiveModelContext` - wraps mechanical state

**Algorithm Classes (Stateless):**
- `Kinetics`, `Mechanics`, `ConstitutiveModel` never store histories
- Operate only on context objects
- Compute and return values
- Results stored back in owner (Constituent/Layer)

### Example Pattern

```python
# Data stored at owner
self.rhoR_alpha_history = []  # Constituent owns this

# Access via context
context = ConstituentKineticsContext(self)

# Algorithm operates on context (stateless)
k_alpha = self.kinetics.compute_k_alpha(context, target_timestep)

# Result stored back at owner
self.k_alpha_history.append(k_alpha)