# svGrowth Project Context

## Project Overview
**svGrowth** is a Python framework for simulating vascular Growth & Remodeling (G&R) using constrained mixture theory with heredity integrals.

## Required Reading at Session Start

### Architecture Documentation (`docs/`)
Always review these design documents first:
- **ARCHITECTURE.md** - Core data ownership principles, separation of concerns
- **FIXED_POINT_ITERATION_CLASS.md** - Coupling algorithm between mass and geometry
- **STRESS_STATE_DESIGN.md** - How stress computations work
- **MECHANICS_INTERFACE_DESIGN_DISCUSSION.md** - Mechanics abstraction layer
- **PRE_ALLOCATION_MIGRATION.md** - Memory management patterns
- **TESTING_FRAMEWORK_DESIGN.md** - Testing philosophy

### API Documentation (`documentation/`)
Reference these for implementation details:
- **configuration.txt** - Configuration class API
- **constituent.txt** - Constituent class API
- **layer.txt** - Layer class API
- **kinetics.txt** - Kinetics computation API
- **io_handler.txt** - File I/O patterns
- **custom_logging.txt** - Logging configuration

## Core Architectural Principles

### 1. Data Ownership Separation
- **Owners**: `Layer` owns geometry/loading histories, `Constituent` owns mass/kinetics histories
- **Algorithms**: `Kinetics`, `Mechanics`, `ConstitutiveModel` are stateless
- **Access**: Context objects provide read-only access without transferring ownership
- **Hierarchy**: `Configuration` → `Layer` → `Constituent`

### 2. No String Literals - Use Enum Registry Pattern
**CRITICAL**: Never use scattered string literals for configuration options, coordinate systems, or method names.

✅ **CORRECT** (Enum with Registry):
```python
class CoordinateSystemType(Enum):
    CYLINDRICAL = 'cylindrical'
    CARTESIAN = 'cartesian'
    SPHERICAL = 'spherical'
```

❌ **WRONG** (String literals):
```python
if coord_system == "cylindrical":  # DON'T DO THIS
    ...
```

**Reference Implementation**: See `src/coordinate_systems.py` for the canonical example of Enum + Registry pattern.

### 3. Best Practices (Always Follow)
- **Type hints**: All functions must have complete type annotations
- **Docstrings**: Google-style docstrings for all public methods
- **Immutability**: Use `@dataclass(frozen=True)` for value objects
- **Pure functions**: Prefer pure functions over stateful methods
- **Explicit over implicit**: No magic values or hidden state
- **DRY principle**: Use registries and factories to eliminate duplication

### 4. C++ Translation Readiness (Bonus)
While not required, prefer patterns that translate well to C++:
- Use explicit data structures over dicts where possible
- Avoid dynamic typing tricks
- Prefer composition over inheritance
- Use dataclasses/named tuples over anonymous dicts
- Keep class hierarchies shallow

## Code Organization
```
src/              # Source code
tests/            # Unit, integration, e2e tests
examples/         # Library of example parameter files and notebooks
benchmarks/       # Performance profiling
docs/             # Architecture documentation
documentation/    # API documentation
```

## Development Workflow
1. **Check documentation first** - Understand existing patterns before suggesting changes
2. **Maintain separation of concerns** - Data owners vs. algorithms 
3. **Follow Enum+Registry pattern** - No string literals for options/types
4. **Type safety** - Complete type hints and enum-based type checking

## When Suggesting Changes
- Maintain consistency with existing patterns
- Suggest enum additions rather than string literal usage
- Consider C++ translation if the pattern is complex
- Provide type hints for all new code
