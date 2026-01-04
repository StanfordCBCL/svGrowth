User: Take a look at the docs folder to familiriaze with the svGrowth architecture. I would like to create a testing framework, comprising units tests, integration tests and end-to-end tests for this codebase. I have created a tests folder for that. Let's brainstorn on what the architecture of the testing framework would look like, some points to think about:
1) what should be the kinds of inputs and outputs strcuture, where and how will they be stored? 
2) how are tests going to be integrated wtih Github actions for conntinuous integration
3) how are tests going to be interacting with the exceptions.py and custom_loging.py calsses? should they have their own exceptions and logging or use the same classes as the rest of the code?
4) should the tests have their own outputs or output format? how to structure it in a modular clean way?

In terms of content here are some examples: 
1) for units tests I would like to test integration methods from integrators.py (naive trapezoidal, backward trapezoidal, naive simpson, backward simpson).
2) for integation tests, I want to test computing survival function from kientics.py or compute_rhoR_alpha functions in constituent.
3) for end-to-end test, we can use the input file for latore2018.yaml where all simulation outputs should be matching known solutions.

Make sure to also give a breif explanation of the purpose of these tests (for example what is the different between integration and unit tests?).

GitHub Copilot: I'll review the documentation to understand the svGrowth architecture and then help you design a comprehensive testing framework. 

Read [](file:///Users/lazaros/workspace/github/svGrowth/docs/ARCHITECTURE.md)

Read [](file:///Users/lazaros/workspace/github/svGrowth/docs/STRESS_STATE_DESIGN.md#1-1), lines 1 to 100

Read [](file:///Users/lazaros/workspace/github/svGrowth/src/exceptions.py#1-1), lines 1 to 50

Read [](file:///Users/lazaros/workspace/github/svGrowth/src/custom_logging.py#1-1), lines 1 to 50

Read [](file:///Users/lazaros/workspace/github/svGrowth/docs/ARCHITECTURE.md#48-48), lines 48 to 48

Read [](file:///Users/lazaros/workspace/github/svGrowth/docs/FIXED_POINT_ITERATION_CLASS.md#1-1), lines 1 to 100

Read [](file:///Users/lazaros/workspace/github/svGrowth/src/integrators.py#1-1), lines 1 to 50

Now let me look at more of the codebase to understand better: 

Read [](file:///Users/lazaros/workspace/github/svGrowth/src/kinetics.py#1-1), lines 1 to 100

Read [](file:///Users/lazaros/workspace/github/svGrowth/src/constituent.py#1-1), lines 1 to 100

Read [](file:///Users/lazaros/workspace/github/svGrowth/tests)

Perfect! Now I have a good understanding of the architecture. Let me create a comprehensive testing framework design document for you.

Based on my review of the svGrowth architecture, I'll help you design a comprehensive testing framework. Let me create a detailed plan:

## **Testing Framework Architecture for svGrowth**

### **Test Types Overview**

**1. Unit Tests** - Test individual functions/methods in isolation
- **Purpose**: Verify that smallest code units (functions, methods) work correctly
- **Examples**: Integration methods, helper functions, single calculations
- **Speed**: Fast (milliseconds)
- **Dependencies**: Minimal - mock external dependencies

**2. Integration Tests** - Test interactions between multiple components
- **Purpose**: Verify that different modules work together correctly
- **Examples**: Kinetics + Integrators, Constituent methods using multiple classes
- **Speed**: Medium (seconds)
- **Dependencies**: Real classes working together, but not full simulation

**3. End-to-End (E2E) Tests** - Test complete simulation workflows
- **Purpose**: Verify entire system produces correct outputs for known inputs
- **Examples**: Full Latorre2018 simulation with validated results
- **Speed**: Slow (minutes)
- **Dependencies**: Everything - full simulation stack

---

## **Proposed Testing Framework Architecture**

```
tests/
├── __init__.py
├── conftest.py                      # Pytest fixtures and configuration
├── test_config.yaml                 # Testing-specific configuration
│
├── unit/                            # Unit tests
│   ├── __init__.py
│   ├── test_integrators.py         # Test integration methods
│   ├── test_tensor_operations.py   # Test tensor utilities
│   ├── test_constitutive_laws.py   # Test stress-strain laws
│   └── test_helper_functions.py    # Test utility functions
│
├── integration/                     # Integration tests
│   ├── __init__.py
│   ├── test_kinetics.py            # Test survival/production computations
│   ├── test_constituent.py         # Test compute_rhoR_alpha workflows
│   ├── test_mechanics.py           # Test stress computation workflows
│   └── test_solver_components.py   # Test solver iterations
│
├── e2e/                            # End-to-end tests
│   ├── __init__.py
│   ├── test_latorre2018.py        # Full Latorre2018 validation
│   └── test_perturbations.py      # Test perturbation scenarios
│
├── fixtures/                       # Test data and configurations
│   ├── inputs/                    # Input YAML files for tests
│   │   ├── latorre2018_test.yaml
│   │   ├── simple_cylinder.yaml
│   │   └── unit_test_configs/
│   ├── expected_outputs/          # Known-good reference outputs
│   │   ├── latorre2018_validated/
│   │   │   ├── rhoR_alpha.npy
│   │   │   ├── geometry.npy
│   │   │   └── metadata.json
│   │   └── integration_test_data/
│   └── test_data/                 # Synthetic data for unit tests
│       ├── integrator_test_data.json
│       └── survival_function_data.json
│
├── utils/                         # Testing utilities
│   ├── __init__.py
│   ├── assertions.py             # Custom assertion helpers
│   ├── comparisons.py            # Numerical comparison utilities
│   ├── fixtures_generators.py    # Generate test data
│   └── mock_objects.py           # Mock classes for testing
│
└── outputs/                       # Test run outputs (gitignored)
    ├── unit/
    ├── integration/
    └── e2e/
```

---

## **Detailed Answers to Your Questions**

### **1. Input/Output Structure and Storage**

**Input Structure:**
```
tests/fixtures/inputs/
├── unit_test_configs/          # Minimal configs for unit tests
│   ├── simple_integrator.yaml
│   └── single_constituent.yaml
├── integration_configs/        # Partial simulation configs
│   ├── kinetics_only.yaml
│   └── two_constituent.yaml
└── e2e_configs/               # Full simulation configs
    ├── latorre2018_test.yaml
    └── validation_cases/
```

**Expected Output Structure:**
```
tests/fixtures/expected_outputs/
├── latorre2018_validated/
│   ├── metadata.json          # Simulation parameters, tolerances
│   ├── geometry.npy           # [timesteps × 3] (a, h, λz)
│   ├── rhoR_alpha.npy         # [timesteps × n_constituents]
│   └── stress.npy             # [timesteps × 3] (σr, σθ, σz)
└── integration_test_data/
    ├── survival_function_expected.json
    └── rhoR_alpha_computation_expected.json
```

**Test Output Structure:**
```
tests/outputs/                  # Gitignored - runtime outputs
├── test_run_<timestamp>/
│   ├── logs/
│   │   ├── test_integration.log
│   │   └── test_e2e.log
│   ├── failures/              # Only written on failure
│   │   ├── test_latorre2018_diff.png
│   │   └── actual_vs_expected.json
│   └── coverage/
│       └── htmlcov/
```

---

### **2. GitHub Actions CI/CD Integration**

**Proposed `.github/workflows/tests.yml`:**
```yaml
name: svGrowth Tests

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.9', '3.10', '3.11']
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install -r tests/requirements_test.txt  # pytest, coverage, etc.
    
    - name: Run unit tests
      run: |
        pytest tests/unit/ -v --cov=src --cov-report=xml
    
    - name: Run integration tests
      run: |
        pytest tests/integration/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Run E2E tests
      run: |
        pytest tests/e2e/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        files: ./coverage.xml
        flags: unittests
        name: codecov-umbrella

  test-fast:
    # Quick smoke test for rapid feedback
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    - run: |
        pip install -r requirements.txt -r tests/requirements_test.txt
        pytest tests/unit/ -v --maxfail=1  # Fail fast
```

**Test Requirements File (`tests/requirements_test.txt`):**
```
pytest>=7.4.0
pytest-cov>=4.1.0
pytest-xdist>=3.3.0      # Parallel test execution
pytest-timeout>=2.1.0     # Timeout protection
hypothesis>=6.82.0        # Property-based testing
numpy-testing>=0.1.0      # Numerical assertions
```

---

### **3. Exceptions and Logging Integration**

**Recommendation: Reuse with Test-Specific Extensions**

**Strategy:**
- **Reuse exceptions.py** - Tests should raise and catch the same exceptions
- **Extend logging with test mode** - Add test-specific handlers

**Implementation:**

Create `tests/utils/test_logging.py`:
```python
"""Test-specific logging configuration."""
from custom_logging import get_logger
import logging
from pathlib import Path

def setup_test_logging(test_name: str, output_dir: Path):
    """Configure logging for test execution.
    
    Args:
        test_name: Name of test for log file
        output_dir: Directory for log files
    """
    log_file = output_dir / f"{test_name}.log"
    
    # Create test-specific handler
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    )
    
    # Add to root logger
    logging.root.addHandler(handler)
    logging.root.setLevel(logging.DEBUG)  # Capture all for tests
    
    return log_file
```

**In test files:**
```python
from src.exceptions import ConvergenceError, PhysicsViolationError
from src.custom_logging import get_logger

logger = get_logger(__name__)

def test_solver_convergence_failure():
    """Test that solver raises ConvergenceError appropriately."""
    with pytest.raises(ConvergenceError) as exc_info:
        # ... code that should fail to converge
    
    assert "max iterations" in str(exc_info.value).lower()
    logger.info(f"Caught expected error: {exc_info.value}")
```

**Benefits:**
- **Consistency**: Same error types across codebase and tests
- **Debugging**: Test logs use same format as production
- **Isolation**: Test logs written separately, don't pollute main logs

---

### **4. Test Output Structure - Modular and Clean**

**Design Principles:**
1. **Gitignore test outputs** - Don't commit runtime artifacts
2. **Structured storage** - Organize by test type and timestamp
3. **Failure-only details** - Only save detailed outputs on failure
4. **Easy comparison** - Store in formats easy to diff (JSON, NPY)

**Implementation:**

Create `tests/utils/output_manager.py`:
```python
"""Manage test outputs in a clean, modular way."""
from pathlib import Path
from datetime import datetime
import json
import numpy as np
import shutil

class TestOutputManager:
    """Handles test output storage and cleanup."""
    
    def __init__(self, test_name: str, test_type: str):
        """
        Args:
            test_name: Name of test (e.g., 'test_latorre2018')
            test_type: 'unit', 'integration', or 'e2e'
        """
        self.test_name = test_name
        self.test_type = test_type
        
        # Create timestamped output directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(__file__).parent.parent / "outputs" / test_type / f"{test_name}_{timestamp}"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def save_array(self, name: str, array: np.ndarray):
        """Save numpy array."""
        np.save(self.output_dir / f"{name}.npy", array)
    
    def save_dict(self, name: str, data: dict):
        """Save dictionary as JSON."""
        with open(self.output_dir / f"{name}.json", 'w') as f:
            json.dump(data, f, indent=2)
    
    def save_failure_comparison(self, actual, expected, name: str):
        """Save detailed comparison for failed test."""
        failure_dir = self.output_dir / "failures"
        failure_dir.mkdir(exist_ok=True)
        
        # Save both arrays
        np.save(failure_dir / f"{name}_actual.npy", actual)
        np.save(failure_dir / f"{name}_expected.npy", expected)
        
        # Save difference
        diff = actual - expected
        np.save(failure_dir / f"{name}_diff.npy", diff)
        
        # Save metadata
        metadata = {
            "max_absolute_error": float(np.max(np.abs(diff))),
            "mean_absolute_error": float(np.mean(np.abs(diff))),
            "relative_error": float(np.max(np.abs(diff / expected))),
        }
        self.save_dict(f"failures/{name}_metrics", metadata)
    
    def cleanup_on_success(self):
        """Remove output directory if test passed."""
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
```

**Usage in tests:**
```python
def test_latorre2018_full_simulation():
    output_mgr = TestOutputManager('latorre2018', 'e2e')
    
    # Run simulation
    results = run_simulation('tests/fixtures/inputs/latorre2018_test.yaml')
    
    # Load expected
    expected = np.load('tests/fixtures/expected_outputs/latorre2018_validated/geometry.npy')
    
    try:
        np.testing.assert_allclose(results.geometry, expected, rtol=1e-5)
        output_mgr.cleanup_on_success()  # Clean up if passed
    except AssertionError:
        output_mgr.save_failure_comparison(results.geometry, expected, 'geometry')
        raise  # Re-raise for pytest
```

---

## **Example Test Files**

### **Unit Test Example: `tests/unit/test_integrators.py`**

```python
"""Unit tests for numerical integrators."""
import pytest
import numpy as np
from src.integrators import (
    TrapezoidIntegrator, 
    BackwardTrapezoidIntegrator,
    SimpsonIntegrator,
    BackwardSimpsonIntegrator
)

class TestTrapezoidIntegrator:
    """Test naive trapezoidal integration."""
    
    def test_constant_function(self):
        """Integral of constant function f(x)=c should be c*(b-a)."""
        dt = 0.1
        f = [5.0] * 10  # Constant function f=5
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=9)
        
        result = integrator.integrate(f)
        expected = 5.0 * 0.9  # c * (9-0) * dt
        
        assert abs(result - expected) < 1e-10
    
    def test_linear_function(self):
        """Integral of f(x)=x should be x²/2."""
        dt = 1.0
        f = [0.0, 1.0, 2.0, 3.0, 4.0]  # f(x) = x
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=4)
        
        result = integrator.integrate(f)
        expected = 8.0  # ∫₀⁴ x dx = 4²/2 = 8
        
        assert abs(result - expected) < 1e-6
    
    @pytest.mark.parametrize("n_points,expected_error", [
        (10, 1e-3),
        (100, 1e-5),
        (1000, 1e-7),
    ])
    def test_convergence_rate(self, n_points, expected_error):
        """Test that error decreases with finer discretization."""
        # Integrate sin(x) from 0 to π
        x = np.linspace(0, np.pi, n_points)
        f = np.sin(x)
        dt = x[1] - x[0]
        
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=n_points-1)
        result = integrator.integrate(f.tolist())
        
        expected = 2.0  # ∫₀^π sin(x) dx = 2
        error = abs(result - expected)
        
        assert error < expected_error
```

### **Integration Test Example: `tests/integration/test_kinetics.py`**

```python
"""Integration tests for kinetics computations."""
import pytest
import numpy as np
from src.kinetics import Kinetics
from src.constituent import Constituent
from src.kinetics_interface import ConstituentKineticsContext

class TestSurvivalFunctionComputation:
    """Test survival function q(s,tau) computation."""
    
    @pytest.fixture
    def simple_constituent(self):
        """Create a minimal constituent for testing."""
        # Use test utility to create mock constituent
        from tests.utils.mock_objects import create_test_constituent
        return create_test_constituent(
            k_alpha_values=[0.01] * 10,  # Constant degradation
            timesteps=10
        )
    
    def test_survival_decreases_monotonically(self, simple_constituent):
        """Survival function should decrease with cohort age."""
        context = ConstituentKineticsContext(simple_constituent)
        kinetics = simple_constituent.kinetics
        
        current_timestep = 5
        q_values = kinetics.compute_survival_history(context, current_timestep)
        
        # q(s,s) = 1, and q should decrease for older cohorts
        assert q_values[-1] == pytest.approx(1.0)  # Most recent cohort
        assert all(q_values[i] <= q_values[i+1] for i in range(len(q_values)-1))
    
    def test_survival_with_known_solution(self):
        """Test against analytical solution for constant k_alpha."""
        k_const = 0.1  # 1/day
        dt = 1.0  # day
        s = 10  # Current time
        tau = 5  # Deposition time
        
        # Analytical: q(s,tau) = exp(-k * (s-tau))
        expected = np.exp(-k_const * (s - tau))
        
        # Create constituent with constant k_alpha
        from tests.utils.fixtures_generators import create_constant_degradation_constituent
        constituent = create_constant_degradation_constituent(k_const, s+1, dt)
        
        context = ConstituentKineticsContext(constituent)
        kinetics = constituent.kinetics
        q_history = kinetics.compute_survival_history(context, s)
        
        actual = q_history[tau]
        assert actual == pytest.approx(expected, rel=1e-5)
```

### **E2E Test Example: `tests/e2e/test_latorre2018.py`**

```python
"""End-to-end test for Latorre2018 validation case."""
import pytest
import numpy as np
from pathlib import Path
from src.simulation import Simulation
from tests.utils.output_manager import TestOutputManager
from tests.utils.comparisons import assert_simulation_matches

class TestLatorre2018Validation:
    """Validate against published Latorre et al. 2018 results."""
    
    @pytest.fixture
    def config_path(self):
        return Path(__file__).parent.parent / "fixtures" / "inputs" / "latorre2018_test.yaml"
    
    @pytest.fixture
    def expected_outputs(self):
        base_path = Path(__file__).parent.parent / "fixtures" / "expected_outputs" / "latorre2018_validated"
        return {
            'geometry': np.load(base_path / 'geometry.npy'),
            'rhoR_alpha': np.load(base_path / 'rhoR_alpha.npy'),
            'stress': np.load(base_path / 'stress.npy'),
        }
    
    def test_full_simulation_matches_validation(self, config_path, expected_outputs):
        """Run full Latorre2018 and compare to validated outputs."""
        output_mgr = TestOutputManager('latorre2018_full', 'e2e')
        
        # Run simulation
        sim = Simulation.from_yaml(config_path)
        results = sim.run()
        
        # Extract key outputs
        actual_geometry = np.array([[
            layer.inner_radius_history[i],
            layer.thickness_history[i],
            layer.lambda_z_history[i]
        ] for i in range(len(layer.inner_radius_history))])
        
        try:
            # Compare with tolerances
            assert_simulation_matches(
                actual_geometry, 
                expected_outputs['geometry'],
                rtol=1e-4,
                atol=1e-8,
                output_mgr=output_mgr,
                name='geometry'
            )
            output_mgr.cleanup_on_success()
        except AssertionError:
            # Outputs saved automatically by assert_simulation_matches
            raise
    
    @pytest.mark.slow
    def test_convergence_tolerance_sensitivity(self, config_path):
        """Verify results are insensitive to solver tolerance (within reason)."""
        tolerances = [1e-4, 1e-6, 1e-8]
        results = []
        
        for tol in tolerances:
            sim = Simulation.from_yaml(config_path)
            sim.solver.tolerance = tol
            result = sim.run()
            results.append(result.geometry_final)
        
        # Results should be nearly identical
        for i in range(len(results)-1):
            np.testing.assert_allclose(results[i], results[i+1], rtol=1e-3)
```

---

## **Next Steps**

Would you like me to:

1. **Create the full testing directory structure** with template files?
2. **Implement specific test files** (e.g., `test_integrators.py` with all 4 integrator types)?
3. **Set up the GitHub Actions workflow** file?
4. **Create test utilities** (`assertions.py`, `comparisons.py`, `mock_objects.py`)?
5. **Generate fixture data** for your specific test cases?

Let me know which parts you'd like me to implement first!

User: Before procceding, can you give me some context on how pytest works here? for example what are the pytest.mark.slow or pytest.fixture. Can you also explain how github will now what part of the code is covered? Also, should we have an assertions.py for the main code (to be reused for tests)? What is the difference between exceptions and assertios here? Also what would comparison.py and mock_object.py look like? where would they live?

GitHub Copilot: Excellent questions! Let me explain these concepts in detail before we implement anything.

## **Pytest Fundamentals**

### **Pytest Fixtures (`@pytest.fixture`)**

**What they are:** Reusable setup/teardown code that provides data or objects to tests.

**Think of them as:** "Test ingredients" - you define them once, and pytest automatically provides them to any test that needs them.

```python
import pytest

# Define a fixture
@pytest.fixture
def sample_constituent():
    """Creates a constituent for testing - runs BEFORE each test."""
    constituent = Constituent(name="elastin")
    constituent.rhoR_alpha_history = [1.0, 1.1, 1.2]
    constituent.k_alpha_history = [0.01, 0.01, 0.01]
    
    yield constituent  # Provide to test
    
    # Code after yield runs AFTER test (cleanup)
    constituent.cleanup()  # If needed

# Use the fixture by adding it as a parameter
def test_get_density(sample_constituent):
    # pytest automatically creates and passes sample_constituent
    density = sample_constituent.get_rhoR_alpha(timestep=0)
    assert density == 1.0
```

**Fixture Scopes:**
```python
@pytest.fixture(scope="function")  # Default: new instance per test
def per_test_fixture():
    return "recreated each test"

@pytest.fixture(scope="class")  # Shared across test class
def per_class_fixture():
    return "recreated per test class"

@pytest.fixture(scope="module")  # Shared across entire file
def expensive_setup():
    # Create expensive objects once per file
    return load_large_dataset()

@pytest.fixture(scope="session")  # Once per entire test run
def database_connection():
    # Connect once for all tests
    db = create_connection()
    yield db
    db.close()
```

### **Pytest Marks (`@pytest.mark.*`)**

**What they are:** Labels/tags you attach to tests to categorize or modify behavior.

**Common marks:**

```python
# Mark slow tests
@pytest.mark.slow
def test_full_simulation():
    """Takes 5 minutes to run."""
    pass

# Mark tests that need external resources
@pytest.mark.requires_internet
def test_download_data():
    pass

# Mark expected failures
@pytest.mark.xfail(reason="Known bug in upstream library")
def test_with_known_issue():
    pass

# Skip tests conditionally
@pytest.mark.skipif(sys.version_info < (3, 9), reason="Requires Python 3.9+")
def test_new_feature():
    pass

# Parametrize - run same test with different inputs
@pytest.mark.parametrize("input,expected", [
    (1, 2),
    (2, 4),
    (3, 6),
])
def test_doubling(input, expected):
    assert input * 2 == expected
```

**Running tests selectively:**
```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit/

# Skip slow tests
pytest -m "not slow"

# Run ONLY slow tests
pytest -m "slow"

# Run tests matching a pattern
pytest -k "integrator"

# Verbose output
pytest -v

# Stop at first failure
pytest -x

# Run in parallel (requires pytest-xdist)
pytest -n auto
```

---

## **Code Coverage - How GitHub Knows What's Covered**

### **How Coverage Works**

1. **Instrumentation**: `pytest-cov` (built on `coverage.py`) instruments your code to track which lines execute
2. **Tracking**: During test execution, it records which lines run
3. **Reporting**: Generates reports showing covered/uncovered lines

**Visual representation:**
```python
# src/integrators.py
def integrate(self, f: List[float]) -> float:
    """Compute integral."""
    if len(f) < 2:           # ✅ Covered by test_empty_list
        raise ValueError()    # ✅ Covered by test_empty_list
    
    result = 0.0             # ✅ Covered by test_constant_function
    for i in range(len(f)):  # ✅ Covered by test_constant_function
        result += f[i]       # ✅ Covered by test_constant_function
    
    return result * self.dt  # ❌ NOT COVERED - missing return test!
```

### **GitHub Integration Flow**

```
┌─────────────────────────────────────────────┐
│ 1. Developer pushes code to GitHub         │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 2. GitHub Actions triggers workflow         │
│    (.github/workflows/tests.yml)            │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 3. Workflow runs:                           │
│    pytest --cov=src --cov-report=xml        │
│                                             │
│    Creates: coverage.xml                    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 4. Upload to Codecov (or similar):          │
│    codecov/codecov-action@v3                │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 5. Codecov analyzes coverage.xml            │
│    - Tracks line-by-line coverage           │
│    - Compares to previous commits           │
│    - Comments on PR with coverage change    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 6. Badge shows on README:                   │
│    ![Coverage](codecov.io/badge.svg)        │
│    "92% coverage"                           │
└─────────────────────────────────────────────┘
```

**The coverage.xml file** contains:
```xml
<coverage>
  <packages>
    <package name="src">
      <classes>
        <class name="integrators.py">
          <lines>
            <line number="45" hits="12"/>  <!-- Executed 12 times -->
            <line number="46" hits="0"/>   <!-- Never executed -->
          </lines>
        </class>
      </classes>
    </package>
  </packages>
</coverage>
```

**Example GitHub Actions workflow excerpt:**
```yaml
- name: Run tests with coverage
  run: |
    pytest tests/ \
      --cov=src \                    # Track coverage for src/
      --cov-report=xml \             # Generate XML for Codecov
      --cov-report=html \            # Generate HTML for local viewing
      --cov-report=term-missing      # Show missing lines in terminal

- name: Upload coverage to Codecov
  uses: codecov/codecov-action@v3
  with:
    files: ./coverage.xml
    flags: unittests
    fail_ci_if_error: true           # Fail if upload fails
```

---

## **Assertions vs Exceptions - Key Differences**

### **Exceptions (exceptions.py)**

**Purpose:** Handle runtime errors in **production code**

**When to use:**
- Invalid user input
- Physical violations during simulation
- Solver convergence failures
- Configuration errors

**Example in production code:**
```python
# src/solver.py
from exceptions import ConvergenceError, PhysicsViolationError

def solve_equilibrium(self, max_iter=100):
    """Solve for equilibrium state."""
    for i in range(max_iter):
        residual = self._compute_residual()
        
        if residual < self.tolerance:
            return self.state  # Success
    
    # Raise exception if fails to converge
    raise ConvergenceError(
        f"Failed to converge after {max_iter} iterations. "
        f"Final residual: {residual:.2e}"
    )

def compute_stress(self, stretch):
    """Compute stress from stretch."""
    if stretch < 0:
        # This violates physics - raise exception
        raise PhysicsViolationError(
            f"Negative stretch detected: {stretch}"
        )
    
    return self.material_model(stretch)
```

**Exception hierarchy in your code:**
```python
# src/exceptions.py
class GrowthRemodelingError(Exception):
    """Base class - catch all G&R errors."""
    pass

class SolverError(GrowthRemodelingError):
    """Problems during solving."""
    pass

class ConvergenceError(SolverError):
    """Specific: didn't converge."""
    pass

class PhysicsViolationError(GrowthRemodelingError):
    """Specific: violated physical constraints."""
    pass
```

### **Assertions**

**Purpose:** Verify conditions in **tests** and **development** (debug checks)

**Two types:**

#### **1. Python's `assert` statement** (for development/debugging)
```python
# In production code (optional, disabled with python -O)
def compute_survival(self, k_alpha, dt):
    assert k_alpha > 0, "Degradation rate must be positive"
    assert dt > 0, "Time step must be positive"
    return math.exp(-k_alpha * dt)

# These are REMOVED when running python with optimization flag
# So don't rely on them for critical checks!
```

#### **2. Pytest assertions** (for tests)
```python
# In tests/unit/test_kinetics.py
def test_survival_function():
    result = compute_survival(k_alpha=0.1, dt=1.0)
    
    # Basic assertion
    assert result > 0
    assert result < 1
    
    # With helpful message
    assert result == pytest.approx(0.9048, abs=1e-4), \
        f"Expected ~0.9048, got {result}"
    
    # Numpy assertions (better error messages)
    import numpy.testing as npt
    npt.assert_allclose(result, 0.9048, rtol=1e-5)
```

### **Should You Have `assertions.py` in Main Code?**

**Recommendation: NO - Keep assertions ONLY in tests/**

**Why:**
1. **Production code should use exceptions** for error handling
2. **Assertions can be disabled** (`python -O` removes them)
3. **Tests need different assertion helpers** than production code

**Better structure:**

```
src/
├── exceptions.py          # Production error handling
└── validation.py          # Production input validation (raises exceptions)

tests/
└── utils/
    └── assertions.py      # Test-specific assertion helpers
```

**Example of what SHOULD go in each:**

```python
# src/validation.py (production code)
from exceptions import GeometryError

def validate_thin_wall_assumption(radius, thickness):
    """Validate thin-wall assumption.
    
    Raises:
        GeometryError: If h/a >= 0.2
    """
    ratio = thickness / radius
    if ratio >= 0.2:
        raise GeometryError(
            f"Thin-wall assumption violated: h/a = {ratio:.3f} >= 0.2"
        )

# tests/utils/assertions.py (test helpers)
import numpy as np

def assert_thin_wall_valid(radius, thickness):
    """Test helper: assert thin-wall assumption holds."""
    ratio = thickness / radius
    assert ratio < 0.2, \
        f"Test expected thin-wall (h/a < 0.2), got {ratio:.3f}"

def assert_arrays_close_with_summary(actual, expected, rtol=1e-5):
    """Assert arrays close with detailed error message."""
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol)
    except AssertionError as e:
        # Enhanced error message
        diff = actual - expected
        max_err = np.max(np.abs(diff))
        max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
        
        raise AssertionError(
            f"Arrays not close:\n"
            f"  Max error: {max_err:.2e} at index {max_idx}\n"
            f"  Expected: {expected[max_idx]}\n"
            f"  Actual:   {actual[max_idx]}\n"
            f"  Original error: {e}"
        )
```

---

## **What Would `tests/utils/` Files Look Like?**

### **1. `tests/utils/comparisons.py`**

**Purpose:** Advanced comparison utilities for simulation outputs

```python
"""Numerical comparison utilities for tests."""
import numpy as np
from pathlib import Path
import json

def assert_simulation_matches(actual, expected, rtol=1e-5, atol=1e-8, 
                               output_mgr=None, name='output'):
    """Compare simulation outputs with detailed error reporting.
    
    Args:
        actual: Numpy array of actual results
        expected: Numpy array of expected results
        rtol: Relative tolerance
        atol: Absolute tolerance
        output_mgr: TestOutputManager for saving failures
        name: Name for output files
    
    Raises:
        AssertionError: If arrays don't match within tolerance
    """
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)
    except AssertionError as e:
        if output_mgr:
            # Save detailed comparison
            output_mgr.save_failure_comparison(actual, expected, name)
        
        # Enhanced error message
        diff = np.abs(actual - expected)
        rel_err = diff / (np.abs(expected) + 1e-16)
        
        error_msg = (
            f"\nSimulation output mismatch for '{name}':\n"
            f"  Max absolute error: {np.max(diff):.2e}\n"
            f"  Max relative error: {np.max(rel_err):.2e}\n"
            f"  Mean absolute error: {np.mean(diff):.2e}\n"
            f"  Tolerances: rtol={rtol}, atol={atol}\n"
        )
        
        if output_mgr:
            error_msg += f"  Detailed comparison saved to: {output_mgr.output_dir}\n"
        
        error_msg += f"\nOriginal error:\n{e}"
        raise AssertionError(error_msg)


def compare_history_convergence(history, expected_final, 
                                 convergence_timestep=None):
    """Verify history converges to expected final value.
    
    Args:
        history: List or array of values over time
        expected_final: Expected final converged value
        convergence_timestep: Timestep by which it should converge (optional)
    
    Returns:
        bool: True if converged appropriately
    """
    history = np.array(history)
    
    # Check final value
    final_value = history[-1]
    is_close = np.isclose(final_value, expected_final, rtol=1e-4)
    
    if not is_close:
        raise AssertionError(
            f"Final value {final_value:.4e} doesn't match "
            f"expected {expected_final:.4e}"
        )
    
    # Check convergence timing if specified
    if convergence_timestep is not None:
        values_after_convergence = history[convergence_timestep:]
        if not np.allclose(values_after_convergence, expected_final, rtol=1e-3):
            raise AssertionError(
                f"Values after timestep {convergence_timestep} not converged"
            )
    
    return True


def compare_constituent_histories(actual_constituent, expected_constituent, 
                                   check_fields=None):
    """Compare all history fields between two constituents.
    
    Args:
        actual_constituent: Constituent from test
        expected_constituent: Reference constituent
        check_fields: List of field names to check (default: all histories)
    """
    if check_fields is None:
        check_fields = [
            'rhoR_alpha_history',
            'k_alpha_history',
            'mR_alpha_history',
            'survival_history',
        ]
    
    errors = []
    
    for field in check_fields:
        actual = getattr(actual_constituent, field)
        expected = getattr(expected_constituent, field)
        
        try:
            if isinstance(actual[0], list):  # Nested lists (e.g., survival_history)
                for i, (a, e) in enumerate(zip(actual, expected)):
                    np.testing.assert_allclose(a, e, rtol=1e-5,
                        err_msg=f"{field}[{i}]")
            else:
                np.testing.assert_allclose(actual, expected, rtol=1e-5,
                    err_msg=field)
        except AssertionError as e:
            errors.append(str(e))
    
    if errors:
        raise AssertionError(
            f"Constituent histories don't match:\n" + "\n".join(errors)
        )
```

### **2. `tests/utils/mock_objects.py`**

**Purpose:** Create fake/simplified objects for isolated testing

```python
"""Mock objects for testing without full simulation setup."""
import numpy as np
from typing import List, Optional
from src.constituent import Constituent
from src.kinetics import Kinetics, ConstantDegradationRate
from src.layer import Layer

class MockLayer:
    """Minimal Layer for testing without full geometry."""
    
    def __init__(self, timesteps=10):
        self.inner_radius_history = [1.0] * timesteps
        self.thickness_history = [0.1] * timesteps
        self.lambda_z_history = [1.0] * timesteps
        self.intramural_stress_history = [100.0] * timesteps
        
        # Kinematics mock
        self.kinematics = MockKinematics()
    
    def get_inner_radius(self, timestep):
        return self.inner_radius_history[timestep]


class MockKinematics:
    """Minimal kinematics for testing stress components."""
    
    def get_component_names(self):
        return ['r', 'theta', 'z']
    
    def get_component_index(self, name):
        mapping = {'r': 0, 'theta': 1, 'z': 2}
        return mapping[name]
    
    def get_intramural_component_name(self):
        return 'theta'


def create_test_constituent(name="test_constituent", 
                            k_alpha_values=None, 
                            timesteps=10,
                            homeostatic_density=1.0):
    """Create a constituent with preset histories for testing.
    
    Args:
        name: Constituent name
        k_alpha_values: List of degradation rates (default: constant 0.01)
        timesteps: Number of timesteps to initialize
        homeostatic_density: Initial density
    
    Returns:
        Constituent: Ready for testing
    """
    from src.constituent import SingleConstituent
    
    # Create minimal constituent
    constituent = SingleConstituent(name=name)
    constituent.homeostatic_referential_density = homeostatic_density
    constituent.layer = MockLayer(timesteps)
    
    # Initialize kinetics
    if k_alpha_values is None:
        k_alpha_values = [0.01] * timesteps
    
    constituent.kinetics = Kinetics(
        degradation_model=ConstantDegradationRate(k_alpha_h=0.01),
        production_model=None,
        integrator_type='trapezoid',
        dt=1.0
    )
    
    # Populate histories
    constituent.k_alpha_history = k_alpha_values
    constituent.rhoR_alpha_history = [homeostatic_density] * timesteps
    constituent.survival_history = [[1.0]] * timesteps
    
    return constituent


def create_constant_degradation_constituent(k_const, timesteps, dt):
    """Create constituent with constant degradation for analytical testing.
    
    Useful for tests where analytical solution exists.
    
    Args:
        k_const: Constant degradation rate (1/day)
        timesteps: Number of timesteps
        dt: Time step size
    
    Returns:
        Constituent: With constant k_alpha throughout history
    """
    constituent = create_test_constituent(
        name="analytical_test",
        k_alpha_values=[k_const] * timesteps,
        timesteps=timesteps
    )
    constituent.kinetics.dt = dt
    
    return constituent


class MockIntegrator:
    """Simple mock integrator for testing without real integration."""
    
    def __init__(self, return_value=1.0):
        self.return_value = return_value
        self.calls = []  # Track calls for verification
    
    def integrate(self, f):
        """Record call and return preset value."""
        self.calls.append(f)
        return self.return_value
    
    def was_called_with(self, expected_function):
        """Check if integrator was called with expected values."""
        return any(
            np.allclose(call, expected_function, rtol=1e-6) 
            for call in self.calls
        )
```

### **3. `tests/utils/fixtures_generators.py`**

**Purpose:** Generate test data programmatically

```python
"""Generate test fixtures programmatically."""
import numpy as np
import json
from pathlib import Path

def generate_survival_test_data(output_path: Path):
    """Generate reference data for survival function tests.
    
    Creates analytical test cases where exact solution is known.
    """
    test_cases = []
    
    # Case 1: Constant degradation
    k_const = 0.1
    dt = 1.0
    timesteps = 10
    
    for s in range(timesteps):
        q_values = []
        for tau in range(s + 1):
            # Analytical: q(s,τ) = exp(-k*(s-τ)*dt)
            q_analytical = np.exp(-k_const * (s - tau) * dt)
            q_values.append(q_analytical)
        
        test_cases.append({
            'current_timestep': s,
            'k_alpha_history': [k_const] * (s + 1),
            'dt': dt,
            'expected_q': q_values
        })
    
    # Save to JSON
    with open(output_path, 'w') as f:
        json.dump(test_cases, f, indent=2)
    
    print(f"Generated {len(test_cases)} test cases in {output_path}")


def generate_integration_test_data():
    """Generate test data for numerical integrators.
    
    Returns functions with known analytical integrals.
    """
    test_functions = []
    
    # Constant function: ∫c dx = c*x
    test_functions.append({
        'name': 'constant',
        'f': lambda x: 5.0,
        'integral': lambda a, b: 5.0 * (b - a),
        'description': 'f(x) = 5'
    })
    
    # Linear: ∫x dx = x²/2
    test_functions.append({
        'name': 'linear',
        'f': lambda x: x,
        'integral': lambda a, b: (b**2 - a**2) / 2,
        'description': 'f(x) = x'
    })
    
    # Quadratic: ∫x² dx = x³/3
    test_functions.append({
        'name': 'quadratic',
        'f': lambda x: x**2,
        'integral': lambda a, b: (b**3 - a**3) / 3,
        'description': 'f(x) = x²'
    })
    
    # Sin: ∫sin(x) dx = -cos(x)
    test_functions.append({
        'name': 'sine',
        'f': lambda x: np.sin(x),
        'integral': lambda a, b: -np.cos(b) + np.cos(a),
        'description': 'f(x) = sin(x)'
    })
    
    return test_functions
```

### **Where These Files Live**

```
tests/
├── __init__.py
├── conftest.py                 # Pytest configuration (fixtures available to ALL tests)
│
├── unit/
│   ├── __init__.py
│   └── test_integrators.py     # Uses: fixtures, comparisons
│
├── integration/
│   ├── __init__.py
│   └── test_kinetics.py        # Uses: mock_objects, comparisons
│
├── e2e/
│   ├── __init__.py
│   └── test_latorre2018.py     # Uses: all utils
│
└── utils/                      # Shared testing utilities
    ├── __init__.py
    ├── assertions.py           # Test-specific assertion helpers
    ├── comparisons.py          # Numerical comparison functions
    ├── mock_objects.py         # Fake objects for isolated testing
    ├── fixtures_generators.py  # Generate test data
    └── output_manager.py       # Manage test outputs
```

**Import patterns:**
```python
# In any test file
from tests.utils.comparisons import assert_simulation_matches
from tests.utils.mock_objects import create_test_constituent, MockLayer
from tests.utils.assertions import assert_arrays_close_with_summary
```

---

## **Summary Table**

| Concept | Location | Purpose | Used By |
|---------|----------|---------|---------|
| **Exceptions** | exceptions.py | Production error handling | Production code + tests (to verify errors) |
| **Validation** | `src/validation.py` | Input validation (raises exceptions) | Production code |
| **Assertions (pytest)** | Test files | Verify test expectations | Tests only |
| **assertions.py** | `tests/utils/` | Custom test assertion helpers | Tests only |
| **comparisons.py** | `tests/utils/` | Numerical comparison utilities | Tests only |
| **mock_objects.py** | `tests/utils/` | Fake objects for isolated testing | Tests (especially unit/integration) |
| **fixtures_generators.py** | `tests/utils/` | Generate test data | Setup scripts, tests |

Does this clarify the testing concepts? Would you like me to now proceed with creating the actual testing framework structure?

User: Before we proceed, I want you to create a .md document with the contents of that explanation. In addition, add the preliminary plan for the testing framework based on the previous replies. The goal is to reference this document in the future, to give you an idea of what are the design goals, principles and key concepts for the testing framework implementation.

GitHub Copilot: Read [](file:///Users/lazaros/workspace/github/svGrowth/docs/ARCHITECTURE.md)

Now I'll create a comprehensive testing framework design document:

```markdown
# svGrowth Testing Framework Design

**Document Purpose:** Reference guide for implementing the testing framework for the svGrowth codebase. This document outlines design goals, principles, key concepts, and the complete testing architecture.

**Created:** January 3, 2026  
**Status:** Design Phase

---

## Table of Contents

1. [Testing Philosophy & Goals](#testing-philosophy--goals)
2. [Test Types Overview](#test-types-overview)
3. [Pytest Fundamentals](#pytest-fundamentals)
4. [Code Coverage Strategy](#code-coverage-strategy)
5. [Exceptions vs Assertions](#exceptions-vs-assertions)
6. [Testing Framework Architecture](#testing-framework-architecture)
7. [Test Data Management](#test-data-management)
8. [GitHub Actions CI/CD](#github-actions-cicd)
9. [Utility Modules](#utility-modules)
10. [Implementation Roadmap](#implementation-roadmap)

---

## Testing Philosophy & Goals

### Design Principles

1. **Comprehensive Coverage:** Test at multiple granularities (unit → integration → end-to-end)
2. **Isolation:** Unit tests should not depend on full simulation setup
3. **Reproducibility:** All tests must produce deterministic results
4. **Fast Feedback:** Unit tests run in milliseconds, integration in seconds
5. **Clear Failure Messages:** When tests fail, immediately understand why
6. **Maintainability:** Tests should be easy to update as code evolves

### Success Criteria

- **Coverage Target:** ≥85% code coverage across all modules
- **Speed:** Unit test suite completes in <10 seconds
- **Reliability:** Zero flaky tests (non-deterministic failures)
- **Documentation:** Every test clearly documents what it verifies

---

## Test Types Overview

### 1. Unit Tests
**Purpose:** Verify that smallest code units (functions, methods) work correctly in isolation

**Characteristics:**
- **Speed:** Fast (milliseconds per test)
- **Dependencies:** Minimal - use mocks for external dependencies
- **Scope:** Single function or method
- **Location:** `tests/unit/`

**Examples for svGrowth:**
- Integration methods in `integrators.py` (trapezoidal, Simpson's rule)
- Tensor operations in `tensor_operations.py`
- Individual stress computation functions in `constitutive_laws.py`
- Helper functions and utilities

**Example:**
```python
def test_trapezoid_integrator_constant_function():
    """Integral of f(x)=c should equal c*(b-a)."""
    dt = 0.1
    f = [5.0] * 10  # Constant function
    integrator = TrapezoidIntegrator(dt=dt, start=0, stop=9)
    
    result = integrator.integrate(f)
    expected = 5.0 * 0.9  # c * (9-0) * dt
    
    assert abs(result - expected) < 1e-10
```

### 2. Integration Tests
**Purpose:** Verify that different modules work together correctly

**Characteristics:**
- **Speed:** Medium (seconds per test)
- **Dependencies:** Real classes working together, but not full simulation
- **Scope:** Multiple interacting components
- **Location:** `tests/integration/`

**Examples for svGrowth:**
- Survival function computation (kinetics + integrators)
- `compute_rhoR_alpha()` workflow (constituent + kinetics + layer data)
- Stress computation workflow (mechanics + constitutive laws + kinematics)
- Solver iteration logic (solver + mechanics + layer)

**Example:**
```python
def test_survival_function_with_constant_degradation():
    """Test survival computation with analytical solution."""
    k_const = 0.1  # Constant degradation rate
    constituent = create_test_constituent(k_alpha=[k_const] * 10)
    
    context = ConstituentKineticsContext(constituent)
    q_history = constituent.kinetics.compute_survival_history(context, timestep=5)
    
    # Compare to analytical: q(s,τ) = exp(-k*(s-τ))
    for tau, q_value in enumerate(q_history):
        expected = np.exp(-k_const * (5 - tau))
        assert abs(q_value - expected) < 1e-6
```

### 3. End-to-End (E2E) Tests
**Purpose:** Verify entire system produces correct outputs for known inputs

**Characteristics:**
- **Speed:** Slow (minutes per test)
- **Dependencies:** Full simulation stack
- **Scope:** Complete workflows from input YAML to output files
- **Location:** `tests/e2e/`

**Examples for svGrowth:**
- Full Latorre2018 simulation matching published results
- Pressure perturbation scenarios
- Axial stretch perturbation validation
- Convergence to homeostasis verification

**Example:**
```python
def test_latorre2018_full_validation():
    """Verify full simulation matches Latorre et al. 2018 results."""
    sim = Simulation.from_yaml('tests/fixtures/inputs/latorre2018_test.yaml')
    results = sim.run()
    
    # Load validated reference data
    expected = np.load('tests/fixtures/expected_outputs/latorre2018/geometry.npy')
    
    assert_simulation_matches(results.geometry, expected, rtol=1e-4)
```

---

## Pytest Fundamentals

### Pytest Fixtures

**What are they?** Reusable setup/teardown code that provides data or objects to tests.

**Think of them as:** "Test ingredients" - define once, pytest automatically provides to any test that needs them.

**Basic Example:**
```python
import pytest

@pytest.fixture
def sample_constituent():
    """Creates a constituent for testing - runs BEFORE each test."""
    constituent = Constituent(name="elastin")
    constituent.rhoR_alpha_history = [1.0, 1.1, 1.2]
    constituent.k_alpha_history = [0.01, 0.01, 0.01]
    
    yield constituent  # Provide to test
    
    # Code after yield runs AFTER test (cleanup)
    constituent.cleanup()  # If needed

# Use fixture by adding it as parameter
def test_get_density(sample_constituent):
    # pytest automatically creates and passes sample_constituent
    density = sample_constituent.get_rhoR_alpha(timestep=0)
    assert density == 1.0
```

**Fixture Scopes:**
```python
@pytest.fixture(scope="function")  # Default: new instance per test
def per_test_fixture():
    return "recreated each test"

@pytest.fixture(scope="class")  # Shared across test class
def per_class_fixture():
    return "recreated per test class"

@pytest.fixture(scope="module")  # Shared across entire file
def expensive_setup():
    return load_large_dataset()

@pytest.fixture(scope="session")  # Once per entire test run
def database_connection():
    db = create_connection()
    yield db
    db.close()
```

### Pytest Marks

**What are they?** Labels/tags to categorize or modify test behavior.

**Common Marks:**
```python
# Mark slow tests
@pytest.mark.slow
def test_full_simulation():
    """Takes 5 minutes to run."""
    pass

# Mark tests needing external resources
@pytest.mark.requires_internet
def test_download_data():
    pass

# Mark expected failures
@pytest.mark.xfail(reason="Known bug in upstream library")
def test_with_known_issue():
    pass

# Skip tests conditionally
@pytest.mark.skipif(sys.version_info < (3, 9), reason="Requires Python 3.9+")
def test_new_feature():
    pass

# Parametrize - run same test with different inputs
@pytest.mark.parametrize("input,expected", [
    (1, 2),
    (2, 4),
    (3, 6),
])
def test_doubling(input, expected):
    assert input * 2 == expected
```

**Running Tests Selectively:**
```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit/

# Skip slow tests
pytest -m "not slow"

# Run ONLY slow tests
pytest -m "slow"

# Run tests matching pattern
pytest -k "integrator"

# Verbose output
pytest -v

# Stop at first failure
pytest -x

# Run in parallel (requires pytest-xdist)
pytest -n auto
```

---

## Code Coverage Strategy

### How Coverage Works

1. **Instrumentation:** `pytest-cov` (built on `coverage.py`) instruments code to track line execution
2. **Tracking:** During test execution, records which lines run
3. **Reporting:** Generates reports showing covered/uncovered lines

**Visual Example:**
```python
# src/integrators.py
def integrate(self, f: List[float]) -> float:
    """Compute integral."""
    if len(f) < 2:           # ✅ Covered by test_empty_list
        raise ValueError()    # ✅ Covered by test_empty_list
    
    result = 0.0             # ✅ Covered by test_constant_function
    for i in range(len(f)):  # ✅ Covered by test_constant_function
        result += f[i]       # ✅ Covered by test_constant_function
    
    return result * self.dt  # ❌ NOT COVERED - missing return test!
```

### GitHub Integration Flow

```
┌─────────────────────────────────────────────┐
│ 1. Developer pushes code to GitHub         │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 2. GitHub Actions triggers workflow         │
│    (.github/workflows/tests.yml)            │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 3. Workflow runs:                           │
│    pytest --cov=src --cov-report=xml        │
│                                             │
│    Creates: coverage.xml                    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 4. Upload to Codecov:                       │
│    codecov/codecov-action@v3                │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 5. Codecov analyzes coverage.xml            │
│    - Tracks line-by-line coverage           │
│    - Compares to previous commits           │
│    - Comments on PR with coverage change    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 6. Badge shows on README:                   │
│    ![Coverage](codecov.io/badge.svg)        │
│    "92% coverage"                           │
└─────────────────────────────────────────────┘
```

**Coverage XML Structure:**
```xml
<coverage>
  <packages>
    <package name="src">
      <classes>
        <class name="integrators.py">
          <lines>
            <line number="45" hits="12"/>  <!-- Executed 12 times -->
            <line number="46" hits="0"/>   <!-- Never executed -->
          </lines>
        </class>
      </classes>
    </package>
  </packages>
</coverage>
```

### Running Coverage Locally

```bash
# Run tests with coverage
pytest tests/ --cov=src --cov-report=html

# View results
open htmlcov/index.html

# Terminal output with missing lines
pytest tests/ --cov=src --cov-report=term-missing

# Example output:
# Name                    Stmts   Miss  Cover   Missing
# -----------------------------------------------------
# src/integrators.py         45      3    93%   102-104
# src/kinetics.py           120      8    93%   45, 67-73
```

---

## Exceptions vs Assertions

### Exceptions (Production Code)

**Purpose:** Handle runtime errors in production code

**Location:** exceptions.py

**When to Use:**
- Invalid user input
- Physical violations during simulation
- Solver convergence failures
- Configuration errors

**Example:**
```python
# src/solver.py
from exceptions import ConvergenceError, PhysicsViolationError

def solve_equilibrium(self, max_iter=100):
    """Solve for equilibrium state."""
    for i in range(max_iter):
        residual = self._compute_residual()
        
        if residual < self.tolerance:
            return self.state  # Success
    
    # Raise exception if fails to converge
    raise ConvergenceError(
        f"Failed to converge after {max_iter} iterations. "
        f"Final residual: {residual:.2e}"
    )

def compute_stress(self, stretch):
    """Compute stress from stretch."""
    if stretch < 0:
        raise PhysicsViolationError(
            f"Negative stretch detected: {stretch}"
        )
    return self.material_model(stretch)
```

**Exception Hierarchy:**
```python
# src/exceptions.py
class GrowthRemodelingError(Exception):
    """Base class - catch all G&R errors."""
    pass

class SolverError(GrowthRemodelingError):
    """Problems during solving."""
    pass

class ConvergenceError(SolverError):
    """Specific: didn't converge."""
    pass

class PhysicsViolationError(GrowthRemodelingError):
    """Specific: violated physical constraints."""
    pass
```

### Assertions (Test Code)

**Purpose:** Verify conditions in tests

**Location:** Test files (tests)

**Two Types:**

#### 1. Python's `assert` Statement (Development/Debug)
```python
# In production code (optional, disabled with python -O)
def compute_survival(self, k_alpha, dt):
    assert k_alpha > 0, "Degradation rate must be positive"
    assert dt > 0, "Time step must be positive"
    return math.exp(-k_alpha * dt)

# These are REMOVED when running python with optimization flag!
# Don't rely on them for critical checks in production!
```

#### 2. Pytest Assertions (Testing)
```python
# In tests/unit/test_kinetics.py
def test_survival_function():
    result = compute_survival(k_alpha=0.1, dt=1.0)
    
    # Basic assertion
    assert result > 0
    assert result < 1
    
    # With helpful message
    assert result == pytest.approx(0.9048, abs=1e-4), \
        f"Expected ~0.9048, got {result}"
    
    # Numpy assertions (better error messages)
    import numpy.testing as npt
    npt.assert_allclose(result, 0.9048, rtol=1e-5)
```

### Design Decision: No `assertions.py` in Main Code

**Recommendation:** Keep assertions ONLY in `tests/utils/assertions.py`

**Rationale:**
1. Production code should use exceptions for error handling
2. Assertions can be disabled (`python -O` removes them)
3. Tests need different assertion helpers than production code

**File Organization:**

| File | Purpose | Used By |
|------|---------|---------|
| exceptions.py | Production error handling | Production code + tests |
| `src/validation.py` | Input validation (raises exceptions) | Production code |
| `tests/utils/assertions.py` | Test-specific assertion helpers | Tests only |

**Examples:**

```python
# src/validation.py (production code)
from exceptions import GeometryError

def validate_thin_wall_assumption(radius, thickness):
    """Validate thin-wall assumption.
    
    Raises:
        GeometryError: If h/a >= 0.2
    """
    ratio = thickness / radius
    if ratio >= 0.2:
        raise GeometryError(
            f"Thin-wall assumption violated: h/a = {ratio:.3f} >= 0.2"
        )

# tests/utils/assertions.py (test helpers)
import numpy as np

def assert_thin_wall_valid(radius, thickness):
    """Test helper: assert thin-wall assumption holds."""
    ratio = thickness / radius
    assert ratio < 0.2, \
        f"Test expected thin-wall (h/a < 0.2), got {ratio:.3f}"

def assert_arrays_close_with_summary(actual, expected, rtol=1e-5):
    """Assert arrays close with detailed error message."""
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol)
    except AssertionError as e:
        diff = actual - expected
        max_err = np.max(np.abs(diff))
        max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
        
        raise AssertionError(
            f"Arrays not close:\n"
            f"  Max error: {max_err:.2e} at index {max_idx}\n"
            f"  Expected: {expected[max_idx]}\n"
            f"  Actual:   {actual[max_idx]}\n"
            f"  Original error: {e}"
        )
```

---

## Testing Framework Architecture

### Directory Structure

```
tests/
├── __init__.py
├── conftest.py                      # Pytest fixtures and configuration
├── test_config.yaml                 # Testing-specific configuration
├── requirements_test.txt            # Test dependencies
│
├── unit/                            # Unit tests
│   ├── __init__.py
│   ├── test_integrators.py         # Test integration methods
│   ├── test_tensor_operations.py   # Test tensor utilities
│   ├── test_constitutive_laws.py   # Test stress-strain laws
│   ├── test_deformation_kinematics.py
│   └── test_helper_functions.py    # Test utility functions
│
├── integration/                     # Integration tests
│   ├── __init__.py
│   ├── test_kinetics.py            # Test survival/production computations
│   ├── test_constituent.py         # Test compute_rhoR_alpha workflows
│   ├── test_mechanics.py           # Test stress computation workflows
│   └── test_solver_components.py   # Test solver iterations
│
├── e2e/                            # End-to-end tests
│   ├── __init__.py
│   ├── test_latorre2018.py        # Full Latorre2018 validation
│   └── test_perturbations.py      # Test perturbation scenarios
│
├── fixtures/                       # Test data and configurations
│   ├── inputs/                    # Input YAML files for tests
│   │   ├── latorre2018_test.yaml
│   │   ├── simple_cylinder.yaml
│   │   ├── unit_test_configs/
│   │   │   ├── simple_integrator.yaml
│   │   │   └── single_constituent.yaml
│   │   └── integration_configs/
│   │       ├── kinetics_only.yaml
│   │       └── two_constituent.yaml
│   ├── expected_outputs/          # Known-good reference outputs
│   │   ├── latorre2018_validated/
│   │   │   ├── metadata.json
│   │   │   ├── geometry.npy
│   │   │   ├── rhoR_alpha.npy
│   │   │   └── stress.npy
│   │   └── integration_test_data/
│   │       ├── survival_function_expected.json
│   │       └── rhoR_alpha_computation_expected.json
│   └── test_data/                 # Synthetic data for unit tests
│       ├── integrator_test_data.json
│       └── survival_function_data.json
│
├── utils/                         # Testing utilities
│   ├── __init__.py
│   ├── assertions.py             # Custom assertion helpers
│   ├── comparisons.py            # Numerical comparison utilities
│   ├── fixtures_generators.py    # Generate test data
│   ├── mock_objects.py           # Mock classes for testing
│   └── output_manager.py         # Manage test outputs
│
└── outputs/                       # Test run outputs (gitignored)
    ├── .gitkeep
    ├── unit/
    ├── integration/
    └── e2e/
```

### Key Files

#### `tests/conftest.py`
Central pytest configuration with shared fixtures available to all tests.

```python
"""Pytest configuration and shared fixtures."""
import pytest
from pathlib import Path

@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "fixtures" / "test_data"

@pytest.fixture(scope="session")
def expected_outputs_dir():
    """Path to expected outputs directory."""
    return Path(__file__).parent / "fixtures" / "expected_outputs"

@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for test outputs."""
    return tmp_path / "test_outputs"
```

#### `tests/requirements_test.txt`
Testing-specific dependencies (not needed for production).

```
pytest>=7.4.0
pytest-cov>=4.1.0
pytest-xdist>=3.3.0      # Parallel test execution
pytest-timeout>=2.1.0     # Timeout protection
hypothesis>=6.82.0        # Property-based testing
numpy-testing>=0.1.0      # Numerical assertions
```

#### `tests/.gitignore`
Prevent committing test outputs.

```
# Test outputs
outputs/
*.pyc
__pycache__/
.pytest_cache/
htmlcov/
.coverage
coverage.xml
```

---

## Test Data Management

### Input Data Structure

```
tests/fixtures/inputs/
├── unit_test_configs/          # Minimal configs for unit tests
│   ├── simple_integrator.yaml
│   └── single_constituent.yaml
├── integration_configs/        # Partial simulation configs
│   ├── kinetics_only.yaml
│   └── two_constituent.yaml
└── e2e_configs/               # Full simulation configs
    ├── latorre2018_test.yaml
    └── validation_cases/
        ├── pressure_perturbation.yaml
        └── axial_stretch_perturbation.yaml
```

### Expected Output Structure

```
tests/fixtures/expected_outputs/
├── latorre2018_validated/
│   ├── metadata.json          # Simulation parameters, tolerances
│   ├── geometry.npy           # [timesteps × 3] (a, h, λz)
│   ├── rhoR_alpha.npy         # [timesteps × n_constituents]
│   └── stress.npy             # [timesteps × 3] (σr, σθ, σz)
└── integration_test_data/
    ├── survival_function_expected.json
    └── rhoR_alpha_computation_expected.json
```

**`metadata.json` Structure:**
```json
{
  "description": "Latorre et al. 2018 validation case",
  "source": "Published paper results",
  "timesteps": 100,
  "dt": 1.0,
  "tolerances": {
    "geometry": {"rtol": 1e-4, "atol": 1e-8},
    "stress": {"rtol": 1e-3, "atol": 1e-6}
  },
  "generated_date": "2026-01-03",
  "generated_by": "validation_plots.py"
}
```

### Test Output Structure (Gitignored)

```
tests/outputs/                  # Gitignored - runtime outputs
├── test_run_<timestamp>/
│   ├── logs/
│   │   ├── test_integration.log
│   │   └── test_e2e.log
│   ├── failures/              # Only written on failure
│   │   ├── test_latorre2018_geometry_diff.png
│   │   ├── test_latorre2018_actual_vs_expected.json
│   │   └── test_survival_function_metrics.json
│   └── coverage/
│       └── htmlcov/
```

### Data Management Principles

1. **Version Control:**
   - ✅ Commit: Input YAMLs, expected outputs, test data
   - ❌ Don't commit: Test run outputs, logs, temporary files

2. **Data Formats:**
   - Numpy arrays: `.npy` (binary, fast, exact)
   - Metadata: `.json` (human-readable, version-controllable)
   - Logs: `.log` (text, easy to grep)

3. **Storage Organization:**
   - Group by test type (unit/integration/e2e)
   - Use timestamps for run-specific outputs
   - Save detailed outputs only on failure

---

## GitHub Actions CI/CD

### Workflow File: `.github/workflows/tests.yml`

```yaml
name: svGrowth Tests

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.9', '3.10', '3.11']
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install -r tests/requirements_test.txt
    
    - name: Run unit tests
      run: |
        pytest tests/unit/ -v --cov=src --cov-report=xml --cov-report=term-missing
    
    - name: Run integration tests
      run: |
        pytest tests/integration/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Run E2E tests
      run: |
        pytest tests/e2e/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        files: ./coverage.xml
        flags: unittests
        name: codecov-umbrella
        fail_ci_if_error: true

  test-fast:
    # Quick smoke test for rapid feedback
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    
    - uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    
    - name: Install dependencies
      run: |
        pip install -r requirements.txt -r tests/requirements_test.txt
    
    - name: Run unit tests only (fast)
      run: |
        pytest tests/unit/ -v --maxfail=1  # Fail fast
```

### Workflow Triggers

- **Push to main/develop:** Run full test suite
- **Pull requests:** Run full test suite + coverage comparison
- **Scheduled (optional):** Daily full test run to catch environment issues

### Integration with Codecov

1. Sign up at [codecov.io](https://codecov.io)
2. Add repository to Codecov
3. Get upload token (stored in GitHub Secrets)
4. Codecov automatically:
   - Comments on PRs with coverage changes
   - Shows line-by-line coverage in GitHub UI
   - Generates coverage badges

**Badge for README:**
```markdown
[![codecov](https://codecov.io/gh/username/svGrowth/branch/main/graph/badge.svg)](https://codecov.io/gh/username/svGrowth)
```

---

## Utility Modules

### `tests/utils/comparisons.py`

**Purpose:** Advanced numerical comparison utilities for simulation outputs.

```python
"""Numerical comparison utilities for tests."""
import numpy as np
from pathlib import Path
import json

def assert_simulation_matches(actual, expected, rtol=1e-5, atol=1e-8, 
                               output_mgr=None, name='output'):
    """Compare simulation outputs with detailed error reporting.
    
    Args:
        actual: Numpy array of actual results
        expected: Numpy array of expected results
        rtol: Relative tolerance
        atol: Absolute tolerance
        output_mgr: TestOutputManager for saving failures
        name: Name for output files
    
    Raises:
        AssertionError: If arrays don't match within tolerance
    """
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)
    except AssertionError as e:
        if output_mgr:
            output_mgr.save_failure_comparison(actual, expected, name)
        
        # Enhanced error message
        diff = np.abs(actual - expected)
        rel_err = diff / (np.abs(expected) + 1e-16)
        
        error_msg = (
            f"\nSimulation output mismatch for '{name}':\n"
            f"  Max absolute error: {np.max(diff):.2e}\n"
            f"  Max relative error: {np.max(rel_err):.2e}\n"
            f"  Mean absolute error: {np.mean(diff):.2e}\n"
            f"  Tolerances: rtol={rtol}, atol={atol}\n"
        )
        
        if output_mgr:
            error_msg += f"  Detailed comparison saved to: {output_mgr.output_dir}\n"
        
        error_msg += f"\nOriginal error:\n{e}"
        raise AssertionError(error_msg)


def compare_history_convergence(history, expected_final, 
                                 convergence_timestep=None):
    """Verify history converges to expected final value.
    
    Args:
        history: List or array of values over time
        expected_final: Expected final converged value
        convergence_timestep: Timestep by which it should converge (optional)
    
    Returns:
        bool: True if converged appropriately
    """
    history = np.array(history)
    
    # Check final value
    final_value = history[-1]
    is_close = np.isclose(final_value, expected_final, rtol=1e-4)
    
    if not is_close:
        raise AssertionError(
            f"Final value {final_value:.4e} doesn't match "
            f"expected {expected_final:.4e}"
        )
    
    # Check convergence timing if specified
    if convergence_timestep is not None:
        values_after_convergence = history[convergence_timestep:]
        if not np.allclose(values_after_convergence, expected_final, rtol=1e-3):
            raise AssertionError(
                f"Values after timestep {convergence_timestep} not converged"
            )
    
    return True


def compare_constituent_histories(actual_constituent, expected_constituent, 
                                   check_fields=None):
    """Compare all history fields between two constituents.
    
    Args:
        actual_constituent: Constituent from test
        expected_constituent: Reference constituent
        check_fields: List of field names to check (default: all histories)
    """
    if check_fields is None:
        check_fields = [
            'rhoR_alpha_history',
            'k_alpha_history',
            'mR_alpha_history',
            'survival_history',
        ]
    
    errors = []
    
    for field in check_fields:
        actual = getattr(actual_constituent, field)
        expected = getattr(expected_constituent, field)
        
        try:
            if isinstance(actual[0], list):  # Nested lists (e.g., survival_history)
                for i, (a, e) in enumerate(zip(actual, expected)):
                    np.testing.assert_allclose(a, e, rtol=1e-5,
                        err_msg=f"{field}[{i}]")
            else:
                np.testing.assert_allclose(actual, expected, rtol=1e-5,
                    err_msg=field)
        except AssertionError as e:
            errors.append(str(e))
    
    if errors:
        raise AssertionError(
            f"Constituent histories don't match:\n" + "\n".join(errors)
        )
```

### `tests/utils/mock_objects.py`

**Purpose:** Create fake/simplified objects for isolated testing without full simulation setup.

```python
"""Mock objects for testing without full simulation setup."""
import numpy as np
from typing import List, Optional
from src.constituent import Constituent
from src.kinetics import Kinetics, ConstantDegradationRate
from src.layer import Layer

class MockLayer:
    """Minimal Layer for testing without full geometry."""
    
    def __init__(self, timesteps=10):
        self.inner_radius_history = [1.0] * timesteps
        self.thickness_history = [0.1] * timesteps
        self.lambda_z_history = [1.0] * timesteps
        self.intramural_stress_history = [100.0] * timesteps
        
        # Kinematics mock
        self.kinematics = MockKinematics()
    
    def get_inner_radius(self, timestep):
        return self.inner_radius_history[timestep]


class MockKinematics:
    """Minimal kinematics for testing stress components."""
    
    def get_component_names(self):
        return ['r', 'theta', 'z']
    
    def get_component_index(self, name):
        mapping = {'r': 0, 'theta': 1, 'z': 2}
        return mapping[name]
    
    def get_intramural_component_name(self):
        return 'theta'


def create_test_constituent(name="test_constituent", 
                            k_alpha_values=None, 
                            timesteps=10,
                            homeostatic_density=1.0):
    """Create a constituent with preset histories for testing.
    
    Args:
        name: Constituent name
        k_alpha_values: List of degradation rates (default: constant 0.01)
        timesteps: Number of timesteps to initialize
        homeostatic_density: Initial density
    
    Returns:
        Constituent: Ready for testing
    """
    from src.constituent import SingleConstituent
    
    # Create minimal constituent
    constituent = SingleConstituent(name=name)
    constituent.homeostatic_referential_density = homeostatic_density
    constituent.layer = MockLayer(timesteps)
    
    # Initialize kinetics
    if k_alpha_values is None:
        k_alpha_values = [0.01] * timesteps
    
    constituent.kinetics = Kinetics(
        degradation_model=ConstantDegradationRate(k_alpha_h=0.01),
        production_model=None,
        integrator_type='trapezoid',
        dt=1.0
    )
    
    # Populate histories
    constituent.k_alpha_history = k_alpha_values
    constituent.rhoR_alpha_history = [homeostatic_density] * timesteps
    constituent.survival_history = [[1.0]] * timesteps
    
    return constituent


def create_constant_degradation_constituent(k_const, timesteps, dt):
    """Create constituent with constant degradation for analytical testing.
    
    Useful for tests where analytical solution exists.
    
    Args:
        k_const: Constant degradation rate (1/day)
        timesteps: Number of timesteps
        dt: Time step size
    
    Returns:
        Constituent: With constant k_alpha throughout history
    """
    constituent = create_test_constituent(
        name="analytical_test",
        k_alpha_values=[k_const] * timesteps,
        timesteps=timesteps
    )
    constituent.kinetics.dt = dt
    
    return constituent


class MockIntegrator:
    """Simple mock integrator for testing without real integration."""
    
    def __init__(self, return_value=1.0):
        self.return_value = return_value
        self.calls = []  # Track calls for verification
    
    def integrate(self, f):
        """Record call and return preset value."""
        self.calls.append(f)
        return self.return_value
    
    def was_called_with(self, expected_function):
        """Check if integrator was called with expected values."""
        return any(
            np.allclose(call, expected_function, rtol=1e-6) 
            for call in self.calls
        )
```

### `tests/utils/output_manager.py`

**Purpose:** Manage test outputs in a clean, modular way.

```python
"""Manage test outputs in a clean, modular way."""
from pathlib import Path
from datetime import datetime
import json
import numpy as np
import shutil

class TestOutputManager:
    """Handles test output storage and cleanup."""
    
    def __init__(self, test_name: str, test_type: str):
        """
        Args:
            test_name: Name of test (e.g., 'test_latorre2018')
            test_type: 'unit', 'integration', or 'e2e'
        """
        self.test_name = test_name
        self.test_type = test_type
        
        # Create timestamped output directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(__file__).parent.parent / "outputs" / test_type / f"{test_name}_{timestamp}"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def save_array(self, name: str, array: np.ndarray):
        """Save numpy array."""
        np.save(self.output_dir / f"{name}.npy", array)
    
    def save_dict(self, name: str, data: dict):
        """Save dictionary as JSON."""
        with open(self.output_dir / f"{name}.json", 'w') as f:
            json.dump(data, f, indent=2)
    
    def save_failure_comparison(self, actual, expected, name: str):
        """Save detailed comparison for failed test."""
        failure_dir = self.output_dir / "failures"
        failure_dir.mkdir(exist_ok=True)
        
        # Save both arrays
        np.save(failure_dir / f"{name}_actual.npy", actual)
        np.save(failure_dir / f"{name}_expected.npy", expected)
        
        # Save difference
        diff = actual - expected
        np.save(failure_dir / f"{name}_diff.npy", diff)
        
        # Save metadata
        metadata = {
            "max_absolute_error": float(np.max(np.abs(diff))),
            "mean_absolute_error": float(np.mean(np.abs(diff))),
            "relative_error": float(np.max(np.abs(diff / (expected + 1e-16)))),
        }
        self.save_dict(f"failures/{name}_metrics", metadata)
    
    def cleanup_on_success(self):
        """Remove output directory if test passed."""
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
```

### `tests/utils/fixtures_generators.py`

**Purpose:** Generate test data programmatically.

```python
"""Generate test fixtures programmatically."""
import numpy as np
import json
from pathlib import Path

def generate_survival_test_data(output_path: Path):
    """Generate reference data for survival function tests.
    
    Creates analytical test cases where exact solution is known.
    """
    test_cases = []
    
    # Case 1: Constant degradation
    k_const = 0.1
    dt = 1.0
    timesteps = 10
    
    for s in range(timesteps):
        q_values = []
        for tau in range(s + 1):
            # Analytical: q(s,τ) = exp(-k*(s-τ)*dt)
            q_analytical = np.exp(-k_const * (s - tau) * dt)
            q_values.append(q_analytical)
        
        test_cases.append({
            'current_timestep': s,
            'k_alpha_history': [k_const] * (s + 1),
            'dt': dt,
            'expected_q': q_values
        })
    
    # Save to JSON
    with open(output_path, 'w') as f:
        json.dump(test_cases, f, indent=2)
    
    print(f"Generated {len(test_cases)} test cases in {output_path}")


def generate_integration_test_data():
    """Generate test data for numerical integrators.
    
    Returns functions with known analytical integrals.
    """
    test_functions = []
    
    # Constant function: ∫c dx = c*x
    test_functions.append({
        'name': 'constant',
        'f': lambda x: 5.0,
        'integral': lambda a, b: 5.0 * (b - a),
        'description': 'f(x) = 5'
    })
    
    # Linear: ∫x dx = x²/2
    test_functions.append({
        'name': 'linear',
        'f': lambda x: x,
        'integral': lambda a, b: (b**2 - a**2) / 2,
        'description': 'f(x) = x'
    })
    
    # Quadratic: ∫x² dx = x³/3
    test_functions.append({
        'name': 'quadratic',
        'f': lambda x: x**2,
        'integral': lambda a, b: (b**3 - a**3) / 3,
        'description': 'f(x) = x²'
    })
    
    # Sin: ∫sin(x) dx = -cos(x)
    test_functions.append({
        'name': 'sine',
        'f': lambda x: np.sin(x),
        'integral': lambda a, b: -np.cos(b) + np.cos(a),
        'description': 'f(x) = sin(x)'
    })
    
    return test_functions
```

### `tests/utils/assertions.py`

**Purpose:** Custom test-specific assertion helpers.

```python
"""Custom assertion helpers for tests."""
import numpy as np

def assert_thin_wall_valid(radius, thickness):
    """Test helper: assert thin-wall assumption holds."""
    ratio = thickness / radius
    assert ratio < 0.2, \
        f"Test expected thin-wall (h/a < 0.2), got {ratio:.3f}"


def assert_arrays_close_with_summary(actual, expected, rtol=1e-5):
    """Assert arrays close with detailed error message."""
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol)
    except AssertionError as e:
        diff = actual - expected
        max_err = np.max(np.abs(diff))
        max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
        
        raise AssertionError(
            f"Arrays not close:\n"
            f"  Max error: {max_err:.2e} at index {max_idx}\n"
            f"  Expected: {expected[max_idx]}\n"
            f"  Actual:   {actual[max_idx]}\n"
            f"  Original error: {e}"
        )


def assert_physical_constraints(geometry_dict):
    """Assert physical constraints on geometry."""
    assert geometry_dict['radius'] > 0, "Radius must be positive"
    assert geometry_dict['thickness'] > 0, "Thickness must be positive"
    assert geometry_dict['lambda_z'] > 0, "Axial stretch must be positive"
```

---

## Implementation Roadmap

### Phase 1: Foundation (Week 1)

**Goal:** Set up testing infrastructure

- [ ] Create directory structure (tests, subdirectories)
- [ ] Create `conftest.py` with basic fixtures
- [ ] Create `requirements_test.txt`
- [ ] Set up .gitignore for test outputs
- [ ] Create utility modules (stubs):
  - `output_manager.py`
  - `mock_objects.py`
  - `comparisons.py`
  - `assertions.py`
  - `fixtures_generators.py`

### Phase 2: Unit Tests (Week 2)

**Goal:** Test individual functions in isolation

- [ ] `test_integrators.py`:
  - Naive trapezoidal integration
  - Backward trapezoidal integration
  - Naive Simpson's rule
  - Backward Simpson's rule
  - Convergence rate tests
- [ ] `test_tensor_operations.py`:
  - Matrix operations
  - Tensor transformations
- [ ] `test_constitutive_laws.py`:
  - Fung exponential model
  - Neo-Hookean model
  - Stress computation functions

**Target:** 50+ unit tests, >80% coverage of tested modules

### Phase 3: Integration Tests (Week 3)

**Goal:** Test component interactions

- [ ] `test_kinetics.py`:
  - Survival function computation
  - Production function computation
  - Integration with degradation models
- [ ] `test_constituent.py`:
  - `compute_rhoR_alpha()` workflow
  - History management
  - Context usage
- [ ] `test_mechanics.py`:
  - Stress computation workflow
  - Equilibrium solving
  - Kinematics integration

**Target:** 30+ integration tests

### Phase 4: End-to-End Tests (Week 4)

**Goal:** Validate complete workflows

- [ ] Generate expected outputs:
  - Run validated Latorre2018 simulation
  - Save geometry, stress, density histories
  - Document tolerances and metadata
- [ ] `test_latorre2018.py`:
  - Full simulation validation
  - Convergence to homeostasis
  - Sensitivity to solver tolerance
- [ ] `test_perturbations.py`:
  - Pressure perturbation
  - Axial stretch perturbation

**Target:** 10+ E2E tests, full simulation coverage

### Phase 5: CI/CD Integration (Week 5)

**Goal:** Automate testing

- [ ] Create `.github/workflows/tests.yml`
- [ ] Set up Codecov integration
- [ ] Add coverage badge to README
- [ ] Configure test reports in GitHub Actions
- [ ] Set up branch protection (require tests to pass)

### Phase 6: Documentation & Refinement (Week 6)

**Goal:** Finalize and document

- [ ] Add docstrings to all test files
- [ ] Create `tests/README.md` with:
  - How to run tests
  - How to add new tests
  - Interpreting test failures
- [ ] Review and refactor test utilities
- [ ] Add property-based tests (hypothesis) for edge cases
- [ ] Performance benchmarking tests (optional)

---

## Best Practices Summary

### Writing Tests

1. **One concept per test:** Each test should verify one specific behavior
2. **Descriptive names:** `test_survival_decreases_with_cohort_age()` not `test_survival()`
3. **Arrange-Act-Assert pattern:**
   ```python
   def test_example():
       # Arrange: Set up test data
       constituent = create_test_constituent()
       
       # Act: Execute the code under test
       result = constituent.compute_rhoR_alpha(timestep=5)
       
       # Assert: Verify the result
       assert result > 0
   ```
4. **Use fixtures for setup:** Avoid duplication
5. **Test edge cases:** Empty lists, zero values, boundary conditions
6. **Test error handling:** Verify exceptions are raised appropriately

### Maintaining Tests

1. **Keep tests independent:** No test should depend on another test's state
2. **Fast feedback:** Unit tests should run quickly
3. **Clear failure messages:** Use custom messages in assertions
4. **Regular updates:** Update tests when code changes
5. **Review coverage:** Aim for high coverage but prioritize meaningful tests

### Testing Numerical Code

1. **Use appropriate tolerances:** `np.testing.assert_allclose(rtol=1e-5)`
2. **Test against analytical solutions** when possible
3. **Test convergence properties:** Verify error decreases with finer discretization
4. **Test physical constraints:** Non-negativity, conservation laws
5. **Document expected behavior:** Why should this value be what it is?

---

## References & Resources

### Pytest Documentation
- [Pytest Official Docs](https://docs.pytest.org/)
- [Pytest Fixtures](https://docs.pytest.org/en/latest/fixture.html)
- [Parametrize](https://docs.pytest.org/en/latest/parametrize.html)

### Coverage Tools
- [Coverage.py](https://coverage.readthedocs.io/)
- [Pytest-cov](https://pytest-cov.readthedocs.io/)
- [Codecov](https://docs.codecov.io/)

### Testing Best Practices
- [Testing Scientific Code](https://www.software.ac.uk/guide/testing-your-software)
- [NumPy Testing Guidelines](https://numpy.org/doc/stable/reference/testing.html)

---

## Appendix: Example Test Files

### Example 1: Unit Test for Integrators

```python
"""Unit tests for numerical integrators."""
import pytest
import numpy as np
from src.integrators import (
    TrapezoidIntegrator, 
    BackwardTrapezoidIntegrator,
    SimpsonIntegrator,
    BackwardSimpsonIntegrator
)

class TestTrapezoidIntegrator:
    """Test naive trapezoidal integration."""
    
    def test_constant_function(self):
        """Integral of constant function f(x)=c should be c*(b-a)."""
        dt = 0.1
        f = [5.0] * 10  # Constant function f=5
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=9)
        
        result = integrator.integrate(f)
        expected = 5.0 * 0.9  # c * (9-0) * dt
        
        assert abs(result - expected) < 1e-10
    
    def test_linear_function(self):
        """Integral of f(x)=x should be x²/2."""
        dt = 1.0
        f = [0.0, 1.0, 2.0, 3.0, 4.0]  # f(x) = x
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=4)
        
        result = integrator.integrate(f)
        expected = 8.0  # ∫₀⁴ x dx = 4²/2 = 8
        
        assert abs(result - expected) < 1e-6
    
    @pytest.mark.parametrize("n_points,expected_error", [
        (10, 1e-3),
        (100, 1e-5),
        (1000, 1e-7),
    ])
    def test_convergence_rate(self, n_points, expected_error):
        """Test that error decreases with finer discretization."""
        # Integrate sin(x) from 0 to π
        x = np.linspace(0, np.pi, n_points)
        f = np.sin(x)
        dt = x[1] - x[0]
        
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=n_points-1)
        result = integrator.integrate(f.tolist())
        
        expected = 2.0  # ∫₀^π sin(x) dx = 2
        error = abs(result - expected)
        
        assert error < expected_error
```

### Example 2: Integration Test for Kinetics

```python
"""Integration tests for kinetics computations."""
import pytest
import numpy as np
from src.kinetics import Kinetics
from src.constituent import Constituent
from src.kinetics_interface import ConstituentKineticsContext
from tests.utils.mock_objects import create_test_constituent

class TestSurvivalFunctionComputation:
    """Test survival function q(s,tau) computation."""
    
    @pytest.fixture
    def simple_constituent(self):
        """Create a minimal constituent for testing."""
        return create_test_constituent(
            k_alpha_values=[0.01] * 10,  # Constant degradation
            timesteps=10
        )
    
    def test_survival_decreases_monotonically(self, simple_constituent):
        """Survival function should decrease with cohort age."""
        context = ConstituentKineticsContext(simple_constituent)
        kinetics = simple_constituent.kinetics
        
        current_timestep = 5
        q_values = kinetics.compute_survival_history(context, current_timestep)
        
        # q(s,s) = 1, and q should decrease for older cohorts
        assert q_values[-1] == pytest.approx(1.0)  # Most recent cohort
        assert all(q_values[i] <= q_values[i+1] for i in range(len(q_values)-1))
    
    def test_survival_with_known_solution(self):
        """Test against analytical solution for constant k_alpha."""
        k_const = 0.1  # 1/day
        dt = 1.0  # day
        s = 10  # Current time
        tau = 5  # Deposition time
        
        # Analytical: q(s,tau) = exp(-k * (s-tau))
        expected = np.exp(-k_const * (s - tau))
        
        # Create constituent with constant k_alpha
        from tests.utils.mock_objects import create_constant_degradation_constituent
        constituent = create_constant_degradation_constituent(k_const, s+1, dt)
        
        context = ConstituentKineticsContext(constituent)
        kinetics = constituent.kinetics
        q_history = kinetics.compute_survival_history(context, s)
        
        actual = q_history[tau]
        assert actual == pytest.approx(expected, rel=1e-5)
```

### Example 3: End-to-End Test

```python
"""End-to-end test for Latorre2018 validation case."""
import pytest
import numpy as np
from pathlib import Path
from src.simulation import Simulation
from tests.utils.output_manager import TestOutputManager
from tests.utils.comparisons import assert_simulation_matches

class TestLatorre2018Validation:
    """Validate against published Latorre et al. 2018 results."""
    
    @pytest.fixture
    def config_path(self):
        return Path(__file__).parent.parent / "fixtures" / "inputs" / "latorre2018_test.yaml"
    
    @pytest.fixture
    def expected_outputs(self):
        base_path = Path(__file__).parent.parent / "fixtures" / "expected_outputs" / "latorre2018_validated"
        return {
            'geometry': np.load(base_path / 'geometry.npy'),
            'rhoR_alpha': np.load(base_path / 'rhoR_alpha.npy'),
            'stress': np.load(base_path / 'stress.npy'),
        }
    
    @pytest.mark.slow
    def test_full_simulation_matches_validation(self, config_path, expected_outputs):
        """Run full Latorre2018 and compare to validated outputs."""
        output_mgr = TestOutputManager('latorre2018_full', 'e2e')
        
        # Run simulation
        sim = Simulation.from_yaml(config_path)
        results = sim.run()
        
        # Extract key outputs
        actual_geometry = results.get_geometry_history()
        
        try:
            # Compare with tolerances
            assert_simulation_matches(
                actual_geometry, 
                expected_outputs['geometry'],
                rtol=1e-4,
                atol=1e-8,
                output_mgr=output_mgr,
                name='geometry'
            )
            output_mgr.cleanup_on_success()
        except AssertionError:
            # Outputs saved automatically by assert_simulation_matches
            raise
```

---

**End of Document**

This design document should be referenced throughout the implementation of the testing framework to ensure consistency with the design goals and principles outlined here.
```

I'll save this as a markdown document in the docs folder:

Since I don't have file editing capabilities, here's what you should do:

**Create the file:** `docs/TESTING_FRAMEWORK_DESIGN.md`

**Copy the entire markdown content** I've provided above into that file.

This document now serves as your comprehensive reference for:
- **Testing concepts** (pytest, coverage, fixtures, marks)
- **Design principles** (what goes where and why)
- **Complete architecture** (directory structure, file organization)
- **Implementation roadmap** (6-week phased approach)
- **Best practices** (how to write and maintain tests)
- **Example code** (ready-to-use templates)

When you're ready to start implementing, you can reference specific sections of this document to guide your work. Would you like me to proceed with creating the actual testing framework files, starting with Phase 1 (foundation)?I'll save this as a markdown document in the docs folder:

Since I don't have file editing capabilities, here's what you should do:

**Create the file:** `docs/TESTING_FRAMEWORK_DESIGN.md`

**Copy the entire markdown content** I've provided above into that file.

This document now serves as your comprehensive reference for:
- **Testing concepts** (pytest, coverage, fixtures, marks)
- **Design principles** (what goes where and why)
- **Complete architecture** (directory structure, file organization)
- **Implementation roadmap** (6-week phased approach)
- **Best practices** (how to write and maintain tests)
- **Example code** (ready-to-use templates)

When you're ready to start implementing, you can reference specific sections of this document to guide your work. Would you like me to proceed with creating the actual testing framework files, starting with Phase 1 (foundation)?

User: Can you re-write your very last response in a way to be copy-pasted directly to a .md file? Currently I cannot easily copy-paste the entire reply.

GitHub Copilot: ```markdown
# svGrowth Testing Framework Design

**Document Purpose:** Reference guide for implementing the testing framework for the svGrowth codebase. This document outlines design goals, principles, key concepts, and the complete testing architecture.

**Created:** January 3, 2026  
**Status:** Design Phase

---

## Table of Contents

1. [Testing Philosophy & Goals](#testing-philosophy--goals)
2. [Test Types Overview](#test-types-overview)
3. [Pytest Fundamentals](#pytest-fundamentals)
4. [Code Coverage Strategy](#code-coverage-strategy)
5. [Exceptions vs Assertions](#exceptions-vs-assertions)
6. [Testing Framework Architecture](#testing-framework-architecture)
7. [Test Data Management](#test-data-management)
8. [GitHub Actions CI/CD](#github-actions-cicd)
9. [Utility Modules](#utility-modules)
10. [Implementation Roadmap](#implementation-roadmap)

---

## Testing Philosophy & Goals

### Design Principles

1. **Comprehensive Coverage:** Test at multiple granularities (unit → integration → end-to-end)
2. **Isolation:** Unit tests should not depend on full simulation setup
3. **Reproducibility:** All tests must produce deterministic results
4. **Fast Feedback:** Unit tests run in milliseconds, integration in seconds
5. **Clear Failure Messages:** When tests fail, immediately understand why
6. **Maintainability:** Tests should be easy to update as code evolves

### Success Criteria

- **Coverage Target:** ≥85% code coverage across all modules
- **Speed:** Unit test suite completes in <10 seconds
- **Reliability:** Zero flaky tests (non-deterministic failures)
- **Documentation:** Every test clearly documents what it verifies

---

## Test Types Overview

### 1. Unit Tests
**Purpose:** Verify that smallest code units (functions, methods) work correctly in isolation

**Characteristics:**
- **Speed:** Fast (milliseconds per test)
- **Dependencies:** Minimal - use mocks for external dependencies
- **Scope:** Single function or method
- **Location:** `tests/unit/`

**Examples for svGrowth:**
- Integration methods in `integrators.py` (trapezoidal, Simpson's rule)
- Tensor operations in `tensor_operations.py`
- Individual stress computation functions in `constitutive_laws.py`
- Helper functions and utilities

**Example:**
```python
def test_trapezoid_integrator_constant_function():
    """Integral of f(x)=c should equal c*(b-a)."""
    dt = 0.1
    f = [5.0] * 10  # Constant function
    integrator = TrapezoidIntegrator(dt=dt, start=0, stop=9)
    
    result = integrator.integrate(f)
    expected = 5.0 * 0.9  # c * (9-0) * dt
    
    assert abs(result - expected) < 1e-10
```

### 2. Integration Tests
**Purpose:** Verify that different modules work together correctly

**Characteristics:**
- **Speed:** Medium (seconds per test)
- **Dependencies:** Real classes working together, but not full simulation
- **Scope:** Multiple interacting components
- **Location:** `tests/integration/`

**Examples for svGrowth:**
- Survival function computation (kinetics + integrators)
- `compute_rhoR_alpha()` workflow (constituent + kinetics + layer data)
- Stress computation workflow (mechanics + constitutive laws + kinematics)
- Solver iteration logic (solver + mechanics + layer)

**Example:**
```python
def test_survival_function_with_constant_degradation():
    """Test survival computation with analytical solution."""
    k_const = 0.1  # Constant degradation rate
    constituent = create_test_constituent(k_alpha=[k_const] * 10)
    
    context = ConstituentKineticsContext(constituent)
    q_history = constituent.kinetics.compute_survival_history(context, timestep=5)
    
    # Compare to analytical: q(s,τ) = exp(-k*(s-τ))
    for tau, q_value in enumerate(q_history):
        expected = np.exp(-k_const * (5 - tau))
        assert abs(q_value - expected) < 1e-6
```

### 3. End-to-End (E2E) Tests
**Purpose:** Verify entire system produces correct outputs for known inputs

**Characteristics:**
- **Speed:** Slow (minutes per test)
- **Dependencies:** Full simulation stack
- **Scope:** Complete workflows from input YAML to output files
- **Location:** `tests/e2e/`

**Examples for svGrowth:**
- Full Latorre2018 simulation matching published results
- Pressure perturbation scenarios
- Axial stretch perturbation validation
- Convergence to homeostasis verification

**Example:**
```python
def test_latorre2018_full_validation():
    """Verify full simulation matches Latorre et al. 2018 results."""
    sim = Simulation.from_yaml('tests/fixtures/inputs/latorre2018_test.yaml')
    results = sim.run()
    
    # Load validated reference data
    expected = np.load('tests/fixtures/expected_outputs/latorre2018/geometry.npy')
    
    assert_simulation_matches(results.geometry, expected, rtol=1e-4)
```

---

## Pytest Fundamentals

### Pytest Fixtures

**What are they?** Reusable setup/teardown code that provides data or objects to tests.

**Think of them as:** "Test ingredients" - define once, pytest automatically provides to any test that needs them.

**Basic Example:**
```python
import pytest

@pytest.fixture
def sample_constituent():
    """Creates a constituent for testing - runs BEFORE each test."""
    constituent = Constituent(name="elastin")
    constituent.rhoR_alpha_history = [1.0, 1.1, 1.2]
    constituent.k_alpha_history = [0.01, 0.01, 0.01]
    
    yield constituent  # Provide to test
    
    # Code after yield runs AFTER test (cleanup)
    constituent.cleanup()  # If needed

# Use fixture by adding it as parameter
def test_get_density(sample_constituent):
    # pytest automatically creates and passes sample_constituent
    density = sample_constituent.get_rhoR_alpha(timestep=0)
    assert density == 1.0
```

**Fixture Scopes:**
```python
@pytest.fixture(scope="function")  # Default: new instance per test
def per_test_fixture():
    return "recreated each test"

@pytest.fixture(scope="class")  # Shared across test class
def per_class_fixture():
    return "recreated per test class"

@pytest.fixture(scope="module")  # Shared across entire file
def expensive_setup():
    return load_large_dataset()

@pytest.fixture(scope="session")  # Once per entire test run
def database_connection():
    db = create_connection()
    yield db
    db.close()
```

### Pytest Marks

**What are they?** Labels/tags to categorize or modify test behavior.

**Common Marks:**
```python
# Mark slow tests
@pytest.mark.slow
def test_full_simulation():
    """Takes 5 minutes to run."""
    pass

# Mark tests needing external resources
@pytest.mark.requires_internet
def test_download_data():
    pass

# Mark expected failures
@pytest.mark.xfail(reason="Known bug in upstream library")
def test_with_known_issue():
    pass

# Skip tests conditionally
@pytest.mark.skipif(sys.version_info < (3, 9), reason="Requires Python 3.9+")
def test_new_feature():
    pass

# Parametrize - run same test with different inputs
@pytest.mark.parametrize("input,expected", [
    (1, 2),
    (2, 4),
    (3, 6),
])
def test_doubling(input, expected):
    assert input * 2 == expected
```

**Running Tests Selectively:**
```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit/

# Skip slow tests
pytest -m "not slow"

# Run ONLY slow tests
pytest -m "slow"

# Run tests matching pattern
pytest -k "integrator"

# Verbose output
pytest -v

# Stop at first failure
pytest -x

# Run in parallel (requires pytest-xdist)
pytest -n auto
```

---

## Code Coverage Strategy

### How Coverage Works

1. **Instrumentation:** `pytest-cov` (built on `coverage.py`) instruments code to track line execution
2. **Tracking:** During test execution, records which lines run
3. **Reporting:** Generates reports showing covered/uncovered lines

**Visual Example:**
```python
# src/integrators.py
def integrate(self, f: List[float]) -> float:
    """Compute integral."""
    if len(f) < 2:           # ✅ Covered by test_empty_list
        raise ValueError()    # ✅ Covered by test_empty_list
    
    result = 0.0             # ✅ Covered by test_constant_function
    for i in range(len(f)):  # ✅ Covered by test_constant_function
        result += f[i]       # ✅ Covered by test_constant_function
    
    return result * self.dt  # ❌ NOT COVERED - missing return test!
```

### GitHub Integration Flow

```
┌─────────────────────────────────────────────┐
│ 1. Developer pushes code to GitHub         │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 2. GitHub Actions triggers workflow         │
│    (.github/workflows/tests.yml)            │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 3. Workflow runs:                           │
│    pytest --cov=src --cov-report=xml        │
│                                             │
│    Creates: coverage.xml                    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 4. Upload to Codecov:                       │
│    codecov/codecov-action@v3                │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 5. Codecov analyzes coverage.xml            │
│    - Tracks line-by-line coverage           │
│    - Compares to previous commits           │
│    - Comments on PR with coverage change    │
└──────────────────┬──────────────────────────┘
                   ↓
┌─────────────────────────────────────────────┐
│ 6. Badge shows on README:                   │
│    ![Coverage](codecov.io/badge.svg)        │
│    "92% coverage"                           │
└─────────────────────────────────────────────┘
```

**Coverage XML Structure:**
```xml
<coverage>
  <packages>
    <package name="src">
      <classes>
        <class name="integrators.py">
          <lines>
            <line number="45" hits="12"/>  <!-- Executed 12 times -->
            <line number="46" hits="0"/>   <!-- Never executed -->
          </lines>
        </class>
      </classes>
    </package>
  </packages>
</coverage>
```

### Running Coverage Locally

```bash
# Run tests with coverage
pytest tests/ --cov=src --cov-report=html

# View results
open htmlcov/index.html

# Terminal output with missing lines
pytest tests/ --cov=src --cov-report=term-missing

# Example output:
# Name                    Stmts   Miss  Cover   Missing
# -----------------------------------------------------
# src/integrators.py         45      3    93%   102-104
# src/kinetics.py           120      8    93%   45, 67-73
```

---

## Exceptions vs Assertions

### Exceptions (Production Code)

**Purpose:** Handle runtime errors in production code

**Location:** exceptions.py

**When to Use:**
- Invalid user input
- Physical violations during simulation
- Solver convergence failures
- Configuration errors

**Example:**
```python
# src/solver.py
from exceptions import ConvergenceError, PhysicsViolationError

def solve_equilibrium(self, max_iter=100):
    """Solve for equilibrium state."""
    for i in range(max_iter):
        residual = self._compute_residual()
        
        if residual < self.tolerance:
            return self.state  # Success
    
    # Raise exception if fails to converge
    raise ConvergenceError(
        f"Failed to converge after {max_iter} iterations. "
        f"Final residual: {residual:.2e}"
    )

def compute_stress(self, stretch):
    """Compute stress from stretch."""
    if stretch < 0:
        raise PhysicsViolationError(
            f"Negative stretch detected: {stretch}"
        )
    return self.material_model(stretch)
```

**Exception Hierarchy:**
```python
# src/exceptions.py
class GrowthRemodelingError(Exception):
    """Base class - catch all G&R errors."""
    pass

class SolverError(GrowthRemodelingError):
    """Problems during solving."""
    pass

class ConvergenceError(SolverError):
    """Specific: didn't converge."""
    pass

class PhysicsViolationError(GrowthRemodelingError):
    """Specific: violated physical constraints."""
    pass
```

### Assertions (Test Code)

**Purpose:** Verify conditions in tests

**Location:** Test files (tests)

**Two Types:**

#### 1. Python's `assert` Statement (Development/Debug)
```python
# In production code (optional, disabled with python -O)
def compute_survival(self, k_alpha, dt):
    assert k_alpha > 0, "Degradation rate must be positive"
    assert dt > 0, "Time step must be positive"
    return math.exp(-k_alpha * dt)

# These are REMOVED when running python with optimization flag!
# Don't rely on them for critical checks in production!
```

#### 2. Pytest Assertions (Testing)
```python
# In tests/unit/test_kinetics.py
def test_survival_function():
    result = compute_survival(k_alpha=0.1, dt=1.0)
    
    # Basic assertion
    assert result > 0
    assert result < 1
    
    # With helpful message
    assert result == pytest.approx(0.9048, abs=1e-4), \
        f"Expected ~0.9048, got {result}"
    
    # Numpy assertions (better error messages)
    import numpy.testing as npt
    npt.assert_allclose(result, 0.9048, rtol=1e-5)
```

### Design Decision: No `assertions.py` in Main Code

**Recommendation:** Keep assertions ONLY in `tests/utils/assertions.py`

**Rationale:**
1. Production code should use exceptions for error handling
2. Assertions can be disabled (`python -O` removes them)
3. Tests need different assertion helpers than production code

**File Organization:**

| File | Purpose | Used By |
|------|---------|---------|
| exceptions.py | Production error handling | Production code + tests |
| `src/validation.py` | Input validation (raises exceptions) | Production code |
| `tests/utils/assertions.py` | Test-specific assertion helpers | Tests only |

**Examples:**

```python
# src/validation.py (production code)
from exceptions import GeometryError

def validate_thin_wall_assumption(radius, thickness):
    """Validate thin-wall assumption.
    
    Raises:
        GeometryError: If h/a >= 0.2
    """
    ratio = thickness / radius
    if ratio >= 0.2:
        raise GeometryError(
            f"Thin-wall assumption violated: h/a = {ratio:.3f} >= 0.2"
        )

# tests/utils/assertions.py (test helpers)
import numpy as np

def assert_thin_wall_valid(radius, thickness):
    """Test helper: assert thin-wall assumption holds."""
    ratio = thickness / radius
    assert ratio < 0.2, \
        f"Test expected thin-wall (h/a < 0.2), got {ratio:.3f}"

def assert_arrays_close_with_summary(actual, expected, rtol=1e-5):
    """Assert arrays close with detailed error message."""
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol)
    except AssertionError as e:
        diff = actual - expected
        max_err = np.max(np.abs(diff))
        max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
        
        raise AssertionError(
            f"Arrays not close:\n"
            f"  Max error: {max_err:.2e} at index {max_idx}\n"
            f"  Expected: {expected[max_idx]}\n"
            f"  Actual:   {actual[max_idx]}\n"
            f"  Original error: {e}"
        )
```

---

## Testing Framework Architecture

### Directory Structure

```
tests/
├── __init__.py
├── conftest.py                      # Pytest fixtures and configuration
├── test_config.yaml                 # Testing-specific configuration
├── requirements_test.txt            # Test dependencies
│
├── unit/                            # Unit tests
│   ├── __init__.py
│   ├── test_integrators.py         # Test integration methods
│   ├── test_tensor_operations.py   # Test tensor utilities
│   ├── test_constitutive_laws.py   # Test stress-strain laws
│   ├── test_deformation_kinematics.py
│   └── test_helper_functions.py    # Test utility functions
│
├── integration/                     # Integration tests
│   ├── __init__.py
│   ├── test_kinetics.py            # Test survival/production computations
│   ├── test_constituent.py         # Test compute_rhoR_alpha workflows
│   ├── test_mechanics.py           # Test stress computation workflows
│   └── test_solver_components.py   # Test solver iterations
│
├── e2e/                            # End-to-end tests
│   ├── __init__.py
│   ├── test_latorre2018.py        # Full Latorre2018 validation
│   └── test_perturbations.py      # Test perturbation scenarios
│
├── fixtures/                       # Test data and configurations
│   ├── inputs/                    # Input YAML files for tests
│   │   ├── latorre2018_test.yaml
│   │   ├── simple_cylinder.yaml
│   │   ├── unit_test_configs/
│   │   │   ├── simple_integrator.yaml
│   │   │   └── single_constituent.yaml
│   │   └── integration_configs/
│   │       ├── kinetics_only.yaml
│   │       └── two_constituent.yaml
│   ├── expected_outputs/          # Known-good reference outputs
│   │   ├── latorre2018_validated/
│   │   │   ├── metadata.json
│   │   │   ├── geometry.npy
│   │   │   ├── rhoR_alpha.npy
│   │   │   └── stress.npy
│   │   └── integration_test_data/
│   │       ├── survival_function_expected.json
│   │       └── rhoR_alpha_computation_expected.json
│   └── test_data/                 # Synthetic data for unit tests
│       ├── integrator_test_data.json
│       └── survival_function_data.json
│
├── utils/                         # Testing utilities
│   ├── __init__.py
│   ├── assertions.py             # Custom assertion helpers
│   ├── comparisons.py            # Numerical comparison utilities
│   ├── fixtures_generators.py    # Generate test data
│   ├── mock_objects.py           # Mock classes for testing
│   └── output_manager.py         # Manage test outputs
│
└── outputs/                       # Test run outputs (gitignored)
    ├── .gitkeep
    ├── unit/
    ├── integration/
    └── e2e/
```

### Key Files

#### `tests/conftest.py`
Central pytest configuration with shared fixtures available to all tests.

```python
"""Pytest configuration and shared fixtures."""
import pytest
from pathlib import Path

@pytest.fixture(scope="session")
def test_data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "fixtures" / "test_data"

@pytest.fixture(scope="session")
def expected_outputs_dir():
    """Path to expected outputs directory."""
    return Path(__file__).parent / "fixtures" / "expected_outputs"

@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for test outputs."""
    return tmp_path / "test_outputs"
```

#### `tests/requirements_test.txt`
Testing-specific dependencies (not needed for production).

```
pytest>=7.4.0
pytest-cov>=4.1.0
pytest-xdist>=3.3.0      # Parallel test execution
pytest-timeout>=2.1.0     # Timeout protection
hypothesis>=6.82.0        # Property-based testing
numpy-testing>=0.1.0      # Numerical assertions
```

#### `tests/.gitignore`
Prevent committing test outputs.

```
# Test outputs
outputs/
*.pyc
__pycache__/
.pytest_cache/
htmlcov/
.coverage
coverage.xml
```

---

## Test Data Management

### Input Data Structure

```
tests/fixtures/inputs/
├── unit_test_configs/          # Minimal configs for unit tests
│   ├── simple_integrator.yaml
│   └── single_constituent.yaml
├── integration_configs/        # Partial simulation configs
│   ├── kinetics_only.yaml
│   └── two_constituent.yaml
└── e2e_configs/               # Full simulation configs
    ├── latorre2018_test.yaml
    └── validation_cases/
        ├── pressure_perturbation.yaml
        └── axial_stretch_perturbation.yaml
```

### Expected Output Structure

```
tests/fixtures/expected_outputs/
├── latorre2018_validated/
│   ├── metadata.json          # Simulation parameters, tolerances
│   ├── geometry.npy           # [timesteps × 3] (a, h, λz)
│   ├── rhoR_alpha.npy         # [timesteps × n_constituents]
│   └── stress.npy             # [timesteps × 3] (σr, σθ, σz)
└── integration_test_data/
    ├── survival_function_expected.json
    └── rhoR_alpha_computation_expected.json
```

**`metadata.json` Structure:**
```json
{
  "description": "Latorre et al. 2018 validation case",
  "source": "Published paper results",
  "timesteps": 100,
  "dt": 1.0,
  "tolerances": {
    "geometry": {"rtol": 1e-4, "atol": 1e-8},
    "stress": {"rtol": 1e-3, "atol": 1e-6}
  },
  "generated_date": "2026-01-03",
  "generated_by": "validation_plots.py"
}
```

### Test Output Structure (Gitignored)

```
tests/outputs/                  # Gitignored - runtime outputs
├── test_run_<timestamp>/
│   ├── logs/
│   │   ├── test_integration.log
│   │   └── test_e2e.log
│   ├── failures/              # Only written on failure
│   │   ├── test_latorre2018_geometry_diff.png
│   │   ├── test_latorre2018_actual_vs_expected.json
│   │   └── test_survival_function_metrics.json
│   └── coverage/
│       └── htmlcov/
```

### Data Management Principles

1. **Version Control:**
   - ✅ Commit: Input YAMLs, expected outputs, test data
   - ❌ Don't commit: Test run outputs, logs, temporary files

2. **Data Formats:**
   - Numpy arrays: `.npy` (binary, fast, exact)
   - Metadata: `.json` (human-readable, version-controllable)
   - Logs: `.log` (text, easy to grep)

3. **Storage Organization:**
   - Group by test type (unit/integration/e2e)
   - Use timestamps for run-specific outputs
   - Save detailed outputs only on failure

---

## GitHub Actions CI/CD

### Workflow File: `.github/workflows/tests.yml`

```yaml
name: svGrowth Tests

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.9', '3.10', '3.11']
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install -r tests/requirements_test.txt
    
    - name: Run unit tests
      run: |
        pytest tests/unit/ -v --cov=src --cov-report=xml --cov-report=term-missing
    
    - name: Run integration tests
      run: |
        pytest tests/integration/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Run E2E tests
      run: |
        pytest tests/e2e/ -v --cov=src --cov-append --cov-report=xml
    
    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        files: ./coverage.xml
        flags: unittests
        name: codecov-umbrella
        fail_ci_if_error: true

  test-fast:
    # Quick smoke test for rapid feedback
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    
    - uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    
    - name: Install dependencies
      run: |
        pip install -r requirements.txt -r tests/requirements_test.txt
    
    - name: Run unit tests only (fast)
      run: |
        pytest tests/unit/ -v --maxfail=1  # Fail fast
```

### Workflow Triggers

- **Push to main/develop:** Run full test suite
- **Pull requests:** Run full test suite + coverage comparison
- **Scheduled (optional):** Daily full test run to catch environment issues

### Integration with Codecov

1. Sign up at [codecov.io](https://codecov.io)
2. Add repository to Codecov
3. Get upload token (stored in GitHub Secrets)
4. Codecov automatically:
   - Comments on PRs with coverage changes
   - Shows line-by-line coverage in GitHub UI
   - Generates coverage badges

**Badge for README:**
```markdown
[![codecov](https://codecov.io/gh/username/svGrowth/branch/main/graph/badge.svg)](https://codecov.io/gh/username/svGrowth)
```

---

## Utility Modules

### `tests/utils/comparisons.py`

**Purpose:** Advanced numerical comparison utilities for simulation outputs.

```python
"""Numerical comparison utilities for tests."""
import numpy as np
from pathlib import Path
import json

def assert_simulation_matches(actual, expected, rtol=1e-5, atol=1e-8, 
                               output_mgr=None, name='output'):
    """Compare simulation outputs with detailed error reporting.
    
    Args:
        actual: Numpy array of actual results
        expected: Numpy array of expected results
        rtol: Relative tolerance
        atol: Absolute tolerance
        output_mgr: TestOutputManager for saving failures
        name: Name for output files
    
    Raises:
        AssertionError: If arrays don't match within tolerance
    """
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)
    except AssertionError as e:
        if output_mgr:
            output_mgr.save_failure_comparison(actual, expected, name)
        
        # Enhanced error message
        diff = np.abs(actual - expected)
        rel_err = diff / (np.abs(expected) + 1e-16)
        
        error_msg = (
            f"\nSimulation output mismatch for '{name}':\n"
            f"  Max absolute error: {np.max(diff):.2e}\n"
            f"  Max relative error: {np.max(rel_err):.2e}\n"
            f"  Mean absolute error: {np.mean(diff):.2e}\n"
            f"  Tolerances: rtol={rtol}, atol={atol}\n"
        )
        
        if output_mgr:
            error_msg += f"  Detailed comparison saved to: {output_mgr.output_dir}\n"
        
        error_msg += f"\nOriginal error:\n{e}"
        raise AssertionError(error_msg)


def compare_history_convergence(history, expected_final, 
                                 convergence_timestep=None):
    """Verify history converges to expected final value.
    
    Args:
        history: List or array of values over time
        expected_final: Expected final converged value
        convergence_timestep: Timestep by which it should converge (optional)
    
    Returns:
        bool: True if converged appropriately
    """
    history = np.array(history)
    
    # Check final value
    final_value = history[-1]
    is_close = np.isclose(final_value, expected_final, rtol=1e-4)
    
    if not is_close:
        raise AssertionError(
            f"Final value {final_value:.4e} doesn't match "
            f"expected {expected_final:.4e}"
        )
    
    # Check convergence timing if specified
    if convergence_timestep is not None:
        values_after_convergence = history[convergence_timestep:]
        if not np.allclose(values_after_convergence, expected_final, rtol=1e-3):
            raise AssertionError(
                f"Values after timestep {convergence_timestep} not converged"
            )
    
    return True


def compare_constituent_histories(actual_constituent, expected_constituent, 
                                   check_fields=None):
    """Compare all history fields between two constituents.
    
    Args:
        actual_constituent: Constituent from test
        expected_constituent: Reference constituent
        check_fields: List of field names to check (default: all histories)
    """
    if check_fields is None:
        check_fields = [
            'rhoR_alpha_history',
            'k_alpha_history',
            'mR_alpha_history',
            'survival_history',
        ]
    
    errors = []
    
    for field in check_fields:
        actual = getattr(actual_constituent, field)
        expected = getattr(expected_constituent, field)
        
        try:
            if isinstance(actual[0], list):  # Nested lists (e.g., survival_history)
                for i, (a, e) in enumerate(zip(actual, expected)):
                    np.testing.assert_allclose(a, e, rtol=1e-5,
                        err_msg=f"{field}[{i}]")
            else:
                np.testing.assert_allclose(actual, expected, rtol=1e-5,
                    err_msg=field)
        except AssertionError as e:
            errors.append(str(e))
    
    if errors:
        raise AssertionError(
            f"Constituent histories don't match:\n" + "\n".join(errors)
        )
```

### `tests/utils/mock_objects.py`

**Purpose:** Create fake/simplified objects for isolated testing without full simulation setup.

```python
"""Mock objects for testing without full simulation setup."""
import numpy as np
from typing import List, Optional
from src.constituent import Constituent
from src.kinetics import Kinetics, ConstantDegradationRate
from src.layer import Layer

class MockLayer:
    """Minimal Layer for testing without full geometry."""
    
    def __init__(self, timesteps=10):
        self.inner_radius_history = [1.0] * timesteps
        self.thickness_history = [0.1] * timesteps
        self.lambda_z_history = [1.0] * timesteps
        self.intramural_stress_history = [100.0] * timesteps
        
        # Kinematics mock
        self.kinematics = MockKinematics()
    
    def get_inner_radius(self, timestep):
        return self.inner_radius_history[timestep]


class MockKinematics:
    """Minimal kinematics for testing stress components."""
    
    def get_component_names(self):
        return ['r', 'theta', 'z']
    
    def get_component_index(self, name):
        mapping = {'r': 0, 'theta': 1, 'z': 2}
        return mapping[name]
    
    def get_intramural_component_name(self):
        return 'theta'


def create_test_constituent(name="test_constituent", 
                            k_alpha_values=None, 
                            timesteps=10,
                            homeostatic_density=1.0):
    """Create a constituent with preset histories for testing.
    
    Args:
        name: Constituent name
        k_alpha_values: List of degradation rates (default: constant 0.01)
        timesteps: Number of timesteps to initialize
        homeostatic_density: Initial density
    
    Returns:
        Constituent: Ready for testing
    """
    from src.constituent import SingleConstituent
    
    # Create minimal constituent
    constituent = SingleConstituent(name=name)
    constituent.homeostatic_referential_density = homeostatic_density
    constituent.layer = MockLayer(timesteps)
    
    # Initialize kinetics
    if k_alpha_values is None:
        k_alpha_values = [0.01] * timesteps
    
    constituent.kinetics = Kinetics(
        degradation_model=ConstantDegradationRate(k_alpha_h=0.01),
        production_model=None,
        integrator_type='trapezoid',
        dt=1.0
    )
    
    # Populate histories
    constituent.k_alpha_history = k_alpha_values
    constituent.rhoR_alpha_history = [homeostatic_density] * timesteps
    constituent.survival_history = [[1.0]] * timesteps
    
    return constituent


def create_constant_degradation_constituent(k_const, timesteps, dt):
    """Create constituent with constant degradation for analytical testing.
    
    Useful for tests where analytical solution exists.
    
    Args:
        k_const: Constant degradation rate (1/day)
        timesteps: Number of timesteps
        dt: Time step size
    
    Returns:
        Constituent: With constant k_alpha throughout history
    """
    constituent = create_test_constituent(
        name="analytical_test",
        k_alpha_values=[k_const] * timesteps,
        timesteps=timesteps
    )
    constituent.kinetics.dt = dt
    
    return constituent


class MockIntegrator:
    """Simple mock integrator for testing without real integration."""
    
    def __init__(self, return_value=1.0):
        self.return_value = return_value
        self.calls = []  # Track calls for verification
    
    def integrate(self, f):
        """Record call and return preset value."""
        self.calls.append(f)
        return self.return_value
    
    def was_called_with(self, expected_function):
        """Check if integrator was called with expected values."""
        return any(
            np.allclose(call, expected_function, rtol=1e-6) 
            for call in self.calls
        )
```

### `tests/utils/output_manager.py`

**Purpose:** Manage test outputs in a clean, modular way.

```python
"""Manage test outputs in a clean, modular way."""
from pathlib import Path
from datetime import datetime
import json
import numpy as np
import shutil

class TestOutputManager:
    """Handles test output storage and cleanup."""
    
    def __init__(self, test_name: str, test_type: str):
        """
        Args:
            test_name: Name of test (e.g., 'test_latorre2018')
            test_type: 'unit', 'integration', or 'e2e'
        """
        self.test_name = test_name
        self.test_type = test_type
        
        # Create timestamped output directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(__file__).parent.parent / "outputs" / test_type / f"{test_name}_{timestamp}"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def save_array(self, name: str, array: np.ndarray):
        """Save numpy array."""
        np.save(self.output_dir / f"{name}.npy", array)
    
    def save_dict(self, name: str, data: dict):
        """Save dictionary as JSON."""
        with open(self.output_dir / f"{name}.json", 'w') as f:
            json.dump(data, f, indent=2)
    
    def save_failure_comparison(self, actual, expected, name: str):
        """Save detailed comparison for failed test."""
        failure_dir = self.output_dir / "failures"
        failure_dir.mkdir(exist_ok=True)
        
        # Save both arrays
        np.save(failure_dir / f"{name}_actual.npy", actual)
        np.save(failure_dir / f"{name}_expected.npy", expected)
        
        # Save difference
        diff = actual - expected
        np.save(failure_dir / f"{name}_diff.npy", diff)
        
        # Save metadata
        metadata = {
            "max_absolute_error": float(np.max(np.abs(diff))),
            "mean_absolute_error": float(np.mean(np.abs(diff))),
            "relative_error": float(np.max(np.abs(diff / (expected + 1e-16)))),
        }
        self.save_dict(f"failures/{name}_metrics", metadata)
    
    def cleanup_on_success(self):
        """Remove output directory if test passed."""
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
```

### `tests/utils/fixtures_generators.py`

**Purpose:** Generate test data programmatically.

```python
"""Generate test fixtures programmatically."""
import numpy as np
import json
from pathlib import Path

def generate_survival_test_data(output_path: Path):
    """Generate reference data for survival function tests.
    
    Creates analytical test cases where exact solution is known.
    """
    test_cases = []
    
    # Case 1: Constant degradation
    k_const = 0.1
    dt = 1.0
    timesteps = 10
    
    for s in range(timesteps):
        q_values = []
        for tau in range(s + 1):
            # Analytical: q(s,τ) = exp(-k*(s-τ)*dt)
            q_analytical = np.exp(-k_const * (s - tau) * dt)
            q_values.append(q_analytical)
        
        test_cases.append({
            'current_timestep': s,
            'k_alpha_history': [k_const] * (s + 1),
            'dt': dt,
            'expected_q': q_values
        })
    
    # Save to JSON
    with open(output_path, 'w') as f:
        json.dump(test_cases, f, indent=2)
    
    print(f"Generated {len(test_cases)} test cases in {output_path}")


def generate_integration_test_data():
    """Generate test data for numerical integrators.
    
    Returns functions with known analytical integrals.
    """
    test_functions = []
    
    # Constant function: ∫c dx = c*x
    test_functions.append({
        'name': 'constant',
        'f': lambda x: 5.0,
        'integral': lambda a, b: 5.0 * (b - a),
        'description': 'f(x) = 5'
    })
    
    # Linear: ∫x dx = x²/2
    test_functions.append({
        'name': 'linear',
        'f': lambda x: x,
        'integral': lambda a, b: (b**2 - a**2) / 2,
        'description': 'f(x) = x'
    })
    
    # Quadratic: ∫x² dx = x³/3
    test_functions.append({
        'name': 'quadratic',
        'f': lambda x: x**2,
        'integral': lambda a, b: (b**3 - a**3) / 3,
        'description': 'f(x) = x²'
    })
    
    # Sin: ∫sin(x) dx = -cos(x)
    test_functions.append({
        'name': 'sine',
        'f': lambda x: np.sin(x),
        'integral': lambda a, b: -np.cos(b) + np.cos(a),
        'description': 'f(x) = sin(x)'
    })
    
    return test_functions
```

### `tests/utils/assertions.py`

**Purpose:** Custom test-specific assertion helpers.

```python
"""Custom assertion helpers for tests."""
import numpy as np

def assert_thin_wall_valid(radius, thickness):
    """Test helper: assert thin-wall assumption holds."""
    ratio = thickness / radius
    assert ratio < 0.2, \
        f"Test expected thin-wall (h/a < 0.2), got {ratio:.3f}"


def assert_arrays_close_with_summary(actual, expected, rtol=1e-5):
    """Assert arrays close with detailed error message."""
    try:
        np.testing.assert_allclose(actual, expected, rtol=rtol)
    except AssertionError as e:
        diff = actual - expected
        max_err = np.max(np.abs(diff))
        max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
        
        raise AssertionError(
            f"Arrays not close:\n"
            f"  Max error: {max_err:.2e} at index {max_idx}\n"
            f"  Expected: {expected[max_idx]}\n"
            f"  Actual:   {actual[max_idx]}\n"
            f"  Original error: {e}"
        )


def assert_physical_constraints(geometry_dict):
    """Assert physical constraints on geometry."""
    assert geometry_dict['radius'] > 0, "Radius must be positive"
    assert geometry_dict['thickness'] > 0, "Thickness must be positive"
    assert geometry_dict['lambda_z'] > 0, "Axial stretch must be positive"
```

---

## Implementation Roadmap

### Phase 1: Foundation (Week 1)

**Goal:** Set up testing infrastructure

- [ ] Create directory structure (tests, subdirectories)
- [ ] Create `conftest.py` with basic fixtures
- [ ] Create `requirements_test.txt`
- [ ] Set up .gitignore for test outputs
- [ ] Create utility modules (stubs):
  - `output_manager.py`
  - `mock_objects.py`
  - `comparisons.py`
  - `assertions.py`
  - `fixtures_generators.py`

### Phase 2: Unit Tests (Week 2)

**Goal:** Test individual functions in isolation

- [ ] `test_integrators.py`:
  - Naive trapezoidal integration
  - Backward trapezoidal integration
  - Naive Simpson's rule
  - Backward Simpson's rule
  - Convergence rate tests
- [ ] `test_tensor_operations.py`:
  - Matrix operations
  - Tensor transformations
- [ ] `test_constitutive_laws.py`:
  - Fung exponential model
  - Neo-Hookean model
  - Stress computation functions

**Target:** 50+ unit tests, >80% coverage of tested modules

### Phase 3: Integration Tests (Week 3)

**Goal:** Test component interactions

- [ ] `test_kinetics.py`:
  - Survival function computation
  - Production function computation
  - Integration with degradation models
- [ ] `test_constituent.py`:
  - `compute_rhoR_alpha()` workflow
  - History management
  - Context usage
- [ ] `test_mechanics.py`:
  - Stress computation workflow
  - Equilibrium solving
  - Kinematics integration

**Target:** 30+ integration tests

### Phase 4: End-to-End Tests (Week 4)

**Goal:** Validate complete workflows

- [ ] Generate expected outputs:
  - Run validated Latorre2018 simulation
  - Save geometry, stress, density histories
  - Document tolerances and metadata
- [ ] `test_latorre2018.py`:
  - Full simulation validation
  - Convergence to homeostasis
  - Sensitivity to solver tolerance
- [ ] `test_perturbations.py`:
  - Pressure perturbation
  - Axial stretch perturbation

**Target:** 10+ E2E tests, full simulation coverage

### Phase 5: CI/CD Integration (Week 5)

**Goal:** Automate testing

- [ ] Create `.github/workflows/tests.yml`
- [ ] Set up Codecov integration
- [ ] Add coverage badge to README
- [ ] Configure test reports in GitHub Actions
- [ ] Set up branch protection (require tests to pass)

### Phase 6: Documentation & Refinement (Week 6)

**Goal:** Finalize and document

- [ ] Add docstrings to all test files
- [ ] Create `tests/README.md` with:
  - How to run tests
  - How to add new tests
  - Interpreting test failures
- [ ] Review and refactor test utilities
- [ ] Add property-based tests (hypothesis) for edge cases
- [ ] Performance benchmarking tests (optional)

---

## Best Practices Summary

### Writing Tests

1. **One concept per test:** Each test should verify one specific behavior
2. **Descriptive names:** `test_survival_decreases_with_cohort_age()` not `test_survival()`
3. **Arrange-Act-Assert pattern:**
   ```python
   def test_example():
       # Arrange: Set up test data
       constituent = create_test_constituent()
       
       # Act: Execute the code under test
       result = constituent.compute_rhoR_alpha(timestep=5)
       
       # Assert: Verify the result
       assert result > 0
   ```
4. **Use fixtures for setup:** Avoid duplication
5. **Test edge cases:** Empty lists, zero values, boundary conditions
6. **Test error handling:** Verify exceptions are raised appropriately

### Maintaining Tests

1. **Keep tests independent:** No test should depend on another test's state
2. **Fast feedback:** Unit tests should run quickly
3. **Clear failure messages:** Use custom messages in assertions
4. **Regular updates:** Update tests when code changes
5. **Review coverage:** Aim for high coverage but prioritize meaningful tests

### Testing Numerical Code

1. **Use appropriate tolerances:** `np.testing.assert_allclose(rtol=1e-5)`
2. **Test against analytical solutions** when possible
3. **Test convergence properties:** Verify error decreases with finer discretization
4. **Test physical constraints:** Non-negativity, conservation laws
5. **Document expected behavior:** Why should this value be what it is?

---

## References & Resources

### Pytest Documentation
- [Pytest Official Docs](https://docs.pytest.org/)
- [Pytest Fixtures](https://docs.pytest.org/en/latest/fixture.html)
- [Parametrize](https://docs.pytest.org/en/latest/parametrize.html)

### Coverage Tools
- [Coverage.py](https://coverage.readthedocs.io/)
- [Pytest-cov](https://pytest-cov.readthedocs.io/)
- [Codecov](https://docs.codecov.io/)

### Testing Best Practices
- [Testing Scientific Code](https://www.software.ac.uk/guide/testing-your-software)
- [NumPy Testing Guidelines](https://numpy.org/doc/stable/reference/testing.html)

---

## Appendix: Example Test Files

### Example 1: Unit Test for Integrators

```python
"""Unit tests for numerical integrators."""
import pytest
import numpy as np
from src.integrators import (
    TrapezoidIntegrator, 
    BackwardTrapezoidIntegrator,
    SimpsonIntegrator,
    BackwardSimpsonIntegrator
)

class TestTrapezoidIntegrator:
    """Test naive trapezoidal integration."""
    
    def test_constant_function(self):
        """Integral of constant function f(x)=c should be c*(b-a)."""
        dt = 0.1
        f = [5.0] * 10  # Constant function f=5
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=9)
        
        result = integrator.integrate(f)
        expected = 5.0 * 0.9  # c * (9-0) * dt
        
        assert abs(result - expected) < 1e-10
    
    def test_linear_function(self):
        """Integral of f(x)=x should be x²/2."""
        dt = 1.0
        f = [0.0, 1.0, 2.0, 3.0, 4.0]  # f(x) = x
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=4)
        
        result = integrator.integrate(f)
        expected = 8.0  # ∫₀⁴ x dx = 4²/2 = 8
        
        assert abs(result - expected) < 1e-6
    
    @pytest.mark.parametrize("n_points,expected_error", [
        (10, 1e-3),
        (100, 1e-5),
        (1000, 1e-7),
    ])
    def test_convergence_rate(self, n_points, expected_error):
        """Test that error decreases with finer discretization."""
        # Integrate sin(x) from 0 to π
        x = np.linspace(0, np.pi, n_points)
        f = np.sin(x)
        dt = x[1] - x[0]
        
        integrator = TrapezoidIntegrator(dt=dt, start=0, stop=n_points-1)
        result = integrator.integrate(f.tolist())
        
        expected = 2.0  # ∫₀^π sin(x) dx = 2
        error = abs(result - expected)
        
        assert error < expected_error
```

### Example 2: Integration Test for Kinetics

```python
"""Integration tests for kinetics computations."""
import pytest
import numpy as np
from src.kinetics import Kinetics
from src.constituent import Constituent
from src.kinetics_interface import ConstituentKineticsContext
from tests.utils.mock_objects import create_test_constituent

class TestSurvivalFunctionComputation:
    """Test survival function q(s,tau) computation."""
    
    @pytest.fixture
    def simple_constituent(self):
        """Create a minimal constituent for testing."""
        return create_test_constituent(
            k_alpha_values=[0.01] * 10,  # Constant degradation
            timesteps=10
        )
    
    def test_survival_decreases_monotonically(self, simple_constituent):
        """Survival function should decrease with cohort age."""
        context = ConstituentKineticsContext(simple_constituent)
        kinetics = simple_constituent.kinetics
        
        current_timestep = 5
        q_values = kinetics.compute_survival_history(context, current_timestep)
        
        # q(s,s) = 1, and q should decrease for older cohorts
        assert q_values[-1] == pytest.approx(1.0)  # Most recent cohort
        assert all(q_values[i] <= q_values[i+1] for i in range(len(q_values)-1))
    
    def test_survival_with_known_solution(self):
        """Test against analytical solution for constant k_alpha."""
        k_const = 0.1  # 1/day
        dt = 1.0  # day
        s = 10  # Current time
        tau = 5  # Deposition time
        
        # Analytical: q(s,tau) = exp(-k * (s-tau))
        expected = np.exp(-k_const * (s - tau))
        
        # Create constituent with constant k_alpha
        from tests.utils.mock