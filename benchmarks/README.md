# Benchmarks & Profiling

This directory contains tools for performance analysis and profiling of svGrowth simulations.

## Quick Start

```bash
# Profile a complete simulation
python benchmarks/profile_simulation.py

# Profile and open results in browser
python benchmarks/profile_simulation.py --browser

# Profile specific configuration
python benchmarks/profile_simulation.py --config latorre2018_test.yaml

# Profile constituent-level operations
python benchmarks/profile_constituents.py
```

## Available Tools

### `profile_simulation.py`
Profiles a complete simulation run using PyInstrument.

**Features:**
- Interactive HTML flame graphs
- Text-based console output
- Automatic result archiving with timestamps

**Output:**
- `results/profile_<config>_<timestamp>.html` - Interactive visualization
- `results/profile_<config>_<timestamp>.txt` - Text report

### `profile_constituents.py`
Deep-dive profiling of constituent-level calculations.

**Focuses on:**
- Stress computation (`compute_sigma_alpha`)
- Mass density calculation (`compute_rhoR_alpha`)
- Survival function evaluation

## Installation

Install profiling dependencies:

```bash
pip install -r requirements-dev.txt
```

Or manually:
```bash
pip install pyinstrument line-profiler memory-profiler snakeviz
```

## Interpreting Results

### PyInstrument Output

```
  _     ._   __/__   _ _  _  _ _/_   Recorded: 14:32:18  Duration: 45.234s
 /_//_/// /_\ / //_// / //_'/ //    
/   _/                      v4.6.0

45.234 main  src/main.py:1
└─ 42.891 Simulation.run  simulation.py:45
   ├─ 28.456 Layer.solve_equilibrium  layer.py:234
   │  └─ 25.123 SingleConstituent.compute_sigma_alpha  constituent.py:234  ← HOTSPOT
   └─ 8.901 integrate  integrators.py:145
```

**Key columns:**
- **Left number**: Total time spent in function (seconds)
- **Function name**: Which function is running
- **File:line**: Where it's defined
- **Tree structure**: Call hierarchy (indent = nested call)

**Look for:**
- 🔴 Functions taking >30% of total time
- 🟡 Unexpected function calls (shouldn't be there)
- 🟢 Deep nesting (opportunity to optimize call stack)

### Performance Optimization Workflow

1. **Run profile_simulation.py** → Identify top 5 slowest functions
2. **Run profile_constituents.py** → Deep-dive into slow functions
3. **Optimize** → Vectorize, cache, or use Numba
4. **Re-profile** → Verify improvement

## Results Directory

All profiling results are saved to `results/` with timestamps:

```
results/
├── profile_latorre2018_20260108_143015.html
├── profile_latorre2018_20260108_143015.txt
├── profile_constituents_20260108_144523.html
└── ...
```

**Note:** Results are gitignored by default (large files, local performance data).

## Advanced Usage

### Profile with Custom Config

```bash
python benchmarks/profile_simulation.py --config my_custom_config.yaml
```

### Console-Only Output

```bash
python benchmarks/profile_simulation.py --no-html
```

### Automated Benchmarking

See `benchmark_suite.py` for tracking performance over time.

## Troubleshooting

### ImportError: No module named 'pyinstrument'

```bash
pip install pyinstrument
```

### Profiling Hangs or Crashes

- Try profiling smaller simulations first
- Check memory usage (`top` or Activity Monitor)
- Reduce number of timesteps in config file

### No Output / Empty Results

- Ensure `src/main.py` runs successfully without profiling
- Check that simulation actually runs (not just imports)

## Next Steps

After profiling:

1. **Identify bottlenecks** (functions taking >30% time)
2. **Optimize hot functions**:
   - Vectorize with NumPy
   - Add Numba `@jit` decorators
   - Cache repeated calculations
3. **Re-profile and compare**
4. **Track improvements over time** with benchmark suite
