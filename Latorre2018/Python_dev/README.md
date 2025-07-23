
---

# Python G&R Simulation for Vessels

This repository contains a Python implementation of a **Growth and Remodeling (G&R)** model for blood vessels. It is adapted from a set of C++ routines (including `main()`, `update_time_step`, `find_iv_geom`, `find_tf_geom`, `find_equil_geom`, and various objective functions) that were originally written to study arterial adaptation under altered mechanical loading (pressure, flow, axial stretch).

## Table of Contents
1. [Overview](#overview)  
2. [Key Features](#key-features)  
3. [File/Function Descriptions](#filefunction-descriptions)  
4. [Usage](#usage)  
5. [Assumptions](#assumptions)  
6. [Dependencies](#dependencies)  
7. [Example References and Literature](#example-references-and-literature)

---

## Overview

Blood vessels undergo continuous remodeling in response to changes in mechanical 
loads (e.g., increased pressure or flow). This simulation implements a 
**mixture-based G&R model** in Python, capturing how wall geometry 
(radius, thickness) and constituent mass (e.g., elastin, collagen, smooth muscle) 
adapt over time.

Key functionalities include:

- **Initialization** of vessel geometry and material parameters.
- **Time-stepping** routines that update mass production/degradation, 
   solve for equilibrium geometry, and track changes in wall stress and 
   shear over many days of simulated adaptation.
- **Traction-free geometry** and **loaded geometry** calculations to isolate 
    unloaded vs. loaded configurations.
- **Equilibrium** solver for final mechanobiological steady-state.

---

## Key Features

1. **Evolution Loop**: The code increments time and recalculates 
    how vessel radius and thickness evolve (via `update_time_step_py`)
    for the time-dependent (transient) formulation of the G&R model.
2. **Multi-Constituent Mixture**:
   - Elastin, Muscle, Collagen families (3 constituents).
   - Tracks separate mass densities, production rates, and stress contributions.
3. **Root-Finding**:
   - **`brentq`** or **`toms748`** methods (1D) for traction-free or partial solves.
   - Multi-dimensional solvers (`method='hybr'`) for equilibrium geometry.
4. **Stress Analyses**:
   - Intramural stress (circumferential, axial).
   - Wall shear stress from simplified Poiseuille approximation.
5. **Pythonic**: High readability, straightforward function calls, 
     easily modified for play-ground testing of G&R formulations. 

---

## File/Function Descriptions

Below is a short overview of the key Python functions:

1. **`main()`**  
   - Sets up the `Vessel` object with initial geometry and model parameters.
   - Runs time-stepping for a specified number of days (`n_days`), logs geometry and mass changes to an output file.
   - Example usage:  
     ```bash
     python main.py
     ```
   - **Purpose**: Orchestrates the entire G&R simulation.

2. **`update_kinetics(curr_vessel)`**  
   - Computes updated mass production/degradation rates for each 
     constituent based on current stress and shear state.
   - Called each time step.

3. **`update_time_step_py(curr_vessel)`**  
   - Mimics the C++ `update_time_step`.  
   - Copies previous geometry and mass densities, updates kinetics, 
     solves for in vivo geometry repeatedly until mass changes converge.

4. **`update_sigma(curr_vessel)`**  
   - Calculates the intramural stresses given the current geometry and 
     mass distribution.  
   - In the mixture approach, each constituent’s stress contribution is summed.

5. **`find_iv_geom(curr_vessel)`**  
   - Finds the loaded (in vivo) geometry (e.g., radius) that satisfies 
     equilibrium given current mass and pressure/flow loads.  
   - Uses a root-finding method (`brentq` or another bracketed approach).

6. **`find_tf_geom(curr_vessel)`**  
   - Solves for the traction-free (unloaded) geometry (e.g., removing 
     external loads to see how the vessel shape changes).  
   - Typically uses a multi-dimensional root solver (`method='hybr'` or similar).

7. **`find_equil_geom(curr_vessel)`**  
   - Solves for a “mechanobiological” equilibrium given changes in pressure, flow, 
     etc.  
   - Uses `equil_obj_f_py` as an objective function to set up 4 equations for 
     final geometry and mass states.

8. **Objective Functions**  
   - **`tf_obj_f_py(vars, curr_vessel)`**: 2D system for traction-free geometry.  
   - **`equil_obj_f_py(vars, curr_vessel)`**: 4D system for final equilibrium 
     geometry and mass distribution.  

---

## Usage

1. **Clone or Download** this repository.
2. **Install Dependencies** (see [Dependencies](#dependencies)).
3. **Run the main code**:
   ```bash
   python main.py
   ```

4. **Check Outputs**:
   - Typical outputs might go to **`GnR_out.txt`**, containing lines of radius, 
     thickness, mass fractions, and a measure of equilibrium at each time step.
   - The code also prints basic console messages about solver convergence.

5. **Customizing**:
   - Modify `main()` to change `n_days`, or to apply different mechanical 
     perturbations (pressure ramp, flow step changes).
   - Add/alter constituents or the functional forms in `update_kinetics` if
     you have different G&R laws.

---

## Assumptions

1. **Homogeneous, Axisymmetric Vessel**:
   - The code currently assumes a cylindrical vessel with uniform properties
     in the circumferential and axial directions (aside from the specified 
     angles for collagen, muscle, etc.).
2. **Mixture Theory**:
   - Each constituent (elastin, muscle, collagen) is modeled with an 
     exponential-type stress-strain law. Summation leads to total stress.
3. **Limited Biological Complexity**:
   - No explicit modeling of endothelial cell function, microporosity, 
     or advanced biochemical signaling. Instead, stress/shear triggers 
     mass changes in muscle and collagen.
4. **Time Step**:
   - One day per step. If your application needs finer resolution, 
     adjust `curr_vessel.dt`.

---

## Dependencies

- **Python 3.7+** (should run on any modern Python).
- **NumPy** for array operations.
- **SciPy** for root-finding functions (`scipy.optimize.brentq`, 
  `scipy.optimize.root`, `scipy.optimize.toms748`).
- **math** (standard library) for exponentials, pi, etc.

Optionally:
- **matplotlib** if you want to add plotting of the results 
  (not strictly required by the code above, but useful for visualizing
  G&R behavior).

---

## Example References and Literature

1. **Latorre M, Humphrey JD.**  
   “A mechanobiologically equilibrated constrained mixture model for growth and remodeling
    of soft tissues.”  
   *ZAMM*, 98:2048-2071, 2018.

2. **Humphrey JD, Rajagopal KR.**  
   “A constrained mixture model for growth and remodeling of soft tissues.”  
   *Mathematical Models and Methods in Applied Sciences*, 12(03):407–30, 2002.

3. **Valentin A, Cardamone L, Baek S, Humphrey JD.**
    "Complementary vasoactivity and matrix remodelling in arterial adaptations to
    altered flow and pressure."
   *J. R. Soc. Interface*, 6:293-306, 2009.

This code recreates the figure 2 from reference 1. For greater background on the
transient formulation refer to reference 2 & 3. 

---
