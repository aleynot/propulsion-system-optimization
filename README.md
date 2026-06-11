# Propulsion System Performance Optimization

MATLAB implementations of thermodynamic cycle analysis and metaheuristic
(Particle Swarm Optimization) optimization for air-breathing propulsion
systems, including turbojet, three-shaft turbofan, and reverse-Brayton
(precooled) engine cycles.

These scripts were developed as part of coursework in propulsion system
design and performance optimization, applying ideal-cycle thermodynamics
combined with PSO-based multi-variable, multi-regime optimization.

## Repository Structure

```
propulsion-system-optimization/
├── README.md
├── turbojet/
│   ├── turbojet_op.m          # Parametric sweep + PSO (OPR, TIT)
│   └── engine_perf_4.m        # Turbojet vs Ramjet vs Reverse-Brayton comparison
├── reverse_brayton/
│   └── ters_brayton.m         # 8-variable PSO, takeoff + hypersonic cruise
├── turbofan_3shaft/
│   ├── turbofan3.m            # 6-variable PSO, takeoff + cruise design
│   └── turbofan_opt.m         # 5-variable PSO, SFC minimization
└── controls/
    └── backstepping_output_feedback.m  # Nonlinear output-feedback controller
```

## Contents

### 1. Turbojet Cycle Optimization (`turbojet/turbojet_op.m`)

Ideal turbojet performance model as a function of overall pressure ratio
(OPR) and turbine inlet temperature (TIT), evaluated at M = 0.9, 9 km
altitude.

- **Method:** Full parametric sweep (OPR ∈ [2, 50], TIT ∈ [1500, 2100] K)
  combined with PSO for direct comparison.
- **Objectives:** Minimize thrust-specific fuel consumption (TSFC) and
  maximize net thrust under various thrust/TSFC constraints.
- **Outputs:** Performance maps (thrust and TSFC contours over the
  OPR–TIT design space), optimum design points for six scenarios (A1–A3,
  B1–B3), and a side-by-side comparison of parametric vs. PSO results.

### 2. Engine Cycle Comparison (`turbojet/engine_perf_4.m`)

Design-point performance comparison of three engine architectures —
turbojet, ramjet, and reverse-Brayton cycle — at sea-level static and
across a Mach sweep (0–8).

- Computes specific impulse and TSFC for each cycle as a function of
  flight Mach number.
- Includes a sensitivity study of reverse-Brayton takeoff thrust versus
  the cooling-section temperature drop.

### 3. Reverse-Brayton Cycle PSO Optimization (`reverse_brayton/ters_brayton.m`)

Multi-regime (takeoff + hypersonic cruise, M = 5.5, 25 km) optimization of
a reverse-Brayton (precooled) engine cycle using PSO.

- **Design variables (8):** preburner and afterburner exit temperatures,
  cooling-section ΔT, and turbine pressure ratio, each defined separately
  for the takeoff and cruise operating points.
- **Objective:** Minimize cruise specific fuel consumption (SFC) subject
  to penalty-based constraints on takeoff thrust (≥ 120 kN), cruise
  thrust (≥ 300 kN), and component temperature/pressure feasibility
  limits.
- **Outputs:** PSO convergence history for all 8 design variables, optimum
  design point, and a thrust comparison bar chart for the two operating
  regimes.

### 4. Three-Shaft Turbofan Optimization (`turbofan_3shaft/`)

Two complementary PSO-based studies of a three-shaft turbofan engine
(fan + LPC/IPC + HPC, each on its own spool), modeled with separate cold-
and hot-section gas properties.

- **`turbofan3.m`** — 6 design variables (fan polytropic efficiency,
  fan/LPC/HPC pressure ratios, takeoff TIT, bypass ratio), evaluated at
  takeoff (SL, M = 0) and cruise (11 km, M = 0.85), including a simplified
  engine weight model.
- **`turbofan_opt.m`** — 5 design variables (LPC, HPC, and fan pressure
  ratios, T4, bypass ratio), minimizing SFC at cruise (10 km, M = 0.85)
  subject to a net thrust constraint of ≥ 100 kN.

### 5. Nonlinear Output-Feedback Control (`controls/backstepping_output_feedback.m`)

A backstepping-based output-feedback controller for a third-order
nonlinear system, using a "dirty derivative" filter to estimate an
unmeasured state. Included as a supporting example of nonlinear control
design and numerical simulation (`ode23t`) used in related coursework.

## Methods Summary

| Script | Design Variables | Optimization Method | Objective |
|---|---|---|---|
| `turbojet_op.m` | 2 (OPR, TIT) | Parametric sweep + PSO | Min TSFC / Max thrust |
| `ters_brayton.m` | 8 (takeoff + cruise) | PSO | Min cruise SFC |
| `turbofan3.m` | 6 | PSO | Min SFC, weight-aware |
| `turbofan_opt.m` | 5 | PSO | Min SFC, thrust-constrained |

All cycle models use ideal-gas thermodynamic relations with cold/hot
section gas properties (γ, c_p) where applicable, and PSO implementations
use penalty-function methods to enforce thrust, temperature, and pressure
feasibility constraints.

## Requirements

- MATLAB (tested with base MATLAB; no additional toolboxes required —
  PSO is implemented from scratch)

## Usage

Each script is self-contained. Run directly in MATLAB:

```matlab
run('turbojet/turbojet_op.m')
```

Each script prints a summary table of optimization results to the console
and generates the corresponding performance/convergence plots.

## Author

Aleyna Demirci
Aeronautical and Space Engineer
