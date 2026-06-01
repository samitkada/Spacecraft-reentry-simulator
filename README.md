# Atmospheric Re-Entry Trajectory Simulator

A high-fidelity MATLAB simulation of an Orion-class spacecraft performing atmospheric re-entry from low Earth orbit. The model predicts trajectory evolution, aerodynamic heating, dynamic pressure, and G-loading while visualizing the complete re-entry sequence through animated 3D flight paths.

## Overview

This project models the physics of atmospheric entry using coupled nonlinear equations of motion and aerodynamic heating correlations commonly used in preliminary spacecraft design studies.

The simulator tracks a capsule from **120 km entry interface altitude at Mach 32** through atmospheric deceleration and parachute deployment conditions while computing:

* Re-entry trajectory and downrange distance
* Velocity and altitude profiles
* Dynamic pressure history
* Aerodynamic deceleration (G-load)
* Stagnation-point heating rates
* Peak thermal loading conditions
* 3D animated Earth-relative flight visualization

The project was developed to strengthen aerospace engineering skills in:

* Flight mechanics
* Numerical simulation
* Orbital and atmospheric dynamics
* Aerothermodynamics
* MATLAB modeling and visualization

---

## Physics & Engineering Models

### Atmospheric Model

An exponential atmosphere approximation is used:

[
\rho = \rho_0 e^{-h/H}
]

where atmospheric density decreases with altitude according to a scale height model.

### Vehicle Dynamics

The simulation numerically integrates the full nonlinear re-entry equations of motion:

[
\dot r = V\sin\gamma
]

[
\dot\theta = \frac{V\cos\gamma}{r}
]

[
\dot V = -\frac{D}{m}-g\sin\gamma
]

[
\dot\gamma = \frac{L}{mV}+\left(\frac{V}{r}-\frac{g}{V}\right)\cos\gamma
]

including:

* Variable gravity
* Aerodynamic drag
* Lift generation
* Earth curvature effects

### Aerodynamic Heating

Convective heating is estimated using the Sutton-Graves stagnation point correlation:

[
\dot q = k\sqrt{\frac{\rho}{R_n}}V^3
]

allowing prediction of peak thermal loads experienced by the heat shield during entry.

---

## Example Results

For an Orion-like capsule:

| Parameter         | Value     |
| ----------------- | --------- |
| Entry Altitude    | 120 km    |
| Entry Velocity    | 10.8 km/s |
| Entry Mach Number | ~32       |
| Flight Path Angle | -6.5°     |
| Vehicle Mass      | 9,300 kg  |
| Diameter          | 5.0 m     |

The simulation outputs:

* Peak heating rate location
* Maximum deceleration loads
* Total downrange distance
* Flight time to parachute deployment altitude
* Temperature-colored re-entry trajectory

---

## Visualization Features

### Mission Analysis Dashboard

The simulator automatically generates a six-panel engineering dashboard displaying:

* Entry corridor
* Heating rate history
* G-load history
* Altitude profile
* Dynamic pressure profile
* Velocity profile

### 3D Re-Entry Animation

A real-time animated visualization includes:

* Earth rendering
* Atmospheric shell visualization
* Heating-colored trajectory path
* Dynamic spacecraft glow intensity based on heating rate
* Live flight telemetry display

This animation can be exported directly to MP4 for engineering portfolios, presentations, or LinkedIn demonstrations.

---

## Technical Highlights

✔ Nonlinear atmospheric entry simulation

✔ ODE45 numerical integration

✔ Event-based parachute deployment termination

✔ Sutton-Graves heating model

✔ Dynamic pressure and G-load analysis

✔ Publication-quality MATLAB visualizations

✔ 3D mission animation and export capability

---

## Future Improvements

Planned upgrades include:

* Standard Atmosphere 1976 implementation
* Six-degree-of-freedom vehicle dynamics
* Bank-angle guidance algorithms
* Skip-entry trajectories
* Ablative heat shield mass loss modeling
* Monte Carlo uncertainty analysis
* Lunar-return and Mars-return mission profiles

---

### Engineering Relevance

Atmospheric entry remains one of the most challenging phases of spacecraft design due to the simultaneous interaction of extreme aerodynamic, thermal, and structural loads. This project demonstrates the application of numerical methods and aerospace engineering principles to analyze those conditions and evaluate vehicle performance during hypersonic flight.

