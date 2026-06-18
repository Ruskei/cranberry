# Cranberry

Cranberry is a from-scratch C++ implementation of a 3D electromagnetic particle-in-cell simulation for studying plasma dynamics and numerical electrodynamics. It implements Yee-grid FDTD field evolution, relativistic particle pushing, Esirkepov current deposition, electrostatic initialization with a multigrid-preconditioned conjugate gradient solver, absorbing boundary conditions, threaded updates, and VTK output for Paraview visualization.

The main demonstration is of beam-drivne plasma wakefield acceleration (PWFA), where a relativistic driver bunch of electrons excites a wake in a plasma and a witness bunch is placed in the accelerating phase of the resulting longitudinal electric field.

<img width="991" height="732" alt="pwfa_magnitude_collage" src="https://github.com/user-attachments/assets/49d36e9e-3b65-4018-9b46-576778dbe595" />

(Magnitude of E-field at 3 snapshots of PWFA simulation)

### Features
- 3D Yee-grid FDTD solver for electromagnetic field updates
- Relativistic particle pusher
- Esirkepov charge-conserving current deposition
- Higher order particle shapes and smoothing to reduce noise
- Electrostatic initialization using PCG Poisson solve
- Absorbing boundary conditions
- VTK output for Paraview visualization
- Parallelized field update, particle puhs, and current deposition with std::thread

### Main simulation: beam-driven PWFA

The demo initializes a neutral plasma with mobile electrons and the much heavier protons, then injects relativistic driver and witness bunches of electrons with a Lorentz factor of 100.

As the driver bunch moves through the plasma, it pushes plasma electrons out of its path, leaving a wake of alternating field behind it. The witness bunch is initialized behind the driver such that it's accelerated by the wake.

<img width="991" height="732" alt="pwfa_y_collage" src="https://github.com/user-attachments/assets/a62f71c2-800e-44dd-ab79-8b07c5314545" />

(Longitudinal E-field of PWFA simulation)

### Gauss' law verification

Esirkepov current deposition satisfies discrete charge continuity, and Yee FDTD preserves Gauss' law in its domain. However, the absorbing boundary conditions are not fully constraint-preserving, causing drift:

<img width="1440" height="1120" alt="gauss_errors" src="https://github.com/user-attachments/assets/edbe300d-6bab-4f1f-aa0a-53eccf37bdda" />

(Data taken from half the simulation timestep of the PWFA simulation shown initially). However this error stays relatively low (<5%), so its impact on the shown simulation is likely not significant.

### Numerical Cherenkov Instability

The most challening problem was numerical Cherenkov instability (NCI). When highly relativistic particles move through the grid, they are near the information transfer speed of FDTD (since it's local), causing an unphysical shockwave to form.

This was mitigated by using higher-order current and charge deposition, aka using a cubic shape for the particles and binomial smoothing for the fields. This leaves invariants intact, but significantly reduced NCI. Still, NCI is visible in the above snapshots, but to a reduced extent.

### Learnings

I worked on this project parallel to working through _Griffith's Introduction to Electrodynamics_ after wanting to visualize the fields I was solving problems about.

The FDTD algorithm presented here in large part came from _Understanding the Finite-Difference Time-Domain Method by John B. Schneider_ and various research papers for a deeper understanding of EMPIC code.

- Understanding how Maxwell's equations can be discretized and computed
- Writing an efficient electrostatic initializer
- How to match artifacts in visualizations to issues in implementation
- Applied a subset of C++ to writing a memory-intensive simulation, and how to use low level parallelism constructs

### Building

This was only tested on my MacOS M1 machine, but the results should be reproducible with:
- C++20 compiler
- CMake >= 3.16
- VTK

```bash
cmake -B build -S .
cmake --build build
./build/cranberry
```

results are written to out/ and can be previewed in the FOSS Paraview.

### Limitations

I wrote this project because I wanted to gain a deeper understanding of electrodynamics and vector calculus, and thus there are some limitations:

- Boundary conditions are simple and create unphysical artifacts near edges
- Performance and accuracy is worse than that of established PIC codes
- Inflexible code design
- No collisions or advanced corrections
