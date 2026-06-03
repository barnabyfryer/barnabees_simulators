# barnabees_simulators

A collection of fluid flow, thermal, mechanical, and earthquake simulators implemented in MATLAB and Python.

The purpose of this repository is to develop a progression of numerical simulators, beginning with simple one-dimensional flow problems and advancing toward fully coupled Thermo-Hydro-Mechanical (THM) reservoir and geomechanical models.

---

## Repository Contents

This repository contains:

- MATLAB implementations
- Python implementations
- Mathematical formulations
- Verification and benchmark examples
- Progressive development from simple hydraulic models to coupled THM simulators

---

## Naming Convention

The simulator names follow the convention:

- **H** = Hydraulic (fluid flow)
- **HM** = Hydro-Mechanical
- **HT** = Hydro-Thermal
- **THM** = Thermo-Hydro-Mechanical

Additional descriptors indicate:

- Numerical method (e.g. FVM, FEM)
- Spatial dimension (e.g. 1D, 2D, 3D)
- Physical assumptions (e.g. incompressible, heterogeneous)
- Solution type (e.g. steady-state, transient)

Example:

`01_H_Simulator_FVM_1D_Incompressible_Homogeneous`

represents a:

- Hydraulic simulator
- Finite Volume Method
- One-dimensional domain
- Incompressible fluid
- Homogeneous permeability field

---

## Current Simulators

| Simulator | Physics | Numerical Method | MATLAB | Python |
|------------|----------|------------------|---------|---------|
| 01_H_Simulator_FVM_1D_Incompressible_Homogeneous | Hydraulic | Finite Volume Method | ✓ | ✓ |
| 02_H_Simulator_FVM_1D_Incompressible_Heterogeneous | Hydraulic | Finite Volume Method | ✓ | ✓ |

## Recommended Learning Order

1. 01_H_FVM_1D_Incompressible_Homogeneous
2. 02_H_FVM_1D_Incompressible_Heterogeneous

---

# Simulator Descriptions

## 01_H_Simulator_FVM_1D_Incompressible_Homogeneous

**Directory**

`H_Simulator_FVM_1D/01_H_Simulator_FVM_1D_Incompressible_Homogeneous`

**Description**

A one-dimensional steady-state hydraulic flow simulator based on the Finite Volume Method (FVM).

**Assumptions**

- Single-phase fluid
- Incompressible fluid
- Homogeneous permeability
- Constant permeability
- No gravity
- Fixed-pressure (Dirichlet) boundary conditions
- Steady-state formulation

---

## 02_H_Simulator_FVM_1D_Incompressible_Heterogeneous

**Directory**

`H_Simulator_FVM_1D/02_H_Simulator_FVM_1D_Incompressible_Heterogeneous`

**Description**

A one-dimensional steady-state hydraulic flow simulator based on the Finite Volume Method (FVM) with spatially varying permeability.

**Assumptions**

- Single-phase fluid
- Incompressible fluid
- Heterogeneous permeability
- Constant permeability within each cell
- No gravity
- Fixed-pressure (Dirichlet) boundary conditions
- Steady-state formulation

---

## Planned Future Simulators

### Hydraulic (H)

- 1D transient flow
- Compressible fluid flow
- Variable viscosity
- Gravity-driven flow
- Multiphase flow

### Hydro-Mechanical (HM)

- 1D poroelasticity
- Biot consolidation
- Coupled pressure-displacement formulations

### Hydro-Thermal (HT)

- Heat transport
- Advection-diffusion
- Thermally coupled flow

### Thermo-Hydro-Mechanical (THM)

- Fully coupled THM reservoir simulators
- Geothermal reservoir models
- Fault slip and earthquake simulations

---

## Author

**Barnaby Fryer**  
2026