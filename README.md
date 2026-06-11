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

| Simulator | Physics | MATLAB | Python | Verification | Formulation |
|------------|----------|---------|---------|------------|------------|
| 01_H_Simulator_FVM_1D_Incompressible_Homogeneous | Hydraulic | ✓ | ✓ | ✓ | ✓ |
| 02_H_Simulator_FVM_1D_Incompressible_Heterogeneous | Hydraulic | ✓ | ✓ | ✓ | ✓ |
| 03_H_Simulator_FVM_1D_Compressible_Homogeneous | Hydraulic | ✓ | ✓ | ✓ | ✓ |
| 04_H_Simulator_FVM_1D_Compressible_Heterogeneous | Hydraulic | ✓ | ✓ |  | ✓ |
| 05_H_Simulator_FVM_1D_Neumann_Wells | Hydraulic | ✓ | ✓ |  | ✓ |
| 06_H_Simulator_FVM_1D_CompressibleRock | Hydraulic |  |  |  | ✓ |
| 07_H_Simulator_FVM_1D_2Phase_Immiscible | Hydraulic |  |  |  |  |

## Recommended Learning Order

1. 01_H_FVM_1D_Incompressible_Homogeneous
2. 02_H_FVM_1D_Incompressible_Heterogeneous
3. 03_H_FVM_1D_Compressible_Homogeneous
4. 04_H_FVM_1D_Compressible_Heterogeneous
5. 05_H_FVM_1D_Neumann_Wells
6. 06_H_Simulator_FVM_1D_CompressibleRock
7. 07_H_Simulator_FVM_1D_2Phase_Immiscible

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

## 03_H_Simulator_FVM_1D_Compressible_Homogeneous

**Directory**

`H_Simulator_FVM_1D/03_H_Simulator_FVM_1D_Compressible_Homogeneous`

**Description**

A one-dimensional hydraulic flow simulator based on the Finite Volume Method (FVM).

**Assumptions**

- Single-phase fluid
- Slightly compressible fluid
- Homogeneous permeability
- Constant permeability
- No gravity
- Fixed-pressure (Dirichlet) boundary conditions

---

## 04_H_Simulator_FVM_1D_Compressible_Heterogeneous

**Directory**

`H_Simulator_FVM_1D/04_H_Simulator_FVM_1D_Compressible_Heterogeneous`

**Description**

A one-dimensional hydraulic flow simulator based on the Finite Volume Method (FVM) with spatially varying permeability and porosity.

**Assumptions**

- Single-phase fluid
- Slightly compressible fluid
- Heterogeneous permeability and porosity
- Constant permeability within each cell
- No gravity
- Fixed-pressure (Dirichlet) boundary conditions

---

## 05_H_Simulator_FVM_1D_Neumann_Wells

**Directory**

`H_Simulator_FVM_1D/05_H_Simulator_FVM_1D_Neumann_Wells`

**Description**

A one-dimensional hydraulic flow simulator based on the Finite Volume Method (FVM) with spatially varying permeability and porosity and including source terms (wells).

**Assumptions**

- Single-phase fluid
- Slightly compressible fluid
- Heterogeneous permeability and porosity
- Constant permeability within each cell
- No gravity
- No-flow (Neumann) boundary conditions
- Source terms in domain (wells)

---

## 06_H_Simulator_FVM_1D_CompressibleRock

**Directory**

`H_Simulator_FVM_1D/06_H_Simulator_FVM_1D_CompressibleRock`

**Description**

A one-dimensional hydraulic flow simulator based on the Finite Volume Method (FVM) with spatially varying permeability and porosity which depend on fluid pressure and including source terms (wells).

**Assumptions**

- Single-phase fluid
- Slightly compressible fluid
- Heterogeneous permeability and porosity which depends on fluid pressure
- Constant permeability within each cell
- No gravity
- No-flow (Neumann) boundary conditions
- Source terms in domain (wells)

---

## 07_H_Simulator_FVM_1D_2Phase_Immiscible

**Directory**

`H_Simulator_FVM_1D/07_H_Simulator_FVM_1D_2Phase_Immiscible`

**Description**

A one-dimensional two-phase hydraulic flow simulator based on the Finite Volume Method (FVM) with spatially varying permeability and porosity which depend on fluid pressure and including source terms (wells).

**Assumptions**

- Two-phase immiscible fluids
- Slightly compressible fluids
- Heterogeneous permeability and porosity which depends on fluid pressure
- Constant permeability within each cell
- No gravity
- No-flow (Neumann) boundary conditions
- Source terms in domain (wells)

---

## Planned Future Simulators

### Hydraulic (H)

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