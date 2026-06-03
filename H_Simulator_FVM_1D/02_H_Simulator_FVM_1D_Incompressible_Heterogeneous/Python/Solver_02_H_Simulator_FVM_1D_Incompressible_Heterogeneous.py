# - About
# This reservoir simulator uses a FVM formulation to solve the continuity of mass
# balance equation in 1-D. The boundary conditions are fixed pressure
# at the edges (Dirichlet). It uses a constant, heterogeneous permeability.
# There is no gravity and the simulator is single phase. The fluid is
# considered to be incompressible. The formulation is therefore steady
# state.

# Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Inputs
# =============================================================================

# Number of cells for discretization [-]
Nx = 100
# Viscosity [Pa sec]
mu = 1
# Permeability, of outer edges [m^2]
k_out = 1
# Permeability, of inner zone [m^2]
k_in = 5
# Permeability edge bands length [m]
L_k = 0.2
# Length [m]
Lx = 1
#Fixed boundary pressure on left [Pa]
PL = 10e6
# Fixed boundary pressure on right [Pa]
PR = 0

# =============================================================================
# Preparation
# =============================================================================

# Predefine memory for transmissivity at each cell edge
HTx = np.zeros(Nx + 1)
# Predefine memory for Jacobian
A = np.zeros((Nx, Nx))
# Define cell size
dx = Lx/Nx
# Define grid position
x = np.linspace(dx/2, Lx - dx/2, Nx)
# Define edge positions
x_edge = np.linspace(0, Lx, Nx+1)
# Define permeability of each cell
k = np.ones(Nx)*k_in
k[x < L_k] = k_out
k[x > Lx-L_k] = k_out
# Predefine memory for flow through each cell
q = np.zeros(Nx)
# Predefine memory for velocity across cell boundary
ux = np.zeros(Nx + 1)

# =============================================================================
# Build transmissivities
# =============================================================================

# Find harmonic average between cells for permeability
kHx = 2*k[:-1]*k[1:]/(k[:-1]+k[1:])
# Calculate transmissivity
HTx[1:-1] = kHx/(mu*dx)
# Deal with edge cell boundaries
HTx[0] = k[0]/(mu*dx/2)
HTx[Nx] = k[Nx-1]/(mu*dx/2)

# =============================================================================
# Build Jacobian
# =============================================================================

for j in range(Nx):
    # Add cell transmissivity
    if j > 0:
        A[j,j] += HTx[j]
        A[j,j-1] -= HTx[j]
    # Deal with first cell
    else:
        A[j,j] += HTx[j]
        q[j] += PL*HTx[j]
    # Add cell transmissivity
    if j < Nx-1:
        A[j, j] += HTx[j+1]
        A[j, j+1] -= HTx[j+1]
    # Deal with last cell
    else:
        A[j, j] += HTx[j+1]
        q[j] += PR*HTx[j+1]

# =============================================================================
# Solve for presure in each cell
# =============================================================================

P = np.linalg.solve(A, q)

# =============================================================================
# Find velocity between each cell
# =============================================================================

# Velocity in middle interfaces
ux[1:-1] = (P[:-1] - P[1:]) * HTx[1:-1]
# Find velocity at edges
ux[0] = (PL - P[0]) * HTx[0]
ux[-1] = (P[-1] - PR) * HTx[-1]

# =============================================================================
# Plotting pressure
# =============================================================================

fig, ax = plt.subplots()
# Figure size
fig.set_size_inches(6, 4.5)
# Plot
ax.plot(x, P/1e6, 'k-', linewidth=1)
# Labels
ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
ax.set_ylabel(r'Pressure, $P$ [MPa]', fontsize=10)
# Font size
ax.tick_params(labelsize=7)
# Tick direction
ax.tick_params(direction='out')
# Box off
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

# =============================================================================
# Plotting permeability
# =============================================================================

fig, ax = plt.subplots()
# Figure size
fig.set_size_inches(6, 4.5)
# Plot
ax.plot(x, k, 'k-', linewidth=1)
# Labels
ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
ax.set_ylabel(r'Permeability, $k$ [m$^2$]', fontsize=10)
# Font size
ax.tick_params(labelsize=7)
# Tick direction
ax.tick_params(direction='out')
# Box off
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()

# =============================================================================
# Plotting velocity
# =============================================================================

fig, ax = plt.subplots()
# Figure size
fig.set_size_inches(6, 4.5)
# Plot
ax.plot(x_edge, ux, 'k-', linewidth=1)
# Labels
ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
ax.set_ylabel(r'Darcy Flux [m/s]', fontsize=10)
# Font size
ax.tick_params(labelsize=7)
# Tick direction
ax.tick_params(direction='out')
# Box off
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.show()