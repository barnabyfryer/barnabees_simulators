# - About
# Elastic FEM Mechanics Simulator For Quadrilatoral Elements Only 2-D

# Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np

from src.build_force import build_force
from src.input_data import input_data
from src.global_stiffness_builder import global_stiffness_builder
from src.fix_x import fix_x
from src.fix_y import fix_y
from src.solve_stresses import solve_stresses
from src.plotting_sim import plotting_sim

# =============================================================================
# Inputs
# =============================================================================

Gen, Plotting, Pos = input_data()

# =============================================================================
# Force vector, determine where to apply forces here
# =============================================================================

f = build_force(Gen)

# =============================================================================
# Build stiffness matrix
# =============================================================================

k_global = global_stiffness_builder(Gen)

# =============================================================================
# Boundary conditions
# =============================================================================

#Apply fixes to global matrix
f, k_global = fix_x(f, k_global, Gen["nodes_left"])
f, k_global = fix_y(f, k_global, Gen["nodes_bottom"])

# =============================================================================
# Solve
# =============================================================================

Pos["dx"] = np.linalg.solve(k_global, f)

# =============================================================================
# Process displacements
# =============================================================================

# Displacement in x-direction (u)
Pos["u"] = Pos["dx"][0::2]

# Displacement in y-direction (v)
Pos["v"] = Pos["dx"][1::2]

# Reshape to grid form
Pos["du"] = Pos["u"].reshape((Gen["Ny"] + 1, Gen["Nx"] + 1))
Pos["dv"] = Pos["v"].reshape((Gen["Ny"] + 1, Gen["Nx"] + 1))

# New node locations
Pos["x_new"] = Pos["x"] + Pos["du"]
Pos["y_new"] = Pos["y"] + Pos["dv"]

# =============================================================================
# Solve for stresses and volumetric strain
# =============================================================================

Sigma, e_vol = solve_stresses(Gen, Pos)

# =============================================================================
# Plotting
# =============================================================================

plotting_sim(f, Gen, Plotting, Pos, Sigma, e_vol)

