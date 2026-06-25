# - About
# Elastic FEM Mechanics Simulator For Quadrilatoral Elements Only 2-D

# Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np

from src.src_mech.build_force import build_force
from src.src_mech.global_stiffness_builder import global_stiffness_builder
from src.src_mech.fix_x import fix_x
from src.src_mech.fix_y import fix_y
from src.src_mech.solve_stresses import solve_stresses
from src.src_mech.add_pressure import add_pressure

def M_Simulator_FEM_FVM_2D(Gen,Pos,State):

    # =============================================================================
    # Force vector, determine where to apply body forces here
    # =============================================================================

    f = build_force(Gen)

    # =============================================================================
    # Force vector, determine where to pore-pressure-induced forces
    # =============================================================================

    f = add_pressure(f, Gen, State)

    # =============================================================================
    # Build stiffness matrix
    # =============================================================================

    k_global = global_stiffness_builder(Gen)

    # =============================================================================
    # Boundary conditions
    # =============================================================================
    x_fixes = np.concatenate((Gen["nodes_left"],Gen["nodes_right"]),axis=0)
    y_fixes = np.concatenate((Gen["nodes_bottom"], Gen["nodes_top"]),axis=0)
    #Apply fixes to global matrix
    f, k_global = fix_x(f, k_global, x_fixes)
    f, k_global = fix_y(f, k_global, y_fixes)

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
    Pos["du"] = Pos["u"].reshape((Gen["Nx"] + 1, Gen["Ny"] + 1))
    Pos["dv"] = Pos["v"].reshape((Gen["Nx"] + 1, Gen["Ny"] + 1))

    # New node locations
    Pos["x_new"] = Pos["x"] + Pos["du"]
    Pos["y_new"] = Pos["y"] + Pos["dv"]

    # =============================================================================
    # Solve for stresses and volumetric strain
    # =============================================================================

    Sigma, State["e_vol"] = solve_stresses(Gen, Pos, State)

    #Reshape to vector
    State["e_vol"] = State["e_vol"].reshape(-1, order="F")

    return State

