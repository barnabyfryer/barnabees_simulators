# - About
# This reservoir simulator uses a FVM formulation to solve the continuity of mass
# balance equation in 2-D. The boundary conditions are fixed pressure
# at the edges (Dirichlet). It uses a pressure dependent, heterogeneous, anisotropic permeability.
# There is no gravity and the simulator is single phase. The fluid is
# considered to be slightly compressible.
# That flow model is coupled to a FEM mechanics model which is elastic and assumes small strains.
# The coupling is two-way and iterative.

#Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np
import copy

from src.Plotting_file import Plotting_file
from src.src_mech.input_mech_data import input_mech_data
from src.src_flow.input_flow_data import input_flow_data
from src.src_flow.FIMPressure_2D_1Phase import FIMPressure_2D_1Phase
from src.src_mech.M_Simulator_FEM_2D import M_Simulator_FEM_FVM_2D

# =============================================================================
# Inputs
# =============================================================================

Flow, Gen, State, Storage, Wells = input_flow_data()
Gen, Plotting, Pos, State = input_mech_data(Gen,State)

# =============================================================================
# Run simulation
# =============================================================================

while State["t"] < Gen["tf"]:

    err = 1
    while err > Gen["tol_all"]:
        #Save state before iteration
        P0 = State["P"].copy()
        e0 = State["e_vol"].copy()
        #Solve for new pressure
        State = FIMPressure_2D_1Phase(Flow,Gen,State,Wells)
        #Solve for stresses and volumetric strain
        State = M_Simulator_FEM_FVM_2D(Gen, Pos, State)
        #Check relative change in pressure
        errP = np.max(np.abs(State["P"] - P0) / np.maximum(np.abs(State["P"]), 1.0) )
        #Check relative change in volumetric strain
        erre = np.max(np.abs(State["e_vol"] - e0) / np.maximum(np.abs(State["e_vol"]), 1e-12) )
        #Find max error to check for convergence of coupled system
        err = max(errP, erre)
        State["step"] = 1
        Storage["P"][State["step"], :] = State["P"]
        Storage["phi"][State["step"], :] = State["phi"]
        Storage["kx"][State["step"], :] = State["kx"]
        Storage["ky"][State["step"], :] = State["ky"]
        Storage["flux"][State["step"], :] = State["flux"]
        Storage["e_vol"][State["step"], :] = State["e_vol"]
        Plotting_file(Flow, Gen, Storage)


    #Update to new time
    State["t"] += Gen["tstep"]

    #Store results
    if np.any(np.isclose(State["t"], Storage["TStorage"])):
        State["step"] += 1
        Storage["P"][State["step"], :] = State["P"]
        Storage["phi"][State["step"], :] = State["phi"]
        Storage["kx"][State["step"], :] = State["kx"]
        Storage["ky"][State["step"], :] = State["ky"]
        Storage["flux"][State["step"], :] = State["flux"]
        Storage["e_vol"][State["step"], :] = State["e_vol"]

# =============================================================================
# Plotting
# =============================================================================

Plotting_file(Flow,Gen,Storage)