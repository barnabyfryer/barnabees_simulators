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
from src.Plotting_file import Plotting_file
from src.src_mech.input_mech_data import input_mech_data
from src.src_flow.input_flow_data import input_flow_data
from src.src_flow.FIMPressure_2D_1Phase import FIMPressure_2D_1Phase

# =============================================================================
# Inputs
# =============================================================================

Flow, Gen, State, Storage, Wells = input_flow_data()
#Gen, Pos = input_mech_data(Gen)


# =============================================================================
# Run simulation
# =============================================================================

while State["t"] < Gen["tf"]:
    #Solve for new pressure
    State = FIMPressure_2D_1Phase(Flow,Gen,State,Wells)
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

# =============================================================================
# Plotting
# =============================================================================

Plotting_file(Flow,Gen,Storage)