# - About
# This reservoir simulator uses a FVM formulation to solve the continuity of mass
# balance equation in 1-D. The boundary conditions are fixed pressure
# at the edges (Dirichlet). It uses a constant, heterogeneous permeability.
# There is no gravity and the simulator is single phase. The fluid is
# considered to be slightly compressible.

# Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np

from src.FIMPressure_1D_1Phase import FIMPressure_1D_1Phase
from src.input_data import input_data
from src.Plotting_file import Plotting_file

# =============================================================================
# Inputs
# =============================================================================

Flow, Gen, State, Storage, Wells = input_data()

# =============================================================================
# Run simulation
# =============================================================================

while State["t"] < Gen["tf"]:
    #Solve for new pressure
    State = FIMPressure_1D_1Phase(Flow,Gen,State,Wells)
    #Update to new time
    State["t"] += Gen["tstep"]

    #Store results
    if np.any(np.isclose(State["t"], Storage["TStorage"])):
        State["step"] += 1
        Storage["P"][State["step"], :] = State["P"]
        Storage["phi"][State["step"], :] = State["phi"]
        Storage["kx"][State["step"], :] = State["kx"]
        Storage["flux"][State["step"], :] = State["flux"]



# =============================================================================
# Plotting
# =============================================================================

Plotting_file(Flow,Gen,Storage)


