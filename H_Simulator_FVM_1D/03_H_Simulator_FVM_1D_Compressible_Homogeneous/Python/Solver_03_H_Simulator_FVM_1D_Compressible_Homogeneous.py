# - About
# This reservoir simulator uses a FVM formulation to solve the continuity of mass
# balance equation in 1-D. The boundary conditions are fixed pressure
# at the edges (Dirichlet). It uses a constant, homogeneous permeability.
# There is no gravity and the simulator is single phase. The fluid is
# considered to be slightly compressible.

# Barnaby Fryer, 2026

# =============================================================================
# Import libraries
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

from src.FIMPressure_1D_1Phase import FIMPressure_1D_1Phase
from src.input_data import input_data

# =============================================================================
# Inputs
# =============================================================================

Flow, Gen, State, Storage = input_data()

# =============================================================================
# Run simulation
# =============================================================================

while State["t"] < Gen["tf"]:

    State = FIMPressure_1D_1Phase(Flow,Gen,State)


    State["t"] += Gen["tstep"]

# =============================================================================
# Plotting
# =============================================================================

# =============================================================================
# Plotting pressure
# =============================================================================

fig, ax = plt.subplots()
# Figure size
fig.set_size_inches(6, 4.5)
# Plot
ax.plot(Gen["x"], State["P"]/1e6, 'k-', linewidth=1, label='Simulation')
#ax.plot(x, P_an/1e6, 'r--', linewidth=1, label='Analytical Soln.')
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
lgd = ax.legend(fontsize=7)
lgd.set_frame_on(False)
#fig.savefig('../Verification/Pp_Python.jpg',
            #dpi=300,
            #bbox_inches='tight')
plt.show()

