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
import matplotlib.pyplot as plt

from src.Plotting_file import Plotting_file
from src.src_mech.input_mech_data import input_mech_data
from src.src_flow.input_flow_data import input_flow_data
from src.src_flow.FIMPressure_2D_1Phase import FIMPressure_2D_1Phase
from src.src_mech.M_Simulator_FEM_2D import M_Simulator_FEM_2D
from src.Mandel.Mandel_comp import Mandel_comp

# =============================================================================
# Inputs
# =============================================================================

Flow, Gen, Plotting, State, Storage, Wells = input_flow_data()
Gen, Pos, State, Storage = input_mech_data(Flow, Gen,State, Storage)

# =============================================================================
# Run simulation
# =============================================================================

while State["t"] < Gen["tf"]:

    if State["t"] >= Gen["t_applied"]:
        Gen["F"] = Gen["F_2"]

    #Initialize error
    err = 1
    #Initialize number of iterations
    it = 0
    # Store the state variable from previous time step for use in residual
    State0 = {"P": State["P"].copy(),"e_vol": State["e_vol"].copy()}
    # Store the state variable from previous iteration for density update
    State_phi = State.copy()
    while err > Gen["tol_all"]:
        #Save state before iteration
        P0 = State["P"].copy()
        e0 = State["e_vol"].copy()
        #Solve for new pressure
        State = FIMPressure_2D_1Phase(Flow,Gen,State,State0,State_phi,Wells)
        # Store the state variable from previous iteration for density update
        State_phi = State.copy()
        #Solve for stresses and volumetric strain
        State = M_Simulator_FEM_2D(Gen, Pos, State)
        #Check relative change in pressure
        dP = State["P"] - P0
        State["errP"] = np.abs(dP) / np.maximum(np.abs(State["P"]), 1.0)
        #Check relative change in volumetric strain
        State["erre"] = np.abs(State["e_vol"] - e0) / np.maximum(np.abs(State["e_vol"]), 1e-12)
        #Find max error to check for convergence of coupled system
        err = max(max(State["errP"]), max(State["erre"]))
        it += 1

    #Update to new time
    State["t"] += Gen["tstep"]
    #Print time
    print(f"\rIterations = {it} | t = {State['t']:.3f} s", flush=True)

    #Store results
    if np.any(np.isclose(State["t"], Storage["TStorage"])):
        State["step"] += 1
        Storage["P"][State["step"], :] = State["P"]
        Storage["phi"][State["step"], :] = State["phi"]
        Storage["kx"][State["step"], :] = State["kx"]
        Storage["ky"][State["step"], :] = State["ky"]
        Storage["flux"][State["step"], :] = State["flux"]
        Storage["e_vol"][State["step"], :] = State["e_vol"]
        Storage["errP"][State["step"], :] = State["errP"]
        Storage["erre"][State["step"], :] = State["erre"]
        Storage["fx"][State["step"], :] = State["fx"]
        Storage["fy"][State["step"], :] = State["fy"]
        Storage["s_xx"][State["step"], :] = State["s_xx"]
        Storage["s_yy"][State["step"], :] = State["s_yy"]
        Storage["s_xy"][State["step"], :] = State["s_xy"]
        Storage["e_xx"][State["step"], :] = State["e_xx"]
        Storage["e_yy"][State["step"], :] = State["e_yy"]
        Storage["gamma_xy"][State["step"], :] = State["gamma_xy"]

# =============================================================================
# Plotting
# =============================================================================

Plotting_file(Gen,Pos,Storage)

# =============================================================================
# Plotting Mandel
# =============================================================================
if Gen["Mandel"] == 1:
    #Drained Young's Modulus
    E = Gen["E"]
    #Drained poisson's ratio
    vd = Gen["nu"]
    #Size of quarter domain
    a = Gen["Lx"]
    b = Gen["Ly"]
    #Fluid compressibility
    cf = Flow["cf"]
    #Permeability
    k = Flow["kx0"][0]
    #Fluid viscosity
    muf = Flow["muf"]
    #Porosity
    phi = Flow["phi0"][0]
    #Grain bulk modulus
    Ks = Flow["ks"]
    #Load applied [Pa]
    F = Gen["F"]

    #Dimensionless x for simulator
    x_d = Storage["x"]/a

    # index of the row closest to the domain centre
    I = np.argmin(np.abs(Storage["y"] - Gen["Ly"] / 2))
    # corresponding y-coordinate
    y_mid = Storage["y"][I]
    # all cells lying on that y-coordinate
    idx = np.where(Storage["y"] == y_mid)[0]


    # =============================================================================
    # Plotting pressure Mandel
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    for i in range(7, len(Storage["TStorage"]), int(round(len(Storage["TStorage"])/7))):

        t = Storage["TStorage"][i] - Gen["t_applied"]
        Mandel = Mandel_comp(a,b,cf,E,F,k,Ks,muf,phi,t,vd)

        ax.plot(x_d[idx], (Storage["P"][i, idx] - Storage["P"][0, idx])/Mandel["p0"], 'k-', linewidth=1, label='Simulation' if i == 7 else None)

        ax.plot(Mandel["x"], Mandel["Pd"], 'r--', linewidth=1, label='Analytical Soln.' if i == 7 else None)

    # Labels
    ax.set_xlabel(r'Distance, $x/a$', fontsize=10)
    ax.set_ylabel(r'Pressure, $(P-P(t=0))/P_0$', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/Mandel_python.jpg',
    dpi=300,
    bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting pressure in time
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot((Storage["TStorage"] - Gen["t_applied"])/Mandel["tnorm"], Storage["P"][:, idx[0]] / (Mandel["p0"] + Gen["Pi"]), 'k-', linewidth=1)
    # Labels
    ax.set_xlabel(r'Time, $t \cdot c_v/a^2$', fontsize=10)
    ax.set_ylabel(r'Pressure, $P/(P_0+P(t=0))$', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, 1)
    fig.savefig('../Verification/Mandel_pressure_python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Show match to Mandel solution, displacement x
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    for i in range(7, len(Storage["TStorage"]), int(round(len(Storage["TStorage"]) / 7))):
        t = Storage["TStorage"][i] - Gen["t_applied"]
        Mandel = Mandel_comp(a, b, cf, E, F, k, Ks, muf, phi, t, vd)

        ax.plot(x_d[idx], np.cumsum(Storage["e_xx"][i,idx] * Gen["dx"]), 'k-', linewidth=1,
                label='Simulation' if i == 7 else None)

        ax.plot(Mandel["x"], Mandel["ux"], 'r--', linewidth=1, label='Analytical Soln.' if i == 7 else None)

    # Labels
    ax.set_xlabel(r'Distance, $x/a$', fontsize=10)
    ax.set_ylabel(r'x-displacement [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/Mandel_disp_x_python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Show match to Mandel solution, displacement y
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    for i in range(7, len(Storage["TStorage"]), int(round(len(Storage["TStorage"]) / 7))):
        t = Storage["TStorage"][i] - Gen["t_applied"]
        Mandel = Mandel_comp(a, b, cf, E, F, k, Ks, muf, phi, t, vd)

        ax.plot(x_d[idx], np.cumsum(Storage["e_yy"][i, idx] * Gen["dx"]), 'k-', linewidth=1,
                label='Simulation' if i == 7 else None)

        ax.plot(Mandel["x"], Mandel["uy"], 'r--', linewidth=1, label='Analytical Soln.' if i == 7 else None)

    # Labels
    ax.set_xlabel(r'Distance, $x/a$', fontsize=10)
    ax.set_ylabel(r'y-displacement [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/Mandel_disp_y_python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()