import matplotlib.pyplot as plt
from scipy.special import erfc
import numpy as np
from src.perm import perm
from src.phiCalc import phiCalc
from src.density import density

def Plotting_file(Flow,Gen,Storage,Wells):

    # =============================================================================
    # Validation
    # =============================================================================
    #Here we find an "analytical" solution by reaching steady state with one
    #fixed pressure well and one fixed injection mass rate. Then, knowing the
    #mass rate everywhere, we find the pressure at the next cell. This results
    #in an update to permeability, which in turn updates pressure again. We
    #keep iterating until convergence on that cell. Then we move to the next
    #cell.
    Q = - Wells["Q"][0]
    tol = 1e-6
    P = np.zeros(len(Storage["x"])+1) + Wells["P"][0]
    for i in range(len(Storage["x"])):
        err = 1.0

        # Initial guess
        P[i + 1] = P[i]

        while err > tol:
            P_old = P[i + 1]

            State = {}
            State["P"] = P_old

            # Update permeability
            kx,_ = perm(Flow, State)
            k = kx[0]

            # Darcy equation
            P[i+1] = (P[i] - Q * Flow["muf"] * Gen["dx"] / (Flow["Rho0"] * k) * (Gen["Ly"] * Gen["Lz"]))

            err = abs(P[i + 1] - P_old)




    # Find permeability
    State = {}
    State["P"] = Storage["P"][1, :]
    kx, dkxdP = perm(Flow, State)
    phi, _ = phiCalc(Flow, State)

    # =============================================================================
    # Plotting pressure
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], Storage["P"][1, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], P[:-1]/1e6, 'r--', linewidth=1, label='Analytical Soln.')
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
    # fig.savefig('../Verification/Pp_Python.jpg',
    # dpi=300,
    # bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting velocity
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], Storage["flux"][1, :], 'k-', linewidth=1, label='Simulation')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Mass rate [kg/s]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/Flux_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting permeability
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], kx, 'k-', linewidth=1, label='Simulation')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Permeability [m^2]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yscale('log')
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    plt.show()

    # =============================================================================
    # Plotting porosity
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], phi, 'k-', linewidth=1, label='Simulation')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Porosity [-]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/phi_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting mass
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    #Find analytical injected mass
    t = np.linspace(0, np.max(Storage["TStorage"]), 100)
    ax.plot(t, t*Wells["Q"][0], 'r--', linewidth=1, label='Analytical Soln.')
    #Find simulator total mass
    mass = np.zeros(len(Storage["TStorage"]))
    V = Gen["dx"] * Gen["Ly"] * Gen["Lz"]
    for i in range(len(Storage["TStorage"])):
        rho,_ = density(Flow, Storage["P"][i, :])
        State = {}
        State["P"]= Storage["P"][i, :]
        phi,_ = phiCalc(Flow, State)
        mass[i] = np.sum(rho * phi * V)
        # Plot
        ax.plot(Storage["TStorage"][i], mass[i]-mass[0], 'ko', linewidth=1, label='Simulation' if i == 0 else None)

    # Labels
    ax.set_xlabel(r'Time, $t$ [s]', fontsize=10)
    ax.set_ylabel(r'Mass change [kg]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lgd = ax.legend(fontsize=7)
    lgd.set_frame_on(False)
    fig.savefig('../Verification/mass_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    return