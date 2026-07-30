import matplotlib.pyplot as plt
from scipy.special import erfc
import numpy as np

def Plotting_file(Flow,Gen,Storage, Wells):

    # =============================================================================
    # Validation
    # =============================================================================
    #Steady state analytcal solution for different fixed pressure wells
    # idx2 = Gen["x"] > Wells["xP"][1]
    # P_an = Wells["P"][0] - (Wells["P"][0] - Wells["P"][1]) * Storage["x"]/Wells["xP"][1]
    # P_an[idx2] = Wells["P"][1] - (Wells["P"][1] - Wells["P"][2]) * (Storage["x"][idx2] - Wells["xP"][1])/(Wells["xP"][2] - Wells["xP"][1])

    #Steady state for the constant flow rate example given in verification file
    P_an = 100e6 - 10e6*Gen["x"]


    #Error calculation
    #err = np.abs((Storage["P"][1, :] - P_an) / (Gen["PL"] - Gen["PR"]))
    #Ep = np.max(err)
    #print(f'Maximum relative pressure error     = {Ep:.6e}')

    # =============================================================================
    # Plotting pressure
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], Storage["P"][1, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], P_an/1e6, 'r--', linewidth=1, label='Analytical Soln.')
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
    fig.savefig('../Verification/Pp_Python_q.jpg',
    dpi=300,
    bbox_inches='tight')
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
    fig.savefig('../Verification/Flux_Python_q.jpg',
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
    ax.plot(Gen["x"], Flow["kx"], 'k-', linewidth=1, label='Simulation')
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
    ax.plot(Gen["x"], Flow["phi"], 'k-', linewidth=1, label='Simulation')
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
    plt.show()

    return