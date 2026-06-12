import matplotlib.pyplot as plt
from scipy.special import erfc
import numpy as np

def Plotting_file(Flow,Gen,Storage):

    # =============================================================================
    # Validation
    # =============================================================================
    #Semi infinite solution, valid at early time
    alpha = Flow["kx"] / (Flow["phi"] * Flow["muf"] * Flow["cf"])
    P_an = 1e5 + (Gen["PL"] - 1e5) * erfc(Storage["x"] / (2 * np.sqrt(alpha * Storage["TStorage"][1])))
    #Error calculation
    err = np.abs((Storage["P"][1, :] - P_an) / (Gen["PL"] - Gen["PR"]))
    Ep = np.max(err)
    print(f'Maximum relative pressure error     = {Ep:.6e}')

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
    fig.savefig('../Verification/Pp_Python.jpg',
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
    ax.set_ylabel(r'Darcy Flux [m/s]', fontsize=10)
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
    fig.savefig('../Verification/Flux_Python.jpg',
                dpi=300,
                bbox_inches='tight')
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
    fig.savefig('../Verification/Flux_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    return