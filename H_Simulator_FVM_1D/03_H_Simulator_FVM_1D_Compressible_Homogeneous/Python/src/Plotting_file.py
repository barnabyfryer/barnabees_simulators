import matplotlib.pyplot as plt

def Plotting_file(Gen,Storage):

    # =============================================================================
    # Plotting pressure
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(Gen["x"], Storage["P"][0, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["P"][1, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["P"][2, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["P"][3, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["P"][4, :] / 1e6, 'k-', linewidth=1, label='Simulation')
    # ax.plot(x, P_an/1e6, 'r--', linewidth=1, label='Analytical Soln.')
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
    ax.plot(Gen["x"], Storage["flux"][0, :], 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["flux"][1, :], 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["flux"][2, :], 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["flux"][3, :], 'k-', linewidth=1, label='Simulation')
    ax.plot(Gen["x"], Storage["flux"][4, :], 'k-', linewidth=1, label='Simulation')
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


    return