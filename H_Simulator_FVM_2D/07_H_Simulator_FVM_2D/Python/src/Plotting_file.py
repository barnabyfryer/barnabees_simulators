import matplotlib.pyplot as plt
from scipy.special import erfc
import numpy as np
from src.perm import perm
from src.phiCalc import phiCalc
from scipy.special import exp1

def Plotting_file(Flow,Gen,Storage,Wells):

    # =============================================================================
    # Basic calculations
    # =============================================================================

    X = Storage["x"].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    Y = Storage["y"].reshape((Gen["Nx"], Gen["Ny"]), order='F').T

    P2D = Storage["P"][1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    flux_2D = Storage["flux"][1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    kx_2D = Storage["kx"][1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    ky_2D = Storage["ky"][1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    phi_2D = Storage["phi"][1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T


    # =============================================================================
    # Plotting pressure
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, P2D/1e6, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    #Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label('Pressure [MPa]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/Pp_2D_Python.jpg',
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
    pcm = ax.pcolormesh(X, Y, flux_2D, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label('Mass rate [kg/s]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # fig.savefig('../Verification/Flux_Python.jpg',
    #             dpi=300,
    #             bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting permeability
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, kx_2D, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label('x-Permeability [m$^2$]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # fig.savefig('../Verification/kx_Python.jpg',
    #             dpi=300,
    #             bbox_inches='tight')
    plt.show()


    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, ky_2D, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label('y-Permeability [m$^2$]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # fig.savefig('../Verification/ky_Python.jpg',
    #             dpi=300,
    #             bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting porosity
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, phi_2D, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label('Porosity [-]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # fig.savefig('../Verification/phi_Python.jpg',
    #             dpi=300,
    #             bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting permeability
    # =============================================================================

    dist = np.sqrt((Wells["xQ"][0] - Storage["x"])**2 + (Wells["yQ"][0] - Storage["y"])**2)

    dist_max = np.max(dist)
    dist_dummy = np.linspace(0, dist_max, num=1000)
    t_eval = Storage["TStorage"][-1]
    ct = Flow["cf"]
    alpha = Flow["kx0"][0]/(Flow["phi0"][0]*Flow["muf"]*ct)
    u = dist_dummy ** 2 / (4 * alpha * t_eval)
    P_Theis = Storage["P"][0,0] + Flow["muf"]*Wells["Q"][0]/(4*np.pi*Flow["kx0"][0]*Gen["Lz"]*Flow["Rho0"])*exp1(u)


    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    ax.plot(dist, Storage["P"][-1, :] / 1e6, 'k.', linewidth=1, label='Simulation')
    ax.plot(dist_dummy, P_Theis / 1e6, 'r--', linewidth=1, label='Analytical Soln.')
    # Labels
    ax.set_xlabel(r'Radial distance, $r$ [m]', fontsize=10)
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
    fig.savefig('../Verification/Pp_Theis_Python.jpg',
    dpi=300,
    bbox_inches='tight')
    plt.show()

    return