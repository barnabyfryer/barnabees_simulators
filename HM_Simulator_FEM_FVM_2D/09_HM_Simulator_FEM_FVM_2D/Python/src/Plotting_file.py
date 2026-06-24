import matplotlib.pyplot as plt
from scipy.special import erfc
import numpy as np
from src.src_flow.perm import perm
from src.src_flow.phiCalc import phiCalc

def Plotting_file(Flow,Gen,Storage):

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
    cbar.set_label('Mass flux [kg/s]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
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
    fig.savefig('../Verification/Pp_Python.jpg',
                dpi=300,
                bbox_inches='tight')
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
    fig.savefig('../Verification/Pp_Python.jpg',
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
    fig.savefig('../Verification/Pp_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    return