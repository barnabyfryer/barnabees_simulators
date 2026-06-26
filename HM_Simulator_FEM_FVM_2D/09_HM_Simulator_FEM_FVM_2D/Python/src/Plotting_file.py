import matplotlib.pyplot as plt

def Plotting_file(Gen,Pos,Storage):

    # =============================================================================
    # Basic calculations
    # =============================================================================

    X = Storage["x"].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    Y = Storage["y"].reshape((Gen["Nx"], Gen["Ny"]), order='F').T

    P2D = Storage["P"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    flux_2D = Storage["flux"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    kx_2D = Storage["kx"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    ky_2D = Storage["ky"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    phi_2D = Storage["phi"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    e_vol = Storage["e_vol"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T
    fx = Storage["fx"][-1, :].reshape((Gen["Nx"]+1, Gen["Ny"]+1), order='F').T
    fy = Storage["fy"][-1, :].reshape((Gen["Nx"]+1, Gen["Ny"]+1), order='F').T
    s_xx = Storage["s_xx"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T + P2D
    s_yy = Storage["s_yy"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T + P2D
    s_xy = Storage["s_xy"][-1, :].reshape((Gen["Nx"], Gen["Ny"]), order='F').T

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

    # =============================================================================
    # Plotting volumetric strain
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, e_vol, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'$\epsilon_{v}$ [-]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/e_vol_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting total stress change, sxx
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, s_xx/1e6, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'Total stress change, $\Delta S_{xx}$ [MPa]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/s_xx_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting total stress change, syy
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, s_yy / 1e6, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'Total stress change, $\Delta S_{yy}$ [MPa]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/s_yy_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting shear stress change, sxy
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(X, Y, s_xy / 1e6, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'Stress change, $\Delta S_{xy}$ [MPa]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/s_xy_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting force in x direction
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(Pos["x"], Pos["y"], fx, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'$f_{x}$ [N/m]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/fx_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    # =============================================================================
    # Plotting force in x direction
    # =============================================================================

    fig, ax = plt.subplots()
    # Figure size
    fig.set_size_inches(6, 4.5)
    # Plot
    pcm = ax.pcolormesh(Pos["x"], Pos["y"], fy, shading='nearest')
    # Labels
    ax.set_xlabel(r'Position, $x$ [m]', fontsize=10)
    ax.set_ylabel(r'Position, $y$ [m]', fontsize=10)
    # Font size
    ax.tick_params(labelsize=7)
    # Tick direction
    ax.tick_params(direction='out')
    # Colorbar
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(r'$f_{y}$ [N/m]')
    # Box off
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.savefig('../Verification/fy_Python.jpg',
                dpi=300,
                bbox_inches='tight')
    plt.show()

    return