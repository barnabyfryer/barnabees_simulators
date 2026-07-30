import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, BoundaryNorm

def plotting_sim(f, Gen, Plotting, Pos, Sig, Strain, e_vol):

    plt.close('all')

    # =============================================================================
    # Unpack
    # =============================================================================

    Sigx = Sig[:, 0].reshape(Gen["Ny"], Gen["Nx"])
    Sigy = Sig[:, 1].reshape(Gen["Ny"], Gen["Nx"])
    Sigxy = Sig[:, 2].reshape(Gen["Ny"], Gen["Nx"])

    e_xx = Strain[:, 0].reshape(Gen["Ny"], Gen["Nx"])
    e_yy = Strain[:, 1].reshape(Gen["Ny"], Gen["Nx"])
    e_xy = Strain[:, 2].reshape(Gen["Ny"], Gen["Nx"])

    x = Pos["x"]
    y = Pos["y"]

    x_new = Pos["x_new"]
    y_new = Pos["y_new"]

    du = Pos["du"]
    dv = Pos["dv"]

    nodesx = Gen["nodes_left"]
    nodesy = Gen["nodes_bottom"]

    nodes_fixed = x.reshape(Gen["Nn"],1)*0
    nodes_fixed[nodesx] += 1
    nodes_fixed[nodesy] += 2
    nodes_fixed = nodes_fixed.reshape(Gen["Ny"] + 1,Gen["Nx"] +  1)

    lwidth_1col = Plotting["lwidth_1col"]
    Position_1col_matrix = Plotting["Position_1col_matrix"]
    fsize_1col = Plotting["fsize_1col"]

    # =============================================================================
    # Cell centers
    # =============================================================================

    xc = np.zeros((Gen["Ny"], Gen["Nx"]))
    yc = np.zeros((Gen["Ny"], Gen["Nx"]))

    # =============================================================================
    # Forces
    # =============================================================================

    Fx = f[0::2].reshape(Gen["Ny"] + 1, Gen["Nx"] + 1)
    Fy = f[1::2].reshape(Gen["Ny"] + 1, Gen["Nx"] + 1)

    # =============================================================================
    # Cell centre coordinates
    # =============================================================================

    for j in range(Gen["Ny"]):
        for i in range(Gen["Nx"]):
            xc[j, i] = (
                               Pos["x_new"][j, i] +
                               Pos["x_new"][j, i + 1] +
                               Pos["x_new"][j + 1, i] +
                               Pos["x_new"][j + 1, i + 1]
                       ) / 4

            yc[j, i] = (
                               Pos["y_new"][j, i] +
                               Pos["y_new"][j, i + 1] +
                               Pos["y_new"][j + 1, i] +
                               Pos["y_new"][j + 1, i + 1]
                       ) / 4

    # =============================================================================
    # Plot deformed + undeformed grid
    # =============================================================================

    fig, ax = plt.subplots()

    # grid styling equivalent
    for j in range(Gen["Ny"] + 1):
        ax.plot(x[j, :], y[j, :], "b-")
        ax.plot(x[j, :], y[j, :], "b*")

        ax.plot(x_new[j, :], y_new[j, :], "r-")
        ax.plot(x_new[j, :], y_new[j, :], "r*")

    for k in range(Gen["Nx"] + 1):
        ax.plot(x[:, k], y[:, k], "b-")
        ax.plot(x_new[:, k], y_new[:, k], "r-")

    ax.set_xlabel("x-Location [m]")
    ax.set_ylabel("y-Location [m]")

    ax.set_aspect("equal")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.savefig('Plots/Grid.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot x displacements
    # =============================================================================

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    im2 = ax2.imshow(
        du,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax2.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax2.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax2.tick_params(direction='out', labelsize=fsize_1col)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    cbar2 = fig2.colorbar(im2)
    cbar2.set_label(r'x-displacement [m]', fontsize=fsize_1col)
    cbar2.ax.tick_params(labelsize=fsize_1col)

    fig2.patch.set_facecolor('white')

    fig2.savefig(
        'Plots/du_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot y displacements
    # =============================================================================

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)

    im3 = ax3.imshow(
        dv,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax3.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax3.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax3.tick_params(direction='out', labelsize=fsize_1col)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    cbar3 = fig3.colorbar(im3)
    cbar3.set_label(r'y-displacement [m]', fontsize=fsize_1col)
    cbar3.ax.tick_params(labelsize=fsize_1col)

    fig3.patch.set_facecolor('white')

    fig3.savefig(
        'Plots/dv_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot x forces
    # =============================================================================

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)

    im4 = ax4.imshow(
        Fx * Gen["Ny"] / Gen["Ly"],
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax4.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax4.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax4.tick_params(direction='out', labelsize=fsize_1col)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    cbar4 = fig4.colorbar(im4)
    cbar4.set_label(r'$S_{x}$ [Pa]', fontsize=fsize_1col)
    cbar4.ax.tick_params(labelsize=fsize_1col)

    fig4.patch.set_facecolor('white')

    fig4.savefig(
        'Plots/Fx_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot y forces
    # =============================================================================

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111)

    im5 = ax5.imshow(
        Fy * Gen["Nx"] / Gen["Lx"],
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax5.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax5.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax5.tick_params(direction='out', labelsize=fsize_1col)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)

    cbar5 = fig5.colorbar(im5)
    cbar5.set_label(r'$S_{y}$ [Pa]', fontsize=fsize_1col)
    cbar5.ax.tick_params(labelsize=fsize_1col)

    fig5.patch.set_facecolor('white')

    fig5.savefig(
        'Plots/Fy_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot x stresses
    # =============================================================================

    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111)

    im6 = ax6.imshow(
        Sigx,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax6.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax6.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax6.tick_params(direction='out', labelsize=fsize_1col)
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)

    cbar6 = fig6.colorbar(im6)
    cbar6.set_label(r'$\sigma_{xx}$ [Pa]', fontsize=fsize_1col)
    cbar6.ax.tick_params(labelsize=fsize_1col)

    fig6.patch.set_facecolor('white')

    fig6.savefig(
        'Plots/Sigxx_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot y stresses
    # =============================================================================

    fig7 = plt.figure()
    ax7 = fig7.add_subplot(111)

    im7 = ax7.imshow(
        Sigy,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax7.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax7.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax7.tick_params(direction='out', labelsize=fsize_1col)
    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)

    cbar7 = fig7.colorbar(im7)
    cbar7.set_label(r'$\sigma_{yy}$ [Pa]', fontsize=fsize_1col)
    cbar7.ax.tick_params(labelsize=fsize_1col)

    fig7.patch.set_facecolor('white')

    fig7.savefig(
        'Plots/Sigyy_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot xy stresses
    # =============================================================================

    fig8 = plt.figure()
    ax8 = fig8.add_subplot(111)

    im8 = ax8.imshow(
        Sigxy,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax8.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax8.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax8.tick_params(direction='out', labelsize=fsize_1col)
    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)

    cbar8 = fig8.colorbar(im8)
    cbar8.set_label(r'$\sigma_{xy}$ [Pa]', fontsize=fsize_1col)
    cbar8.ax.tick_params(labelsize=fsize_1col)

    fig8.patch.set_facecolor('white')

    fig8.savefig(
        'Plots/Sigxy_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot volumetric strain
    # =============================================================================

    fig9 = plt.figure()
    ax9 = fig9.add_subplot(111)

    im9 = ax9.imshow(
        e_vol,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax9.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax9.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax9.tick_params(direction='out', labelsize=fsize_1col)
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)

    cbar9 = fig9.colorbar(im9)
    cbar9.set_label(r'$\epsilon_{v}$', fontsize=fsize_1col)
    cbar9.ax.tick_params(labelsize=fsize_1col)

    fig9.patch.set_facecolor('white')

    fig9.savefig(
        'Plots/e_vol_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot fixed nodes
    # =============================================================================

    cmap = ListedColormap(["purple", "blue", "green", "yellow"])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, cmap.N)

    facevals = np.maximum.reduce([
        nodes_fixed[:-1, :-1],
        nodes_fixed[1:, :-1],
        nodes_fixed[:-1, 1:],
        nodes_fixed[1:, 1:]
    ])
    facecolors = cmap(norm(facevals))

    fig10 = plt.figure()
    ax10 = fig10.add_subplot(111, projection='3d')
    pcm = ax10.plot_surface(x, y, nodes_fixed, facecolors=facecolors, norm=norm, shade=False)

    ax10.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax10.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    fig10.patch.set_facecolor('white')

    # Colorbar
    cbar = fig10.colorbar(
        ScalarMappable(cmap=cmap, norm=norm),
        ax=ax9,
        ticks=[0, 1, 2, 3]
    )
    cbar.ax.set_yticklabels([
        "Free node",
        "x-fixed",
        "y-fixed",
        "x & y fixed"
    ])

    fig10.savefig('Plots/Fixed_nodes.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot x-direction strain
    # =============================================================================

    fig11 = plt.figure()
    ax11 = fig11.add_subplot(111)

    im11 = ax11.imshow(
        e_xx,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax11.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax11.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax11.tick_params(direction='out', labelsize=fsize_1col)
    ax11.spines['top'].set_visible(False)
    ax11.spines['right'].set_visible(False)

    cbar11 = fig11.colorbar(im11)
    cbar11.set_label(r'$\epsilon_{xx}$', fontsize=fsize_1col)
    cbar11.ax.tick_params(labelsize=fsize_1col)

    fig11.patch.set_facecolor('white')

    fig11.savefig(
        'Plots/e_xx_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot y-direction strain
    # =============================================================================

    fig12 = plt.figure()
    ax12 = fig12.add_subplot(111)

    im12 = ax12.imshow(
        e_yy,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax12.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax12.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax12.tick_params(direction='out', labelsize=fsize_1col)
    ax12.spines['top'].set_visible(False)
    ax12.spines['right'].set_visible(False)

    cbar12 = fig12.colorbar(im12)
    cbar12.set_label(r'$\epsilon_{yy}$', fontsize=fsize_1col)
    cbar12.ax.tick_params(labelsize=fsize_1col)

    fig12.patch.set_facecolor('white')

    fig12.savefig(
        'Plots/e_yy_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    # =============================================================================
    # Plot engineering shear strain
    # =============================================================================

    fig13 = plt.figure()
    ax13 = fig13.add_subplot(111)

    im13 = ax13.imshow(
        e_xy,
        origin='lower',
        extent=[xc.min(), xc.max(), yc.min(), yc.max()],
        aspect='auto',
        cmap='viridis'
    )

    ax13.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax13.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    ax13.tick_params(direction='out', labelsize=fsize_1col)
    ax13.spines['top'].set_visible(False)
    ax13.spines['right'].set_visible(False)

    cbar13 = fig13.colorbar(im13)
    cbar13.set_label(r'$\gamma_{xy}$', fontsize=fsize_1col)
    cbar13.ax.tick_params(labelsize=fsize_1col)

    fig13.patch.set_facecolor('white')

    fig13.savefig(
        'Plots/gamma_xy_python.pdf',
        bbox_inches='tight'
    )

    plt.show()

    return