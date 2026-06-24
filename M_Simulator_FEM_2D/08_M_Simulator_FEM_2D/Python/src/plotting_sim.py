import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, BoundaryNorm

def plotting_sim(f, Gen, Plotting, Pos, Sig):

    plt.close('all')

    # =============================================================================
    # Unpack
    # =============================================================================

    Sigx = Sig[:, 0].reshape(Gen["Ny"], Gen["Nx"])
    Sigy = Sig[:, 1].reshape(Gen["Ny"], Gen["Nx"])
    Sigxy = Sig[:, 2].reshape(Gen["Ny"], Gen["Nx"])

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
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.plot_surface(x, y, du, cmap='viridis')

    ax2.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax2.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax2.set_zlabel('x-displacement [m]', fontsize=fsize_1col)

    fig2.patch.set_facecolor('white')

    fig2.savefig('Plots/du.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot y displacements
    # =============================================================================

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.plot_surface(x, y, dv, cmap='viridis')

    ax3.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax3.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax3.set_zlabel('y-displacement [m]', fontsize=fsize_1col)

    fig3.patch.set_facecolor('white')

    fig3.savefig('Plots/dv.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot x forces
    # =============================================================================

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111, projection='3d')
    ax4.plot_surface(x, y, Fx, cmap='viridis')

    ax4.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax4.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax4.set_zlabel('x-direction force [N]', fontsize=fsize_1col)

    fig4.patch.set_facecolor('white')

    fig4.savefig('Plots/Fx.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot y forces
    # =============================================================================

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111, projection='3d')
    ax5.plot_surface(x, y, Fy, cmap='viridis')

    ax5.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax5.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax5.set_zlabel('y-direction force [N]', fontsize=fsize_1col)

    fig5.patch.set_facecolor('white')

    fig5.savefig('Plots/Fy.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot x stresses
    # =============================================================================

    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111, projection='3d')
    ax6.plot_surface(xc, yc, Sigx, cmap='viridis')

    ax6.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax6.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax6.set_zlabel(r'$\sigma_{xx}$ [Pa]', fontsize=fsize_1col)

    fig6.patch.set_facecolor('white')

    fig6.savefig('Plots/Sigxx.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot y stresses
    # =============================================================================

    fig7 = plt.figure()
    ax7 = fig7.add_subplot(111, projection='3d')
    ax7.plot_surface(xc, yc, Sigy, cmap='viridis')

    ax7.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax7.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax7.set_zlabel(r'$\sigma_{yy}$ [Pa]', fontsize=fsize_1col)

    fig7.patch.set_facecolor('white')

    fig7.savefig('Plots/Sigyy.jpg',
                dpi=300,
                bbox_inches='tight')

    plt.show()

    # =============================================================================
    # Plot xy stresses
    # =============================================================================

    fig8 = plt.figure()
    ax8 = fig8.add_subplot(111, projection='3d')
    ax8.plot_surface(xc, yc, Sigxy, cmap='viridis')

    ax8.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax8.set_ylabel('y-Location [m]', fontsize=fsize_1col)
    ax8.set_zlabel(r'$\sigma_{xy}$ [Pa]', fontsize=fsize_1col)

    fig8.patch.set_facecolor('white')

    fig8.savefig('Plots/Sigxy.jpg',
                 dpi=300,
                 bbox_inches='tight')

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

    fig9 = plt.figure()
    ax9 = fig9.add_subplot(111, projection='3d')
    pcm = ax9.plot_surface(x, y, nodes_fixed, facecolors=facecolors, norm=norm, shade=False)

    ax9.set_xlabel('x-Location [m]', fontsize=fsize_1col)
    ax9.set_ylabel('y-Location [m]', fontsize=fsize_1col)

    fig9.patch.set_facecolor('white')

    # Colorbar
    cbar = fig9.colorbar(
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

    fig9.savefig('Plots/Fixed_nodes.jpg',
                 dpi=300,
                 bbox_inches='tight')

    plt.show()

    return