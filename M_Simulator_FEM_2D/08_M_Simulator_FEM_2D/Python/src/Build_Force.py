import numpy as np

def Build_Force(Gen):


    # =============================================================================
    # Find where to apply force
    # =============================================================================

    left_nodes = np.arange(0, Gen["Nn"], Gen["Nx"] + 1)
    right_nodes = np.arange(Gen["Nx"], Gen["Nn"], Gen["Nx"] + 1)
    bottom_nodes = np.arange(0, Gen["Nx"] + 1)
    top_nodes = np.arange(Gen["Nn"] - Gen["Nx"] - 1, Gen["Nn"])

    # DOFs
    left_x = 2 * left_nodes
    right_x = 2 * right_nodes
    bottom_y = 2 * bottom_nodes + 1
    top_y = 2 * top_nodes + 1


    # =============================================================================
    # Make force vector
    # =============================================================================

    F = np.zeros(2*Gen["Nn"])

    #Forces to be applied (per unit of out of plane thickness) [N/m]
    fy = -10e7
    fx = -8.5e7

    #To which nodes to we apply it
    #On top nodes
    F[top_y[1:-1]] = fy
    #On right nodes
    F[right_x[1:-1]] = fx
    #On bottom nodes
    # F[bottom_y[1:-1]] = -fy
    #On left nodes
    # F[left_x[1:-1]] = -fx

    #Only apply half the force on the ends
    #Top left corner
    F[top_y[0]] = fy / 2
    #Top right corner
    F[top_y[-1]] = fy / 2

    #Right bottom corner
    F[right_x[0]] = fx/2
    #Right top corner
    F[right_x[-1]] = fx/2

    #Bottom left corner
    # F[bottom_y[0]] = -fy/2
    #Bottom right corner
    # F[bottom_y[-1]] = -fy/2

    #Left bottom corner
    # F[left_x[0]] = -fx/2
    #Left top corner
    # F[left_x[-1]] = -fx/2

    return F