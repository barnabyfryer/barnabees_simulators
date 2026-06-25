import numpy as np

def build_force(Gen):


    # =============================================================================
    # Find where to apply force
    # =============================================================================

    left_nodes = Gen["nodes_left"]
    right_nodes = Gen["nodes_right"]
    bottom_nodes = Gen["nodes_bottom"]
    top_nodes = Gen["nodes_top"]

    # =============================================================================
    # Make force vector
    # =============================================================================

    F = np.zeros(2*Gen["Nn"])

    #Forces to be applied (per unit of out of plane thickness - note mesh dependent as written here) [N/m]
    fy = 0
    fx = 0

    #To which nodes to we apply it
    #On top nodes
    F[top_nodes[1:-1] * 2 + 1] = fy
    #On right nodes
    F[right_nodes[1:-1]*2] = fx
    #On bottom nodes
    # F[bottom_nodes[1:-1]*2 + 1] = -fy
    #On left nodes
    # F[left_nodes[1:-1]*2] = -fx

    #Only apply half the force on the ends
    #Top left corner
    F[2 * top_nodes[0] + 1] = fy / 2
    #Top right corner
    F[2*top_nodes[-1] + 1] = fy / 2

    #Right bottom corner
    F[2*right_nodes[0]] = fx / 2
    #Right top corner
    F[2*right_nodes[-1]] = fx / 2

    #Bottom left corner
    # F[2*bottom_nodes[0] + 1] = -fy/2
    #Bottom right corner
    # F[2*bottom_nodes[-1] + 1] = -fy/2

    #Left bottom corner
    # F[2*left_nodes[0]] = -fx / 2
    #Left top corner
    # F[2*left_nodes[-1]] = -fx / 2

    return F