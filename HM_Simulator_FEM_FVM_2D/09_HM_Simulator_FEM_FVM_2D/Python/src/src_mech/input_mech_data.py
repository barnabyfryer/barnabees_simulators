import numpy as np

def input_mech_data(Gen, State):

    # =============================================================================
    # General inputs
    # =============================================================================

    #Young's modulus [Pa]
    Gen["E"] = 1e9
    #Poisson's ratio
    Gen["nu"] = 0.3
    #Choose plane strain ("plane") or plane stress ("stress")
    Gen["plane"] = "strain"

    # =============================================================================
    # Basic calculations
    # =============================================================================

    #Number of elements
    Gen["Ne"] = Gen["Nx"] * Gen["Ny"]
    #Number of nodes
    Gen["Nn"] = (Gen["Nx"]+1) * (Gen["Ny"]+1)
    #Matrix with the nodes making up each element
    Gen["Ref"] = build_connectivity(Gen["Nx"], Gen["Ny"])

    #Original node locations
    x = np.linspace(0, Gen["Lx"], Gen["Nx"] + 1)
    y = np.linspace(0, Gen["Ly"], Gen["Ny"] + 1)

    Pos = {}
    Pos["x"], Pos["y"] = np.meshgrid(x, y, indexing="ij")

    # =============================================================================
    # Fixed displacement nodes
    # =============================================================================

    #Get nodes for fixing later
    #Left nodes
    Gen["nodes_left"] = np.arange(0, Gen["Nn"], Gen["Nx"] + 1)
    #Right nodes
    Gen["nodes_right"] = np.arange(Gen["Nx"], Gen["Nn"], Gen["Nx"] + 1)
    #Top nodes
    Gen["nodes_top"] = np.arange(Gen["Ny"]*(Gen["Nx"]+1), (Gen["Ny"]+1)*(Gen["Nx"]+1))
    #Bottom nodes
    Gen["nodes_bottom"] = np.arange(0, Gen["Nx"] + 1, 1)

    # =============================================================================
    # Plotting parameters
    # =============================================================================

    Plotting = {}
    Plotting["lwidth_1col"] = 0.75
    Plotting["Position_1col_matrix"] = [2.2, 1.8, 6, 4.5]
    Plotting["fsize_1col"] = 7

    return Gen, Plotting, Pos, State

# =============================================================================
# Connectivity function
# =============================================================================

def build_connectivity(Nx, Ny):
    elements = []

    for j in range(Ny):
        for i in range(Nx):
            n1 = j * (Nx + 1) + i  # bottom-left
            n2 = n1 + 1  # bottom-right
            n3 = n1 + (Nx + 1)  # top-left
            n4 = n3 + 1  # top-right

            elements.append([n1, n2, n3, n4])

    return np.array(elements)