import numpy as np

def input_data():

    # =============================================================================
    # General inputs
    # =============================================================================

    Gen = {}
    #Elements in x direction
    Gen["Nx"] = 21
    #Elements in y direction
    Gen["Ny"] = 21
    #Reservoir length [m]
    Gen["Lx"] = 30
    Gen["Ly"] = 10
    #Young's modulus [Pa]
    Gen["E"] = 1e9
    #Poisson's ratio
    Gen["v"] = 0.3

    # =============================================================================
    # Basic calculations
    # =============================================================================

    #Element edge lengths
    Gen["dx"] = Gen["Lx"] / Gen["Nx"]
    Gen["dy"] = Gen["Ly"] / Gen["Ny"]

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
    Pos["x"], Pos["y"] = np.meshgrid(x, y, indexing="xy")

    # =============================================================================
    # Fixed displacement nodes
    # =============================================================================

    #Fix x-direction
    #Fix left nodes
    Gen["nodesx"] = np.arange(0, Gen["Nn"], Gen["Nx"] + 1)
    #Fix right nodes
    # Gen["nodesx"] = np.arange(Gen["Nx"], Gen["Nn"], Gen["Nx"] + 1)

    #Fix y-direction
    #Fix bottom nodes
    Gen["nodesy"] = np.arange(0, Gen["Nx"] + 1, 1)
    #Fix top nodes
    Gen["nodesy"] = np.arange(Gen["Nn"] - Gen["Nx"] - 1, Gen["Nn"], 1)

    return Gen, Pos

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