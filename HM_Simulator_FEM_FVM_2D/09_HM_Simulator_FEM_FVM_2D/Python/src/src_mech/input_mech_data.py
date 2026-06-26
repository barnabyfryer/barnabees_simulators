import numpy as np
from src.src_mech.M_Simulator_FEM_2D import M_Simulator_FEM_2D
from src.src_flow.perm import perm
from src.src_flow.phiCalc import phiCalc

def input_mech_data(Flow, Gen, State, Storage):

    # =============================================================================
    # General inputs
    # =============================================================================

    #Young's modulus [Pa]
    Gen["E"] = 50e7
    #Poisson's ratio
    Gen["nu"] = 0.3
    #Choose plane strain ("plane") or plane stress ("stress")
    Gen["plane"] = "strain"
    #Biot coefficient
    Gen["biot"] = 1

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
    Pos["x"], Pos["y"] = np.meshgrid(x, y, indexing="xy")

    # =============================================================================
    # Fixed displacement nodes
    # =============================================================================

    #Get nodes for fixing later
    #Left nodes
    Gen["nodes_left"] = np.arange(0, Gen["Nn"], Gen["Nx"] + 1)
    #Right nodes
    Gen["nodes_right"] = np.arange(Gen["Nx"], Gen["Nn"], Gen["Nx"] + 1)
    #Top nodes
    Gen["nodes_top"] = np.arange(Gen["Ny"]*(Gen["Nx"]+1), Gen["Nn"])
    #Bottom nodes
    Gen["nodes_bottom"] = np.arange(0, Gen["Nx"] + 1, 1)

    # ------------------------------------------------------------------
    # Initialize phi and e_vol
    # ------------------------------------------------------------------

    #Initialize e_vol
    State = M_Simulator_FEM_2D(Gen, Pos, State)
    #Initialize permeability
    State["kx"], State["ky"], _, _, _, _ = perm(Flow, State)
    #Initialize porosity
    State["phi"], _ = phiCalc(Flow, State)

    # ------------------------------------------------------------------
    # Storage Matrices
    # ------------------------------------------------------------------

    TStore = min(20, int(np.floor(Gen["tf"] / Gen["tstep"])))

    Storage["TStorage"] = np.arange(
        0,
        Gen["tf"] + Gen["tf"] / TStore,
        Gen["tf"] / TStore
    )

    Storage["TStorage"] = (
            np.floor(Storage["TStorage"] / Gen["tstep"])
            * Gen["tstep"]
    )

    # Store the initial pressure and prepare storage space
    Storage["P"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["P"][0, :] = State["P"]

    # Store the initial porosity
    Storage["phi"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["phi"][0, :] = State["phi"]

    # Store the initial permeability
    Storage["kx"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["kx"][0, :] = State["kx"]
    Storage["ky"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["ky"][0, :] = State["ky"]

    # Store the initial flux (assumed zero) and prepare storage space
    Storage["flux"] = np.zeros((TStore + 1, Gen["Ne"]))

    # Store the initial volumetric strain (assumed zero) and prepare storage space
    Storage["e_vol"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["e_vol"][0, :] = State["e_vol"]

    #Initialize errors
    Storage["errP"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["erre"] = np.zeros((TStore + 1, Gen["Ne"]))

    # Initialize forces
    Storage["fx"] = np.zeros((TStore + 1, Gen["Nn"]))
    Storage["fy"] = np.zeros((TStore + 1, Gen["Nn"]))

    # Initialize effective stresses
    Storage["s_xx"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["s_yy"] = np.zeros((TStore + 1, Gen["Ne"]))
    Storage["s_xy"] = np.zeros((TStore + 1, Gen["Ne"]))

    return Gen, Pos, State, Storage

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