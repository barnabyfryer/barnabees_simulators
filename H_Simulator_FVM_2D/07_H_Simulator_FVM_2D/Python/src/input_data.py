import numpy as np

from src.phiCalc import phiCalc
from src.perm import perm

def input_data():

    # ------------------------------------------------------------------
    # General Parameters
    # ------------------------------------------------------------------

    Gen = {}

    Gen["tf"] = 1000.0  # Final time [sec]
    Gen["tstep"] = 10.0  # Time step [sec]
    Gen["tol"] = 1e-5  # Tolerance [-]

    Gen["Nx"] = 31  # Number of cells in x direction [-]
    Gen["Ny"] = 31  # Number of cells in y direction [-]

    Gen["Lx"] = 10.0  # Reservoir length x [m]
    Gen["Ly"] = 10.0  # Reservoir length y [m]
    Gen["Lz"] = 1.0  # Reservoir height [m]

    # ------------------------------------------------------------------
    # Basic Calculations
    # ------------------------------------------------------------------

    Gen["dx"] = Gen["Lx"] / Gen["Nx"]
    Gen["dy"] = Gen["Ly"] / Gen["Ny"]
    Gen["dz"] = Gen["Lz"]

    # Define grid position, cell centers
    xc = np.linspace(Gen["dx"] / 2, Gen["Lx"] - Gen["dx"] / 2, Gen["Nx"])
    yc = np.linspace(Gen["dy"] / 2, Gen["Ly"] - Gen["dy"] / 2, Gen["Ny"])
    # Mesh of cell centers
    xc_grid, yc_grid = np.meshgrid(xc, yc)

    Storage = {}

    Storage["x"] = xc_grid.T.reshape(-1)
    Storage["y"] = yc_grid.T.reshape(-1)

    # ------------------------------------------------------------------
    # Flow Model
    # ------------------------------------------------------------------

    Flow = {}

    # Permeability in x direction [m^2]
    #Permeability on left side
    kx_left = 1e-12
    #Permeability on right side
    kx_right = 1e-13
    #Permeability band size on left side
    Lx_k = 0
    #Reference permeability
    Flow["kx0"] = np.zeros(Gen["Nx"]*Gen["Ny"]) + kx_right  # Permeability [m²]
    Flow["kx0"][Storage["x"] < Lx_k ] = kx_left

    #Permeability in y direction [m^2]
    # Permeability on top side
    ky_top = 1e-12
    # Permeability on bot side
    ky_bot = 1e-13
    # Permeability band size on left side
    Ly_k = 0
    # Reference permeability
    Flow["ky0"] = np.zeros(Gen["Nx"]*Gen["Ny"]) + ky_bot  # Permeability [m²]
    Flow["ky0"][Storage["y"] < Ly_k] = ky_top

    #"Compressibility" of permeability
    Flow["ck"] = 1e-8
    #Reference pressure [Pa]
    Flow["kP0"] = 1e5

    #Porosity on the left side
    phi_left = 0.3
    #Porosity on the right side
    phi_right = 0.2
    #Porosity band size on left side
    L_phi = 0
    #Reference porosity
    Flow["phi0"] = np.zeros(Gen["Nx"]*Gen["Ny"]) + phi_right  # Porosity [-]
    Flow["phi0"][Storage["x"] < L_phi] = phi_left
    #"Compressibility" of permeability
    Flow["cphi"] = 1e-9
    #Reference pressure [Pa]
    Flow["phiP0"] = 1e5

    Flow["cf"] = 1e-8  # Fluid compressibility [1/Pa]
    Flow["muf"] = 0.1  # Fluid viscosity [Pa.s]

    Flow["Rho0"] = 1000.0  # Reference density [kg/m³]
    Flow["RhoP"] = 1e5  # Reference pressure [Pa]

    # ------------------------------------------------------------------
    # Wells
    # ------------------------------------------------------------------

    #Constant pressure wells
    Wells = {}
    Wells["P"] = [5e7, 1e7] #Row vector of well pressurese [Pa]
    Wells["WI"] = [0, 0]    #Well indices [m]
    Wells["xP"] = [10, 0]   #Well locations in x [m]
    Wells["yP"] = [10, 0]  #Well locations in y [m]

    #Find the cells for these wells
    Wells["Loc_P"] = np.zeros(np.size(Wells["P"]), dtype=int)
    for i in range(len(Wells["xP"])):
        Wells["Loc_P"][i] = np.argmin((Storage["x"] - Wells["xP"][i])**2 + (Storage["y"] - Wells["yP"][i])**2)

    # Constant rate wells
    Wells["Q"] = [1]  # Row vector of well rates [kg/sec]
    Wells["xQ"] = [5]  # Well locations in x [m]
    Wells["yQ"] = [5]  # Well locations in y [m]

    # Find the cells for these wells
    Wells["Loc_Q"] = np.zeros(np.size(Wells["Q"]), dtype=int)
    for i in range(len(Wells["xQ"])):
        Wells["Loc_Q"][i] = np.argmin((Storage["x"] - Wells["xQ"][i])**2 + (Storage["y"] - Wells["yQ"][i])**2)


    # ------------------------------------------------------------------
    # Initialize State
    # ------------------------------------------------------------------

    State = {}

    State["t"] = 0.0

    State["P"] = np.full(Gen["Nx"]*Gen["Ny"],1e5)

    State["step"] = 0

    State["kx"],State["ky"],_,_,_,_ = perm(Flow,State)

    State["phi"],_ = phiCalc(Flow,State)

    # ------------------------------------------------------------------
    # Storage Matrices
    # ------------------------------------------------------------------

    TStore = 1

    Storage["TStorage"] = np.arange(
        0,
        Gen["tf"] + Gen["tf"] / TStore,
        Gen["tf"] / TStore
    )

    Storage["TStorage"] = (
            np.floor(Storage["TStorage"] / Gen["tstep"])
            * Gen["tstep"]
    )


    #Store the initial pressure and prepare storage space
    Storage["P"] = np.zeros((TStore+1, Gen["Nx"]*Gen["Ny"]))
    Storage["P"][0, :] = State["P"]

    # Store the initial porosity
    Storage["phi"] = np.zeros((TStore + 1, Gen["Nx"]*Gen["Ny"]))
    Storage["phi"][0, :] = State["phi"]

    # Store the initial permeability
    Storage["kx"] = np.zeros((TStore + 1, Gen["Nx"]*Gen["Ny"]))
    Storage["kx"][0, :] = State["kx"]
    Storage["ky"] = np.zeros((TStore + 1, Gen["Nx"] * Gen["Ny"]))
    Storage["ky"][0, :] = State["ky"]

    #Store the initial flux (assumed zero) and prepare storage space
    Storage["flux"] = np.zeros((TStore+1, Gen["Nx"]*Gen["Ny"]))

    return Flow, Gen, State, Storage, Wells