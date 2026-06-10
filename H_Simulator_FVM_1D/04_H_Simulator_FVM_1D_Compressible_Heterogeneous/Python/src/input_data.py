import numpy as np

def input_data():

    # ------------------------------------------------------------------
    # General Parameters
    # ------------------------------------------------------------------

    Gen = {}

    Gen["tf"] = 10000.0  # Final time [sec]
    Gen["tstep"] = 100.0  # Time step [sec]
    Gen["tol"] = 1e-6  # Tolerance [-]

    Gen["Nx"] = 5000  # Number of cells [-]

    Gen["Lx"] = 10.0  # Reservoir length x [m]
    Gen["Ly"] = 1.0  # Reservoir length y [m]
    Gen["Lz"] = 1.0  # Reservoir height [m]

    Gen["PL"] = 10e6  # Left pressure [Pa]
    Gen["PR"] = 1e6  # Right pressure [Pa]

    # ------------------------------------------------------------------
    # Basic Calculations
    # ------------------------------------------------------------------

    Gen["dx"] = Gen["Lx"] / Gen["Nx"]

    # Define grid position
    Gen["x"] = np.linspace(Gen["dx"] / 2, Gen["Lx"] - Gen["dx"] / 2, Gen["Nx"])
    # Define edge positions
    Gen["x_edge"] = np.linspace(0, Gen["Lx"], Gen["Nx"] + 1)

    Storage = {}

    Storage["x"] = np.linspace(
        Gen["dx"] / 2,
        Gen["Lx"] - Gen["dx"] / 2,
        Gen["Nx"]
    )

    # ------------------------------------------------------------------
    # Flow Model
    # ------------------------------------------------------------------

    Flow = {}

    #Permeability on left side
    k_left = 1e-12
    #Permeability on right side
    k_right = 1e-13
    #Permeability band size on left side
    L_k = 2

    Flow["kx"] = np.zeros(Gen["Nx"]) + k_right  # Permeability [m²]
    Flow["kx"][Storage["x"] < L_k ] = k_left

    #Porosity on the left side
    phi_left = 0.3
    #Porosity on the right side
    phi_right = 0.2
    #Porosity band size on left side
    L_phi = 2

    Flow["phi"] = np.zeros(Gen["Nx"]) + phi_right  # Porosity [-]
    Flow["phi"][Storage["x"] < L_phi] = phi_left

    Flow["cf"] = 1e-8  # Fluid compressibility [1/Pa]
    Flow["muf"] = 0.1  # Fluid viscosity [Pa.s]

    Flow["Rho0"] = 1000.0  # Reference density [kg/m³]
    Flow["RhoP"] = 1e5  # Reference pressure [Pa]

    # ------------------------------------------------------------------
    # Initialize State
    # ------------------------------------------------------------------

    State = {}

    State["t"] = 0.0

    State["P"] = np.full(
        Gen["Nx"],
        1e5
    )

    State["step"] = 0

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
    Storage["P"] = np.zeros(
        (TStore+1, Gen["Nx"])
    )

    Storage["P"][0, :] = State["P"]

    #Store the initial flux (assumed zero) and prepare storage space
    Storage["flux"] = np.zeros(
        (TStore+1, Gen["Nx"])
    )

    return Flow, Gen, State, Storage