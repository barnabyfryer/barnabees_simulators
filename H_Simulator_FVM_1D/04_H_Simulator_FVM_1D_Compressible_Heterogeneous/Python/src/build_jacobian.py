from src.density import density
from src.upwind import upwind
from src.trans_fluid import trans_fluid
from src.add_bc import add_bc
from scipy.sparse import diags

def build_jacobian(Flow,Gen,State,trans):


    #Basic grid parameters
    # Length of cell in x-direction [m]
    dx = Gen["dx"]
    # Length of cell in y-direction [m]
    dy = Gen["Ly"] / 1
    # Length of cell in z-direction [m]
    dz = Gen["Lz"] / 1
    # Cell volume
    V = dx * dy * dz
    #Number of cells
    N = Gen["Nx"] * 1

    #Find density and derivative
    Rho, dRhodP = density(Flow,State["P"])
    #Find upwinding
    A = upwind(Gen,State,trans)

    #Returns the transmissibilities between the two cells
    FTrans, dFTrans = trans_fluid(A,Flow,Gen,State)

    #The sum of all transmissibilities is also located in the main diagonal basically this is so you can multiply
    #FTrans*P (where P is a vector of pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc It also
    #returns the derivative of each cell's transmissibility wrt to the pressure of every cell, multiplied by the
    #difference in pressure of those two cells. This is effectively the effect compressibility has

    # =============================================================================
    # Derivative of accumulation terms
    # =============================================================================

    Acc = dRhodP * Flow["phi"] * V / Gen["tstep"]

    # =============================================================================
    # Derivative of convection terms
    # =============================================================================

    #Due to change in pressure
    Conv = FTrans
    #Due to compressibility
    Comp = dFTrans

    # =============================================================================
    # Add boundary conditions
    # =============================================================================

    _, dBC = add_bc(dRhodP, Gen, Rho, State, trans)

    # =============================================================================
    # Combine terms
    # =============================================================================

    #Diagonal terms
    DiagVecs = Acc + dBC  # (N,)
    Diag = diags(DiagVecs, offsets=0, shape=(N, N), format="csr")

    # --- Combine full Jacobian ---
    jac = Diag + Conv + Comp


    return jac