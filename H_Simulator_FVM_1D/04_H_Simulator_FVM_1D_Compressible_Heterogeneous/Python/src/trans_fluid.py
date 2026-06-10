import numpy as np
from scipy.sparse import diags

def trans_fluid(A,dRhodP,Flow,Gen,Rho,Trans):

    #Returns the transmissibilities between the two cells
    #The sum of all transmissibilities is also located in the main diagonal
    #basically this is so you can multiply FTrans*P (where P is a vector of
    #pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc
    #It also returns the derivative of each cell's transmissibility wrt to the
    #pressure of every cell, multiplied by the difference in pressure of those
    #two cells. This is effectively the effect compressibility has.


    #Basic grid properties
    Nx = Gen["Nx"]
    Ny = 1
    N = Nx*Ny

    #Predefine memory
    Tx = np.zeros((Nx+1,Ny))

    # =============================================================================
    # Transmissibility
    # =============================================================================

    #Apply upwind operator to find fluid mobility
    #This takes the upwind fluid mobility for the right interface of cell i
    Mupx = (A["x"] @ (Rho / Flow["muf"]))
    #Reshape these
    Mupx = np.asarray(Mupx).reshape(Nx, Ny)

    #Update transmissibility to include fluid mobility
    #Applies the upwind fluid mobility to the right interface
    Tx[1:Nx,0] = Trans["x"][1:-1] * Mupx[0:Nx-1,0]

    #Construct matrix
    #Grabs left interface trans (inc mobility) of a cell
    x1 = Tx[0:Nx].reshape(N,1)
    #Right interface trans of a cell
    x2 = Tx[1:Nx+1].reshape(N,1)

    #The values are the transmissibilities between the two cells
    FTrans = diags(
        [-x2.flatten(), (x1 + x2).flatten(), -x1[1:].flatten()],
        offsets=[-1, 0, 1],
        shape=(N, N)
    ).tocsr()

    #The sum of all transmissibilities is also located in the main diagonal basically this is so you can multiply
    # FTrans*P (where P is a vector of pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc

    # =============================================================================
    # dTransdP
    # =============================================================================

    #Fluid dependence on pressure at right x-interfaces
    dMupx = (A["x"] @ (dRhodP / Flow["muf"]))
    dMupx = np.asarray(dMupx).reshape(Nx, Ny)

    #Multiply the fluid derivative by the appropriate rock transmissibility
    dTxdp = np.zeros((Nx + 1, Ny))
    dTxdp[1:Nx, 0] = Trans["x"][1:-1] * dMupx[:-1, 0]

    #Finds the diagonals that will go into the dTransmissibility matrix
    L = dTxdp[:Nx, 0].reshape(N, 1)
    R = dTxdp[1:Nx + 1, 0].reshape(N, 1)

    #Create a matrix of the dependencies in the correct row
    dFTrans = diags(
        [-R.flatten(), (L + R).flatten(), -L[1:].flatten()],
        offsets=[-1, 0, 1],
        shape=(N, N)
    ).tocsr()

    return FTrans, dFTrans