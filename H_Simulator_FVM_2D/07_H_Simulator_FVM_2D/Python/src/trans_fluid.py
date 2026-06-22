import numpy as np
from src.density import density
from src.trans_rock import trans_rock
from scipy.sparse import diags

def trans_fluid(A,Flow,Gen,State):

    #Returns the transmissibilities between the two cells
    #The sum of all transmissibilities is also located in the main diagonal
    #basically this is so you can multiply FTrans*P (where P is a vector of
    #pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc
    #It also returns the derivative of each cell's transmissibility wrt to the
    #pressure of every cell, multiplied by the difference in pressure of those
    #two cells. This is effectively the effect compressibility has.


    #Basic grid properties
    Nx = Gen["Nx"]
    Ny = Gen["Ny"]
    N = Nx*Ny

    #Grab pressure
    P = np.reshape(State["P"],(Nx,Ny))
    #Get density
    Rho,dRhodP = density(Flow,State["P"])
    #Get rock transmissibility
    Trans, dTransdP_LB, dTransdP_RT = trans_rock(Flow,Gen,State)

    #Predefine memory
    Tx = np.zeros((Nx+1,Ny))
    Ty = np.zeros((Nx, Ny+1))

    # =============================================================================
    # Transmissibility
    # =============================================================================

    #Apply upwind operator to find fluid mobility
    #This takes the upwind fluid mobility for the right interface of cell i
    Mupx = (A["x"] @ (Rho / Flow["muf"]))
    Mupy = (A["y"] @ (Rho / Flow["muf"]))
    #Reshape these
    Mupx = np.asarray(Mupx).reshape(Nx, Ny, order='F')
    Mupy = np.asarray(Mupy).reshape(Nx, Ny, order='F')

    #Update transmissibility to include fluid mobility
    #Applies the upwind fluid mobility to the right interface
    Tx[1:Nx,:] = Trans["x"][1:-1,:] * Mupx[0:Nx-1,:]
    #Do same thing for top interface
    Ty[:, 1:Ny] = Trans["y"][:,1:-1] * Mupy[:,0:Ny - 1]

    #Construct matrix
    #Grabs left interface trans (inc mobility) of a cell
    x1 = Tx[0:Nx,:].reshape(N,1, order='F')
    #Right interface trans of a cell
    x2 = Tx[1:Nx+1,:].reshape(N,1, order='F')
    #Bottom
    y1 = Ty[:,0:Ny].reshape(N,1, order='F')
    #Top
    y2 = Ty[:,1:Ny+1].reshape(N,1, order='F')


    #The values are the transmissibilities between the two cells
    FTrans = diags(
        [-y2.flatten(),-x2.flatten(), (y2 + x1 + x2 + y1).flatten(), -x1[1:].flatten(),-y1[Nx:].flatten()],
        offsets=[-Nx,-1, 0, 1,Nx],
        shape=(N, N)
    ).tocsr()

    #The sum of all transmissibilities is also located in the main diagonal basically this is so you can multiply
    # FTrans*P (where P is a vector of pressure) and get the equations T_R * (P-P_R) + T_L * (P-P_L) etc

    # =============================================================================
    # dTransdP, part 1 (density dependence on pressure)
    # =============================================================================

    #Fluid dependence on pressure at right x-interfaces
    dMupx = (A["x"] @ (dRhodP / Flow["muf"]))
    dMupx = np.asarray(dMupx).reshape(Nx, Ny, order='F')

    dMupy = (A["y"] @ (dRhodP / Flow["muf"]))
    dMupy = np.asarray(dMupy).reshape(Nx, Ny, order='F')

    #Multiply the fluid derivative by the appropriate rock transmissibility
    dTxdp_L = np.zeros((Nx + 1, Ny))
    dTxdp_R = np.zeros((Nx + 1, Ny))
    dTxdp_L[1:Nx, :] = Trans["x"][1:-1,:] * dMupx[1:, :] * (P[1:,:] - P[:-1,:])
    dTxdp_R[1:Nx, :] = Trans["x"][1:-1,:] * dMupx[:-1, :] * (P[1:,:] - P[:-1,:])
    #Same for top and bottom
    dTydp_B = np.zeros((Nx, Ny+1))
    dTydp_T = np.zeros((Nx, Ny+1))
    dTydp_B[:,1:Ny] = Trans["y"][:,1:-1] * dMupy[:,1:] * (P[:,1:] - P[:,:-1])
    dTydp_T[:,1:Ny] = Trans["y"][:,1:-1] * dMupy[:,:-1] * (P[:,1:] - P[:,:-1])

    #Finds the diagonals that will go into the dTransmissibility matrix
    L = dTxdp_L[:Nx, :].reshape(N, 1, order='F')
    R = -dTxdp_R[1:Nx + 1, :].reshape(N, 1, order='F')
    B = dTydp_B[:,:Ny].reshape(N, 1, order='F')
    T = -dTydp_T[:,1:Ny + 1].reshape(N, 1, order='F')

    #Create a matrix of the dependencies in the correct row
    dFTrans_1 = diags(
        [-T.flatten(),-R.flatten(), (T+L + R+B).flatten(), -L[1:].flatten(),-B[Nx:].flatten()],
        offsets=[-Nx,-1, 0, 1,Nx],
        shape=(N, N)
    ).tocsr()

    # print("dFTrans=", dFTrans_1.toarray())

    # =============================================================================
    # dTransdP, part 2 (permeability dependence on pressure)
    # =============================================================================

    # Fluid dependence on pressure at right x-interfaces
    Mupx = (A["x"] @ (Rho / Flow["muf"]))
    Mupx = np.asarray(Mupx).reshape(Nx, Ny, order='F')

    Mupy = (A["y"] @ (Rho / Flow["muf"]))
    Mupy = np.asarray(Mupy).reshape(Nx, Ny, order='F')

    # Multiply the fluid derivative by the appropriate rock transmissibility
    dTxdp_L = np.zeros((Nx + 1, Ny))
    dTxdp_R = np.zeros((Nx + 1, Ny))
    dTxdp_L[1:Nx, :] = dTransdP_LB["x"][1:-1,:] * Mupx[1:, :] * (P[1:, :] - P[:-1, :])
    dTxdp_R[1:Nx, :] = dTransdP_RT["x"][1:-1,:] * Mupx[:-1, :] * (P[1:, :] - P[:-1, :])

    dTydp_B = np.zeros((Nx, Ny+1))
    dTydp_T = np.zeros((Nx, Ny+1))
    dTydp_B[:,1:Ny] = dTransdP_LB["y"][:,1:-1] * Mupy[:,1:] * (P[:,1:] - P[:,:-1])
    dTydp_T[:,1:Ny] = dTransdP_RT["y"][:,1:-1] * Mupy[:,:-1] * (P[:,1:] - P[:,:-1])

    # Finds the diagonals that will go into the dTransmissibility matrix
    L = dTxdp_L[:Nx, :].reshape(N, 1, order='F')
    R = -dTxdp_R[1:Nx + 1, :].reshape(N, 1, order='F')
    B = dTydp_B[:,:Ny].reshape(N, 1, order='F')
    T = -dTydp_T[:,1:Ny + 1].reshape(N, 1, order='F')

    # Create a matrix of the dependencies in the correct row
    dFTrans = diags(
        [-T.flatten(),-R.flatten(), (T+L + R+B).flatten(), -L[1:].flatten(), -B[Nx:].flatten()],
        offsets=[-Nx,-1, 0, 1,Nx],
        shape=(N, N)
    ).tocsr() + dFTrans_1



    return FTrans, dFTrans