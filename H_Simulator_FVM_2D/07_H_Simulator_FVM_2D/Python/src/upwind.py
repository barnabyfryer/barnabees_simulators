import numpy as np
from scipy.sparse import diags

def upwind(Gen,State,Trans):
    #Upwind calculator
    #Returns structure of 2 matrices containing whether a cell (represented by the row) is the upstream direction for flow
    #between it and the cell (column) on the right (or top) of row cell

    #Basic definitions
    Nx = Gen["Nx"]
    Ny = Gen["Ny"]
    N = Nx * Ny
    P = np.asarray(State["P"]).reshape(Nx,Ny, order='F')

    #Find velocities at edges. These velocities are still missing all fluid properties which need to be upwinded (kr, rho, mu)
    Ux = np.zeros((Nx+1,Ny))
    Uy = np.zeros((Nx, Ny+1))
    #Velocity at cell edge in x-direction. If positive flow to right
    Ux[1:Nx,:] = (P[0:Nx-1,:] - P[1:Nx,:]) * Trans["x"][1:Nx,:]
    # Velocity at cell edge in y-direction. If positive flow to top
    Uy[:,1:Ny] = (P[:,0:Ny - 1] - P[:,1:Ny]) * Trans["y"][:,1:Ny]

    #Use velocity to build upwind operator
    #The U.x >0 returns a 1 for all locations with non-negative velocity
    R = (Ux[1:Nx+1,:] >= 0).astype(np.float64).reshape(-1, order='F')
    #Then the vector is reshaped to be the size of the reservoir
    L = (Ux[0:Nx,:] < 0).astype(np.float64).reshape(-1, order='F')
    T = (Uy[:,1:Ny + 1] >= 0).astype(np.float64).reshape(-1, order='F')
    B = (Uy[:,0:Ny] < 0).astype(np.float64).reshape(-1, order='F')

    #Stores the vectors containing if the cell is upstream on its right and left interfaces
    A = {}
    A["x"] = diags([R, L[1:]],[0, 1],shape=(N, N)).tocsr()
    A["y"] = diags([T, B[Nx:]], [0, Nx], shape=(N, N)).tocsr()

    # print("Ax=", A["x"].toarray())
    # print("Ay=", A["y"].toarray())

    return A