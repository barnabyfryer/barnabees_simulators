import numpy as np

def trans_rock(Flow,Gen):

    #Permeability
    kx = Flow["kx"]

    #Length of cell in x-direction [m]
    dx = Gen["dx"]
    # Length of cell in y-direction [m]
    dy = Gen["Ly"]/1
    # Length of cell in z-direction [m]
    dz = Gen["Lz"] / 1

    #Cross-sectional area in x-direction [m^2]
    Ax = dy*dz
    
    #Rock transmissibilities
    trans = {}
    #Initialize x-direction trans (of left of cell edge)
    trans["x"] = np.zeros(Gen["Nx"]+1)
    #Use harmonic average to find transmissibilities
    trans["x"][1:-1] = 2*kx[:-1] * kx[1:] * Ax / ((kx[:-1] + kx[1:]) * dx)
    #Do boundaries
    trans["x"][0] = kx[0]*Ax / (dx/2)
    trans["x"][-1] = kx[-1] * Ax / (dx/2)



    return trans