import numpy as np

from src.perm import perm

def trans_rock(Flow,Gen,State):

    #Permeability
    kx,dkxdP = perm(Flow,State)
    kx = np.reshape(kx,(Gen["Nx"],1))
    kx = kx.squeeze()
    dkxdP = np.reshape(dkxdP, (Gen["Nx"], 1))
    dkxdP = dkxdP.squeeze()

    #Length of cell in x-direction [m]
    dx = Gen["dx"]
    # Length of cell in y-direction [m]
    dy = Gen["Ly"]/1
    # Length of cell in z-direction [m]
    dz = Gen["Lz"] / 1

    #Cross-sectional area in x-direction [m^2]
    Ax = dy*dz
    
    #Rock transmissibilities and derivatives
    trans = {}
    dTransdP_LB = {}
    dTransdP_RT = {}

    # =============================================================================
    # Find transmissibilities
    # =============================================================================

    #Initialize x-direction trans (of left of cell edge)
    trans["x"] = np.zeros(Gen["Nx"]+1)
    dTransdP_LB["x"] = np.zeros(Gen["Nx"] + 1)
    dTransdP_RT["x"] = np.zeros(Gen["Nx"] + 1)
    #Use harmonic average to find transmissibilities
    trans["x"][1:-1] = 2*kx[:-1] * kx[1:] * Ax / ((kx[:-1] + kx[1:]) * dx)
    #Do boundaries
    trans["x"][0] = kx[0]*Ax / (dx/2)
    trans["x"][-1] = kx[-1] * Ax / (dx/2)

    # =============================================================================
    # Find derivative of transmissibilities wrt pressure
    # =============================================================================

    #Derivative wrt to left or bottom cell
    dTransdP_LB["x"][1:Gen["Nx"]] = (2*Ax/dx) * kx[1:]**2 / (kx[:Gen["Nx"]-1] + kx[1:])**2 * dkxdP[:Gen["Nx"]-1]
    #Right-most interface
    dTransdP_LB["x"][-1] = (2*Ax/dx) * dkxdP[-1]

    # Derivative wrt to right or top cell
    dTransdP_RT["x"][1:Gen["Nx"]] = (2 * Ax / dx) * kx[0:Gen["Nx"]-1] ** 2 / (kx[:Gen["Nx"]-1] + kx[1:]) ** 2 * dkxdP[1:]
    # Left-most interface
    dTransdP_RT["x"][0] = (2 * Ax / dx) * dkxdP[0]


    return trans, dTransdP_LB, dTransdP_RT