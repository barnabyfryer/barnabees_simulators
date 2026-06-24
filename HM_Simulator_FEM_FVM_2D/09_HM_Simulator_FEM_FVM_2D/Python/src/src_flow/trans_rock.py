import numpy as np

from src.src_flow.perm import perm

def trans_rock(Flow,Gen,State):

    #Permeability
    kx,ky,_,dkx_dP,dky_dP,_ = perm(Flow,State)
    kx = np.reshape(kx,(Gen["Nx"],Gen["Ny"]))
    kx = kx.squeeze()
    ky = np.reshape(ky, (Gen["Nx"], Gen["Ny"]))
    ky = ky.squeeze()
    dkx_dP = np.reshape(dkx_dP, (Gen["Nx"], Gen["Ny"]))
    dkx_dP = dkx_dP.squeeze()
    dky_dP = np.reshape(dky_dP, (Gen["Nx"], Gen["Ny"]))
    dky_dP = dky_dP.squeeze()

    #Length of cell in x-direction [m]
    dx = Gen["dx"]
    # Length of cell in y-direction [m]
    dy = Gen["dy"]
    # Length of cell in z-direction [m]
    dz = Gen["Lz"] / 1

    #Cross-sectional area in x-direction [m^2]
    Ax = dy*dz
    #Cross-sectional area in y-direction [m^2]
    Ay = dx*dz
    
    #Rock transmissibilities and derivatives
    trans = {}
    dTransdP_LB = {}
    dTransdP_RT = {}

    # =============================================================================
    # Find transmissibilities
    # =============================================================================

    #Initialize x-direction trans (of left of cell edge)
    trans["x"] = np.zeros((Gen["Nx"]+1,Gen["Ny"]))
    dTransdP_LB["x"] = np.zeros((Gen["Nx"] + 1,Gen["Ny"]))
    dTransdP_RT["x"] = np.zeros((Gen["Nx"] + 1,Gen["Ny"]))
    # Initialize y-direction trans (of bottom of cell edge)
    trans["y"] = np.zeros((Gen["Nx"],Gen["Ny"] + 1))
    dTransdP_LB["y"] = np.zeros((Gen["Nx"],Gen["Ny"] + 1))
    dTransdP_RT["y"] = np.zeros((Gen["Nx"],Gen["Ny"] + 1))
    #Use harmonic average to find transmissibilities
    trans["x"][1:-1,:] = 2*kx[:-1,:] * kx[1:,:] * Ax / ((kx[:-1,:] + kx[1:,:]) * dx)
    trans["y"][:,1:-1] = 2 * ky[:,:-1] * ky[:,1:] * Ay / ((ky[:,:-1] + ky[:,1:]) * dy)
    #Do boundaries
    trans["x"][0,:] = kx[0,:]*Ax / (dx/2)
    trans["x"][-1,:] = kx[-1,:] * Ax / (dx/2)
    trans["y"][:,0] = ky[:,0] * Ay / (dy / 2)
    trans["y"][:,-1] = ky[:,-1] * Ay / (dy / 2)

    # =============================================================================
    # Find derivative of transmissibilities wrt pressure
    # =============================================================================

    #Derivative wrt to left or bottom cell
    dTransdP_LB["x"][1:Gen["Nx"],:] = (2*Ax/dx) * kx[1:,:]**2 / (kx[:Gen["Nx"]-1,:] + kx[1:,:])**2 * dkx_dP[:Gen["Nx"]-1,:]
    #Right-most interface
    dTransdP_LB["x"][-1,:] = (2*Ax/dx) * dkx_dP[-1,:]
    # Derivative wrt to left or bottom cell
    dTransdP_LB["y"][:,1:Gen["Ny"]] = (2 * Ay / dy) * ky[:,1:] ** 2 / (ky[:,:Gen["Ny"] - 1] + ky[:,1:]) ** 2 * dky_dP[:,:Gen["Ny"] - 1]
    #
    dTransdP_LB["y"][:,-1] = (2 * Ay / dy) * dky_dP[:,-1]

    # Derivative wrt to right or top cell
    dTransdP_RT["x"][1:Gen["Nx"],:] = (2 * Ax / dx) * kx[0:Gen["Nx"]-1,:] ** 2 / (kx[:Gen["Nx"]-1,:] + kx[1:,:]) ** 2 * dkx_dP[1:,:]
    # Left-most interface
    dTransdP_RT["x"][0,:] = (2 * Ax / dx) * dkx_dP[0,:]
    # Derivative wrt to right or top cell
    dTransdP_RT["y"][:,1:Gen["Ny"]] = (2 * Ay / dy) * ky[:,0:Gen["Ny"] - 1] ** 2 / (ky[:,:Gen["Ny"] - 1] + ky[:,1:]) ** 2 * dky_dP[:,1:]
    #
    dTransdP_RT["y"][:,0] = (2 * Ay / dy) * dky_dP[:,0]


    return trans, dTransdP_LB, dTransdP_RT