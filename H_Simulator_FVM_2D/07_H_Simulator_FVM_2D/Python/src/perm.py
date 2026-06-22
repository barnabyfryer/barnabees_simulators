import numpy as np

def perm(Flow,State):

    #Find permeability using slightly compressible formulation
    kx = Flow["kx0"] * np.exp(Flow["ck"] * (State["P"] - Flow["kP0"]))
    ky = Flow["ky0"] * np.exp(Flow["ck"] * (State["P"] - Flow["kP0"]))
    #Find an effective permeability of the cell for use in wells
    keff = (kx * ky) ** 0.5
    #Find derivative wrt Pressure
    dkx_dP = Flow["ck"] * kx
    dky_dP = Flow["ck"] * kx
    #Derivative wrt pressure for effective permeability
    dkeff_dP = 0.5 * (kx * ky) ** (-0.5) * (kx * dky_dP + ky * dkx_dP)

    return kx, ky, keff, dkx_dP, dky_dP, dkeff_dP
