import numpy as np

def perm(Flow,State):

    #Find permeability using slightly compressible formulation
    kx = Flow["kx0"] * np.exp(Flow["ck"] * (State["P"] - Flow["kP0"]))
    #Find derivative wrt Pressure
    dkdP = Flow["ck"] * kx

    return kx, dkdP
