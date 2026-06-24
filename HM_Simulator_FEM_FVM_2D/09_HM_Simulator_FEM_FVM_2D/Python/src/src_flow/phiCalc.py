import numpy as np

def phiCalc(Flow,State):

    #Find permeability using slightly compressible formulation
    phi = Flow["phi0"] * np.exp(Flow["cphi"] * (State["P"] - Flow["phiP0"]))
    #Find derivative wrt Pressure
    dphidP = Flow["cphi"] * phi

    return phi, dphidP