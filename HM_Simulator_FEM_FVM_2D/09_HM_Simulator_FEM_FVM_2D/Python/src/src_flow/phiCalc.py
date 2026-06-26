import numpy as np

def phiCalc(Flow,State):


    #Reference porosity
    phi0 = Flow["phi0"]
    #Volumetric strain
    e_vol = State["e_vol"]
    #Pressure
    P = State["P"]

    #Porosity change from deformation
    phi_mech = phi0 - (1 - phi0) * e_vol

    #Find permeability using slightly compressible formulation
    phi = phi_mech * np.exp(Flow["cphi"] * (P - Flow["phiP0"]))
    #Find derivative wrt Pressure
    dphidP = Flow["cphi"] * phi

    return phi, dphidP