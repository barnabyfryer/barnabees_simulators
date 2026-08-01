import numpy as np

def phiCalc(Flow,Gen,State):

    # Reference porosity
    phi0 = Flow["phi0"]
    # Volumetric strain
    e_vol = State["e_vol"]
    # Pressure
    P = State["P"]

    if Gen["Mandel"] == 1:
        #Porosity and its derivative wrt to pressure after Coussy 2007; assume initial e_vol = 0
        phi = -Gen["biot"] * (e_vol + e_vol*0) + P * (Gen["biot"] - phi0) / Flow["ks"] + phi0

        #Derivative wrt pressure
        dphidP = (Gen["biot"] - phi0) / Flow["ks"]

    else:
        #Porosity change from deformation
        phi_mech = phi0 - (1 - phi0) * e_vol

        #Find permeability using slightly compressible formulation
        phi = phi_mech * np.exp(Flow["cphi"] * (P - Flow["phiP0"]))
        #Find derivative wrt Pressure
        dphidP = Flow["cphi"] * phi

    return phi, dphidP