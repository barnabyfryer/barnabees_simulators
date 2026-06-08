import numpy as np

def density(Flow,P):


    cf = Flow["cf"] #Compressibility [1/Pa]
    Rho0 = Flow["Rho0"] #Reference density [kg/m^3]
    RhoP = Flow["RhoP"] #Reference pressure [Pa]

    Rho = Rho0 * np.exp(cf * (P - RhoP))
    dRhodP = Rho * cf #Derivative wrt pressure

    return Rho, dRhodP