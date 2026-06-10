import numpy as np

def add_bc(dRhodP, Gen, Rho, State, trans):

    #Predefine memory
    BC = np.zeros(Gen["Nx"])
    dBC = np.zeros(Gen["Nx"])

    #Find "source" terms from boundary conditions
    BC[0] = (trans["x"][0] / 2) * (Rho[0] * (State["P"][0] - Gen["PL"]))
    BC[-1] = (trans["x"][-1] / 2) * (Rho[-1] * (State["P"][-1] - Gen["PR"]))

    dBC[0] = (trans["x"][0] / 2) * (Rho[0] + dRhodP[0] * (State["P"][0] - Gen["PL"]))
    dBC[-1] = (trans["x"][-1] / 2) * (Rho[-1] + dRhodP[-1] * (State["P"][-1] - Gen["PR"]))

    return BC, dBC