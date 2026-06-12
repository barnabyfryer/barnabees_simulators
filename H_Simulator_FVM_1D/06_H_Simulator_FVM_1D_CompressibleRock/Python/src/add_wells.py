import numpy as np
from src.density import density
from src.perm import perm


def add_wells(Flow,State,Wells):

    N = len(State["P"])

    # =============================================================================
    # Add fixed-pressure wells
    # =============================================================================

    #Find permeability
    kx,dkxdP = perm(Flow,State)

    Pwf = np.zeros(N)  #Wellbore flowing pressure [Pa]
    Pwf[Wells["Loc_P"]] = Wells["P"]        #Only application location with wells

    WI = np.zeros(N)  # Well index [m^3]
    WI[Wells["Loc_P"]] = Wells["WI"]  # Only application location with wells

    Rho = np.zeros(N)  # Density [kg/m^3]
    dRhodP = np.zeros(N)  # Derivative of density wrt pressure [kg/(Pa m^3)]

    #Upstreamed density of fluid being injected
    inj = Pwf >= State["P"]
    if np.any(inj):
        Rho[inj], dRhodP[inj] = density(
            Flow,
            Pwf[inj]
        )

        # Injected fluid density does not depend on reservoir pressure
        dRhodP[inj] = 0.0

    #Upstreamed density of fluid being produced
    prod = Pwf < State["P"]
    if np.any(prod):
        Rho[prod], dRhodP[prod] = density(
            Flow,
            State["P"][prod]
        )

    #Find flow rate and derivatives
    Q = kx * (Pwf - State["P"]) * WI * Rho / Flow["muf"] #Flow rate [kg/sec]
    dQdP = kx * WI * (dRhodP * (Pwf - State["P"]) - Rho) / Flow["muf"] + dkxdP * WI * Rho / Flow["muf"]  #Derivative wrt pressure [kg/(Pa sec)]

    # =============================================================================
    # Add fixed-rate wells
    # =============================================================================

    Q_wells = np.zeros(N)
    Q_wells[Wells["Loc_Q"]] = Wells["Q"]
    #Add this contribution to source terms
    Q += Q_wells


    return Q, dQdP