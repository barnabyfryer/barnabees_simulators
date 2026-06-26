from src.src_flow.phiCalc import phiCalc
from src.src_flow.density import density

def PStrainUpdate(Flow,Gen,State,State_phi):

    #Porosity before mechanics update
    phi_prior, _ = phiCalc(Flow, State_phi)
    #Density before mechanics update
    rho_prior, _ = density(Flow, State_phi["P"])

    # Cell volume
    V = Gen["dx"] * Gen["dy"] * Gen["dz"]
    #Mass prior to update
    m0 = V * rho_prior * phi_prior

    error = 1
    while error > Gen["tol"]:
        #Calculate new fluid densities
        rho_p, drho_dp = density(Flow, State["P"])
        #Calculate new porosity
        phi_p, dphi_dp = phiCalc(Flow, State)
        #Find error in this term that should be equal to m0
        r = rho_p * phi_p * V - m0
        #Jacobian
        J = (drho_dp * phi_p + rho_p * dphi_dp) * V
        #Find linearized derivative wrt pressure (like dx*dr/dx = -dr)
        dP = -r/J
        #Update solution
        State["P"] = State["P"] + dP
        #Check error
        error = max(abs(r) / (abs(m0) + 1e-12))

    return State