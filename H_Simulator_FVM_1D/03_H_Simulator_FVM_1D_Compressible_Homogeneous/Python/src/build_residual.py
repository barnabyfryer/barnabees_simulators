from src.density import density
from src.upwind import upwind
from src.trans_fluid import trans_fluid
from src.add_bc import add_bc

def build_residual(Flow,Gen,State,State0,trans):

    A = upwind(Gen,State,trans)

    #Old time step density
    RhoOld,_ = density(Flow,State0["P"])
    #Current iteration density
    Rho,dRhodP = density(Flow,State["P"])

    #Transmissibility matrix
    FTrans,_ = trans_fluid(A,dRhodP,Flow,Gen,Rho,trans)

    #Cross sectional area in x direction
    Ax = Gen["Ly"]/1 * Gen["Lz"]/1
    #Cell volume
    V = Gen["Ly"]/1 * Gen["Lz"]/1 * Gen["dx"]

    # =============================================================================
    # Accumulation terms
    # =============================================================================

    Acc = (Rho - RhoOld) * Flow["phi"] * (V/Gen["tstep"])

    # =============================================================================
    # Convection terms
    # =============================================================================

    Conv = FTrans @ State["P"]

    # =============================================================================
    # Add boundary conditions
    # =============================================================================

    BC,_ = add_bc(dRhodP,Gen,Rho,State,trans)

    # =============================================================================
    # Find residual [kg/sec]
    # =============================================================================

    res = Acc + Conv + BC

    return res