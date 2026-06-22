from src.density import density
from src.upwind import upwind
from src.add_wells import add_wells
from src.trans_fluid import trans_fluid
from src.phiCalc import phiCalc

def build_residual(Flow,Gen,State,State0,trans,Wells):

    A = upwind(Gen,State,trans)

    #Old time step density
    RhoOld,_ = density(Flow,State0["P"])
    #Current iteration density
    Rho,_ = density(Flow,State["P"])
    # Find porosity and derivative
    phi, _ = phiCalc(Flow, State)

    #Transmissibility matrix
    FTrans,_ = trans_fluid(A,Flow,Gen,State)

    #Cell volume
    V = Gen["dx"] * Gen["dy"] * Gen["dz"]

    # =============================================================================
    # Accumulation terms
    # =============================================================================

    Acc = (Rho - RhoOld) * phi * (V/Gen["tstep"])

    # =============================================================================
    # Convection terms
    # =============================================================================

    Conv = FTrans @ State["P"]

    # =============================================================================
    # Add source terms
    # =============================================================================

    Q,_ = add_wells(Flow,State,Wells)

    # =============================================================================
    # Find residual [kg/sec]
    # =============================================================================

    res = Acc + Conv - Q

    return res