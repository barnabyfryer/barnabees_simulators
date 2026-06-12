import numpy as np

from src.trans_rock import trans_rock
from src.trans_fluid import trans_fluid
from src.upwind import upwind
from src.density import density
from src.build_jacobian import build_jacobian
from src.build_residual import build_residual
from scipy.sparse.linalg import spsolve


def FIMPressure_1D_1Phase(Flow,Gen,State,Wells):

    #Find the non-fluid rock transmissibility
    trans = trans_rock(Flow,Gen)
    #Store the state variable from previous time step for use in residual
    State0 = {"P": State["P"].copy()}

    # =============================================================================
    # Solve flow model
    # =============================================================================
    #Convergence checker
    Converged = 0
    #Number of iterations
    it = 1
    while Converged != 1:

        #Build residual
        res = build_residual(Flow,Gen,State,State0,trans,Wells)
        #Build Jacobian
        jac = build_jacobian(Flow,Gen,State,trans,Wells)

        #Solve
        x = spsolve(jac, -res)

        #Update solution
        State["P"] = State["P"] + x

        # Build residual
        res = build_residual(Flow, Gen, State, State0, trans,Wells)

        #Check convergence
        #Computes infinite norm of residual
        norm = np.linalg.norm(res, ord=np.inf)
        if norm < Gen["tol"]:
            #Set to converged
            Converged = 1
        else:
            #Update iteration tracker
            it += 1

    # =============================================================================
    # Find fluid flux
    # =============================================================================

    #Upwind
    A = upwind(Gen,State,trans)
    #Fluid transmissivity
    Ftrans,_ = trans_fluid(A,Flow,Gen,State)
    #Fluid flux into cell
    State["flux"] = -Ftrans @ State["P"]

    return State