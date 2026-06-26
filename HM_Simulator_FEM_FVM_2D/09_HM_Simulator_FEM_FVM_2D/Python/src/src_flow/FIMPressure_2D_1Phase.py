import numpy as np

from src.src_flow.phiCalc import phiCalc
from src.src_flow.trans_rock import trans_rock
from src.src_flow.trans_fluid import trans_fluid
from src.src_flow.upwind import upwind
from src.src_flow.perm import perm
from src.src_flow.build_jacobian import build_jacobian
from src.src_flow.build_residual import build_residual
from src.src_flow.PStrainUpdate import PStrainUpdate

from scipy.sparse.linalg import spsolve


def FIMPressure_2D_1Phase(Flow,Gen,State,State0,State_phi,Wells):

    # =============================================================================
    # Update pressure based on volumetric strain from mechanics model, conserve mass
    # =============================================================================

    State = PStrainUpdate(Flow,Gen,State,State_phi)

    # =============================================================================
    # Solve flow model
    # =============================================================================
    #Convergence checker
    Converged = 0
    #Number of iterations
    it = 1
    while Converged != 1:

        # Find the non-fluid rock transmissibility
        trans,_,_ = trans_rock(Flow, Gen, State)

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
    #Find permeability
    State["kx"],State["ky"],_,_,_,_ = perm(Flow,State)
    #Find porosity
    State["phi"],_ = phiCalc(Flow,State)

    return State