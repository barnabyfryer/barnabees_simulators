import numpy as np

def solve_stresses(Gen, Pos):

    # =============================================================================
    # Unpack
    # =============================================================================

    #Poisson's ratio
    nu = Gen["nu"]
    #Horizontal and vertical displacements
    u = Pos["u"].reshape(-1)
    v = Pos["v"].reshape(-1)
    #Original node locations but in vector form
    x = Pos["x"].reshape(-1)
    y = Pos["y"].reshape(-1)

    #Predefine for strain
    Eps = np.zeros((Gen["Ne"], 3))

    # =============================================================================
    # Decide plane strain or plane stress
    # =============================================================================

    if Gen["plane"] == "strain":
        #For plane strain conditions
        AA = Gen["E"] / ((1 + nu) * (1 - 2 * nu))
        BB = 1 - nu
        CC = nu
        DD = (1 - 2 * nu) / 2
    else:
        #For plane stress conditions
        AA = Gen["E"] / (1 - nu ** 2)
        BB = 1
        CC = nu
        DD = (1 - nu) / 2

    D = np.array([
        [BB, CC, 0],
        [CC, BB, 0],
        [0, 0, DD]
    ]) * AA

    # =============================================================================
    # Preallocate memory
    # =============================================================================

    Ne = Gen["Ne"]

    E_dx = np.zeros((Ne, 4))
    E_dy = np.zeros((Ne, 4))
    Loc_x = np.zeros((Ne, 4))
    Loc_y = np.zeros((Ne, 4))

    # =============================================================================
    # Element extraction loop
    # =============================================================================

    for i in range(Ne):
        nodes = Gen["Ref"][i]
        BL, BR, TL, TR = nodes

        # displacements
        E_dx[i, 0] = u[BL]
        E_dx[i, 1] = u[BL + 1]
        E_dx[i, 2] = u[TL]
        E_dx[i, 3] = u[TL + 1]

        E_dy[i, 0] = v[BL]
        E_dy[i, 1] = v[BL + 1]
        E_dy[i, 2] = v[TL]
        E_dy[i, 3] = v[TL + 1]

        # coordinates
        Loc_x[i, 0] = x[BL]
        Loc_x[i, 1] = x[BL + 1]
        Loc_x[i, 2] = x[TL]
        Loc_x[i, 3] = x[TL + 1]

        Loc_y[i, 0] = y[BL]
        Loc_y[i, 1] = y[BL + 1]
        Loc_y[i, 2] = y[TL]
        Loc_y[i, 3] = y[TL + 1]

    # =============================================================================
    # Stress loop
    # =============================================================================

    eta = 0.0
    xi = 0.0

    Sigma = np.zeros((Ne, 3))

    for i in range(Ne):

        coord = np.column_stack((Loc_x[i, :], Loc_y[i, :]))

        disp = np.zeros(8)
        disp[0::2] = E_dx[i, :]
        disp[1::2] = E_dy[i, :]

        # shape function derivatives at (xi, eta)
        dN = np.zeros((2, 4))

        dN[0, 0] = -0.25 * (1 - eta)
        dN[1, 0] = -0.25 * (1 - xi)

        dN[0, 1] = 0.25 * (1 - eta)
        dN[1, 1] = -0.25 * (1 + xi)

        dN[0, 2] = -0.25 * (1 + eta)
        dN[1, 2] = 0.25 * (1 - xi)

        dN[0, 3] = 0.25 * (1 + eta)
        dN[1, 3] = 0.25 * (1 + xi)

        Jac = dN @ coord
        deriv = np.linalg.solve(Jac, dN)

        B = np.zeros((3, 8))

        for a in range(4):
            B[0, 2 * a] = deriv[0, a]
            B[1, 2 * a + 1] = deriv[1, a]
            B[2, 2 * a] = deriv[1, a]
            B[2, 2 * a + 1] = deriv[0, a]

        strain = B @ disp
        Eps[i, :] = strain
        Sigma[i, :] = D @ strain

    # =============================================================================
    # Stress loop
    # =============================================================================

    Eps_xx = Eps[:, 0].reshape(Gen["Ny"], Gen["Nx"])
    Eps_yy = Eps[:, 1].reshape(Gen["Ny"], Gen["Nx"])

    eps_vol = Eps_xx + Eps_yy

    return Sigma, eps_vol