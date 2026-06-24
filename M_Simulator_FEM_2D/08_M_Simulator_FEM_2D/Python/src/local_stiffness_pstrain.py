import numpy as np

def local_stiffness_pstrain(Gen):

    #Length of local coordinate
    a = Gen["dx"]
    b = Gen["dy"]

    #Poisson's ratio
    nu = Gen["nu"]

    #Define constants specific to plane strain
    A = Gen["E"]/((1+nu) * (1-2*nu))
    B = 1-nu
    C = nu
    D = (1-2*nu)/2

    #Predefine memory
    k_local = np.zeros((8, 8))

    #Build each stiffness-matrix component
    # Row 1
    k_local[0, 0] = B * b / (3 * a) + D * (a / (3 * b))
    k_local[0, 1] = C * (1 / 4) + D * (1 / 4)
    k_local[0, 2] = -B * b / (3 * a) + D * (a / (6 * b))
    k_local[0, 3] = C * (1 / 4) + D * (-1 / 4)
    k_local[0, 4] = B * b / (6 * a) + D * (-a / (3 * b))
    k_local[0, 5] = C * (-1 / 4) + D * (1 / 4)
    k_local[0, 6] = -B * b / (6 * a) + D * (-a / (6 * b))
    k_local[0, 7] = C * (-1 / 4) + D * (-1 / 4)

    # Row 2
    k_local[1, 0] = C * (1 / 4) + D * (1 / 4)
    k_local[1, 1] = B * a / (3 * b) + D * (b / (3 * a))
    k_local[1, 2] = C * (-1 / 4) + D * (1 / 4)
    k_local[1, 3] = B * a / (6 * b) + D * (-b / (3 * a))
    k_local[1, 4] = C * (1 / 4) + D * (-1 / 4)
    k_local[1, 5] = -B * a / (3 * b) + D * (b / (6 * a))
    k_local[1, 6] = C * (-1 / 4) + D * (-1 / 4)
    k_local[1, 7] = -B * a / (6 * b) + D * (-b / (6 * a))

    # Row 3
    k_local[2, 0] = -B * b / (3 * a) + D * (a / (6 * b))
    k_local[2, 1] = C * (-1 / 4) + D * (1 / 4)
    k_local[2, 2] = B * b / (3 * a) + D * (a / (3 * b))
    k_local[2, 3] = C * (-1 / 4) + D * (-1 / 4)
    k_local[2, 4] = -B * b / (6 * a) + D * (-a / (6 * b))
    k_local[2, 5] = C * (1 / 4) + D * (1 / 4)
    k_local[2, 6] = B * b / (6 * a) + D * (-a / (3 * b))
    k_local[2, 7] = C * (1 / 4) + D * (-1 / 4)

    # Row 4
    k_local[3, 0] = C * (1 / 4) + D * (-1 / 4)
    k_local[3, 1] = B * a / (6 * b) + D * (-b / (3 * a))
    k_local[3, 2] = C * (-1 / 4) + D * (-1 / 4)
    k_local[3, 3] = B * a / (3 * b) + D * (b / (3 * a))
    k_local[3, 4] = C * (1 / 4) + D * (1 / 4)
    k_local[3, 5] = -B * a / (6 * b) + D * (-b / (6 * a))
    k_local[3, 6] = C * (-1 / 4) + D * (1 / 4)
    k_local[3, 7] = -B * a / (3 * b) + D * (b / (6 * a))

    # Row 5
    k_local[4, 0] = B * b / (6 * a) + D * (-a / (3 * b))
    k_local[4, 1] = C * (1 / 4) + D * (-1 / 4)
    k_local[4, 2] = -B * b / (6 * a) + D * (-a / (6 * b))
    k_local[4, 3] = C * (1 / 4) + D * (1 / 4)
    k_local[4, 4] = B * b / (3 * a) + D * (a / (3 * b))
    k_local[4, 5] = C * (-1 / 4) + D * (-1 / 4)
    k_local[4, 6] = -B * b / (3 * a) + D * (a / (6 * b))
    k_local[4, 7] = C * (-1 / 4) + D * (1 / 4)

    # Row 6
    k_local[5, 0] = C * (-1 / 4) + D * (1 / 4)
    k_local[5, 1] = -B * a / (3 * b) + D * (b / (6 * a))
    k_local[5, 2] = C * (1 / 4) + D * (1 / 4)
    k_local[5, 3] = -B * a / (6 * b) + D * (-b / (6 * a))
    k_local[5, 4] = C * (-1 / 4) + D * (-1 / 4)
    k_local[5, 5] = B * a / (3 * b) + D * (b / (3 * a))
    k_local[5, 6] = C * (1 / 4) + D * (-1 / 4)
    k_local[5, 7] = B * a / (6 * b) + D * (-b / (3 * a))

    # Row 7
    k_local[6, 0] = -B * b / (6 * a) + D * (-a / (6 * b))
    k_local[6, 1] = C * (-1 / 4) + D * (-1 / 4)
    k_local[6, 2] = B * b / (6 * a) + D * (-a / (3 * b))
    k_local[6, 3] = C * (-1 / 4) + D * (1 / 4)
    k_local[6, 4] = -B * b / (3 * a) + D * (a / (6 * b))
    k_local[6, 5] = C * (1 / 4) + D * (-1 / 4)
    k_local[6, 6] = B * b / (3 * a) + D * (a / (3 * b))
    k_local[6, 7] = C * (1 / 4) + D * (1 / 4)

    # Row 8
    k_local[7, 0] = C * (-1 / 4) + D * (-1 / 4)
    k_local[7, 1] = -B * a / (6 * b) + D * (-b / (6 * a))
    k_local[7, 2] = C * (1 / 4) + D * (-1 / 4)
    k_local[7, 3] = -B * a / (3 * b) + D * (b / (6 * a))
    k_local[7, 4] = C * (-1 / 4) + D * (1 / 4)
    k_local[7, 5] = B * a / (6 * b) + D * (-b / (3 * a))
    k_local[7, 6] = C * (1 / 4) + D * (1 / 4)
    k_local[7, 7] = B * a / (3 * b) + D * (b / (3 * a))

    # Apply constitutive prefactor
    k_local *= A

    return k_local