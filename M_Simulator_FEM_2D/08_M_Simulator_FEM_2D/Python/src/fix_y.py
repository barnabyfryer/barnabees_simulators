def fix_y(f, k_global,nodesy):

    #Enforce zero y-displacement (v = 0) at given nodes.

    # y-DOFs
    dof_y = 2 * nodesy + 1

    # Zero rows
    k_global[dof_y, :] = 0
    k_global[:, dof_y] = 0

    f[dof_y] = 0

    # Enforce constraint
    k_global[dof_y, dof_y] = 1

    return f, k_global