def fix_x(f, k_global,nodesx):

    #Enforce zero x-displacement (u = 0) at given nodes.

    # x-DOFs (0-based indexing in Python)
    dof_x = 2 * nodesx

    # Zero rows and columns
    k_global[dof_x, :] = 0
    k_global[:, dof_x] = 0

    f[dof_x] = 0

    # Enforce constraint
    k_global[dof_x, dof_x] = 1


    return f, k_global