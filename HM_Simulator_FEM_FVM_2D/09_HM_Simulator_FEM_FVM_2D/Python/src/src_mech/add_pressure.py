import numpy as np

def add_pressure(f,Gen,State):
    P = State["P"]

    for e in range(Gen["Ne"]):
        p = Gen["biot"]*P[e]

        F_app_x = p * Gen["dy"] / 2
        F_app_y = p * Gen["dx"] / 2

        nodes = Gen["Ref"][e]

        BL, BR, TL, TR = nodes


        left_nodes = [BL, TL]
        right_nodes = [BR, TR]
        bottom_nodes = [BL, BR]
        top_nodes = [TL, TR]

        # x-direction
        f[2 * np.array(right_nodes)] += F_app_x
        f[2 * np.array(left_nodes)] -= F_app_x

        # y-direction
        f[2 * np.array(top_nodes) + 1] += F_app_y
        f[2 * np.array(bottom_nodes) + 1] -= F_app_y


    return f