import numpy as np
from src.src_mech.local_stiffness_pstrain import local_stiffness_pstrain
from src.src_mech.local_stiffness_pstress import local_stiffness_pstress

def global_stiffness_builder(Gen):

    #Predefine memory of global stiffness matrix
    k_global = np.zeros((2*Gen["Nn"],2*Gen["Nn"]))

    #Catch ill-defined plane strain or stress
    if Gen["plane"] not in ["strain", "stress"]:
        raise ValueError("Gen['plane'] must be 'strain' or 'stress'")
    #Prepare local stiffness matrix of element (only need to do once since elements same size)
    if Gen["plane"] == "strain":
        #For plane strain conditions
        kl = local_stiffness_pstrain(Gen)
    else:
        #For plane stress conditions
        kl = local_stiffness_pstress(Gen)

    #Build global stiffness matrix
    for i in range(Gen["Ne"]):
        #Find nodes in element
        nodes = Gen["Ref"][i]

        # Global degree-of-freedom indices for this element
        # Ordering: [u1, v1, u2, v2, u3, v3, u4, v4]
        dof = np.array([
            2 * nodes[0], 2 * nodes[0] + 1,
            2 * nodes[1], 2 * nodes[1] + 1,
            2 * nodes[2], 2 * nodes[2] + 1,
            2 * nodes[3], 2 * nodes[3] + 1
        ])

        # Add local stiffness contributions to the global DOFs associated with this element
        k_global[np.ix_(dof, dof)] += kl

    return k_global