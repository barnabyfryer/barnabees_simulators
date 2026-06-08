def trans_rock(Flow,Gen):

    #Permeability
    kx = Flow["kx"]

    #Length of cell in x-direction [m]
    dx = Gen["dx"]
    # Length of cell in y-direction [m]
    dy = Gen["Ly"]/1
    # Length of cell in z-direction [m]
    dz = Gen["Lz"] / 1

    #Cross-sectional area in x-direction [m^2]
    Ax = dy*dz
    
    #Rock transmissibilities
    trans = {}
    trans["x"] = kx*Ax/dx


    return trans