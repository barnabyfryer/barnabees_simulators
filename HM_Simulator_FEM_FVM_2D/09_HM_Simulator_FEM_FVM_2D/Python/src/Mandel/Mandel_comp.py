import numpy as np

def Mandel_comp(a,b,cf,E,F,k,Ks,mu,phi,t,vd):

    #Beta points
    n = 100

    # =============================================================================
    # Basic calculations
    # =============================================================================

    #x-location
    x = np.linspace(0, a, 100)
    #y-location
    y = np.linspace(0, b, 100)

    # =============================================================================
    # Poroelastic calculations
    # =============================================================================

    #Mobility
    kmob = k/mu
    #Shear Modulus(same drained and undrained)
    G = E / (2 * (1 + vd))
    #Unixial Drained bulk modulus
    Kv = 2 * G * (1 - vd) / (1 - 2 * vd)
    #Bulk Modulus
    K = 2 * G * (vd + 1) / (3 * (1 - 2 * vd))
    #Biot coefficient
    alpha = 1 - K / Ks
    #Fluid Bulk modulus
    Kf = 1 / cf
    #Biot's Modulus
    M = ((alpha - phi) / Ks + phi / Kf) ** (-1)
    #Skempton Coefficient
    B = alpha * M / (K + alpha ** 2 * M)
    #Undrained Poisson's Ratio
    vu = (B * alpha * (1 - 2 * vd) + 3 * vd) / (3 - B * alpha * (1 - 2 * vd))
    #Fluid diffusivity coefficient
    cv = kmob * M * Kv / (Kv + alpha ** 2 * M)
    #Auxillary elastic constant
    nu = (1 - vd) / (vu - vd)

    # =============================================================================
    # Find beta
    # =============================================================================

    beta = np.zeros(n) + np.pi/2

    for i in range(n):
        Error = 1
        while Error > 1e-8:
            h = i-1
            betaold = beta[i]
            beta[i] = np.atan(beta[i]*2*nu) + h*np.pi
            Error = np.abs((betaold - beta[i])/beta[i])

    # =============================================================================
    # Solve for pore pressure
    # =============================================================================

    Pd = np.zeros((len(x),1))
    for i in range(1):
        for j in range(len(x)):
            sum = 0
            for kk in range(n):
                sum = (
                        sum + np.sin(beta[kk]) *
                (np.cos(beta[kk]*x[j]/a) - np.cos(beta[kk]))*
                np.exp(-beta[kk]**2*cv*t/a**2)/
                (beta[kk] - np.sin(beta[kk])*np.cos(beta[kk]))
                )
            Pd[j,i] = 2*sum


    # =============================================================================
    # Find horizontal displacement
    # =============================================================================

    ux = np.zeros((len(x),1))
    for i in range(1):
        for j in range(len(x)):
            sum = 0
            sum1 = 0
            for kk in range(n):
                sum = (
                    sum +
                    np.sin(beta[kk]) * np.cos(beta[kk])/
                    (beta[kk] - np.sin(beta[kk])*np.cos(beta[kk])) *
                    np.exp(-beta[kk]**2*cv*t/a**2)
                )
                sum1 = (
                    sum1 +
                    np.sin(beta[kk]*x[j]/a) * np.cos(beta[kk]) /
                    (beta[kk] - np.sin(beta[kk]) * np.cos(beta[kk])) *
                    np.exp(-beta[kk]**2*cv*t/a**2)
                )
            ux[j,i] = (F*vd/(2*G*a) - F*vu/(G*a)*sum)*x[j] + F/G*sum1


    # =============================================================================
    # Find vertical displacement
    # =============================================================================

    uy = np.zeros((len(x),1))
    for i in range(1):
        for j in range(len(x)):
            sum = 0
            for kk in range(n):
                sum = (
                    (sum + np.sin(beta[kk])*np.cos(beta[kk])*np.exp(-beta[kk]**2*cv*t/(a**2))/
                (beta[kk] - np.sin(beta[kk])*np.cos(beta[kk])))
                )
            uy[j,i] = (-F*(1-vd)/(2*G*a) + F*(1-vu)*sum/(G*a))*y[j]

    # =============================================================================
    # Outputs
    # =============================================================================

    Mandel = {}
    Mandel["x"] = x
    Mandel["y"] = y
    Mandel["Pd"] = Pd
    Mandel["t"] = t
    Mandel["td"] = cv*t/(a**2)
    Mandel["p0"] = -B*(1+vu)*F/(3*a)
    Mandel["tnorm"] = cv/(a**2)

    Mandel["ux0"] = F*vu/(2*G)
    Mandel["ux"] = ux
    Mandel["uy0"] = -F*(1-vu)/(2*G)
    Mandel["uy"] = uy

    return Mandel