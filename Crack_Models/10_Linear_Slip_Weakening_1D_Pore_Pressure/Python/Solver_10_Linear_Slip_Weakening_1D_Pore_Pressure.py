# coding: utf-8
########################################################################
#
# Nucleation and propagation on pressurized fault with production phase
#
########################################################################
#
# Author: Mathias Lebihain
# Model used in [Fryer et al., JGR: Solid Earth, 2023]
# doi: 10.1029/2022JB025443

##Scientific packages
import numpy as np
import scipy.optimize as optim
import scipy.integrate as integrate
import scipy.interpolate as interpolate
from scipy.special import erfc
##Files managements packages
import os

#Folder creation
if not(os.path.isdir('data')) :
    os.mkdir('data')

if not(os.path.isdir('figures')) :
    os.mkdir('figures')

########################################################################
## Compute static equilibrium
########################################################################
# Note: Implementation following [Garagash & Germanovitch, JGR, 2012]
# doi: 10.1029/2012JB009209

### Notations
## Material properties
# mu_star = mu/(1-nu) (Mode II rupture) or mu (for Mode III rupture)
# mu is the shear modulus of the host medium
# nu is the Poisson's ratio of the host medium
# alpha is the fluid diffusivity along the fault
# f_p is the peak friction
# f_r is the residual friction
# w is the weakening rate
# delta_w is the weakening slip distance delta_w = f_p / w
# delta_r is the critical slip distance delta_r = (f_p - f_r) / w
##Loading properties
# sigma_n is the normal stress acting on the fault
# p_0 is the initial pore pressure on the fault
# Delta_p_inj is the increase of pore pressure during to injection
# Delta_p_prod is the decrease of pore pressure during to production
# dt_prod is the duration of the production phase
# sigma_0 = sigma_n - p_0 is the initial Terzaghi stress
# tau_p = f_p * sigma_0 is the initial strength of the interface
# tau_r = f_r * sigma_0 is the residual strength of the interface

## Normalization
# S: Shear stress is normalized by tau_p
# P: Pressure is normalized by sigma_0 = sigma_n - p_0
# D: Slip is normalized by delta_w = f_p / w
# X: Distance is normalized by a_w = mu_star / tau_p * delta_w
# T: Time is normalized by alpha / a_w^2

## Information
# We use the Piecewise Constant Slip Method (Appendix A2) to consider the influence of a residual friction

##Production phase
#Decrease of pore pressure due to production
Delta_p_prod = -0.7
#Duration of the production phase
dt_prod = 0.000001
##Injection phase
#Increase of pore pressure due to injection
Delta_p_inj = 0.75

#Pressure profile evolution (results from the fluid problem)
Pi = lambda xi : erfc(xi)
Pi_prime = lambda xi : -2./np.sqrt(np.pi)*np.exp(-xi**2)
Pi_inf = 0.
Pi_prime_inf = 0.
#Residual friction
f_r = 0.6
#Critical slip
D_r = (1.-f_r)

# Computed background stress
#tau_b = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, .95])
#tau_b = np.array([0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975])

tau_b = 0.275

tau_b_min = tau_b
#tau_b_min = min(tau_b)

# Computed crack length
#a = np.linspace(0.,10.,num=251, endpoint=True)[1:]
aa = np.linspace(0.,2,num=3*281, endpoint=True)[1:]
bb = np.linspace(2.,30,num=3*281, endpoint=True)
cc = np.linspace(30,700,num=3*281, endpoint=True)
a = np.concatenate((aa,bb,cc), axis=0)
a = np.linspace(0.,5,num=281, endpoint=True)[1:]
aa = np.linspace(0.,2,num=281, endpoint=True)[1:]
bb = np.linspace(2.,50,num=281, endpoint=True)
a = np.concatenate((aa,bb), axis=0)
#a = np.linspace(0.,2,num=3*281, endpoint=True)[1:]


#Use this line to export the slip, put the crack length you want in the first position, need export at end too
a_target = .73
ab = np.linspace(1,a_target,num=100, endpoint=True)
bc = np.linspace(a_target,300,num=100, endpoint=True)[1:]
#a = np.concatenate((ab,bc), axis=0)
ind_target = np.where(a==a_target)

a_ssy = np.array(a.tolist()+np.linspace(a.max(), 100., num=751, endpoint=True)[1:].tolist()+np.logspace(np.log10(100.), np.log10(1000.), num=251, endpoint=True, base=10.)[1:].tolist())

#Case name
case_name = 'FaultPressurization_Pinj'+str(round(Delta_p_inj,3))+'_Pprod'+str(round(-Delta_p_prod,3))+'_Tprod'+str(round(dt_prod,2))+'_fr'+str(round(f_r,2))+'_Taubmin'+str(round(tau_b_min,3))


#Number of points
N = 101
#N = 301

########################################################################
## Piecewise Constant Slip Method
########################################################################

########################################################################
#Preparation
print("----------- Preparation -----------")
print("Preliminary computations for Piecewise method")

#Indices
i = np.arange(N)
j = np.arange(N+1)
#Step
dX = 1./N
#Points of evaluation Y_i (interior) and Y_j (full)
X_i = i*dX
X_j = j*dX
#Matrices
K_ij = -1./(2.*np.pi*dX*((i[:,np.newaxis]-j[np.newaxis,:])**2-0.25)) -1./(2.*np.pi*dX*((i[:,np.newaxis]+j[np.newaxis,:])**2-0.25))
K_ij[:,0] = -1./(2.*np.pi*dX*(i**2-0.25))

#Vectors - Constant pressure
k_const_j = np.zeros(j.size)
integrand_k_const = lambda X : 1. / np.sqrt(1.-X**2)
for n in range(j.size):
    k_const_j[n] = (1./np.pi) * integrate.quad(integrand_k_const, max(-1., (j[n]-0.5)*dX), min(1., (j[n]+0.5)*dX))[0]
    if (n!= 0):
        k_const_j[n] += (1./np.pi) * integrate.quad(integrand_k_const, max(-1., (-j[n]-0.5)*dX), min(1., (-j[n]+0.5)*dX))[0]

##Vectors - Diffusion from production
#Compute
a_over_t = np.array([0.]+np.logspace(-3., -1., num=10, endpoint=False, base=10.).tolist()+np.logspace(-1., 3., num=130, endpoint=False, base=10.).tolist()+np.logspace(3., 4., num=10, endpoint=True, base=10.).tolist())
k_diff_j_of_aOt = np.zeros((a_over_t.size,j.size,))
dk_diff_j_of_aOt = np.zeros((a_over_t.size,j.size,))
for m in range(a_over_t.size):
    print("Pre-computing k_j for Piecewise method: {:d}%".format(int(100.*m/a_over_t.size)), flush=True, end='\r')
    integrand_k_diff = lambda X : Pi(a_over_t[m]*np.abs(X)) / np.sqrt(1.-X**2)
    integrand_dk_diff = lambda X : Pi_prime(a_over_t[m]*np.abs(X)) * np.abs(X) / np.sqrt(1.-X**2)
    for n in range(j.size):
        if ((m==0) | (k_diff_j_of_aOt[m-1,n]!=0.)):
            k_diff_j_of_aOt[m,n] = (2./np.pi) * integrate.quad(integrand_k_diff, max(0., (j[n]-0.5)*dX), min(1., (j[n]+0.5)*dX))[0]
            dk_diff_j_of_aOt[m,n] = (2./np.pi) * integrate.quad(integrand_dk_diff, max(0., (j[n]-0.5)*dX), min(1., (j[n]+0.5)*dX))[0]

#Value at infinity
k_j_of_inf = 0.
dk_j_of_inf = 0.
#Interpolate
interp_k_diff_j = interpolate.interp1d(a_over_t, k_diff_j_of_aOt, axis=0, kind='cubic', copy=True, fill_value=(0., k_j_of_inf), bounds_error=False)
interp_dk_diff_j = interpolate.interp1d(a_over_t, dk_diff_j_of_aOt, axis=0, kind='cubic', copy=True, fill_value=(0., dk_j_of_inf), bounds_error=False)
print("\x1b[KPre-computing k_j for Piecewise method: SUCCESS", flush=True, end='\n')

########################################################################
#Resolution fixing a and tau_b
##Unknowns: the shape coefficients and the time

def Piecewise_in_t(x, a, tau_b, Delta_p_inj, Delta_p_prod, dt_prod):
    #Decompose x
    D_i = np.asarray(x[:-1])
    D_j = np.insert(D_i, [D_i.size], [0.])
    t = x[-1]
    #Check if t is positif
    check_t_prod = t > - dt_prod
    check_t_inj = t > 0.
    #k_j
    if (check_t_prod):
        k_j = k_const_j - Delta_p_prod * interp_k_diff_j(a/(t**2+dt_prod**2)**0.5)
        dk_j_dt = Delta_p_prod * a*t / (t**2+dt_prod**2)**1.5 * interp_dk_diff_j((t**2+dt_prod**2)**0.5)
    else:
        k_j = k_const_j
        dk_j_dt = np.zeros(k_j.size)
    
    if (check_t_inj):
        k_j += - (Delta_p_inj - Delta_p_prod) * interp_k_diff_j(a/t)
        dk_j_dt += (Delta_p_inj - Delta_p_prod) * a/t**2 * interp_dk_diff_j(a/t)
    
    #Filter for zones where the residual friction is reached
    filter_residual = D_i > D_r
    filter_weakening = np.logical_not(filter_residual)
    #Overpressure
    if (check_t_prod):
        P_j = Delta_p_prod * Pi(a*np.abs(X_j)/(t**2+dt_prod**2)**0.5)
        dP_j_dt =  - Delta_p_prod * a*t / (t**2+dt_prod**2)**1.5 * np.abs(X_j) * Pi_prime(a*np.abs(X_j)/(t**2+dt_prod**2)**0.5)
    else:
        P_j = np.zeros(D_j.size)
        dP_j_dt = np.zeros(D_j.size)
    
    if (check_t_inj):
        P_j += (Delta_p_inj - Delta_p_prod) * Pi(a*np.abs(X_j)/t)
        dP_j_dt += - (Delta_p_inj - Delta_p_prod) * a/t**2 * np.abs(X_j) * Pi_prime(a*np.abs(X_j)/t)
    
    P_i = P_j[:-1]
    dP_i_dt = dP_j_dt[:-1]
    #Elastostatic
    KD_i = K_ij.dot(D_j)
    shear_i = tau_b - (1./a)*KD_i
    #Normal stress
    normal_i = (1. - P_i)
    #Friction
    friction_i = np.maximum(f_r, 1. - D_i)
    friction_j = np.insert(friction_i, [friction_i.size], [1.])
    #Write objectives
    F = np.zeros(D_i.size+1)
    F[:-1] = shear_i - friction_i * normal_i
    F[-1] = tau_b - k_j.dot(friction_j)
    #Write jacobian
    J = np.zeros((D_i.size+1,D_i.size+1))
    J[:-1,:-1] = -(1./a)*K_ij[:,:-1] + np.identity(D_i.size)*filter_weakening[np.newaxis,:]*normal_i[:,np.newaxis]
    J[-1,:-1] = k_j[:-1]*filter_weakening
    J[:-1,-1] = friction_i * dP_i_dt
    J[-1,-1] = - dk_j_dt.dot(friction_j)
    
    return F, J

########################################################################
## Small-scale yielding model
########################################################################

########################################################################
#Preparation
print("Preliminary computations for SSY model")

##Vectors - Diffusion from production
#Compute
a_over_t = np.array([0.]+np.logspace(-3., 2., num=251, endpoint=True, base=10.).tolist())
k_diff_of_aOt = np.zeros(a_over_t.size)
dk_diff_of_aOt = np.zeros(a_over_t.size)
for m in range(a_over_t.size):
    print("Pre-computing k_j for Piecewise method: {:d}%".format(int(100.*m/a_over_t.size)), flush=True, end='\r')
    integrand_k_diff = lambda X : Pi(a_over_t[m]*np.abs(X)) / np.sqrt(1.-X**2)
    integrand_dk_diff = lambda X : Pi_prime(a_over_t[m]*np.abs(X)) * np.abs(X) / np.sqrt(1.-X**2)
    k_diff_of_aOt[m] = (1./np.pi) * integrate.quad(integrand_k_diff, -1., 1.)[0]
    dk_diff_of_aOt[m] = (1./np.pi) * integrate.quad(integrand_dk_diff, -1., 1.)[0]

#Value at infinity
a_over_t_inf = np.logspace(4., 12., num=9, endpoint=True, base=10.)
k_of_inf = lambda xi : np.exp(-1)/xi
dk_of_inf = lambda xi : -np.exp(-1)/xi**2
#Interpolate
interp_k_diff = interpolate.interp1d(a_over_t.tolist()+a_over_t_inf.tolist(), k_diff_of_aOt.tolist()+k_of_inf(a_over_t_inf).tolist(), kind='cubic', copy=True, fill_value=(k_diff_of_aOt[0], 0.), bounds_error=False)
interp_dk_diff = interpolate.interp1d(a_over_t.tolist()+a_over_t_inf.tolist(), dk_diff_of_aOt.tolist()+dk_of_inf(a_over_t_inf).tolist(), kind='cubic', copy=True, fill_value=(dk_diff_of_aOt[0], 0.), bounds_error=False)
print("\x1b[KPre-computing k for SSY model: SUCCESS", flush=True, end='\n')

########################################################################
#Resolution fixing a and tau_b
##Unknowns: t, assumed positive

def SSY_in_t(t, a, tau_b, Delta_p_inj, Delta_p_prod, dt_prod):
    #Overpressure
    P = (Delta_p_inj - Delta_p_prod) * Pi(a/t) \
      + Delta_p_prod * Pi(a/(t**2+dt_prod**2)**0.5)
    dP_dt = - (Delta_p_inj - Delta_p_prod) * a/t**2 * Pi_prime(a/t) \
            - Delta_p_prod * a*t/(t**2+dt_prod**2)**1.5 * Pi_prime(a/(t**2+dt_prod**2)**0.5)
    #Effective crack length
    a_eff = a - (np.pi/2.*0.466**2) / (1.-P)
    da_eff_dt = - (np.pi/2.*0.466**2) * dP_dt/ (1.-P)**2
    ##Mode II SIF
    #Diffusive part
    k_II_diff = (Delta_p_inj-Delta_p_prod) * interp_k_diff(a_eff/t) \
              + Delta_p_prod * interp_k_diff(a_eff/(t**2+dt_prod**2)**0.5)
    dk_II_diff_da_eff = (Delta_p_inj-Delta_p_prod) * 1./t * interp_dk_diff(a_eff/t) \
                      + Delta_p_prod * 1./(t**2+dt_prod**2)**0.5 * interp_dk_diff(a_eff/(t**2+dt_prod**2)**0.5)
    dk_II_diff_dt = - (Delta_p_inj-Delta_p_prod) * a/t**2 * interp_dk_diff(a_eff/t) \
                    - Delta_p_prod * a*t/(t**2+dt_prod**2)**1.5 * interp_dk_diff(a_eff/(t**2+dt_prod**2)**0.5)
    #Total SIF
    K_II = np.sqrt(np.pi*a_eff) * (tau_b - f_r + f_r * k_II_diff)
    dK_II_da_eff = 0.5*K_II/a_eff + np.sqrt(np.pi*a_eff) * f_r * dk_II_diff_da_eff
    dK_II_dt = dK_II_da_eff*da_eff_dt + np.sqrt(np.pi*a_eff) * f_r * dk_II_diff_dt
    #Local toughness KIIc
    K_IIc = (1.-f_r) * (1.-P)**0.5
    dK_IIc_dt = - (1.-f_r) * 0.5 * dP_dt / (1.- P)**0.5
    #Write objectives
    F = K_II - K_IIc
    #Write jacobian
    J = dK_II_dt - dK_IIc_dt
    
    return F, J

########################################################################
## Resolution
########################################################################
print("----------- Resolution ------------")

slip = np.full((a.size,N), np.nan)
time = np.full(a.size, np.nan)
slip_at_fault_center = np.full(a.size, np.nan)
x_bounds = [0.]*(N+1), [np.inf]*(N+1)
for n in range(a.size):
    print("Computing Piecewise method with production: {:d}%".format(int(100.*(n+1)/(time.size))), flush=True, end='\r')
    if (n==0):
        #Loop until sucess
        t_guess = 1E-2*a[0]
        bool_success = False
        while ((not bool_success) & (t_guess <= 1E3)):
            #First guess based on last point computed
            x_guess = ([0.]*N)+[t_guess]
            #Solve static problem
            solve_equilibrium = optim.root(Piecewise_in_t, x_guess, jac=True, args=(a[n], tau_b, Delta_p_inj, Delta_p_prod, dt_prod))
            #Update
            bool_success = solve_equilibrium.success
            t_guess *= 2.
        
        if (not bool_success):
            print('\n The initialization could not work for a={:.3f}'.format(a[n]), flush=True, end='\n')
    else:
        #First guess based on last valid point computed
        if np.isnan(time[n-1]):
            n_guess = int(np.arange(a.size)[np.logical_not(np.isnan(time))].max())
        else : 
            n_guess = n-1
        
        x_guess = slip[n_guess].tolist()+[time[n_guess]]
        #Solve static problem
        solve_equilibrium = optim.root(Piecewise_in_t, x_guess, jac=True, args=(a[n], tau_b, Delta_p_inj, Delta_p_prod, dt_prod))
    
    if (solve_equilibrium.success):
        slip_solution = (solve_equilibrium.x)[:-1]
        #Check if the solution has uniform slip sign
        slip_sign = np.sign(slip_solution)
        valid_sign = np.all(slip_sign == slip_sign[0])
        #Check if the solution has continuously decreasing slip from center
        valid_growth = np.all(np.diff(slip_solution) <= 0.)
        #If so store it
        if (valid_sign & valid_growth) :
            slip[n] = (solve_equilibrium.x)[:-1]
            time[n] = (solve_equilibrium.x)[-1]
            slip_at_fault_center[n] = (solve_equilibrium.x)[0]
        else : 
            slip[n,:] = np.nan
            time[n] = np.nan
            slip_at_fault_center[n] = np.nan
    else:
        slip[n,:] = np.nan
        time[n] = np.nan
        slip_at_fault_center[n] = np.nan

print("\x1b[KComputing Piecewise method with production: SUCCESS", flush=True, end='\n')

time_ssy = np.full(a_ssy.size, np.nan)
t_bounds = [1E-12, 1E12]
n_start = np.argmax(a>=3.)
for n in range(n_start, a_ssy.size):
    print("Computing SSY model with production: {:d}%".format(int(100.*(n+1)/(time_ssy.size))), flush=True, end='\r')
    if (n<a.size):
        #First guess from Piecewise method
        t_guess = time[n]
    else:
        #First guess based on last valid point computed
        if np.isnan(time_ssy[n-1]):
            n_guess = int(np.arange(a_ssy.size)[np.logical_not(np.isnan(time_ssy[m]))].max())
        else : 
            n_guess = n-1
        
        t_guess = time_ssy[n_guess]
    
    #Termination condition
    if ((n>n_start+20) & (np.all(np.isnan(time_ssy[max(n_start, n-10):n])))):
        break
    
    #Solve static problem
    solve_equilibrium = optim.root_scalar(SSY_in_t, args=(a_ssy[n], tau_b, Delta_p_inj, Delta_p_prod, dt_prod), method='newton', bracket=t_bounds, x0=t_guess, fprime=True)
    
    if (solve_equilibrium.converged):
        time_ssy[n] = solve_equilibrium.root
    else:
        time_ssy[n] = np.nan

print("\x1b[KComputing SSY model with production: SUCCESS", flush=True, end='\n')

#Export in .csv
#Data for Zenodo
a_03, tau_b_03 = np.meshgrid(a, tau_b)
export_header = "Scenario: Imposed pore pressure; Background stress tau_b/tau_p = {:.3f}; Injection overpressure Delta_p_inj = {:.3f}; Production overpressure Delta_p_prod = {:.2f}; Production duration sqrt(alpha*dt_prod)/a_w = {:.2f}; Residual friction f_r/f_p = {:.2f}".format(tau_b, Delta_p_inj, Delta_p_prod, dt_prod, f_r)
export_column_names = "Background stress tau_b/tau_p; Crack length a/a_w; Time sqrt(alpha*t)/a_w; Time SSY sqrt(alpha*t)/a_w; Slip at the center delta(0)/delta_w; Table shape ({:d},)".format(a.size)
export_table = np.column_stack([tau_b_03.flatten(), a, time, time_ssy[:a.size], slip_at_fault_center])
np.savetxt("data/"+case_name+"_piecewise.csv", export_table, delimiter=";", header = export_header+'\n'+export_column_names, fmt='%.6e')
export_ssy_column_names = "Crack length a/a_w; Time sqrt(alpha*t)/a_w; Table shape ({:d},)".format(a_ssy.size)
export_ssy_table = np.column_stack([a_ssy[n_start:], time_ssy[n_start:]])
np.savetxt("data/"+case_name+"_ssy.csv", export_ssy_table, delimiter=";", header = export_header+'\n'+export_ssy_column_names, fmt='%.6e')

#delta_large_03, X_i_large_03 = np.meshgrid(slip, X_i)
#export_slip_column_names = "Slip; Xi; Table shape ({:d},{:d})".format(slip.size, X_i.size)
#export_slip_table = np.array([slip[ind_target].flatten(), X_i.flatten()])
#export_slip_table = export_slip_table.transpose()
#np.savetxt("data/"+case_name+"_a_"+str(a[ind_target])+"_slip.csv", export_slip_table, delimiter=";", header = export_header+'\n'+export_slip_column_names, fmt='%.6e')