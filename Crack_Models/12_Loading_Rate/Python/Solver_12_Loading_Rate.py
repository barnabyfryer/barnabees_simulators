# Script solving the EoM of (Garagash, Phil. Trans. Roy. Soc, 2021)
# Authors: Mathias Lebihain, Barnaby Fryer, Dmitry Garagash, François Passelègue
# Contact: mathias.lebihain@enpc.fr
# Dependencies: numpy, scipy, jax, jaxellip, diffrax

import argparse
import numpy as np
import jax as jax
import jax.lax as lax
import jax.numpy as jnp
from jax.scipy.special import spence
from jaxellip import ellipe
import diffrax as dfx
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# use double-precision in JAX
jax.config.update("jax_enable_x64", True)

#Function to help with file naming for outputs
def fmt(x):
    return (f"{x:g}"
            .replace(".", "p")
            .replace("-", "m")
            .replace("+", ""))

########################################################################################
# Unwrap parameters

# Direct-to-state friction parameters ratio
a_over_b = 0.9
# Scaled ambiant sliding velocity
V0_over_Vs = 1e-10

# Scaled overstress
Δf0_over_b_i = 1
#Loading rate [friction/ts or d(df0/b)/d(ts), ts = L/V0]
Loading_rate = 1
# Type of R&S law: "slip" or "aging"
rs_type = "slip"

########################################################################################
# Compute associated parameters

# Scaled ambient rupture velocity
v0_over_cs = V0_over_Vs
# Scaled \bar{v}_0
bar_v0_over_cs = jnp.exp(-Δf0_over_b_i) * v0_over_cs

########################################################################################
# EoM functionals

# Safe ellipe for small rupture velocities - E is not differentiable near 1
def safe_m(v):
    # m = 1 - v^2, but never let it hit exactly 1.0
    return jnp.clip(1.0 - v ** 2, 0.0, 1.0 - 1e-14)

# Function g for mode III
g = lambda vr_over_cs : jnp.sqrt(1 - vr_over_cs**2)
# Function k for mode III
k = lambda vr_over_cs : jnp.sqrt(1 - vr_over_cs**2) / ellipe(safe_m(vr_over_cs))

# 𝓕 of P, G21's Eq. (3.12), assuming P = 0
𝓕 = 1
# Prefactor κ0 (between 0.838 and 0.916 from G21) - taken as average value
if (rs_type=='slip'):
    κ0 = 0.838
elif (rs_type=='aging'):
    κ0 = 0.916 # rough approximation
else:
    raise ValueError('The friction law must be the slip law or the aging law.')
# Define real dilog from spence function
Li2 = lambda z : spence(1-z)
# 𝒱 function of G21's Eq. (2.17)
if (rs_type=='slip'):
    𝒱 = lambda ζ : κ0 * jnp.exp(ζ) / jnp.cbrt(ζ)**2
elif (rs_type=='aging'):
    𝒱 = lambda ζ : κ0 * -jnp.exp(ζ) * jnp.cbrt(1-jnp.exp(-ζ)) * jnp.cbrt(Li2(1-jnp.exp(ζ))) / jnp.cbrt(ζ)**4
else:
    raise ValueError('The friction law must be the slip law or the aging law.')
# χ function of G21's Eq. (2.17)
χ = lambda ζ : 𝒱(ζ) * jnp.exp(-ζ)
# Approximate minimum value of 𝒱 from G21
𝒱_min = 2.138

# Invert Δfp_over_b using G21's method
def Δfp_over_b(vr_over_cs):
    # Transform to array
    vr_over_cs = jnp.asarray(vr_over_cs)
    # Value of 𝒱(ζ)
    V = vr_over_cs / bar_v0_over_cs / g(vr_over_cs)
    # Initialize ζ_0
    ζ = jnp.log(V / κ0)
    # Loop
    n_iter = 2
    for i in range(n_iter):
        ζ = jnp.log(V / χ(ζ))
    return ζ

# Scaled fracture energy, G21's Eq. (2.12)
if (rs_type=='slip'):
    Gc = lambda vr_over_cs : Δfp_over_b(vr_over_cs)
elif (rs_type=='aging'):
    Gc = lambda vr_over_cs : - Li2(1 - jnp.exp(Δfp_over_b(vr_over_cs)))
else:
    raise ValueError('The friction law must be the slip law or the aging law.')
# Scaled toughness
Kc = lambda vr_over_cs : jnp.sqrt(2 * Gc(vr_over_cs) * g(vr_over_cs) / k(vr_over_cs)**2)
# Scaled effective velocity
Veff_over_V0 = lambda l_over_lb, vr_over_cs : k(vr_over_cs)/g(vr_over_cs) * 4*Kc(vr_over_cs) * vr_over_cs/v0_over_cs / jnp.sqrt(jnp.pi * l_over_lb) * 𝓕
# Scaled effective stress
Δτ_eff = lambda l_over_lb, vr_over_cs, Δf0_over_b : Δf0_over_b - (a_over_b - 1) * jnp.log(Veff_over_V0(l_over_lb, vr_over_cs))
# Background SIF
K_Δτ = lambda l_over_lb, vr_over_cs, Δf0_over_b : Δτ_eff(l_over_lb, vr_over_cs, Δf0_over_b) * jnp.sqrt(jnp.pi * l_over_lb)

########################################################################################
# EoM solution

# Residue
EoM_objective = lambda l_over_lb, vr_over_cs, Δf0_over_b : K_Δτ(l_over_lb, vr_over_cs, Δf0_over_b) - Kc(vr_over_cs)
EoM_dl = jax.jit(jax.grad(EoM_objective, argnums=0))
EoM_dv = jax.jit(jax.grad(EoM_objective, argnums=1))

# Find rupture velocity for v_r/c_s for any crack length l/l_b - used only to initialize the ODE
def initialize(l_over_lb, method='brentq'):
    # Bounds on tested velocities
    vmin_over_cs = 𝒱_min * bar_v0_over_cs
    vmax_over_cs = 1-1E-12
    # Initial guess for the rupture velocity
    v_guess_over_cs = 0.001
    #Initialize Δf0_over_b at time zero
    Δf0_over_b = Δf0_over_b_i
    # Find root if possible
    if (EoM_objective(l_over_lb, vmin_over_cs, Δf0_over_b) * EoM_objective(l_over_lb, vmax_over_cs, Δf0_over_b) <= 0):
        # Function
        func = lambda vr_over_cs : EoM_objective(l_over_lb, vr_over_cs, Δf0_over_b)
        # Gradient
        grad = lambda vr_over_cs : EoM_dv(l_over_lb, vr_over_cs, Δf0_over_b)
        # Solve for root of the EoM
        results = root_scalar(func, method=method, fprime=grad, \
                                x0=v_guess_over_cs, bracket=(vmin_over_cs, vmax_over_cs))
        # Check convergence
        if (results.converged):
            print("Initialization success.")
            return results.root
        else:
            raise RuntimeError("The root-finding algorithm in l did not converge.")
    else:
        raise RuntimeError("The root-finding algorithm in l could not start. Decrease the initial crack length l_ini/lb.")


# Find rupture velocity for v_r/c_s for any crack length l/l_b - used only to initialize the ODE
def solve_root_in_l(l_over_lb, v_guess_over_cs, Δf0_over_b):
    # Function
    func = lambda vr_over_cs : EoM_objective(l_over_lb, vr_over_cs, Δf0_over_b)
    # Gradient
    grad = lambda vr_over_cs : EoM_dv(l_over_lb, vr_over_cs, Δf0_over_b)
    # Solve for root of the EoM
    results = root_scalar(func, fprime=grad, x0=v_guess_over_cs)
     
    # Check convergence
    if (results.converged):
        return results.root
    else:
        raise RuntimeError("The root-finding in v algorithm did not converge.")

# Find crack length l/l_b for any rupture velocity for v_r/c_s - used to handle the stiff portion of the ODE
def solve_root_in_v(vr_over_cs, l_guess_over_lb, Δf0_over_b):
    # Function
    func = lambda l_over_lb : EoM_objective(l_over_lb, vr_over_cs, Δf0_over_b)
    # Gradient
    grad = lambda l_over_lb : EoM_dl(l_over_lb, vr_over_cs, Δf0_over_b)
    # Solve for root of the EoM
    results = root_scalar(func, fprime=grad, x0=l_guess_over_lb)
    # Check convergence
    if (results.converged):
        return results.root
    else:
        raise RuntimeError("The root-finding in v algorithm did not converge.")

# Find time t/t_s for different crack length l/l_b - using the idea of Dmitry reformulating the fixed points equation into a second-order ODE
def solve_ode_in_l(l_ini_over_lb, l_fin_over_lb, N_steps=1000, Δl_over_lb = 1E-6, rtol=1E-9, atol=1E-9, max_steps=10000):
    # Right-hand side
    def formulate_ode_in_l(l_over_lb, y, args=None):
        # Time t/t_s and its first derivative with respect to crack length l_over_lb
        t_over_ts, t_over_ts_prime = y
        # Rupture velocity
        vr_over_cs = 1/t_over_ts_prime
        #Update overstress
        Δf0_over_b = Δf0_over_b_i + Loading_rate * t_over_ts * V0_over_Vs

        # Kill higher-order derivatives through EoM_dl / EoM_dv for the stiff version
        num = EoM_dl(l_over_lb, vr_over_cs, Δf0_over_b)
        den = EoM_dv(l_over_lb, vr_over_cs, Δf0_over_b)
        
        # Second derivative from the EoM
        t_over_ts_second = 1/vr_over_cs**2 * num / den
        # Total derivative
        d_y = t_over_ts_prime, t_over_ts_second
        return d_y
    # Terms
    terms = dfx.ODETerm(formulate_ode_in_l)
    
    # Initial condition
    t_ini_over_ts = 0
    v_ini_over_cs = initialize(l_ini_over_lb)
    
    # Start with an explicit solver for the non-stiff part
    solver = dfx.Tsit5()
    # Step size controller
    stepsize_controller = dfx.PIDController(rtol=rtol, atol=atol)
    
    # Saved crack length
    l_save_over_lb = jnp.logspace(np.log10(l_ini_over_lb), np.log10(l_fin_over_lb), endpoint=True, num=N_steps, base=10.)
    saveat = dfx.SaveAt(ts=l_save_over_lb)
    
    # Termination procedure if the rupture velocity drops down below ~ the ambient velocity
    v_thres = 𝒱_min * bar_v0_over_cs
    def terminate_if_too_slow(l_over_lb, y, args, **kwargs):
        return 1/y[1] < v_thres
    event_stop = dfx.Event(terminate_if_too_slow)
    
    # Solve
    sol = dfx.diffeqsolve(
        terms,
        solver,
        l_ini_over_lb,
        l_fin_over_lb,
        Δl_over_lb,
        (t_ini_over_ts, 1/v_ini_over_cs),
        saveat=saveat,
        stepsize_controller = stepsize_controller,
        event = event_stop,
        max_steps = max_steps,
        throw = False
    )
    
    # Get solutions
    l_over_lb = sol.ts
    t_over_ts = sol.ys[0]
    vr_over_cs = 1/sol.ys[1]
    Δf0_over_b = Δf0_over_b_i + Loading_rate * t_over_ts * V0_over_Vs
    
    # Heuristic stopping reasons
    # 0: v below threshold
    if (sol.result == dfx.RESULTS.event_occurred):
        reason = 0
    # 1: reached l_fin_over_lb (within tolerance)
    elif (sol.result == dfx.RESULTS.successful):
        reason = 1
    # 2: solver did not reach l_fin_over_lb (likely difficulty); check acceleration
    else:
        # If solver failed, or just didn't reach the endpoint, we may want to switch.
        if (sol.result == dfx.RESULTS.max_steps_reached):
            reason = 2
            print("Switch to a implicit solver to handle the stiff part.")
            # Indice at failure
            i_restart = np.argmax(np.isinf(l_over_lb)) - 1
            # Get bounds on rupture velocity
            v_min_over_cs = vr_over_cs[i_restart]
            v_max_over_cs = 0.8
            # Log-spaced velocities
            N_cont = 100
            v_cont_over_cs = np.logspace(np.log10(v_min_over_cs), np.log10(v_max_over_cs), num=N_cont, endpoint=True, base=10.)
            #Predefine crack length
            l_cont_over_lb = np.zeros(N_cont, dtype=float)
            l_cont_over_lb[0] = float(l_over_lb[i_restart])
            #Predefine overstress
            Δf0_over_b_cont = np.zeros(N_cont, dtype=float)
            Δf0_over_b_cont[0] = float(Δf0_over_b[i_restart])
            #Predefine memory for time
            t_cont_over_ts = np.zeros(N_cont)
            # Compute them
            #Initial time
            t_cont_over_ts[0] = float(t_over_ts[i_restart])
            for i in range(1, N_cont):
                # Initialize time
                t_i = t_cont_over_ts[i-1]
                Err_l = 1
                while Err_l > 1e-8:
                    #Guess overstress
                    Δf0_over_b_guess = Δf0_over_b_i + Loading_rate * t_i * V0_over_Vs
                    #Solve for crack length
                    l_guess_old = solve_root_in_v(v_cont_over_cs[i], l_cont_over_lb[i-1], Δf0_over_b_guess)
                    #Amount of crack growth
                    dl = l_guess_old - l_cont_over_lb[i-1]
                    #Average velocity of step
                    v_ave = (v_cont_over_cs[i] + v_cont_over_cs[i-1])/2
                    #Time increment
                    dt = dl / v_ave
                    #Update time
                    t_i = t_cont_over_ts[i-1] + dt
                    #Update loading rate
                    Δf0_over_b_guess = Δf0_over_b_i + Loading_rate * t_i * V0_over_Vs
                    #Solve for crack length
                    l_guess_new = solve_root_in_v(v_cont_over_cs[i], l_cont_over_lb[i-1], Δf0_over_b_guess)
                    #Check error
                    Err_l = abs(l_guess_old - l_guess_new)/np.max([np.abs(l_guess_new),1])
                #Save time
                t_cont_over_ts[i] = t_i
                #Save crack length
                l_cont_over_lb[i] = l_guess_new
                #Save overstress
                Δf0_over_b_cont[i] = Δf0_over_b_guess

            # Start again with ODE until maximum crack length
            l_end_over_lb = l_save_over_lb[l_save_over_lb>l_cont_over_lb[-1]]
            # Solve
            saveat = dfx.SaveAt(ts=l_end_over_lb)
            sol = dfx.diffeqsolve(
                terms,
                solver,
                l_cont_over_lb[-1],
                l_end_over_lb[-1],
                Δl_over_lb,
                (t_cont_over_ts[-1], 1/v_cont_over_cs[-1]),
                saveat=saveat,
                stepsize_controller = stepsize_controller,
                event = event_stop,
                max_steps = max_steps,
                throw = False
            )
            # Get solutions
            l_end_over_lb = sol.ts
            t_end_over_ts = sol.ys[0]
            v_end_over_cs = 1/sol.ys[1]
            Δf0_over_b_end =  Δf0_over_b_i + Loading_rate * t_end_over_ts * V0_over_Vs
            # Assemble solutions
            l_over_lb = np.concatenate([l_over_lb[:i_restart+1], l_cont_over_lb, l_end_over_lb])
            t_over_ts = np.concatenate([t_over_ts[:i_restart+1], t_cont_over_ts, t_end_over_ts])
            vr_over_cs = np.concatenate([vr_over_cs[:i_restart+1], v_cont_over_cs, v_end_over_cs])
            Δf0_over_b = np.concatenate([Δf0_over_b[:i_restart+1], Δf0_over_b_cont, Δf0_over_b_end])

        else:
            raise RuntimeError("The solver stopped without handled reasons. Check it manually.")
    return t_over_ts, l_over_lb, vr_over_cs, reason

# Solve
l_ini_over_lb, l_fin_over_lb = 3E-5, 1E4
l_ini_over_lb, l_fin_over_lb = 0.5, 1E4
t_over_ts, l_over_lb, vr_over_cs, reason = solve_ode_in_l(l_ini_over_lb, l_fin_over_lb)

# Export file
filename = (
    f"Output/EoM_{rs_type}"
    f"_Dfi_{fmt(Δf0_over_b_i)}"
    f"_dDfdt_{fmt(Loading_rate)}"
    f"_V0_{fmt(V0_over_Vs)}"
    f"_ab_{fmt(a_over_b)}.csv"
)
header = "t_over_ts; l_over_lb; vr_over_cs"
data = np.column_stack((t_over_ts, l_over_lb, vr_over_cs))
np.savetxt(filename, data, delimiter="; ", header=header, comments='')

# Plot crack velocity versus crack length
plt.plot(l_over_lb, vr_over_cs/v0_over_cs, color='black')
plt.xlabel(r'Crack length $\ell/\ell_b$', fontsize=20)
plt.ylabel(r'Crack velocity $v_r/v_0$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1E-1, 1E4)
plt.ylim(0.5, 2e10)
plt.tight_layout()
plt.show()

# Plot time versus rupture velocity
plt.plot(t_over_ts[1:], vr_over_cs[1:], color='black')
plt.xlabel(r'Time $t/t_s$', fontsize=20)
plt.ylabel(r'Rupture velocity $v_r/c_s$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.ylim(bar_v0_over_cs, 1.2)
plt.tight_layout()
plt.show()

# Plot time versus rupture velocity
plt.plot(7.32489240602825e-03 * t_over_ts[1:], 1345 * vr_over_cs[1:], color='black')
plt.xlabel(r'Time $t$ (s)', fontsize=20)
plt.ylabel(r'Rupture velocity $v_r$ (m/s)', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.03, 20)
plt.ylim(0.01, 1000)
plt.tight_layout()
plt.show()

