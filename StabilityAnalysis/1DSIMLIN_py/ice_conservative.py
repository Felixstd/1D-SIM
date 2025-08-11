"""
Here we want to solve the problem in conservative form

rho(uh)_t + (hu^2 - sigma11)_x = 0
h_t + (uh)_x = 0


Just like Gray, 1999. 

This way we have advection and we are solving the complete equations.
I am curious to see if we get fingers for mu_b > 1. 


"""

import os
import sys

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# import calc_derivatives as deriv
import init_cond as ini
import advection as adv

from scipy.linalg import solveh_banded
from global_constants import *
from numba import jit, prange
from tqdm import tqdm

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


@jit(nopython=True)
def first_derivative(x, deltax, N) : 
    """
    Function to evaluate the first derivative of an array x. 
    
    This is the formula for a second order accuracy. 

    Args:
        x (_type_): _description_
        deltax (_type_): _description_
        N (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    derivative = np.zeros_like(x)
    # print('here')
    
    derivative[1:N-1] = (x[2:] - x[:-2])/(2*deltax)
    
    derivative[0] = (-3*x[0] + 4*x[1] - x[2])/(2*deltax)
    derivative[-1] = (3*x[-1] - 4*x[-2] + x[-3])/(2*deltax)
    
    return derivative

@jit(nopython= True)
def sigma_VP(u, A, zeta_max = zeta_max):
    
    dudx = first_derivative(u, dx, Nx)
    P = Pstar*A*np.exp(-C*(1-A))
    delta = np.abs(dudx)*np.sqrt(1+e**(-2))
    
    zeta_max_P = np.empty(len(P)+1)
    zeta_max_P[0] = zeta_max
    zeta_max_P[1:] = P/(2*delta)

    # zeta_max = np.array([zeta_max])
    zeta = np.nanmin(zeta_max_P)
    
    sigma_11 = zeta*(1+e**-2)*dudx - P/2
    
    return sigma_11

@jit(nopython= True)
def sigma_GC(u, A, zeta_max = zeta_max, mub = mub):
    
    dudx = first_derivative(u, dx, Nx)
    P_x = Pstar*A*np.exp(-C*(1-A))
    I = np.empty(len(P_x))
    mu = np.empty(len(P_x))
    zeta = np.empty(len(P_x))
    eta = np.empty(len(P_x))
    for i, P in enumerate(P_x):
        
        I = np.nanmin([1,dmean*abs(dudx[i])*np.sqrt(rhoice*A[i]/(P+1e-12))])
        mu = mu0 - (deltamu)/(I0/(I+1e-12)+1)
        
        zeta[i] = np.nanmin([zeta_max, mub*P/abs(dudx[i]+1e-12)])
        eta[i] = np.nanmin([zeta_max, mu*P/abs(dudx[i]+1e-12)])
    
    sigma_11 = (zeta+eta/2)*dudx - P
    
    
    return sigma_11

# # 
@jit(nopython= True)
def flux(U, psi_local, sigma11_local, C):
# def flux(U, psi_local, C):
    A = U[0]
    M = U[1]
    u = M / (A) if A > 1e-12 else 0.0
    # print(M)
    if flux_Gray:
        
        if thickness:
            return np.array([
                M,
                M*u- psi_local * A 
            ])
        else:
            # print('here')
            return np.array([
                M,
                M*u- psi_local *A* np.exp(-C*(1-A))
                ])
    else:
        return np.array([
                M,
               M*u - sigma11_local/(rhoice)
                ])
# # def flux(U, psi_local, C)
        
            
# @jit(nopython= True)
# def flux(U, psi_local, C):
#     A = U[0]
#     M = U[1]
#     u = M / (A) if A > 1e-12 else 0.0
#     return np.array([
#         M,
#        M*u- psi_local * A * np.exp(-C*(1-A))# M*u - psi_local * A * np.exp(-C*(1-A))
#     ])        
    
    
@jit(nopython= True)
def flux_limiter(r):
    """
    Roe's Superbee
    
    """
    
    # phi = np.maximum(0, np.minimum(2*theta, 1), np.minimum(theta, 2))
    
    return np.maximum(0, np.minimum(2*r,1), np.minimum(r,2))

@jit(nopython = True)
def limited_slope(U):
    dU = np.zeros_like(U)
    for j in range(U.shape[0]):
        du_fwd = U[j,2:] - U[j,1:-1]
        du_bwd = U[j,1:-1] - U[j,:-2]
        r = np.zeros_like(du_fwd)
        mask = np.abs(du_fwd) > 1e-12
        r[mask] = du_bwd[mask] / du_fwd[mask]
        phi = flux_limiter(r)
        dU[j,1:-1] = 0.5 * phi * (U[j,2:] - U[j,:-2])
    return dU

@jit(nopython = True)
def psi_ViscousPlastic(x):
    
    # print(np.shape(x))
    E_pm = np.zeros_like(x)
    # print(np.shape(E_pm))
    dudx = first_derivative(x, dx, Nx)
    # print(np.shape(E_pm))
    for i in range(1, Nx-1):
        if abs(dudx[i]) < dudx_crit: 
            E_pm[i] = (E_pm_div + E_pm_conv)/2 + ((E_pm_div - E_pm_conv)/2)*dudx[i]/dudx_crit    
                
        else:
                
            if dudx[i] < -dudx_crit:
                E_pm[i] = E_pm_conv
                        
                    
            elif dudx[i] > dudx_crit:
                E_pm[i] = E_pm_div
    # print('here')
    return E_pm*(Pstar)/(rhoice)


@jit(nopython = True)
def psi_GC(u, A):
    dudx = first_derivative(u, dx, Nx)
    
    E_pm = np.zeros_like(u)
    for i in range(1, Nx-1):
        P = Pstar*A[i]*np.exp(-C*(1-A[i]))
        I = np.nanmin([1,dmean*abs(dudx[i])*np.sqrt(rhoice*A[i]/(P+1e-12))])
        mu = mu0 - (deltamu)/(I0/(I+1e-12)+1)
        # mu = mu0
        
        E_pm_div = (mu/2+mub)-1
        E_pm_conv = -(mu/2+mub)-1
        
        if abs(dudx[i]) < dudx_crit: 
            E_pm[i] = (E_pm_div + E_pm_conv)/2 + ((E_pm_div - E_pm_conv)/2)*dudx[i]/dudx_crit    
                
        else:
                
            if dudx[i] < -dudx_crit:
                E_pm[i] = E_pm_conv
                        
                    
            elif dudx[i] > dudx_crit:
                E_pm[i] = E_pm_div
    return E_pm*(Pstar)/(rhoice)
    
    

@jit(nopython= True)
def max_wave_speed(U, psi,sigma, C, eps=1e-6):
# def max_wave_speed(U, psi, C, eps=1e-6):
    n = U.shape[0]
    A = U[0]
    M = U[1]

    # Return zero speed for dry cell
    if A < 1e-8 and np.abs(M) < 1e-8:
        return 0.0

    # Compute baseline flux
    F0 = flux(U, psi, sigma, C)
    # F0 = flux(U, psi, C)
    J = np.zeros((n, n))

    for j in range(n):
        dU = np.zeros(n)
        dU[j] = eps
        F1 = flux(U + dU, psi, sigma, C)
        # F1 = flux(U + dU, psi, C)

        for i in range(n):
            dF = F1[i] - F0[i]
            J[i, j] = dF / eps if np.isfinite(dF) else 0.0

    try:
        eigvals = np.linalg.eigvals(J)
        maxval = 0.0
        for k in range(len(eigvals)):
            val = np.abs(eigvals[k])
            if np.isfinite(val) and val > maxval:
                maxval = val
        return maxval
    except:
        return 0.0

@jit(nopython = True)
def safe_u(A,M): 
    
    return M/(A) if A>1e-12 else 0.0


@jit(nopython= True)
def step(U, dx, dt, psi_profile, sigma_11, C):
# def step(U, dx, dt, psi_profile, C):
    Nx = U.shape[1]

    # Compute limited slopes for MUSCL
    dU = limited_slope(U)

    # Left/right states at interfaces
    UL = U - 0.5*dU
    UR = U + 0.5*dU

    U_new = U.copy()

    for i in range(1, Nx-1):
        psi_L  = psi_profile[i]
        psi_R  = psi_profile[i+1]
        psi_Lm = psi_profile[i-1]

        psi_half = 0.5*(psi_L + psi_R)
        psi_half_m = 0.5*(psi_Lm + psi_L)
        
        sigma_L  = sigma_11[i]
        sigma_R  = sigma_11[i+1]
        sigma_Lm = sigma_11[i-1]

        sigma_half = 0.5*(sigma_L + sigma_R)
        sigma_half_m = 0.5*(sigma_L + sigma_Lm)

        # States
        UL_i = UR[:, i-1]
        UR_i = UL[:, i]

        UL_ip = UR[:, i]
        UR_ip = UL[:, i+1]

        # Safe velocities

        A_L = UL_i[0]; u_L = safe_u(A_L, UL_i[1])
        A_R = UR_i[0]; u_R = safe_u(A_R, UR_i[1])
        A_Lp = UL_ip[0]; u_Lp = safe_u(A_Lp, UL_ip[1])
        A_Rp = UR_ip[0]; u_Rp = safe_u(A_Rp, UR_ip[1])

        FL = flux(UL_i, psi_half_m, sigma_half_m, C)
        FR = flux(UR_i, psi_half_m, sigma_half_m, C)
        FLp = flux(UL_ip, psi_half, sigma_half, C)
        FRp = flux(UR_ip, psi_half, sigma_half, C)
        
        # FL = flux(UL_i, psi_half_m, C)
        # FR = flux(UR_i, psi_half_m, C)
        # FLp = flux(UL_ip, psi_half, C)
        # FRp = flux(UR_ip, psi_half, C)

        # # Wave speeds
        # alpha_m = max(max_eigen(A_L, u_L, psi_half_m, C),
        #               max_eigen(A_R, u_R, psi_half_m, C))
        # alpha_p = max(max_eigen(A_Lp, u_Lp, psi_half, C),
        #               max_eigen(A_Rp, u_Rp, psi_half, C))
        
        alpha_m = max(
            max_wave_speed(UL_i, psi_half_m, sigma_half_m, C),
            max_wave_speed(UR_i, psi_half_m, sigma_half_m,C)
            )
        alpha_p = max(
                 max_wave_speed(UL_ip, psi_half, sigma_half, C),
                 max_wave_speed(UR_ip, psi_half, sigma_half, C)
                )
        
        # alpha_m = max(
        #     max_wave_speed(UL_i, psi_half_m, C),
        #     max_wave_speed(UR_i, psi_half_m, C)
        #     )
        # alpha_p = max(
        #          max_wave_speed(UL_ip, psi_half, C),
        #          max_wave_speed(UR_ip, psi_half,  C)
        #         )

        # Numerical fluxes (TVD Lax–Wendroff form)
        F_half_m = 0.5*(FL + FR) - 0.5*alpha_m*(UR_i - UL_i)
        F_half_p = 0.5*(FLp + FRp) - 0.5*alpha_p*(UR_ip - UL_ip)

        U_new[:, i] = U[:, i] - dt/dx * (F_half_p - F_half_m)

    # Dirichlet BCs
    U_new[:, 0] = 0.0
    U_new[:, -1] = 0.0
    
    # U_new[:, 1] = 0.0
    # U_new[:, -2] = 0.0

    if not thickness:
    #     # Cap+A <= 1
        A_cap = 1.0
        A_min = 0
        Avals = U_new[0]
        Mvals = U_new[1]
        over = Avals > A_cap
        # Mvals[over] *= A_cap / Avals[over]
        Avals[over] = A_cap
        
        under = Avals < A_min
        # Mvals[under] *= A_min / (Avals[under] + 1e-12)  # +eps to avoid 0/0
        Avals[under] = A_min

    return U_new


# @jit(nopython = True)
def solve_ice(U, dt, dx, Nx, Nt, Nt_saves, CFL = False):
    
    # time_save = np.linspace(0, Nt*dt, Nt_saves)#, dtype=np.float32)
    time_tot = np.linspace(0, Nt*dt, Nt)#, dtype=np.float32)
    idx_save = np.linspace(0, Nt - 1, 10, dtype = int)
    time_save = time_tot[idx_save]
    print(time_save)
    print(time_tot)
    U_saves = np.zeros((2, Nx, len(time_save)))
    id_saved =0
   
    for t, time in enumerate(time_tot):
        # print(t, time, time_save)
        if (t%1000 == 0):
            print('Solving: ',  t)
        

        A = U[0]
        M = U[1]
        u = np.zeros_like(A)
        mask = A > 1e-12
        u[mask] = M[mask]/(A[mask])
            
        # psi_profile = psi_ViscousPlastic(u)
        # psi_profile = psi_GC(u, A)
        # sigma_11 = sigma_VP(u, A)
        # sigma_11 = sigma_GC(u, A)
        
        if VP:
            psi_profile = psi_ViscousPlastic(u)
            sigma_11 = sigma_VP(u, A)
            # print('here')
        else:
            sigma_11 = sigma_GC(u, A)
            psi_profile = psi_GC(u, A)
            
        if CFL:
            lamb = 0.0
            for i in range(Nx):
                lamb = max(lamb, max_wave_speed(U[:, i], psi_profile[i], sigma_11[i],C))
            CFL = dt* lamb/dx
            
            CFL_adv = ((np.sqrt(1+e**(-2))-1)/2 * Pstar/rhoice)**(1/2) * dt/dx
            print(CFL, CFL_adv)

        # Unew = step(U, dx, dt, psi_profile,  C)
        Unew = step(U, dx, dt, psi_profile, sigma_11, C)
        U = Unew

        if np.any(np.isclose(time, time_save, atol=1e-5)):
        # print('here')
            print('Saving State: ', t+1)
            U_saves[:,:, id_saved] = U
            id_saved += 1

    return U_saves, time_save
   
  # M = A*u
# print()
U_solved, time_save = solve_ice(U, dt, dx, Nx, Nt, Nt_saves, CFL = CFL)


norm = colors.Normalize(vmin=0, vmax=time_save[-1]/(60))  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = cm.ScalarMappable(cmap=cmap, norm=norm)

fig, axs= plt.subplots(2, 1, figsize = (13, 6), sharex = True)
# axs = [ax1, ax2, ax3, ax4, ax5, ax6]
ax1, ax2 = axs.flatten()


for t, time in enumerate(time_save):
    
    A_ice = U_solved[0,:, t]
    u_ice = np.zeros_like(A_ice)
    mask = A_ice > 1e-12
    u_ice[mask] = U_solved[1,:, t][mask]/(A_ice[mask])
    # with np.printoptions(threshold=np.inf):
    #     print(u_ice)
    color = cmap(t/10)
    # print(np.where(u_VP[:, t] > 10))
    #-- VP --#
    ax1.plot(X/1e3, u_ice, color = color)
    ax2.plot(X/1e3, A_ice , color = color)
    
    
for ax in axs.flatten():
    
    ax.grid()
    ax.set_aspect('auto')
    ax.set_xlim(0.3, 0.7)
    
# for ax in axs.flatten()[:3]:
# ax4.set_ylim(-1, 1)
    


fig.align_ylabels()
if not thickness:
    ax2.set_ylabel(r'$A$ (m)')
ax1.set_ylabel(r'$u$ (m/s)')
fig.colorbar(sm, ax=axs, label = 'Time (min)')
fig.supxlabel(r'$x$ (km)', x = 0.43)

if VP:
    ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)
    fig.savefig('conservative_equations_VP_{}.png'.format(dx))

else:
    ax1.set_title(r'GC: $\mu_b = {}$'.format(mub), x = 0.43)
    fig.savefig('conservative_equations_GC_{}.png'.format(dx))
    
plt.show()
plt.close()
    
    
A = U_solved[0,:, -1]
u = np.zeros_like(A)
mask = A > 1e-12
u[mask] = U_solved[1,:, -1][mask]/A[mask]

# ------------------------
# Plot
# ------------------------
plt.figure(figsize=(8,4))
plt.subplot(211)
plt.plot(X/1e3, A0, 'k--', label='Initial A')
plt.plot(X/1e3, A, 'b', label='Final A')
# plt.xlim(300, 700)
plt.legend(); plt.ylabel("A")

plt.subplot(212)
plt.plot(X/1e3, u, 'r', label='Final u')
plt.plot(X/1e3, u0, 'k--', label='Initial u')
plt.legend(); plt.ylabel("u"); plt.xlabel("x")

plt.tight_layout()
plt.show()

plt.savefig('conservatives_eq_VP_last.png', dpi = 500)