"""
This is a python implementation of the 1D linearized equation for 
the Granular-Compressible Rheology (GC) and the Viscous-Plastic Rheology (VP). 

It solves both the linearized equation:

rhoi ho u_t = gamma_pm P* h_x + Gamma u_xx 
rhoi ho u_t = E_pm P*h_x 
h_t + (uh0)_x = 0

where h and u are primed variables.

This isolates the plastic and validates the linear stability analysis made 
for both models.


This uses an explicit approach. 


Félix St-Denis, 24 July 2025
"""

# %%

print('Importing functions and constants')

import os
import sys
sys.path.insert(1, '/aos/home/fstdenis/1D-SIM/StabilityAnalysis/1DSIMLIN_py/')


import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import calc_derivatives as deriv
import init_cond as ini
import advection as adv

from scipy.linalg import solveh_banded
from global_constants import *
from numba import jit, prange
from tqdm import tqdm

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')
#%%
#----- Global Variables -------#
print('Setting Model Variables and Compiling Functions')

Nx = 400+1
dx = 1e4
L = Nx*dx
X = np.arange(0, L+dx, dx)
x = np.arange(0, Nx+1, 1)
dt = 1e-2
numax = 1e12

Nt = 1000000# Nt = 2

h0 = np.ones(Nx+1)
# A0 = np.ones(Nx)
A0 = 1


#---- Rheological Parameters ----#
e = 2
C = 20
rhoice = 900
# Pstar = 27.5e3
Pstar_VP = 27.5e3*2
Pstar_GC = 27.5e3
Pstarstar_VP = Pstar_VP*np.exp(-C*(1-A0))
Pstarstar_GC = Pstar_GC*np.exp(-C*(1-A0))

mu0 = 0.2
muinf = 0.8
deltamu = muinf-mu0
dmean = 1e3
I0 = 1e-3
dudx_crit = 1e-10

DtoverDx = dt/dx

time = np.arange(0, int(Nt), 1)

# @jit(nopython = True)
def solve_mom_VP(u, h, h0, N, method = 'upwind', concentration = False):
    uout = np.zeros_like(u)
    
    dudx = deriv.first_derivative(u, dx, N)
    dhdx = deriv.first_derivative(h, dx, N)
    d2udx2 = deriv.second_derivative(u, dx, N)
    E_pm_div = (np.sqrt(1+e**(-2))-1)/2
    E_pm_conv = (-np.sqrt(1+e**(-2))-1)/2
    
    E_pm = np.zeros_like(u)
    
    for i in prange(1, N-1):
        if abs(dudx[i]) < dudx_crit: 
            E_pm[i] = (E_pm_div + E_pm_conv)/2 + ((E_pm_div - E_pm_conv)/2)*dudx[i]/dudx_crit    
                
        else:
                
            if dudx[i] < -dudx_crit:
                E_pm[i] = E_pm_conv
                        
                    # elif abs(dudx[i]-1e-12) < 1e-12 :
                    #     E_pm = 0
            elif dudx[i] > dudx_crit:
                E_pm[i] = E_pm_div
            
    if method == 'explicit':
    
        for i in prange(1, N-1):


            if concentration:
                uout[i] = u[i] + dt/(rhoice*h0[i])*(E_pm[i]*Pstar_VP*dhdx[i]*C*np.exp(-C*(1-h[i])) )
            else:
                uout[i] = u[i] + dt/(rhoice*h0[i])*(E_pm[i]*Pstar_VP*dhdx[i] )
    
    if method == 'upwind':

        for i in prange(1, N-1):   
            if concentration:
                flux = adv.calc_flux(E_pm[i+1]*Pstar/(rhoice*h0[i+1]), E_pm[i]*Pstar/(rhoice*h0[i]), h[i-1]*C*np.exp(-C*(1-h[i-1])), h[i]*C*np.exp(-C*(1-h[i])), h[i+1]*C*np.exp(-C*(1-h[i+1])))
            else:
                flux = adv.calc_flux(E_pm[i+1]*Pstar/(rhoice*h0[i+1]), E_pm[i]*Pstar/(rhoice*h0[i]), h[i-1], h[i], h[i+1])
            
            uout[i] = u[i] + DtoverDx*flux
        
        # uout[i] = (u[i+1] + u[i-1])/2 + (DtoverDx/2) * (E_pm*Pstar/(rhoice*h0[i+1])*h[i+1] - E_pm*Pstar/(rhoice*h0[i-1])*h[i-1])
    

    
    if method == 'LAX':

        A = E_pm*Pstar/(rhoice)
        
        Qstar = np.zeros_like(u)
        F = np.zeros_like(u)
        sigma = np.zeros_like(u)
        # print(np.shape(A), np.shape(F), np.shape(qplus), np.shape(sigma))
        sigma[1:-1]  = (A[2:]*h[2:] - A[:-2]*h[:-2])/(2.0*dx)
        qplus  = A[1:-1]*h[1:-1] - sigma[1:-1] * dx/2.0  # q^+_{i-1/2}
        qminus = A[:-2]*h[:-2] + sigma[:-2]  * dx/2.0  # q^-_{i-1/2}
        F[1:-1] = 0.5*(qplus+qminus - dx/dt*(qplus-qminus) )# F_{i-1/2}
        
        Qstar[1:-1] = u[1:-1] + dt/dx*(F[2:]-F[:-2])
        Qstar[-1] = 0
        Qstar[0] = 0
        
        sigma[1:-1]  = (Qstar[2:] - Qstar[:-2])/(2.0*dx)
        qplus  = Qstar[1:-1] - sigma[1:-1] * dx/2.0  # q^+_{i-1/2}
        qminus = Qstar[:-2] + sigma[:-2]  * dx/2.0  # q^-_{i-1/2}
        F[1:-1] = 0.5*(qplus + qminus - dx/dt*(qplus-qminus))# F_{i-1/2}
        
        uout[1:-1] = 0.5*u[1:-1] + 0.5*(Qstar[1:-1] - dt/dx*(F[2:]-F[:-2]))
        
        
        
    uout[0] = 0
    uout[-1] = 0
    
    return uout

@jit(nopython = True)
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    
    Taken from: https://gist.github.com/ofan666/1875903 
    
    Inputs:
        a (array, N-2): lower diagonal
        b (array, N-1): diagonal
        c (array, N-2): upper diagonal
        d (array, N-1): Right hand side
    
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc  = np.copy(a), np.copy(b), np.copy(c), np.copy(d)
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc
        
@jit(nopython=True)
def solve_mom_GC(u, h, h0, N, mub, solver = 'explicit'):
    
    uout = np.zeros_like(u)
    
    dudx = deriv.first_derivative(u, dx, N)
    d2udx2 = deriv.second_derivative(u, dx, N)
    dhdx = deriv.first_derivative(h, dx, N)
    
    gamma_div = mu0/2 + mub - 1
    gamma_conv = -mu0/2 - mub - 1
    
    if solver == 'explicit':
        for i in prange(1, N):
            
            if abs(dudx[i]) < dudx_crit: 
                gamma = (gamma_div + gamma_conv)/2 + \
                    ((gamma_div - gamma_conv)/2)*dudx[i]/dudx_crit    
                # Gamma = 0
            else:
        
                if dudx[i] <= -dudx_crit:
                    gamma = gamma_conv
                
            # elif abs(dudx[i]-1e-12) < 1e-12 :
            #     E_pm = 0
                elif dudx[i] >= dudx_crit:
                    gamma = gamma_div
            
            # if dudx[i] < -1e-12:
            #     gamma = -mu0/2 - mub - 1
            # elif abs(dudx[i] - 1e-12) < 1e-12:
            #     gamma = 0
            # elif dudx[i] > 1e-12:
            #     gamma = mu0/2 + mub - 1

                
            Gamma = h0[i]*np.sqrt(Pstar*rhoice)*(dmean*deltamu/(2*I0))
            
            uout[i] = u[i] + dt/(rhoice*h0[i])*(gamma*Pstar*dhdx[i] + Gamma*d2udx2[i])

    elif solver == 'implicit':

        A = np.zeros(N)
        B = np.zeros(N)
        rhs = np.zeros(N)
    
        # Compute gamma and coefficients A, B
        for i in range(1, N-1):
            if abs(dudx[i]) < dudx_crit: 
                gamma = (gamma_div + gamma_conv)/2 + ((gamma_div - gamma_conv)/2)*dudx[i]/dudx_crit    
                # Gamma = 0
            else:
        
                if dudx[i] <= -dudx_crit:
                    gamma = gamma_conv
                
            # elif abs(dudx[i]-1e-12) < 1e-12 :
            #     E_pm = 0
                elif dudx[i] >= dudx_crit:
                    gamma = gamma_div

            Gamma = h0[i] * np.sqrt(Pstar_GC * rhoice) * (dmean * deltamu / (2 * I0))
            Gamma =0
            A[i] = dt * Gamma / (dx**2 * rhoice * h0[i])
            B[i] = dt * gamma * Pstar_GC / (rhoice * h0[i])
            rhs[i] = u[i] + B[i] * dhdx[i]

        # Tridiagonal coefficients
        lower = -A[1:N-2]
        diag = 1 + 2*A[1:N-1]
        
        # Solve tridiagonal system using Thomas algorithm
        uout[1:N-1] = TDMAsolver(lower, diag, lower, rhs[1:N-1])


    # Boundary conditions
    uout[0] = 0
    uout[-1] = 0

    return uout
    


def kinetic_energy():
    
    return

# %%

# %%
print('Initializing the velocity and Ice Thickness')
#- Granular Compressible Model -#
Bmub = 10
Lmub = 1/2

Nt_saves = 10
save_steps = np.linspace(0, Nt-1, Nt_saves, dtype=int)

u_GC_Bmub_saved = np.zeros((Nx+1, int(Nt_saves)))
u_GC_Lmub_saved = np.zeros((Nx+1, int(Nt_saves)))
h_GC_Bmub_saved = np.zeros((Nx+1, int(Nt_saves)))
h_GC_Lmub_saved = np.zeros((Nx+1, int(Nt_saves)))


#- Viscous-PlasticModel -#
u_VP_saved = np.zeros((Nx+1, int(Nt_saves)))
h_VP_saved = np.zeros((Nx+1, int(Nt_saves)))

k_GC_Lmub_saved = np.zeros((Nx+1, int(Nt_saves)))
k_GC_Bmub_saved = np.zeros((Nx+1, int(Nt_saves)))
k_VP_saved = np.zeros((Nx+1, int(Nt_saves)))


#initial conditions
u_GC_Bmub, h_GC_Bmub = ini.init_condition('GrayDiv', 'Step', Nx)
u_GC_Lmub, h_GC_Lmub = ini.init_condition('GrayDiv', 'Step', Nx)
u_VP, h_VP           = ini.init_condition('GrayDiv', 'Step', Nx)

u_GC_Bmub_saved[:, 0] = u_GC_Bmub
u_GC_Lmub_saved[:, 0] = u_GC_Lmub
h_GC_Bmub_saved[:, 0] = h_GC_Bmub
h_GC_Lmub_saved[:, 0] = h_GC_Lmub

u_VP_saved[:, 0] = u_VP
h_VP_saved[:, 0] = h_VP


id_saved = 1



#main loop
for t in time[:-1]:
     
    if (t%1000) == 0:
        print('Solving for u and h, Time step: ', int(t/Nt*100), '%')
        
    h_GC_Bmub_new = adv.upwind_scheme(u_GC_Bmub, 
                                      h_GC_Bmub, 
                                      Nx)
    u_GC_Bmub_new = solve_mom_GC(u_GC_Bmub, 
                                     h_GC_Bmub, 
                                     h0, 
                                     Nx, 
                                     Bmub, 
                                     solver = 'implicit')

    
    h_GC_Lmub_new = adv.upwind_scheme(u_GC_Lmub, 
                                      h_GC_Lmub, 
                                      Nx)
    u_GC_Lmub_new = solve_mom_GC(u_GC_Lmub, 
                                h_GC_Lmub, 
                                h0, 
                                Nx, 
                                Lmub, 
                                solver = 'implicit')
    
    h_VP_new = adv.upwind_scheme(u_VP, 
                                 h_VP, 
                                 Nx, 
                                 concentration = False)
    u_VP_new = solve_mom_VP(u_VP, 
                            h_VP, 
                            h0,
                            Nx, 
                            concentration=False, method = 'explicit')
    
    if t+1 in save_steps:
        # print('here')
        print('Saving State: ', t+1)
        # print(t+1, save_steps, id_saved)
        u_GC_Bmub_saved[:, id_saved] = u_GC_Bmub_new
        u_GC_Lmub_saved[:, id_saved] = u_GC_Lmub_new
        h_GC_Bmub_saved[:, id_saved] = h_GC_Bmub_new
        h_GC_Lmub_saved[:, id_saved] = h_GC_Lmub_new

        u_VP_saved[:, id_saved] = u_VP_new
        h_VP_saved[:, id_saved] = h_VP_new
        id_saved += 1
    
    u_GC_Bmub, h_GC_Bmub = u_GC_Bmub_new, h_GC_Bmub_new
    u_GC_Lmub, h_GC_Lmub = u_GC_Lmub_new, h_GC_Lmub_new
    u_VP, h_VP = u_VP_new, h_VP_new
    
    
    # print(u_VP[:, 0])
    # k_GC_Lmub[:, t] = 1/2*rhoice*u_GC_Lmub[:, t]**2*h_GC_Lmub[:, t]
    # k_GC_Bmub[:, t] = 1/2*rhoice*u_GC_Bmub[:, t]**2*h_GC_Bmub[:, t]
    # k_VP[:, t] = 1/2*rhoice*u_VP[:, t]**2*h_VP[:, t]



# %%

#------ Figures --------#
# print(u_VP_saved)
#plot 10 time steps
# time_plot = np.linspace(0, (Nt-1)*dt, 10, dtype=int)

time_plot = save_steps

# time_plot = [0]

# time_plot = np.linspace(0, (Nt-1)*dt, 10, dtype=int)
# time_plot = np.linspace(0, Nt-1, 10, dtype=int)
#color bar
norm = colors.Normalize(vmin=0, vmax=time_plot[-1]/(60*60))  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = cm.ScalarMappable(cmap=cmap, norm=norm)


fig, axs= plt.subplots(2, 3, figsize = (13, 6), sharex = True)
# axs = [ax1, ax2, ax3, ax4, ax5, ax6]
ax1, ax2, ax3, ax4, ax5, ax6 = axs.flatten()

for k, t in enumerate(time_plot):
    # u_VP[abs(u_VP_saved[:, t])<1e-12 , t] = 0
    color = cmap(k/10)
    # print(np.where(u_VP[:, t] > 10))
    #-- VP --#
    ax1.plot(X/1e3, h_VP_saved[:, k], color = color)
    ax4.plot(X/1e3, u_VP_saved[:, k] , color = color)
    # print(u_VP[:, t])

    #-- GC Stable --#
    ax2.plot(X/1e3, h_GC_Lmub_saved[:, k], color = color)
    ax5.plot(X/1e3, u_GC_Lmub_saved[:, k], color = color)
    
    #-- GC UnStable --#
    ax3.plot(X/1e3, h_GC_Bmub_saved[:, k], color = color)
    ax6.plot(X/1e3, u_GC_Bmub_saved[:, k], color = color)

for ax in axs.flatten():
    
    ax.grid()
    ax.set_aspect('auto')
    # ax.set_xlim(1500, 2500)
    
# for ax in axs.flatten()[:3]:
# ax4.set_ylim(-1, 1)
    
ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)
ax2.set_title(r'GC: $\mu_b = {}$'.format(Lmub), x = 0.43)
ax3.set_title(r'GC: $\mu_b = {}$'.format(Bmub), x = 0.43)

fig.align_ylabels()
ax1.set_ylabel(r'$h$ (m)')
ax4.set_ylabel(r'$u$ (m/s)')
fig.colorbar(sm, ax=axs, label = 'Time (hr)')
fig.supxlabel(r'$x$ (km)', x = 0.43)

fig.savefig('linearized_equations_VP_GC_{}.png'.format(dx))
plt.show()
plt.close()

# fig2, axs = plt.subplots(1, 3, figsize = (12, 4))
# (ax1, ax2, ax3) = axs.flatten()
# for k, t in enumerate(time_plot[:-1]):
    
#     color = cmap(k/10)
#     ax1.plot(X/1e3, k_VP[:, t]/1e3, color = color)
#     ax2.plot(X/1e3, k_GC_Lmub[:, t]/1e3, color = color)
#     ax3.plot(X/1e3, k_GC_Bmub[:, t]/1e3, color = color)
    
# ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)
# ax2.set_title(r'GC: $\mu_b = {}$'.format(Lmub), x = 0.43)
# ax3.set_title(r'GC: $\mu_b = {}$'.format(Bmub), x = 0.43)

# for ax in axs.flatten():
#     ax.grid()
#     ax.set_aspect('auto')
    


# # plt.grid()
# fig2.supxlabel(r'$x$ (km)')
# ax1.set_ylabel(r'$K$ (kJ)')
# fig2.colorbar(sm, ax=ax3, label = 'Time (hr)')
# plt.savefig('kinetic_energy')



# %%

# %%

# %%
