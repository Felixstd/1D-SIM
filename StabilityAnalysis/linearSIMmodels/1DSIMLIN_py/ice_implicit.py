"""
This a code to solve the implicit equation for VP and GC.


We will implement them as GRay such that


rho h (u^n+1 - u^n)/dt - (psi^n+1P _ i+1 +psi^n+1P _ i-1)/dx = 0


This is a trigiagonal matrix  but we need to handle the diagonal terms which can be 0.

"""

#%%
import os
import sys
import pyamg
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import calc_derivatives as deriv
import init_cond as ini
import advection as adv
import scipy as sp

from scipy.linalg import solveh_banded
from global_constants import *
from numba import jit, prange
from tqdm import tqdm

from scipy.sparse.linalg import gmres, LinearOperator
from scipy.sparse.linalg import spilu

import warnings
warnings.filterwarnings("ignore") 
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

# %%
@jit(nopython = True)
def psi_ViscousPlastic(x):
    
    # print(np.shape(x))
    E_pm = np.zeros_like(x)
    # print(np.shape(E_pm))
    dudx = deriv.first_derivative(x, dx, Nx)
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
    dudx = deriv.first_derivative(u, dx, Nx)
    
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

def uv_linearization(unm1, unm2):
    
    u_k = (unm1 + unm2)/2
    
    return u_k

def solve_velocity(unm1, h, P, dt, dx, eps = 1e-8, Nl_max = 500, VP = True):
    
    Nx = np.shape(unm1)[0]
    un = np.copy(unm1)
    
    moverdt = rhoice*h/dt
    moverdt[moverdt < 1e-12] = 1e-12
    rhs = unm1*moverdt
    res = 0
    unm2 = 0
    
    for k in range(Nl_max):
        if k%100 == 0:
            print('Outer loop iterations k:' , k, res)
            
        if  k > 2:
            un = uv_linearization(unm1, unm2) 
    #Computing the psi for VP
        if VP:
            psi = psi_ViscousPlastic(un)
        else:
            psi = psi_GC(un)
        
        
        diag = moverdt
        upper_diag = -psi[1:]*P[1:]/dx
        lower_diag = psi[:-1]*P[:-1]/dx
    
        # un[1:Nx-1] = TDMAsolver(lower_diag, diag, upper_diag, rhs[1:Nx-1])    
        A = sp.sparse.diags([lower_diag, diag, upper_diag], offsets = [-1, 0, 1], shape = (Nx, Nx))
        # ml = pyamg.smoothed_aggregation_solver(A)
        # M = ml.aspreconditioner()
        # from scipy.sparse.linalg import gmres
        # un, info = gmres(A, b)  
        
        # ilu = spilu(A)  # A must be in CSC format
        # M = LinearOperator(A.shape, matvec=ilu.solve)
        # A = sp.sparse.diags([lower_diag, diag, upper_diag], offset = [-1, 0, 1], shape = (Nx, Nx))
        # M_inv_diag = 1.0 / (A.diagonal()+1e-12)
        # M  = LinearOperator(A.shape, matvec=lambda x: M_inv_diag * x)    
    
        # un, info = gmres(A, rhs, rtol = eps, atol = 1e-5, restart = Nx, callback_type = 'legacy', M = M)
        un, info = pyamg.krylov.fgmres(A, rhs, tol = eps,restart = Nx)#, M = M)# callback_type = 'legacy', M = M)
        # un, info = gmres(A, rhs, rtol=1e-8, restart=Nx)
        # print(info)
        
        if info != 0:
            print('GMRES did not converge')
        #     break
        
        res = np.linalg.norm(un - unm1, np.inf)
    
        if res < eps:
            # print(res, eps)
            print('Converged at {} iterations with a res {}'.format(k,res)) 
    
            return un

        unm1 = un
        unm2 = unm1

    return un
        
def ice_strength(A, h):
    
    return Pstar*h*np.exp(-C*(1-A))


def solve_ice(u0, h0, A0, dx, dt, Nt):
    
    time_save = np.linspace(0, Nt*dt, Nt_saves, dtype=np.float32)
    time_tot = np.linspace(0, Nt*dt, Nt, dtype=np.float32)
    
    print(time_save, time_tot)
    # print(time_tot)
    # print(time_tot)
    U_saves = np.zeros((3, Nx, len(time_save)))
    U_saves[0,:, 0] = u0
    U_saves[1,:, 0] = h0
    U_saves[2,:, 0] = A0
    id_saved =1
    
    unm1 = np.copy(u0)
    hnm1 = np.copy(h0)
    Anm1 = np.copy(A0)
    
    for c, t in enumerate(time_tot[1:]):
        # print(t)
        if c%100 == 0:
            print('Time step: {}'.format(t))
        P = ice_strength(hnm1, Anm1)
        
        un = solve_velocity(unm1, hnm1, P, dt, dx)
        
        # hn = adv.upwind_scheme(un, hnm1, Nx)
        An = adv.upwind_scheme(un, Anm1, Nx, concentration=True)
        
        hnm1 = An
        Anm1 = An
        unm1 = un
        
        # print(un, An, hn)
        # print(t, time_save)
        if np.any(np.isclose(t, time_save, atol=1e-4)):
        # print('here')
            # print(An)
            print('Saving State: ', t)
            U_saves[0,:, id_saved] = un
            U_saves[1,:, id_saved] = An
            U_saves[2,:, id_saved] = An
            id_saved += 1
    
    return U_saves, time_save

U_solved, time_save = solve_ice(u0, h0, A0, dx, dt, Nt)


norm = colors.Normalize(vmin=0, vmax=time_save[-1]/(60))  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = cm.ScalarMappable(cmap=cmap, norm=norm)

fig, axs= plt.subplots(2, 1, figsize = (13, 6), sharex = True)
# axs = [ax1, ax2, ax3, ax4, ax5, ax6]
ax1, ax2 = axs.flatten()
with np.printoptions(threshold=np.inf):
    print(U_solved[2,:, -1])

for t, time in enumerate(time_save):
    
    u_ice = U_solved[0,:, t]
    h_ice = U_solved[1,:, t]
    A_ice = U_solved[2,:, t]
    
    color = cmap(t/10)
    # print(np.where(u_VP[:, t] > 10))
    #-- VP --#
    ax1.plot(X/1e3, u_ice, color = color)
    ax2.plot(X/1e3, A_ice , color = color)
    
    
for ax in axs.flatten():
    
    ax.grid()
    ax.set_aspect('auto')
    # ax.set_xlim(1500, 2500)
    
# for ax in axs.flatten()[:3]:
# ax4.set_ylim(-1, 1)
    
ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)

fig.align_ylabels()
# if not thickness:
ax2.set_ylabel(r'$A$ (m)')
ax1.set_ylabel(r'$u$ (m/s)')
fig.colorbar(sm, ax=axs, label = 'Time (min)')
fig.supxlabel(r'$x$ (km)', x = 0.43)

if VP:
    fig.savefig('implicit_equations_VP_{}.png'.format(dx))
else:
    fig.savefig('implicit_equations_GC_{}.png'.format(dx))
plt.show()
plt.close()