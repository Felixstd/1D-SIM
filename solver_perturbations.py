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
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.linalg import solveh_banded
from numba import jit, prange
from tqdm import tqdm

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')
#%%
#----- Global Variables -------#
print('Setting Model Variables and Compiling Functions')

Nx = 800+1
dx = 1e4
L = Nx*dx
X = np.arange(0, L+dx, dx)
x = np.arange(0, Nx+1, 1)
dt = 1
numax = 1e12

Nt = 1e4 # Nt = 2

h0 = np.ones(Nx)
A0 = np.ones(Nx)


#---- Rheological Parameters ----#
e = 2.0
C = 20
rhoice = 900
Pstar = 27.5e3
Pstarstar = Pstar*np.exp(-C*(1-A0))

mu0 = 0.2
muinf = 0.8
deltamu = muinf-mu0
dmean = 1e3
I0 = 1e-3
dudx_crit = 1e-8

DtoverDx = dt/dx

time = np.arange(0, int(Nt), 1)

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
    
    derivative[1:N] = (x[2:] - x[:-2])/(2*deltax)
    
    derivative[0] = (-3*x[0] + 4*x[1] - x[2])/(2*deltax)
    derivative[-1] = (3*x[-1] - 4*x[-2] + x[-3])/(2*deltax)
    
    return derivative

@jit(nopython=True)
def second_derivative(x, deltax, N):
    """
    Function to evaluate the second derivative of an array x. 
    
    This is the formula for a second order accuracy. 

    Args:
        x (_type_): _description_
        deltax (_type_): _description_
        N (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    sec_derivative = np.zeros_like(x)
    
    sec_derivative[1:N] = (x[2:] - 2*x[1:N] + x[:-2])/(deltax**2)
    sec_derivative[0] = (-x[3] + 4*x[2] - 5*x[1] + 2*x[0])/deltax**2
    sec_derivative[-1] = (-x[-4] + 4*x[-3] - 5*x[-2] + 2*x[-1])/deltax**2
    
    
    return sec_derivative

def init_condition(type_u, type_h):
    """
    Function which treats the initial conditions of the 
    simulations. 

    Args:
        type_u (_type_): _description_
        type_h (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    hinit = np.zeros(Nx+1)
    uinit = np.zeros(Nx+1)
    
    x = np.arange(0, Nx+1,1)
    
    if type_h == 'Step':
        hinit = np.exp(-((x-Nx/2)**2/10000)**5)
        
    elif type_h == 'Gray':
        hinit = (10000 - (x-Nx/2)**2)**(1/20)
        hinit[np.isnan(hinit)] = 0
    
    if type_u == 'Step':
        uinit = np.exp(-((x-Nx/2)**2/10000)**5)
    
    elif type_u == 'Gray':
        uinit[0:int(Nx/4)] = 0
        uinit[int(3*Nx/4)] = 0
        # uinit[int(Nx/4):int(3*Nx/4)] = x[int(Nx/4):int(3*Nx/4)]/1e4
        
        uinit[np.where(hinit>0.2)] = x[np.where(hinit>0.2)]/1e4
        # uinit[int(2*Nx/5):int(3*Nx/5)] = x[int(2*Nx/5):int(3*Nx/5)]/1e4
    elif type_u == 'zero':
        uinit = np.zeros(Nx+1)
        
    if type_h == 'Step':
        hinit = np.exp(-((x-Nx/2)**2/10000)**5)
        
    elif type_h == 'Gray':
        hinit = (10000 - (x-Nx/2)**2)**(1/20)
        hinit[np.isnan(hinit)] = 0

        
    uinit[0] = 0
    uinit[-1] = 0
    
    hinit[0] = 0
    hinit[-1] = 0
    

    return uinit, hinit

@jit(nopython=True)
def calc_flux(ui, uip1, Tim1, Ti, Tip1):
    
    if (ui > 0) :                                                                        
        flux1 = ui*Tim1
    else:
        flux1 = ui*Ti


    if (uip1 < 0):                                                                             
        flux2 = uip1*Ti
    else:
        flux2 = uip1*Tip1

    flux = flux2 - flux1
    
    return flux

@jit(nopython=True)
def upwind_scheme(u, h, N) : 
    
    hout = np.zeros_like(h)
    for i in prange(1, N-1):
        
        flux = calc_flux(u[i], u[i+1], h0[i-1], h0[i], h0[i+1])
        
        hout[i] = h[i] - DtoverDx*flux
        hout[i] = max(hout[i], 0)
    

    hout[0] = 0
    hout[-1] = 0
    
    return hout

@jit(nopython = True)
def solve_mom_VP(u, h, h0, N):
    uout = np.zeros_like(u)
    
    dudx = first_derivative(u, dx, N)
    dhdx = first_derivative(h, dx, N)
    d2udx2 = second_derivative(u, dx, N)
    E_pm_div = (np.sqrt(1+e**(-2))-1)/2
    E_pm_conv = (-np.sqrt(1+e**(-2))-1)/2
    
    for i in prange(1, N-1):
        
        if abs(dudx[i]) < dudx_crit: 
            E_pm = (E_pm_div + E_pm_conv)/2 + ((E_pm_div + E_pm_conv)/2)*dudx[i]/dudx_crit    
        
        else:
        
            if dudx[i] < -dudx_crit:
                E_pm = E_pm_conv
                
            # elif abs(dudx[i]-1e-12) < 1e-12 :
            #     E_pm = 0
            elif dudx[i] > dudx_crit:
                E_pm = E_pm_div

        uout[i] = u[i] + dt/(rhoice*h0[i])*(E_pm*Pstar*dhdx[i] )
        
    uout[0] = 0
    uout[-1] = 0
    
    return uout
        
    
    
# @jit(nopython=True)
def solve_mom_GC(u, h, h0, N, mub, solver = 'explicit'):
    
    uout = np.zeros_like(u)
    
    dudx = first_derivative(u, dx, N)
    d2udx2 = second_derivative(u, dx, N)
    dhdx = first_derivative(h, dx, N)
    
    gamma_div = mu0/2 + mub - 1
    gamma_conv = -mu0/2 - mub - 1
    
    if solver == 'explicit':
        for i in prange(1, N-1):
            
            if abs(dudx[i]) < dudx_crit: 
                gamma = (gamma_div + gamma_conv)/2 + ((gamma_div + gamma_conv)/2)*dudx[i]/dudx_crit    
                Gamma = 0
            else:
        
                if dudx[i] < -dudx_crit:
                    gamma = gamma_conv
                
            # elif abs(dudx[i]-1e-12) < 1e-12 :
            #     E_pm = 0
                elif dudx[i] > dudx_crit:
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
                gamma = (gamma_div + gamma_conv)/2 + ((gamma_div + gamma_conv)/2)*dudx[i]/dudx_crit    
                Gamma = 0
            else:
        
                if dudx[i] < -dudx_crit:
                    gamma = gamma_conv
                
            # elif abs(dudx[i]-1e-12) < 1e-12 :
            #     E_pm = 0
                elif dudx[i] > dudx_crit:
                    gamma = gamma_div

                Gamma = h0[i] * np.sqrt(Pstar * rhoice) * (dmean * deltamu / (2 * I0))
            
            A[i] = dt * Gamma / (dx**2 * rhoice * h0[i])
            B[i] = dt * gamma * Pstar / (rhoice * h0[i])
            rhs[i] = u[i] + B[i] * dhdx[i]

        # Tridiagonal coefficients
        lower = -A[1:N-2]
        diag = 1 + 2*A[1:N-1]

        # Solve tridiagonal system using Thomas algorithm
    
        ab = np.zeros((2, N-2))
        ab[0, :] = diag                    # main diagonal
        ab[1, :-1] = lower                 # lower diagonal (subdiagonal)
    # Solve system
        uout[1:N-1] = solveh_banded(ab, rhs[1:N-1], lower=True)
    

    # Boundary conditions
    uout[0] = 0
    uout[-1] = 0

    return uout
    


def kinetic_energy():
    
    return
    

# %%

print('Initializing the velocity and Ice Thickness')
#- Granular Compressible Model -#
Bmub = 10
Lmub = 1/2


u_GC_Bmub = np.zeros((Nx+1, int(Nt)))
u_GC_Lmub = np.zeros((Nx+1, int(Nt)))
h_GC_Bmub = np.zeros((Nx+1, int(Nt)))
h_GC_Lmub = np.zeros((Nx+1, int(Nt)))

k_GC_Lmub = np.zeros((Nx+1, int(Nt)))
k_GC_Bmub = np.zeros((Nx+1, int(Nt)))
k_VP = np.zeros((Nx+1, int(Nt)))


#- Viscous-PlasticModel -#
u_VP = np.zeros((Nx+1, int(Nt)))
h_VP = np.zeros((Nx+1, int(Nt)))

#initial conditions
u_GC_Bmub[:, 0], h_GC_Bmub[:, 0] = init_condition('Gray', 'Step')
u_GC_Lmub[:, 0], h_GC_Lmub[:, 0] = init_condition('Gray', 'Step')
u_VP[:, 0], h_VP[:, 0] = init_condition('Gray', 'Step')

#main loop
for t in time[:-1]:
    
    if (t%1000) == 0:
        print('Solving for u and h, Time step: ', int(t/Nt*100), '%')
        
    h_GC_Bmub[:, t+1] = upwind_scheme(u_GC_Bmub[:, t], h_GC_Bmub[:, t], Nx)
    u_GC_Bmub[:, t+1] = solve_mom_GC(u_GC_Bmub[:, t], h_GC_Bmub[:, t], h0, Nx, Bmub, solver = 'implicit')

    h_GC_Lmub[:, t+1] = upwind_scheme(u_GC_Lmub[:, t], h_GC_Lmub[:, t], Nx)
    u_GC_Lmub[:, t+1] = solve_mom_GC(u_GC_Lmub[:, t], h_GC_Lmub[:, t], h0, Nx, Lmub, solver = 'implicit')
    
    h_VP[:, t+1] = upwind_scheme(u_VP[:, t], h_VP[:, t], Nx)
    u_VP[:, t+1] = solve_mom_VP(u_VP[:, t], h_VP[:, t], h0, Nx)
    
    k_GC_Lmub[:, t] = 1/2*rhoice*u_GC_Lmub[:, t]**2*h_GC_Lmub[:, t]
    k_GC_Bmub[:, t] = 1/2*rhoice*u_GC_Bmub[:, t]**2*h_GC_Bmub[:, t]
    k_VP[:, t] = 1/2*rhoice*u_VP[:, t]**2*h_VP[:, t]


# %%

#%%

#------ Figures --------#

#plot 10 time steps
time_plot = np.linspace(0, (Nt-1)*dt, 10, dtype=int)
time_plot = np.linspace(0, Nt-1, 10, dtype=int)
#color bar
norm = colors.Normalize(vmin=0, vmax=time_plot[-1]/(60*60))  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = cm.ScalarMappable(cmap=cmap, norm=norm)


fig, axs= plt.subplots(2, 3, figsize = (13, 6), sharex = True)
# axs = [ax1, ax2, ax3, ax4, ax5, ax6]
ax1, ax2, ax3, ax4, ax5, ax6 = axs.flatten()

for k, t in enumerate(time_plot):
    
    color = cmap(k/10)

    #-- VP --#
    ax1.plot(X/1e3, h_VP[:, t], color = color)
    ax4.plot(X/1e3, u_VP[:, t], color = color)
    
    #-- GC Stable --#
    ax2.plot(X/1e3, h_GC_Lmub[:, t], color = color)
    ax5.plot(X/1e3, u_GC_Lmub[:, t], color = color)
    
    #-- GC UnStable --#
    ax3.plot(X/1e3, h_GC_Bmub[:, t], color = color)
    ax6.plot(X/1e3, u_GC_Bmub[:, t], color = color)

for ax in axs.flatten():
    
    ax.grid()
    ax.set_aspect('auto')
    
ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)
ax2.set_title(r'GC: $\mu_b = {}$'.format(Lmub), x = 0.43)
ax3.set_title(r'GC: $\mu_b = {}$'.format(Bmub), x = 0.43)

fig.align_ylabels()
ax1.set_ylabel(r'$h$ (m)')
ax4.set_ylabel(r'$u$ (m/s)')
fig.colorbar(sm, ax=axs, label = 'Time (hr)')
fig.supxlabel(r'$x$ (km)', x = 0.43)

fig.savefig('linearized_equations_VP_GC_{}.png')
plt.show()
plt.close()

fig2, axs = plt.subplots(1, 3, figsize = (12, 4))
(ax1, ax2, ax3) = axs.flatten()
for k, t in enumerate(time_plot[:-1]):
    
    color = cmap(k/10)
    ax1.plot(X/1e3, k_VP[:, t]/1e3, color = color)
    ax2.plot(X/1e3, k_GC_Lmub[:, t]/1e3, color = color)
    ax3.plot(X/1e3, k_GC_Bmub[:, t]/1e3, color = color)
    
ax1.set_title(r'VP: $e = {}$'.format(e), x = 0.43)
ax2.set_title(r'GC: $\mu_b = {}$'.format(Lmub), x = 0.43)
ax3.set_title(r'GC: $\mu_b = {}$'.format(Bmub), x = 0.43)

for ax in axs.flatten():
    ax.grid()
    ax.set_aspect('auto')


# plt.grid()
fig2.supxlabel(r'$x$ (km)')
ax1.set_ylabel(r'$K$ (kJ)')
fig2.colorbar(sm, ax=ax3, label = 'Time (hr)')
plt.savefig('kinetic_energy')



# %%

# %%
