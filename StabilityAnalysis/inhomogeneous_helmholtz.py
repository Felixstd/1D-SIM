"""
Python script to 
solve the inhomogeneous Helmholtz equation. 
"""
import numpy as np
import matplotlib.pyplot as plt  
import seaiceparameters as param
from numba import njit, prange
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


nx = 500
nt = 100000
deltax = 4000
deltax2 = deltax**2
deltat = 1
dtoverdx = deltat/deltax

@njit(inline='always')
def calc_flux(ui, uip1, Tim1, Ti, Tip1):
    flux1 = ui*Tim1   if ui   > 0.0 else ui*Ti
    flux2 = uip1*Ti   if uip1 > 0.0 else uip1*Tip1
    return flux2 - flux1
    
    
    

@njit(parallel=True)
def advection(A, h, u, nx, dtoverdx):
    hout = np.zeros(nx+1)
    Aout = np.zeros(nx+1)
    for i in prange(1, nx-1):
        flux = calc_flux(u[i], u[i+1], h[i-1], h[i], h[i+1])
        hi = h[i] - dtoverdx*flux
        hout[i] = hi if hi > 0.0 else 0.0

        flux = calc_flux(u[i], u[i+1], A[i-1], A[i], A[i+1])
        ai = A[i] - dtoverdx*flux
        ai = min(max(ai, 0.0), 1.0)
        Aout[i] = ai
    return hout, Aout


@njit
def solve_helm(A, h, l, nx, deltax2, Pstar, C_s):
    lower  = np.zeros(nx)
    main_d = np.zeros(nx)
    upper  = np.zeros(nx)
    rhs    = np.zeros(nx)
    P_nl   = np.zeros(nx)
    P_loc  = np.zeros(nx)

    for i in range(nx):
        P_loc[i] = Pstar*h[i]*np.exp(-C_s*(1.0 - A[i]))
        l_scale2 = (l*A[i])**2
        lower[i]  = l_scale2/deltax2
        main_d[i] = -(2.0*l_scale2/deltax2 + 1.0)
        upper[i]  = l_scale2/deltax2
        rhs[i]    = -P_loc[i]

    main_d[0] = 1.0
    upper[0]  = 0.0
    rhs[0]    = P_loc[1]

    main_d[nx-1] = 1.0
    lower[nx-1]  = 0.0
    rhs[nx-1]    = P_loc[nx-1]

    for i in range(1, nx):
        w = lower[i] / main_d[i-1]
        main_d[i] -= w*upper[i-1]
        rhs[i]    -= w*rhs[i-1]

    P_nl[nx-1] = rhs[nx-1] / main_d[nx-1]
    for i in range(nx-2, -1, -1):
        P_nl[i] = (rhs[i] - upper[i]*P_nl[i+1]) / main_d[i]

    return P_nl, P_loc

@njit
def run_simulation(A, h, u, l_c, nx, nt, deltax2, dtoverdx, Pstar, C_s):
    n_lc = l_c.shape[0]
    P_nl_lc = np.zeros((n_lc, nt, nx))
    P_loc_t = np.zeros((nt, nx))

    for t in range(nt-1):
        h[t+1], A[t+1] = advection(A[t], h[t], u, nx, dtoverdx)
        for k in range(n_lc):
            P_nl, P_loc = solve_helm(A[t+1], h[t+1], l_c[k], nx, deltax2, Pstar, C_s)
            P_nl_lc[k, t+1] = P_nl
            P_loc_t[t+1] = P_loc

    return P_nl_lc, P_loc_t, h, A

A = np.ones((nt, nx+1))
h = np.zeros_like(A)
for i in range(nx):
    if i < 125:
        h[0, i] = 5.0 
    elif 125 <= i < 250:
        h[0,i] = (-i/125)+3
    else:
        h[0,i] = 1

u = np.ones(nx+1)
u[0:125] = 0.0001
u[125:-1]= 0.5
u[0], u[nx] = 0.0, 0.0

l_c = np.array([10e3, 50e3, 100e3, 150e3, 200e3])

P_nl_lc, P_loc, h, A= run_simulation(A, h, u, l_c, nx, nt,
                                 deltax2, dtoverdx,
                                 param.Pstar, param.C_s)

#-------

x = np.arange(0, nx+1)*deltax



plt.figure()
for t in range(0,nt,int(nt/5)):

    plt.plot(x/1e3, A[t])
    plt.plot(x/1e3, h[t])
# plt.plot(x/1e3, P_loc[-1][0:-1]/1e3, label = 'local')
plt.xlabel('x (km)')
plt.ylabel('A ')
plt.legend()
plt.savefig('a_nl_heml.png')


plt.figure()
for i, P_nl in enumerate(P_nl_lc):

    plt.plot(x[0:-1]/1e3, P_nl[-1]/1e3, label = r'$l_c = {}$ km'.format(int(l_c[i]/1e3)))
plt.plot(x[0:-1]/1e3, P_loc[-1]/1e3, label = 'local')
plt.xlabel('x (km)')
plt.ylabel('P (kN/m)')
plt.legend()
plt.savefig('p_nl_heml.png')
        