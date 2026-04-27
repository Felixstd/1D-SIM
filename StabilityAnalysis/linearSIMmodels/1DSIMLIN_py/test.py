import numpy as np
import matplotlib.pyplot as plt
from numba import jit
# Parameters
Nx = 400
# Lx = 1e4
dx = 1e4

Lx = dx*Nx
x = np.linspace(0, Lx, Nx)
X = np.arange(0, Nx, 1)
# x = np.linspace()
dt = 0.0005
tmax = 10
# a = c1

# Initial condition
h = np.exp(-((X-Nx/2)**2/5000)**5)
uinit_1 = np.linspace(-2.5,2.5, Nx)
uinit_1[h < 0.2] = 0

width = 0.1
mask_tanh = 0.4 * (1 + np.tanh((h - 0.5) / width))
u = uinit_1 * mask_tanh 

plt.plot(x, h, label='h')
plt.plot(x, u, label='u')
plt.savefig('test1.png')

# HLL flux function
# def flux(q):
#     h = q[0]
#     u = q[1]
#     return np.array([u*h, -a*h])

# def hll(qL, qR):
#     # Roe averages
#     hL,uL = qL
#     hR,uR = qR
#     hbar = 0.5*(hL+hR)
#     ubar = 0.5*(uL+uR)
#     c = np.sqrt(ubar**2 + 4*a*hbar) / 2  # approximate characteristic speed
    
#     # Speeds
#     sL = min(ubar - c, 0.0)
#     sR = max(ubar + c, 0.0)
    
#     FL = flux(qL)
#     FR = flux(qR)
    
#     if sL >= 0: return FL
#     if sR <= 0: return FR
#     return (sR*FL - sL*FR + sL*sR*(qR-qL))/(sR - sL)

# # Time integration (1st order Godunov)
# q = np.vstack([h,u])

# t = 0.0
# plt.figure(figsize=(8,4))
# while t < tmax:
#     # Compute fluxes
#     qL = q[:, :-1]
#     qR = q[:, 1:]
#     F = np.zeros_like(q)
#     for i in range(Nx-1):
#         F[:, i+1] = hll(qL[:,i], qR[:,i])
    
#     # Update
#     q[:,1:-1] -= dt/dx * (F[:,2:] - F[:,1:-1])
    
#     t += dt

# # Extract results
#     h, u = q

#     # if t < 0.001:
#     plt.plot(x, h, label='h')
#     plt.plot(x, u, label='u')
    
# # plt.legend()
# plt.title('Shock capturing solution')
# plt.savefig('test_vel.png')

# import numpy as np
# import matplotlib.pyplot as plt

# # Grid
# Lx = dx*Nx
# x = np.linspace(0, Lx, Nx)
# X = np.arange(0, Nx, 1)
dudx_crit = 1e-10
# # x = np.linspace()
# dt = 1
# tmax = 300000
# # a = c1

# # Initial condition
# h = np.exp(-((X-Nx/2)**2/5000)**5)
e = 2
E_pm_div = (np.sqrt(1+e**(-2))-1)/2
E_pm_conv = (-np.sqrt(1+e**(-2))-1)/2
# # Time
# # dt = 0.00002
# # tmax = 0.05

# # Variable a(x)
# # @jit(nopython = True)
# def a_func(x):
    
#     E_pm = np.zeros_like(u)
#     dudx = np.gradient(x, dx)
#     for i in range(1, Nx-1):
#         if abs(dudx[i]) < dudx_crit: 
#             E_pm[i] = (E_pm_div + E_pm_conv)/2 + ((E_pm_div - E_pm_conv)/2)*dudx[i]/dudx_crit    
                
#         else:
                
#             if dudx[i] < -dudx_crit:
#                 E_pm[i] = E_pm_conv
                        
#                     # elif abs(dudx[i]-1e-12) < 1e-12 :
#                     #     E_pm = 0
#             elif dudx[i] > dudx_crit:
#                 E_pm[i] = E_pm_div
#     return E_pm/(900)

# a = a_func(x)
# da_dx = np.gradient(a, dx)

# # Initial conditions
# # h = np.exp(-((x-0.5)**2)/0.002)
# width = 0.1
# mask_tanh = 0.4 * (1 + np.tanh((h - 0.5) / width))
# uinit_1 = np.linspace(-2.5,2.5, Nx)
# u = uinit_1 * mask_tanh 
# # u = np.zeros_like(x)

# # --- Flux functions ---
# @jit(nopython = True)
# def flux_u(h, a):
#     return -a * h

# @jit(nopython = True)
# def flux_h(u, h):
#     return u * h

# # --- TVD minmod limiter ---
# @jit(nopython = True)
# def minmod(a, b):
#     return np.where(a*b <= 0, 0, np.sign(a)*np.minimum(np.abs(a), np.abs(b)))

# # --- Time integration ---
# t = 0.0
# plt.figure(figsize=(8,4))
# while t < tmax:
#     # Slopes for u and h
#     if t %1000 ==0 : 
#         print(t)
#     du = np.zeros_like(u)
#     dh = np.zeros_like(h)
#     du[1:-1] = minmod(u[1:-1]-u[:-2], u[2:]-u[1:-1])
#     dh[1:-1] = minmod(h[1:-1]-h[:-2], h[2:]-h[1:-1])

#     # Interface states
#     uL = u + 0.5*du
#     uR = u - 0.5*du
#     hL = h + 0.5*dh
#     hR = h - 0.5*dh
    
#     a = a_func(u)

#     # Numerical fluxes (Rusanov)
#     F_u = np.zeros_like(u)
#     F_h = np.zeros_like(h)
#     for i in range(1, Nx-1):
#         # Left/right states at interface i+1/2
#         hL_int = hL[i]
#         hR_int = hR[i+1]
#         uL_int = uL[i]
#         uR_int = uR[i+1]
        
        

#         aL = a[i]
#         aR = a[i+1]

#         # Fluxes left/right
#         Fu_L = flux_u(hL_int, aL)
#         Fu_R = flux_u(hR_int, aR)
#         Fh_L = flux_h(uL_int, hL_int)
#         Fh_R = flux_h(uR_int, hR_int)

#         # Estimate wave speed (local)
#         smax = max(
#         abs(uL_int), abs(uR_int)
#         ) + max(np.sqrt(abs(aL*hL_int)), np.sqrt(abs(aR*hR_int)))

#         # Rusanov flux
#         F_u[i] = 0.5*(Fu_L + Fu_R) - 0.5*smax*(uR_int - uL_int)
#         F_h[i] = 0.5*(Fh_L + Fh_R) - 0.5*smax*(hR_int - hL_int)

#     # Update with source term for u
#     unew = np.zeros_like(u)
#     hnew = np.zeros_like(h)
#     unew[1:-1] = u[1:-1] + dt/dx * (F_u[2:] - F_u[:-2])/2 + dt * (-h[1:-1]*da_dx[1:-1])
#     hnew[1:-1] = h[1:-1] - dt/dx * (F_h[2:] - F_h[:-2])/2

#     # Neumann BC
#     unew[0], unew[-1] = 0, 0
#     hnew[0], hnew[-1] = 0,0
    
#     u = unew
#     h = hnew

#     t += dt
#     plt.plot(x, hnew, label='h(x)')
#     plt.plot(x, unew, label='u(x)')

# # --- Plot ---


# # plt.legend()
# plt.title('Coupled u,h with variable a(x) using TVD MUSCL scheme')
# plt.savefig('test.png')

# Plot result

#     plt.plot(x, h, label='h(x)')
#     plt.plot(x, u, label='u(x)')
# # plt.legend()
# plt.title('Coupled u,h with variable a(x) using Lax–Wendroff + source')
# plt.savefig('test.png')

import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Parameters
# ------------------------
Nx = 1000
# Nx = 400
L =1000.0
dx = L / (Nx-1)

x = np.linspace(0, L, Nx)
X = np.arange(0, Nx, 1)
rho = 900
CFL = 0.1
# psi = 1.0     # can be positive or negative
C = 20.0
Pstar = 27.5e3

# ------------------------
# Initial condition
# ------------------------
A0 = np.exp(-((X-Nx/2)**2/6000)**5)   # Gaussian bump
# u0 = np.linspace(-6, 6, Nx)
u0 = np.zeros_like(A0)
u0[A0 < 0.2] =0

U0= np.zeros((2, Nx))
U0[0] = A0
U0[1] = A0 * u0   # M = A*u
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

# @jit(nopython = True)
@jit(nopython = True)
def psi(x):
    # print(np.shape(x))
    E_pm = np.zeros_like(x)
    dudx = first_derivative(x, dx, Nx)
    # print(np.shape(dudx))
    for i in range(1, Nx-1):
        if abs(dudx[i]) < dudx_crit: 
            E_pm[i] = (E_pm_div + E_pm_conv)/2 + ((E_pm_div - E_pm_conv)/2)*dudx[i]/dudx_crit    
                
        else:
                
            if dudx[i] < -dudx_crit:
                E_pm[i] = E_pm_conv
                        
                    # elif abs(dudx[i]-1e-12) < 1e-12 :
                    #     E_pm = 0
            elif dudx[i] > dudx_crit:
                E_pm[i] = E_pm_div
    return E_pm*(Pstar)/(900)

# ------------------------
# Flux and wave speed
# ------------------------


@jit(nopython = True)
def superbee(r):
    return np.maximum(0, np.minimum(2*r,1), np.minimum(r,2))

@jit(nopython = True)
def minmod(theta):
    return np.maximum(0.0, np.minimum(1.0, theta))
# ------------------------
# Compute limited slope
# ------------------------
@jit(nopython = True)
def limited_slope(U):
    dU = np.zeros_like(U)
    for j in range(U.shape[0]):
        du_fwd = U[j,2:] - U[j,1:-1]
        du_bwd = U[j,1:-1] - U[j,:-2]
        r = np.zeros_like(du_fwd)
        mask = np.abs(du_fwd) > 1e-12
        r[mask] = du_bwd[mask] / du_fwd[mask]
        phi = superbee(r)
        dU[j,1:-1] = 0.5 * phi * (U[j,2:] - U[j,:-2])
    return dU

@jit(nopython= True)
def flux(U, psi_local, C):
    A = U[0]
    M = U[1]
    u = M / (A) if A > 1e-12 else 0.0
    return np.array([
        M,
       M*u- psi_local * A * np.exp(-C*(1-A))# M*u - psi_local * A * np.exp(-C*(1-A))
    ])
@jit(nopython= True)
def max_eigen(A, u, psi_local, C):
    c = np.sqrt(abs(psi_local) * np.exp(-C*(1-A)))
    return abs(u) + c

# ------------------------
# Superbee limiter


# ------------------------
# One step (Lax-Wendroff + Superbee)
# ------------------------
# @jit(nopython= True)
# def step(U, dx, dt, psi_profile, C):
#     Nx = U.shape[1]

#     # Compute limited slopes for MUSCL
#     dU = limited_slope(U)

#     # Left/right states at interfaces
#     UL = U - 0.5*dU
#     UR = U + 0.5*dU

#     U_new = U.copy()

#     for i in range(1, Nx-1):
#         psi_L  = psi_profile[i]
#         psi_R  = psi_profile[i+1]
#         psi_Lm = psi_profile[i-1]

#         psi_half = 0.5*(psi_L + psi_R)
#         psi_half_m = 0.5*(psi_Lm + psi_L)

#         # States
#         UL_i = UR[:, i-1]
#         UR_i = UL[:, i]

#         UL_ip = UR[:, i]
#         UR_ip = UL[:, i+1]

#         # Safe velocities
#         def safe_u(A,M): return M/A if A>1e-12 else 0.0

#         A_L = UL_i[0]; u_L = safe_u(A_L, UL_i[1])
#         A_R = UR_i[0]; u_R = safe_u(A_R, UR_i[1])
#         A_Lp = UL_ip[0]; u_Lp = safe_u(A_Lp, UL_ip[1])
#         A_Rp = UR_ip[0]; u_Rp = safe_u(A_Rp, UR_ip[1])

#         FL = flux(UL_i, psi_half_m, C)
#         FR = flux(UR_i, psi_half_m, C)
#         FLp = flux(UL_ip, psi_half, C)
#         FRp = flux(UR_ip, psi_half, C)

#         # Wave speeds
#         alpha_m = max(max_eigen(A_L, u_L, psi_half_m, C),
#                       max_eigen(A_R, u_R, psi_half_m, C))
#         alpha_p = max(max_eigen(A_Lp, u_Lp, psi_half, C),
#                       max_eigen(A_Rp, u_Rp, psi_half, C))
        
#         # def local_alpha(UL, UR, FL, FR, psi_val):
#         #     dU = UR - UL
#         #     dF = FR - FL
#         #     alpha_local = 0.0
#         #     for k in range(len(U)):
#         #         if abs(dU[k]) > 1e-12:
#         #             alpha_local = max(alpha_local, abs(dF[k]/dU[k]))
#         #     # Safe fallback if dU ~ 0 (u=0 or flat field)
#         #     if alpha_local < 1e-12:
#         #         alpha_local = abs(psi_val)
#         #     return alpha_local

#         # alpha_m = local_alpha(UL_i, UR_i, FL, FR, psi_half)
#         # alpha_p = local_alpha(UL_ip, UR_ip, FLp, FRp, psi_half_m)

#         # Scale by dt/dx to make it a CFL-like number
#         alpha_m *= dt/dx
#         alpha_p *= dt/dx
#         dU_up = U[:, i] - U[:, i-1]
#         dU_dn = U[:, i+1] - U[:, i]

#         theta = np.zeros(2)
#         mask = np.abs(dU_dn) > 1e-12
#         theta[mask] = dU_up[mask] / (dU_dn[mask] + 1e-12)

#         phi = superbee(theta)
#         # phi = minmod(theta)

#         # Split flux into LF (first-order) + LW correction
#         F_LF_m = 0.5*(FL + FR) - 0.5*alpha_m*(UR_i - UL_i)
#         F_LF_p = 0.5*(FLp + FRp) - 0.5*alpha_p*(UR_ip - UL_ip)

#         F_LW_m = 0.5 * (FR - FL)
#         F_LW_p = 0.5 * (FRp - FLp)

#         F_half_m = F_LF_m + phi * F_LW_m
#         F_half_p = F_LF_p + phi * F_LW_p

#         U_new[:, i] = U[:, i] - dt/dx * (F_half_p - F_half_m)
#         # print((F_half_p - F_half_m))
#     # Dirichlet BCs
#     U_new[:, 0] = 0.0
#     U_new[:, -1] = 0.0

#     # Cap A <= 1
#     A_cap = 1.0
#     Avals = U_new[0]
#     Mvals = U_new[1]
#     over = Avals > A_cap
#     Mvals[over] *= A_cap / Avals[over]
#     Avals[over] = A_cap

#     return U_new

# @jit(nopython= True)
# def step(U, dx, dt, psi, C):
#     U_new = U.copy()
#     Nx = U.shape[1]
    
#     print
#     for i in range(1, Nx-1):
#         UL = U[:, i]
#         UR = U[:, i+1]
#         ULm = U[:, i-1]

#         # Safe velocities
#         A_L = UL[0]; u_L = UL[1]/A_L if A_L > 1e-12 else 0.0
#         A_R = UR[0]; u_R = UR[1]/A_R if A_R > 1e-12 else 0.0
#         A_Lm = ULm[0]; u_Lm = ULm[1]/A_Lm if A_Lm > 1e-12 else 0.0

#         FL = flux(UL, psi, C)
#         FR = flux(UR, psi, C)
#         FLm = flux(ULm, psi, C)

#         # Lax-Friedrichs fluxes
#         alpha = max(max_eigen(A_L, u_L, psi, C),
#                     max_eigen(A_R, u_R, psi, C))
#         alphaL = max(max_eigen(A_Lm, u_Lm, psi, C),
#                      max_eigen(A_L, u_L, psi, C))

#         F_half = 0.5*(FL + FR) - 0.5*alpha*(UR - UL)
#         F_half_m = 0.5*(FLm + FL) - 0.5*alphaL*(UL - ULm)

#         U_new[:, i] = U[:, i] - dt/dx * (F_half - F_half_m)

#     # Dirichlet BC: A=0, u=0 at both ends
#     U_new[:, 0] = 0.0
#     U_new[:, -1] = 0.0

#     # Cap A to 1 everywhere and adjust M to keep u consistent
#     A_cap = 1.0
#     Avals = U_new[0]
#     Mvals = U_new[1]
#     over = Avals > A_cap
#     Mvals[over] *= A_cap / Avals[over]
#     Avals[over] = A_cap

#     return U_new
@jit(nopython= True)
def max_wave_speed(U, psi, C, eps=1e-6):
    n = U.shape[0]
    A = U[0]
    M = U[1]

    # Return zero speed for dry cell
    if A < 1e-8 and np.abs(M) < 1e-8:
        return 0.0

    # Compute baseline flux
    F0 = flux(U, psi, C)
    J = np.zeros((n, n))

    for j in range(n):
        dU = np.zeros(n)
        dU[j] = eps
        F1 = flux(U + dU, psi, C)

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
    
@jit(nopython= True)
def step(U, dx, dt, psi_profile, C):
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

        # States
        UL_i = UR[:, i-1]
        UR_i = UL[:, i]

        UL_ip = UR[:, i]
        UR_ip = UL[:, i+1]

        # Safe velocities
        def safe_u(A,M): return M/(A) if A>1e-12 else 0.0

        A_L = UL_i[0]; u_L = safe_u(A_L, UL_i[1])
        A_R = UR_i[0]; u_R = safe_u(A_R, UR_i[1])
        A_Lp = UL_ip[0]; u_Lp = safe_u(A_Lp, UL_ip[1])
        A_Rp = UR_ip[0]; u_Rp = safe_u(A_Rp, UR_ip[1])

        FL = flux(UL_i, psi_half_m, C)
        FR = flux(UR_i, psi_half_m, C)
        FLp = flux(UL_ip, psi_half, C)
        FRp = flux(UR_ip, psi_half, C)

        # Wave speeds
        # alpha_m = max(max_eigen(A_L, u_L, psi_half_m, C),
        #               max_eigen(A_R, u_R, psi_half_m, C))
        # alpha_p = max(max_eigen(A_Lp, u_Lp, psi_half, C),
        #               max_eigen(A_Rp, u_Rp, psi_half, C))
        
        alpha_m = max(
            max_wave_speed(UL_i, psi_half_m, C),
            max_wave_speed(UR_i, psi_half_m, C)
            )
        alpha_p = max(
                 max_wave_speed(UL_ip, psi_half, C),
                 max_wave_speed(UR_ip, psi_half, C)
                )

        # Numerical fluxes (TVD Lax–Wendroff form)
        F_half_m = 0.5*(FL + FR) - 0.5*alpha_m*(UR_i - UL_i)
        F_half_p = 0.5*(FLp + FRp) - 0.5*alpha_p*(UR_ip - UL_ip)

        U_new[:, i] = U[:, i] - dt/dx * (F_half_p - F_half_m)

    # Dirichlet BCs
    U_new[:, 0] = 0.0
    U_new[:, -1] = 0.0

    # Cap A <= 1
    # A_cap = 1.0
    # Avals = U_new[0]
    # Mvals = U_new[1]
    # over = Avals > A_cap
    # Mvals[over] *= A_cap / Avals[over]
    # Avals[over] = A_cap

    return U_new


# Time integration
# @jit(nopython = True)
def solve(U):
    t = 0.0
    c= 0
    dt = 1e-2
    tend = 10
    time = np.arange(0, tend, dt)
    for tmax in time:
        if (c%1000 == 0): print(tmax)
        # print(t)
        A = U[0]
        M = U[1]
        u = np.zeros_like(A)
        mask = A > 1e-12
        u[mask] = M[mask]/(A[mask])
        
        # print(np.shape(u))
        
        psi_profile = psi(u)

        # lamb = 1e-12
        # for i in range(1, Nx-1):
        #     FL = flux(U[:,i-1], psi_profile[i-1], C)
        #     FR = flux(U[:,i], psi_profile[i], C)
        #     dF = FR - FL
        #     dU = U[:,i] - U[:,i-1]
        #     local_alpha = np.max(np.abs(dF / (dU + 1e-12)))
        #     lamb = max(lamb, local_alpha)

        lamb = 0.0
        for i in range(Nx):
            lamb = max(lamb, max_wave_speed(U[:, i], psi_profile[i], C))
        CFL = dt* lamb/dx
        
        
        CFL_adv = ((np.sqrt(1+e**(-2))-1)/2 * Pstar/rho)**(1/2) * dt/dx
        # print(CFL, CFL_adv)
        # if t + dt > tmax:
        #     dt = tmax - t

        Unew = step(U, dx, dt, psi_profile, C)
        # print(Unew[0, 500])
        U = Unew
        
        # print(U)
        t += dt
        c += 1
    return U

# print(np.shape(U0))
Unew = solve(U0)

# ------------------------
# Extract results
# ------------------------
A = Unew[0]
u = np.zeros_like(A)
mask = A > 1e-12
u[mask] = Unew[1][mask]/(A[mask])

# ------------------------
# Plot
# ------------------------
plt.figure(figsize=(8,4))
plt.subplot(211)
plt.plot(x, A0, 'k--', label='Initial A')
plt.plot(x, A, 'b', label='Final A')
# plt.xlim(300, 700)
plt.legend(); plt.ylabel("A")

plt.subplot(212)
plt.plot(x, u, 'r', label='Final u')
plt.plot(x, u0, 'k--', label='Final u')
plt.legend(); plt.ylabel("u"); plt.xlabel("x")

plt.tight_layout()
plt.show()

plt.savefig('test.png', dpi = 500)



# def step(U, dx, dt, psi_profile, C):
#     Nx = U.shape[1]
#     dU = limited_slope(U)

#     UL = U - 0.5*dU
#     UR = U + 0.5*dU

#     U_new = U.copy()

#     for i in range(1, Nx-1):
#         psi_half_m = 0.5*(psi_profile[i-1] + psi_profile[i])
#         psi_half_p = 0.5*(psi_profile[i] + psi_profile[i+1])

#         ULm = UR[:, i-1]
#         URm = UL[:, i]

#         ULp = UR[:, i]
#         URp = UL[:, i+1]

#         FL = flux(ULm, psi_half_m, C)
#         FR = flux(URm, psi_half_m, C)
#         FLp = flux(ULp, psi_half_p, C)
#         FRp = flux(URp, psi_half_p, C)

#         # Discrete Jacobian approximation for alpha
#         dF_m = FR - FL
#         dU_m = URm - ULm
#         alpha_m = np.max(np.abs(dF_m / (dU_m + 1e-12)))
#         # alpha_m = dt/dx*np.linalg.norm(FR - FL, ord=np.inf) / (np.linalg.norm(URm - ULm, ord=np.inf) + 1e-12)

#         dF_p = FRp - FLp
#         dU_p = URp - ULp
#         alpha_p = np.max(np.abs(dF_p / (dU_p + 1e-12)))
#         # alpha_p =  dt//dx*np.linalg.norm(FRp - FLp, ord=np.inf) / (np.linalg.norm(URp - ULp, ord=np.inf) + 1e-12)
#         # Numerical fluxes
#         F_half_m = 0.5*(FL + FR) - 0.5*alpha_m*(URm - ULm)
#         F_half_p = 0.5*(FLp + FRp) - 0.5*alpha_p*(URp - ULp)

#         U_new[:, i] = U[:, i] - dt/dx * (F_half_p - F_half_m)

#     # Dirichlet BCs
#     U_new[:, 0] = 0.0
#     U_new[:, -1] = 0.0

#     # Cap A <= 1
#     A_cap = 1.0
#     Avals = U_new[0]
#     Mvals = U_new[1]
#     over = Avals > A_cap
#     Mvals[over] *= A_cap / Avals[over]
#     Avals[over] = A_cap

#     return U_new

    
# def step(U):
    
#     # Nx = U.shape[1]
    
#     h = U[0]
#     M = U[1]
#     u = np.zeros_like(h)
#     mask = h > 1e-12
#     u[mask] = M[mask]/h[mask]
    
#     psi_VP = psi_ViscousPlastic(u)
#     # print(psi_VP)
#     # dU = limited_slope(U)
    
#     f = flux(U, psi_VP)
#     # print(f)
#     #-- Compute courant local number --#

#     nuL = np.zeros_like(f) #Left
#     nuR = np.zeros_like(f) #Right
#     mask = np.zeros_like(nuL)
    
#     for i in range(f.shape[0]):
        
        
        
#         # for j in range(1, Nx):
            
#             df_fwd = f[i,2:] - f[i,1:-1]
            
#             df_bwd = f[i, 1:-1] - f[i, :-2]
#             dU_fwd = U[i, 2:] - U[i, 1:-1]
#             dU_bwd = U[i,1:-1] - U[i,:-2]
            
            
#             mask_fwd = np.abs(dU_fwd) > 1e-12
#             nuR_i = np.zeros_like(df_fwd)
#             nuR_i[mask_fwd] = np.abs(df_fwd[mask_fwd] / dU_fwd[mask_fwd])

#             # Compute nuL only where dU_bwd is nonzero
#             mask_bwd = np.abs(dU_bwd) > 1e-12
#             nuL_i = np.zeros_like(df_bwd)
#             nuL_i[mask_bwd] = np.abs(df_bwd[mask_bwd] / dU_bwd[mask_bwd])

#             # Assign to global arrays
#             nuR[i, 1:-1] = DtoverDx*nuR_i
#             nuL[i, 1:-1] = DtoverDx*nuL_i
#             # nuR[i, j] = DtoverDx*(f[i, 2:] - f[i, 1:-1])/(U[i, 2:] - U[i, 1:-1] + 1e-12)

    
#     F_low_R = np.zeros_like(nuR)
#     F_low_L = np.zeros_like(nuR)
#     F_High_L = np.zeros_like(nuR)
#     F_High_R = np.zeros_like(nuR)
    
#     theta_R = np.zeros_like(nuR)
#     theta_L = np.zeros_like(nuR)
    
#     for j in range(f.shape[0]):
#         for i in range(1, Nx-1):
            
#             #low order fluxes
#             nu_i_R = nuR[j,  i]
#             nu_i_L = nuL[j,  i]
            
#             if nu_i_R > 0:
                
#                 F_low_R[j, i] = f[j, i]
#                 F_High_R[j, i] = 0.5*(1-nu_i_R)*(f[j, i+1] - f[j, i])
            
#             if nu_i_R < 0:
                
#                 F_low_R[j, i] = f[j, i+1]
#                 F_High_R[j, i] = -0.5*(1-nu_i_R)*(f[j, i+1] - f[j, i])
                
#             if nu_i_L >0:
                
#                 F_low_L[j, i] = f[j, i-1]
#                 F_High_L[j, i] = 0.5*(1-nu_i_L)*(f[j, i] - f[j, i-1])
            
#             if nu_i_L <0:
                
#                 F_low_L[j, i] = f[j, i]
#                 F_High_L[j, i] = -0.5*(1-nu_i_L)*(f[j, i] - f[j, i-1])
                
#             # print(nu_i_R)
#             # idx_prime_R = int(i - np.sign(nu_i_R))
#             theta_R[j, i] = (U[j, i] - U[j, i-1])/(U[j, i+1] - U[j, i]+ 1e-12)
            
#             # idx_prime_R = int(i - np.sign(nu_i_L))
#             theta_L[j, i] = (U[j, i-1] - U[j, i-2])/(U[j, i] - U[j, i-1]+ 1e-12)
    
    
#     phi_R = np.maximum(0, np.minimum(2*theta_R, 1), np.minimum(theta_R, 2))
#     phi_L = np.maximum(0, np.minimum(2*theta_L, 1), np.minimum(theta_L, 2))
    
#     Unew = U - DtoverDx*(F_low_R + F_High_R*phi_R - (F_low_L + F_High_L*phi_L))
#     Unew[:, 0] = 0.0
#     Unew[:, -1] = 0.0
#     Unew[0] = np.clip(Unew[0], 0.0, 1.0)
    
#     # print(Unnew)
    
#     print(DtoverDx*(F_low_R + F_High_R*phi_R - (F_low_L + F_High_L*phi_L)))
    
#     return Unew