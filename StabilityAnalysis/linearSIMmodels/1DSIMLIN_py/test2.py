import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Parameters
# ------------------------
Nx = 200
L = 2.0
dx = L / (Nx-1)
x = np.linspace(0, L, Nx)

CFL = 0.4
C = 2.0
tmax = 1.0

# Spatially varying psi(x), can include negatives
psi_profile = -1.0 + 2.0*np.sin(2*np.pi*x)

# Initial condition
A0 = 1.0 + 0.2*np.exp(-100*(x-1.0)**2)
u0 = np.zeros_like(x)

U = np.zeros((2, Nx))
U[0] = A0
U[1] = A0 * u0

# ------------------------
# Flux and eigenvalues
# ------------------------
def flux(U, psi_local, C):
    A = U[0]
    M = U[1]
    u = M / A if A > 1e-12 else 0.0
    return np.array([
        M,
        M*u - psi_local * A * np.exp(-C*(1-A))
    ])

def max_eigen(A, u, psi_local, C):
    c = np.sqrt(abs(psi_local) * np.exp(-C*(1-A)))
    return abs(u) + c

# ------------------------
# Superbee limiter
# ------------------------
def superbee(r):
    return np.maximum(0, np.maximum(np.minimum(2*r,1), np.minimum(r,2)))

# ------------------------
# Compute limited slope
# ------------------------
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

# ------------------------
# One step (Lax-Wendroff + Superbee)
# ------------------------
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
        def safe_u(A,M): return M/A if A>1e-12 else 0.0

        A_L = UL_i[0]; u_L = safe_u(A_L, UL_i[1])
        A_R = UR_i[0]; u_R = safe_u(A_R, UR_i[1])
        A_Lp = UL_ip[0]; u_Lp = safe_u(A_Lp, UL_ip[1])
        A_Rp = UR_ip[0]; u_Rp = safe_u(A_Rp, UR_ip[1])

        FL = flux(UL_i, psi_half_m, C)
        FR = flux(UR_i, psi_half_m, C)
        FLp = flux(UL_ip, psi_half, C)
        FRp = flux(UR_ip, psi_half, C)

        # Wave speeds
        alpha_m = max(max_eigen(A_L, u_L, psi_half_m, C),
                      max_eigen(A_R, u_R, psi_half_m, C))
        alpha_p = max(max_eigen(A_Lp, u_Lp, psi_half, C),
                      max_eigen(A_Rp, u_Rp, psi_half, C))

        # Numerical fluxes (TVD Lax–Wendroff form)
        F_half_m = 0.5*(FL + FR) - 0.5*alpha_m*(UR_i - UL_i)
        F_half_p = 0.5*(FLp + FRp) - 0.5*alpha_p*(UR_ip - UL_ip)

        U_new[:, i] = U[:, i] - dt/dx * (F_half_p - F_half_m)

    # Dirichlet BCs
    U_new[:, 0] = 0.0
    U_new[:, -1] = 0.0

    # Cap A <= 1
    A_cap = 1.0
    Avals = U_new[0]
    Mvals = U_new[1]
    over = Avals > A_cap
    Mvals[over] *= A_cap / Avals[over]
    Avals[over] = A_cap

    return U_new

# ------------------------
# Time integration
# ------------------------
t = 0.0
while t < tmax:
    A = U[0]
    M = U[1]
    u = np.zeros_like(A)
    mask = A > 1e-12
    u[mask] = M[mask]/A[mask]

    lamb = 0.0
    for i in range(Nx):
        lamb = max(lamb, max_eigen(A[i], u[i], psi_profile[i], C))
    dt = CFL * dx / lamb
    if t + dt > tmax:
        dt = tmax - t

    U = step(U, dx, dt, psi_profile, C)
    t += dt

# ------------------------
# Extract final state
# ------------------------
A = U[0]
u = np.zeros_like(A)
mask = A > 1e-12
u[mask] = U[1][mask]/A[mask]

# ------------------------
# Plot
# ------------------------
plt.figure(figsize=(8,4))
plt.subplot(211)
plt.plot(x, A0, 'k--', label='Initial A')
plt.plot(x, A, 'b', label='Final A')
plt.ylim(0, 1.2)
plt.legend(); plt.ylabel("A")

plt.subplot(212)
plt.plot(x, u, 'r', label='Final u')
plt.legend(); plt.ylabel("u"); plt.xlabel("x")

plt.tight_layout()
plt.show()
plt.savefig('test2.png')