import numpy as np

Nx =1001
L = 1000
dx = L/(Nx-1)
# dx = 1
# L = Nx*dx
X = np.arange(0, L+dx, dx)
# print(np.shape(X))
x = np.arange(0, Nx, 1)
# x = np.linspace(0, L, Nx)
# print(np.shape(x))
dt = 1e-4
numax = 1e12
Nt = 1000000
Nt_saves = 10

flux_Gray = True
thickness = False
VP = False
CFL = True

# Nt = 8e4# Nt = 2

# h0 = np.ones(Nx+1)
# # A0 = np.ones(Nx)
# A0 = 1


#---- Rheological Parameters ----#
e = 2
C = 20
rhoice = 900
Pstar = 27.5e3
# Pstarstar = Pstar*np.exp(-C*(1-A0))
zeta_max = 1e12
E_pm_div = (np.sqrt(1+e**(-2))-1)/2
E_pm_conv = (-np.sqrt(1+e**(-2))-1)/2


mu0 = 0.2
muinf = 0.8
mub = 10
deltamu = muinf-mu0
dmean = 1e3
I0 = 1e-3
dudx_crit = 1e-10

DtoverDx = dt/dx

time = np.arange(0, int(Nt), 1)

A0 = np.exp(-((x-Nx/2)**2/6000)**5)   # Gaussian bump
h0 = A0
# A_time = np.ones_like(A0)
# A0 = np.exp(-((x-Nx/2)**2/1000)**5)   # Gaussian bump
# u0 = np.linspace(-3, 3, Nx)
u0 = np.zeros_like(A0)
u0[A0>0] = (x[A0>0]-Nx/2)
u0[A0 < 1e-10] =0
u0 = np.zeros_like(A0)


U = np.zeros((2, Nx))
U[0] = A0
U[1] = A0 * u0 