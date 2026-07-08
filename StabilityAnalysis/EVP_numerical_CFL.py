
import numpy as np
import matplotlib.pyplot as plt

rho = 900
zeta_max = 1e12


def CFL_EVP(dt, t0,N, dx_i):
    
    dx = np.sqrt(zeta_max/(rho*t0))*np.sqrt(dt)/N
    
    N = np.sqrt(zeta_max/(rho*t0))*np.sqrt(dt)/dx_i
    
    return dx, N


length = 100
dt = np.linspace(1,30,length)
t0 = 0.36
N  = np.linspace(500, 1200, length)
dx = np.linspace(100, 1000, length)

dx_EVP = CFL_EVP(dt, t0, N,dx)

dx_EVP_15 = CFL_EVP(15, t0, N,dx)


print('value for dx with t0 = 0.36: ', CFL_EVP(15, 0.36,900,100))
print('Default value for dx with  t0 = 0.8: ', CFL_EVP(15, 0.8,700,100))
print('Default value for dx with  t0 = 0.36: ', CFL_EVP(1, 0.36,300,100))


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot()


# 2. Generate grid coordinate pairs via meshgrid
X, Y = np.meshgrid(dt, N)
dx_EVP,_ = CFL_EVP(X, t0, Y,dx)
# 3. Plot surface mesh layout and assign a color map
surf = ax.contourf(X, Y, dx_EVP, cmap='viridis')

# 4. Finalize labels and attach a reference scale bar
ax.set_xlabel('dt')
ax.set_ylabel('N')
# ax.set_zlabel('dx')
fig.colorbar(surf,ax = ax, aspect=10, label = 'dx') 

plt.savefig('CFL_evp.png', dpi = 400)


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot()

# 3. Plot surface mesh layout and assign a color map

for t0 in np.linspace(0.01,1,10):
    dx_EVP_15, N_EVP_15 = CFL_EVP(15, t0, N, dx)
    surf = ax.plot(N,dx_EVP_15, label = t0)

# 4. Finalize labels and attach a reference scale bar
ax.set_xlabel('N')
ax.set_ylabel('dx')
plt.legend()
plt.axhline(250)
# ax.set_zlabel('dx')
# fig.colorbar(surf,ax = ax, aspect=10, label = 'dx') 

plt.savefig('CFL_evp_2.png', dpi = 400)

