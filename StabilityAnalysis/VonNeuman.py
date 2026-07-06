import numpy as np
import matplotlib.pyplot as plt

zeta_max = 1e12
dx = 1e3
rho = 900

exp = np.arange(1,5)
k = np.linspace(0,1,1000000)
dxs = 10**exp

plt.figure()
for dx in dxs:
    
    amp_factor = 1/(1+zeta_max/(dx**2*rho)*((2-np.cos(k*dx))))
    plt.plot(k,amp_factor, label = 'dx = '+str(dx/1e3)+' km')

plt.legend()
plt.xlabel('k')
plt.ylabel('G')
plt.xscale('log')
plt.yscale('log')
plt.title('Von Neumann Analysis - Viscous Regime')
plt.savefig('VonNeumann.png')
    
