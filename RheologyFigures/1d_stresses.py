import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')
e = 2
Pstar = 27.5e3
mu_b = 1.2
mu0 = 0.2
muinf = 0.8
delta_mu = muinf-mu0
dmean = 1000
rhoice = 900


def VP_rheology(dudx, dudx_c = 1e-8):
    
    sigma_11 = np.zeros_like(dudx)
    psi_d = 1/2*((1+e**(-2))**(1/2)-1)
    psi_c = 1/2*(-(1+e**(-2))**(1/2)-1)
    
    sigma_11[abs(dudx) < dudx_c] = (psi_d+psi_c)/2 + (psi_d - psi_c)/2*(dudx[abs(dudx) < dudx_c]/dudx_c)
    
    sigma_11[dudx > dudx_c] = psi_d
    sigma_11[dudx < -dudx_c] = psi_c
    
    
    return sigma_11

def GC_rheology(dudx, dudx_c = 1e-8):
    sigma_11 = np.zeros_like(dudx)
    
    I = dmean*abs(dudx)*np.sqrt(rhoice/Pstar)
    mu = mu0 + delta_mu/(1e-3/I+1)
    
    psi_d = ((mu/2 + mu_b)-1)
    psi_c = (-(mu/2 + mu_b)-1)
    
    sigma_11[abs(dudx) < dudx_c] = (psi_d[abs(dudx) < dudx_c]+psi_c[abs(dudx) < dudx_c])/2 + (psi_d[abs(dudx) < dudx_c] - psi_c[abs(dudx) < dudx_c])/2*(dudx[abs(dudx) < dudx_c]/dudx_c)
    
    sigma_11[dudx > dudx_c] = psi_d[dudx > dudx_c]
    sigma_11[dudx < -dudx_c] = psi_c[dudx < -dudx_c]
    
    return sigma_11
    

dudx = np.linspace(-1e-6, 1e-6, 1000000)

sigma_11_VP = VP_rheology(dudx)

sigma_11_GC = GC_rheology(dudx)

plt.figure(figsize= (4, 3))
ax = plt.axes()
# ax.set_xscale('symlog', linthresh=1e-8  )
plt.axhline(0, color = 'lightgray')
# plt.axvline(0, color = 'lightgray')
# plt.text(0.5, 1, 'Divergence')
plt.plot(dudx, sigma_11_VP, label = 'VP', color = 'r')
plt.plot(dudx, sigma_11_GC, label = 'GC', color = 'b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.grid()
plt.legend()

# plt.xscale('log')
plt.xlabel(r'$\frac{\partial{u}}{\partial{x}}$')
plt.ylabel(r'$\sigma_{11}/P^*$')
plt.savefig('sigma_11.png')