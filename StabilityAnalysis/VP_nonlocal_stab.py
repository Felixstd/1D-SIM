import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


A0 = 1
h0 = 2
C = 20
Pstar = 2.75e4
Pstarstar = Pstar*np.exp(-C*(1-A0))
rho = 900
delta_min = 2e-4
e = 2
alpha = np.sqrt(1+e**(-2))
l = 100
lamda_div = (alpha-1)/2


def VPs_growth(k):
    
    
    omega = np.sqrt((Pstarstar*(1+C*A0)/rho)*lamda_div*k**2)
    
    return omega

def VP_nonlocal(k):
    
    
    omega = np.sqrt((Pstarstar*(1+C*A0)/rho)*(lamda_div*k**2-l**2*k**4))+0j
    

    return omega

k = np.linspace(0,1, 100000)


omega_vp_div = VPs_growth(k)
omega_nonlocal_div = VP_nonlocal(k)

print(omega_nonlocal_div)

plt.figure()
plt.plot(k, omega_vp_div, color = 'darkred', label = 'VPs')
plt.plot(k, omega_nonlocal_div, color = 'royalblue', label = 'VPr')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$k$ (1/m)')
plt.ylabel('$\omega_+$ (1/s)')
plt.legend()
plt.savefig('VP_nonlocal.png')