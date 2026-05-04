import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

A0 = 1
h0 = 2
C = 20
Pstar = 2.75e4
rho = 900
delta_min = 2e-4
rhow = 1026
rhoprime = (rhow-rho)/rhow*rho
g = 9.8
e = 2
alpha = np.sqrt(1+e**(-2))

lamda_div = (alpha-1)/2
lamda_conv = (-alpha-1)/2

P0 = rhoprime*g*20
deltaP_star = Pstar-P0

print(Pstar, P0)

k = np.linspace(0, 1e-0, 1000000)


omega_plus_div = 1/(rho*h0)*(-lamda_div*(deltaP_star*alpha/delta_min)*h0*k**2 + \
    np.sqrt((lamda_div*deltaP_star*(alpha/delta_min)*h0*k**2)**2 + 4*rho*h0*lamda_div*P0*h0*(1+C*A0)*k**2))

omega_vp = k*np.sqrt(Pstar*np.exp(-C*(1-A0))*lamda_div/rho)


plt.figure(1)
plt.xscale('log')
plt.yscale('log')
plt.plot(k, omega_plus_div, color = 'royalblue', label = 'new P')
plt.plot(k, omega_vp, color = 'darkred', label = 'standard VP')
plt.legend()
plt.ylim(0,1)
plt.xlabel('$k$ (1/m)')
plt.ylabel(r'Re{$\omega_+$} (1/s)')
plt.savefig('omega_p_div_newP.png')