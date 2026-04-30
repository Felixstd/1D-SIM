import numpy as np
import matplotlib.pyplot as plt
import scienceplots


# plt.style.use('science')
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def icestrength(A, h, C = 20, Pstar = 27.5e3) : 
    
    return Pstar*h*np.exp(-C*(1-A))

rho_ice = 917
rho_water = 1030
g = 9.8
dudx = abs(np.linspace(-1e-6, 1e-6, 1000000))
Pstar_prime = 27.5e3
h=3
P_0 = 2*(rho_water - rho_ice)/(rho_ice)*g*h

Pstar = P_0 + (Pstar_prime-P_0)*np.tanh(dudx/(2e-9))

plt.figure()
plt.plot(dudx, Pstar)
# plt.axhline(P_0, color = 'k')
# plt.axhline(Pstar_prime, color = 'k')
plt.xscale('log')
plt.grid()
plt.yticks([P_0, Pstar_prime,(Pstar_prime + P_0)/2], [r'$2\rho gh$', r'$P^*$', ''])
plt.xlabel(r'$\Delta$ (s$^{-1}$)')
plt.ylabel(r'$P^*$ (Nm$^{-2}$)')
plt.savefig('pstar_log.png')


h = np.linspace(0, 2, 1000)
# h = np.ones(1000)
A1 = np.linspace(0, 1, 1000)
h1 = np.ones_like(A1)


h2 = np.linspace(0, 1, 1000)
A2 = np.ones_like(h2)

P1 = icestrength(A1, h1)
P2 = icestrength(A2, h2)

plt.figure()
plt.plot(A1, P1/1e3, color = 'k')
plt.grid()
plt.gca().invert_xaxis()
plt.title(r'$h = 1$ m')
plt.xlabel('Sea Ice Concentration A')
plt.ylabel('Ice Strength (kN/m)')
plt.savefig('icestrength_hconstant.png')

plt.figure()
plt.plot(h2, P2/1e3, color = 'k')
plt.grid()
plt.gca().invert_xaxis()
plt.title(r'$A = 1$')
plt.xlabel('Ice Thickness $h$ (m)')
plt.ylabel('Ice Strength (kN/m)')
plt.savefig('icestrength_Aconstant.png')

