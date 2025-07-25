import numpy as np
import matplotlib.pyplot as plt
import scienceplots


# plt.style.use('science')
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def icestrength(A, h, C = 20, Pstar = 27.5e3) : 
    
    return Pstar*h*np.exp(-C*(1-A))

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

