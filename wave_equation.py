import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from tqdm import tqdm
import matplotlib.colors as colors
import matplotlib.cm as cm


plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def initial_conditions(x, nx):

    # uinit = np.minimum(np.maximum(((10000-(x-nx/2)**2)**(1/20)), 0), 1)
    # uinit[np.isnan(uinit)] = 0
    
    
    # uinit = 10*np.sin(2*np.pi*x/nx)
    
    
    uinit = 1*np.exp(-((x-nx/2)**2/10000)**5)
    
    print(uinit)
    return uinit

@jit(nopython = True)
def solve_uice(Nt, Nx, uice):
    print('solving for uice')
    for t in range(1, Nt-1):
        
        uice[1:Nx, t+1] = 2*uice[1:Nx, t] - uice[1:Nx, t-1] + \
                C2*(uice[:-2, t] - 2*uice[1:Nx, t] + uice[2:, t])
        # print(uice[t])

    return uice
#---- Constants ----#
dx = 10e3
dt = 1

Nt = 4*24*60*60
Nx = 400

total_time = np.arange(0, Nt, 1, dtype = int)
x = np.arange(0, Nx+1)

#-- ICE --#
Pstar = 27.5e3
h0 = 1
rhoice = 900
m = h0*rhoice

c = np.sqrt(Pstar*h0/m)
C = c*dt/dx
C2 = C**2
print(C)

uice = np.zeros((Nx+1, Nt))

print(np.shape(uice))

#set initial conditions

uice[:, 0] = initial_conditions(x, Nx)

uice[0, :] = 0
uice[-1, :] = 0

#--- First step ---#
uice[1:Nx, 1] = uice[1:Nx, 0] + 0.5 * C2 * (uice[:-2, 0] - 2*uice[1:Nx, 0] + uice[2:, 0])
print(uice[:, 1])
uice = solve_uice(Nt, Nx, uice)
print(uice[:, 1000])

idx = [ 0,  34560,  69120, 103680, 138240, 172800 ,207360 ,241920 ,276480, 311040]
norm = colors.Normalize(vmin=0, vmax=idx[-1]/(60*60))  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    

for k, i in enumerate(idx):
    fig = plt.figure()
    ax = plt.axes()
    color = cmap(k/10)
    plt.plot(x*dx/1e3, uice[:, i], color = color)
    plt.ylim(0, 1)
    fig.colorbar(sm, ax = ax, 
                        label = 'Time (hr)')
    plt.savefig('nonviscid/nonviscid_wave_{}.png'.format(k))

    plt.close()
