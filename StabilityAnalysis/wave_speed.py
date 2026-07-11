import numpy as np 
import matplotlib.pyplot as plt 
from seaiceparameters import *

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

okabe_ito = ['#CC79A7','#0072B2','#009E73','#E69F00','#D55E00']

dx= 4000
dt= 15
k_crit = 2*np.pi/dx

CFL = dx/dt

# CFL = 4000/1
print(CFL)

def wavespeed_VP(A0):
    P0 = Pstar*np.exp(-C_s*(1-A0))
    phase_speed = np.sqrt((-lamda_conv*P0*(1+C_s*A0)/rhoi))    
    
    return phase_speed


def wavespeed_nVP(A0, h0, k, l = 10e3):
    
    l0 =l*A0
    P0 = Pstar*h0*np.exp(-C_s*(1-A0))
    
    # P0 = Ps
    # print(k**2*l0**2)
    phase_speed = np.sqrt((-lamda_conv/(rhoi*h0)*(2*P0+(C_s*A0-1)*Pstarstar*h0)*(1/(k**2*l0**2+1))))    
    group_speed = np.sqrt(-lamda_conv/(rhoi*h0)*(2*P0+(C_s*A0-1)*Pstarstar*h0))*(1/(k**2*l0**2+1)**(3/2))
    return phase_speed, group_speed

n = 1000

A = np.linspace(0.000001, 1, n)
h = np.linspace(0.000001,5, n)
k = np.logspace(-9, 0, n)

ls = [10e3, 50e3, 100e3, 150e3, 200e3]
k_sample = [1e-5,1e-4, 1e-1]
linestyles= ['-', '--','-.']

omega_nVP_A_l = np.zeros([len(ls),len(k_sample), n])
omega_nVP_h_l = np.zeros([len(ls),len(k_sample), n])
omega_nVP_k_l = np.zeros([len(ls),n])


wavespeed_vp = wavespeed_VP(A)

for i, l in enumerate(ls):
    omega_nVP_k_l[i],omegag_nVP_k_l[i] = wavespeed_nVP(1,1, k, l = l)
    for j, ks in enumerate(k_sample):
        omega_nVP_A_l[i][j],omegag_nVP_A_l[i][j] = wavespeed_nVP(A, 1, ks, l = l)
        omega_nVP_h_l[i][j],omegag_nVP_h_l[i][j] = wavespeed_nVP(1, h, ks, l = l)
        
    # print(omega_nVP_k_l[i])
    
    

fig, axs = plt.subplots(3, 1, figsize = (5, 7), constrained_layout = True)

plt.rcParams.update({
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "legend.fontsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        })

for i, l in enumerate(ls):
    axs[2].plot(k, omega_nVP_k_l[i], label = 'l = {}'.format(int(l/1e3)), color = okabe_ito[i])
    axs[2].plot(k, omegag_nVP_k_l[i], label = 'l = {}'.format(int(l/1e3)), color = okabe_ito[i], linestyle = '--')

    for j, ks in enumerate(k_sample):
    
        axs[0].plot(A, omega_nVP_A_l[i][j], label = 'l = {}, k = {}'.format(int(l/1e3), k_sample[j]),linestyle = linestyles[j], color = okabe_ito[i])
        axs[0].plot(A, omegag_nVP_A_l[i][j], label = 'l = {}, k = {}'.format(int(l/1e3), k_sample[j]),linestyle = '--', color = okabe_ito[i])
        axs[1].plot(h, omega_nVP_h_l[i][j], label = 'l = {}'.format(int(l/1e3)),linestyle = linestyles[j], color = okabe_ito[i])
        axs[1].plot(h, omegag_nVP_h_l[i][j], label = 'l = {}, k = {}'.format(int(l/1e3), k_sample[j]),linestyle = '--', color = okabe_ito[i])
for ax in axs: 
    ax.axhline(CFL)
    
axs[0].plot(A, wavespeed_vp, label = 'VP')
    
axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[2].set_yscale('log')
# axs[0].legend(bbox_to_anchor = (1.8, 0.5))
axs[0].set_xlabel(r'$A$')
axs[1].set_xlabel(r'$h$')
axs[2].set_xlabel(r'$k$')
fig.supylabel(r'$\omega/k$')
fig.savefig('wavespeed.png')