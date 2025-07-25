"""
Code made by Felix St-Denis to analyse the Plastic 
growth rates of the GC and VP rheologie.

Reference: St-Denis, F. et Al., In Prep. 

23 July 2025
"""

# %% 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

# %%
#----- Global Constants ------#
h0        = 1
A0        = 1
C         = 20
rhoice    = 900
Pstar     = 27.5e3
Pstarstar = Pstar*np.exp(-C*(1-A0))

#---- GC constants ----#
mu0     = 0.2
muinf   = 0.8
deltamu = muinf-mu0
dmean   = 1e3
I0      = 1e-3


def plastic_VP(k, e):
    """
    This is the growth rate for the VP Rheology assuming that 
    
    u = u'
    h = h_0 + h'
    A = A_0
    
    The difference from convergence and divergence are given the parameter
    E_pm. 


    Args:
        k (float, array): wavenumber
        e (float): elliptical ratio
    """
    
    E_div  = (np.sqrt(1+e**(-2))-1)/2
    E_conv = -(np.sqrt(1+e**(-2))-1)/2
    
    # Divergence 
    w_plus_div = np.sqrt(E_div*Pstarstar/rhoice + 0j)*k
    w_minus_div = -np.sqrt(E_div*Pstarstar/rhoice + 0j)*k 
    
    # Convergence 
    w_plus_conv = np.sqrt(E_conv*Pstarstar/rhoice +  0j)*k
    w_minus_conv = -np.sqrt(E_conv*Pstarstar/rhoice + 0j)*k
    
    return [w_plus_div, w_minus_div], [w_plus_conv, w_minus_conv]



def plastic_GC(k, mub):
    
    """
    This is the growth rate for the GC Rheology assuming that 
    
    u = u'
    h = h_0 + h'
    A = A_0
    
    The difference from convergence and divergence are given the parameter
    gamma_pm. 


    Args:
        k (float, array): wavenumber
        e (float): elliptical ratio
    """
    
    gamma_div  = mu0/2 + mub - 1
    gamma_conv = -mu0/2 - mub - 1
    
    Gamma = h0*np.sqrt(Pstarstar*rhoice)*(dmean*deltamu/(2*I0))
    
    #divergence
    w_plus_div = Gamma*k**2/(2*rhoice*h0)*(-1 + np.sqrt(1+4*rhoice*Pstarstar*gamma_div*h0**2/(Gamma**2*k**2)+ 0j))
    w_minus_div = Gamma*k**2/(2*rhoice*h0)*(-1 - np.sqrt(1+4*rhoice*Pstarstar*gamma_div*h0**2/(Gamma**2*k**2)+ 0j))
    
    #convergence
    w_plus_conv = Gamma*k**2/(2*rhoice*h0)*(-1 + np.sqrt(1+4*rhoice*Pstarstar*gamma_conv*h0**2/(Gamma**2*k**2)+ 0j)) 
    w_minus_conv = Gamma*k**2/(2*rhoice*h0)*(-1 - np.sqrt(1+4*rhoice*Pstarstar*gamma_conv*h0**2/(Gamma**2*k**2)+ 0j)) 
    
    
    return [w_plus_div, w_minus_div], [w_plus_conv, w_minus_conv]

# %%
#---- Computing the growth rates -----#
e = 2 
mub_tot = [1/2, 10]
k = np.linspace(1e-7, 1, 10000000)

w_div_VP, w_conv_VP = plastic_VP(k, e)

w_div_GC_mub_list, w_conv_GC_mub_list = zip(*[plastic_GC(k, mub) for mub in mub_tot])
print(np.shape(w_div_GC_mub_list))
# Convert lists of [w₊, w₋] into 2D numpy arrays
w_div_GC_mub  = np.array(w_div_GC_mub_list)   # shape: (N, 2)
w_conv_GC_mub = np.array(w_conv_GC_mub_list)  # shape: (N, 2)

# %%
#---- Plotting -----#
colors = mpl.colormaps['Dark2'].colors

colors = ['darkred', 'darkorange', 'seagreen']

fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,8), sharex = True, sharey = True)

#-- Real part on the top --#

ax1.plot(k, np.real(w_div_VP[0]), color = colors[0])
ax1.plot(k, np.real(w_conv_VP[0]), linestyle = '--', color = colors[0])


ax1.plot(k, np.real(w_div_GC_mub[0][0]), color = colors[1])
ax1.plot(k, np.real(w_conv_GC_mub[0][0]), linestyle = '--', color = colors[1])

ax1.plot(k, np.real(w_div_GC_mub[1][0]), color = colors[2])
ax1.plot(k, np.real(w_conv_GC_mub[1][0]), linestyle = '--', color = colors[2])

#-- Im part on the bottom --#

ax2.plot(k, np.imag(w_div_VP[0]), color = colors[0])
ax2.plot(k, np.imag(w_conv_VP[0]), linestyle = '--', color = colors[0])


ax2.plot(k, np.imag(w_div_GC_mub[0][0]), color = colors[1])
ax2.plot(k, np.imag(w_conv_GC_mub[0][0]), linestyle = '--', color = colors[1])

ax2.plot(k, np.imag(w_div_GC_mub[1][0]), color = colors[2])
ax2.plot(k, np.imag(w_conv_GC_mub[1][0]), linestyle = '--', color = colors[2])

for ax in [ax1, ax2]:
    ax.grid()
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=1e-5)
    ax.set_ylim(-0.001, 0.001)

labels = []

labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='--',
                        label='Convergence')] + labels
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Divergence')] + labels #+ [mlines.Line2D([], [], color='darkviolet', marker='None', linestyle='None',
labels = [mlines.Line2D([], [], color=colors[1], marker='None', linestyle='None',
                        label=r'GC: $\mu_b = {}$'.format(mub_tot[0]))]+labels
labels = [mlines.Line2D([], [], color=colors[2], marker='None', linestyle='None',
                        label=r'GC: $\mu_b = {}$'.format(mub_tot[1]))]+labels
labels = [mlines.Line2D([], [], color=colors[0], marker='None', linestyle='None',
                        label=r'VP')]+labels
                        # label='VP')]
fig.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig.supxlabel(r'$k$ (1/m)', y = 0.05)
ax1.set_ylabel(r'$\text{Re}\{\omega_+\}$ (s$^{-1}$)', fontsize = 13)
fig.align_ylabels()
ax2.set_ylabel(r'$\text{Im}\{\omega_+\}$ (s$^{-1}$)', fontsize = 13)
plt.savefig('plastic_regime.png')


fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10,8), sharex = True, sharey = True)

#-- Real part on the top --#

ax1.plot(k, np.real(w_div_VP[1]), color = colors[0])
ax1.plot(k, np.real(w_conv_VP[1]), linestyle = '--', color = colors[0])


ax1.plot(k, np.real(w_div_GC_mub[0][1]), color = colors[1])
ax1.plot(k, np.real(w_conv_GC_mub[0][1]), linestyle = '--', color = colors[1])

ax1.plot(k, np.real(w_div_GC_mub[1][1]), color = colors[2])
ax1.plot(k, np.real(w_conv_GC_mub[1][1]), linestyle = '--', color = colors[2])

#-- Im part on the bottom --#

ax2.plot(k, np.imag(w_div_VP[1]), color = colors[0])
ax2.plot(k, np.imag(w_conv_VP[1]), linestyle = '--', color = colors[0])


ax2.plot(k, np.imag(w_div_GC_mub[0][1]), color = colors[1])
ax2.plot(k, np.imag(w_conv_GC_mub[0][1]), linestyle = '--', color = colors[1])

ax2.plot(k, np.imag(w_div_GC_mub[1][1]), color = colors[2])
ax2.plot(k, np.imag(w_conv_GC_mub[1][1]), linestyle = '--', color = colors[2])

for ax in [ax1, ax2]:
    ax.grid()
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=1e-5)
    ax.set_ylim(-0.001, 0.001)

labels = []

labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='--',
                        label='Convergence')] + labels
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Divergence')] + labels #+ [mlines.Line2D([], [], color='darkviolet', marker='None', linestyle='None',
labels = [mlines.Line2D([], [], color=colors[1], marker='None', linestyle='None',
                        label=r'GC: $\mu_b = {}$'.format(mub_tot[0]))]+labels
labels = [mlines.Line2D([], [], color=colors[2], marker='None', linestyle='None',
                        label=r'GC: $\mu_b = {}$'.format(mub_tot[1]))]+labels
labels = [mlines.Line2D([], [], color=colors[0], marker='None', linestyle='None',
                        label=r'VP')]+labels
                        # label='VP')]
fig.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig.supxlabel(r'$k$ (1/m)', y = 0.05)
ax1.set_ylabel(r'$\text{Re}\{\omega_-\}$ (s$^{-1}$)', fontsize = 13)
fig.align_ylabels()
ax2.set_ylabel(r'$\text{Im}\{\omega_-\}$ (s$^{-1}$)', fontsize = 13)
plt.savefig('plastic_regime_minus.png')
# %%
