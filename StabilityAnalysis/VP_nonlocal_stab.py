import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.ticker import FixedLocator
import matplotlib as mpl

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

# Colors for plot
start, end = 0.4, 1
original_cmap = plt.get_cmap('Blues')
# 2. Sample 256 colors from that specific range
colors = original_cmap(np.linspace(start, end, 5))

#----------------------------------------
# Defining constants 
A0 = 1
h0 = 2
C = 20
Pstar = 2.75e4
P0 = Pstar
Pstarstar = Pstar*np.exp(-C*(1-A0))
rho = 900
delta_min = 2e-4
e = 2
alpha = np.sqrt(1+e**(-2))
l = 40e3
lamda_div = (alpha-1)/2
alpha = 2/(l**2)*(P0-Pstarstar)
nu_max = 1e9

def every_two(x, pos):
    exp = int(np.round(np.log10(abs(x))))
    if exp % 2 != 0:
        return r'$10^{{{}}}$'.format(exp)
    return ''

def every_two_symlog(x, pos):
    if x == 0:
        return '0'
    exp = int(np.round(np.log10(abs(x))))
    if exp % 2 != 0:
        sign = '-' if x < 0 else ''
        return r'${}10^{{{}}}$'.format(sign, exp)
    return ''

def VPs_growth(k):
    
    omega = np.sqrt((Pstarstar*(1+C*A0)/rho)*lamda_div*k**2)
    
    return omega

def VP_nonlocal(k, l = 40e3):
    
    omega = np.sqrt(lamda_div/(rho*h0)*((Pstarstar/l**2)*h0*(1+C*A0)+alpha*A0)*k**2/(k**2+1/l**2))
    
    return omega


def VP_nonlocal_viscous(k, l=40e3):
    beta = 4*rho*h0*((Pstarstar*h0/l**2*(1+C*A0))+alpha*A0)/nu_max**2
    arg = 1 - beta/(k**2*(1/l**2+k**2))
    omega = nu_max*k**2/(2*rho*h0)*(-1 + np.sqrt(arg.astype(complex)))
    return omega, beta/(k**2*(1/l**2+k**2))

def VP_local_viscous(k):
    arg = 1 - 4*rho*Pstarstar*h0**2*(1+C*A0)/(nu_max*k)**2
    omega = nu_max*k**2/(2*rho*h0)*(-1 + np.sqrt(arg.astype(complex)))
    return omega, 4*rho*Pstarstar*h0**2*(1+C*A0)/(nu_max*k)**2

k = np.linspace(1e-7,1, 1000000)
l_nl = np.linspace(10, 50, 5)*1e3

omega_vp_div = VPs_growth(k)
omega_nonlocal_div = VP_nonlocal(k)
omega_nVP_viscous,bnloverk = VP_nonlocal_viscous(k)
omega_sVP_viscous,bloverk = VP_local_viscous(k)

#------------------------------------------
# Figures


plt.figure()
plt.plot(k,bnloverk, color = 'r')
plt.plot(k,bloverk, color = 'b')
plt.xscale('log')
plt.yscale('log')
plt.savefig('bnloverk.png')


#-----------------------------------------
#Fig 2 in paper 

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(3.7,4.7),
    sharex = True,
    constrained_layout = True
)

axs = [ax1,ax2]

plt.rcParams.update({
			"font.size": 10,
			"axes.titlesize": 10,
			"axes.labelsize": 10,
			"legend.fontsize": 10,
			"xtick.labelsize": 8,
			"ytick.labelsize": 8,
		})

ax2.plot(k, omega_vp_div, color='goldenrod', lw=1.5, alpha = 0.7)
ax1.plot(k, omega_sVP_viscous, color='goldenrod', lw=1.5, label='VP', alpha = 0.7)

for i, l in enumerate(l_nl):
    ax2.plot(k, VP_nonlocal(k, l), color=colors[i], lw=1.2
             )
    ax1.plot(k, VP_nonlocal_viscous(k, l)[0], color=colors[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))

for ax in axs:
    ax.grid(True, alpha=0.5, lw=0.5)
    ax.tick_params(length=3, width=0.8)
    
labels = ['(a)','(b)']
for i, ax in enumerate(axs):
        ax.grid(alpha = 0.5)
        ax.text(0.035, 0.955, labels[i], transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top', horizontalalignment='left')

# fig.canvas.draw()  # force layout computation
ax_pos = ax2.get_position()
center = (ax_pos.x0 + ax_pos.x1) / 2
print(center)
fig.supxlabel(r'$k$ (m$^{-1}$)', x=center-0.015)
# fig.supxlabel(r'$k$ (m$^{-1}$)',x = 0.465)
fig.supylabel(r'$\text{Re}\{\omega_+\}$ (s$^{-1}$)',y=0.52)
# ax2.legend(fontsize=8, frameon=False)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.set_xscale('log')
ax1.set_title('Viscous')
ax2.set_title('Plastic')
ax1.set_yscale('symlog', linthresh=1e-6)
ax1.yaxis.set_minor_locator(SymmetricalLogLocator(
    base=10, linthresh=1e-6, subs=np.arange(1, 10)
))
ax1.set_ylim(-1e-2, 1e-6)
# ax1.yaxis.set_major_locator(FixedLocator(ticks_neg))
leg = ax1.legend(
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),   # just outside ax1
    frameon=False,
    handlelength=0,
    handletextpad=0
)
for handle, text in zip(leg.legend_handles, leg.get_texts()):
    text.set_color(handle.get_color())
ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two_symlog))


ax2.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10))
ax2.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two))
fig.set_constrained_layout_pads(hspace=0.08)
fig.savefig("dispersion_VP_nonlocal.png", bbox_inches="tight")