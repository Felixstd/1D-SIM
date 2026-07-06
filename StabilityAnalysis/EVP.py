import numpy as np
import matplotlib.pyplot as plt
import VP_nonlocal_stab as omega
import matplotlib as mpl
from numba import jit
from matplotlib.ticker import SymmetricalLogLocator

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


dt = 15
e = 2 
h0 = 1
A0 = 1
C_s = 20
Pstarstar = 27.5e3*np.exp(-C_s*(1-A0))
rhoi = 900
T = 0.36*dt
alpha = 1/((1+e**(-2))*T)
# lamda_p = (np.sqrt(1+e**(-2))-1)/2
lamda_div =(np.sqrt(1+e**(-2))-1)/2
lamda_conv = (-np.sqrt(1+e**(-2))-1)/2


def omega_EVP(lamda, k_tot):
    k_tot = np.asarray(k_tot, dtype=np.float64)
    n = len(k_tot)

    gamma = lamda * Pstarstar * alpha * (1 - C_s * A0) / rhoi
    d = -gamma * k_tot**2          # constant term for each k, coeffs = [1, alpha, 0, d]

    # build a stack of companion matrices, shape (n, 3, 3)
    C = np.zeros((n, 3, 3), dtype=np.float64)
    C[:, 1, 0] = 1.0
    C[:, 2, 1] = 1.0
    C[:, 0, 0] = -alpha           # -a1/a0
    C[:, 0, 1] = 0.0              # -a2/a0
    C[:, 0, 2] = -d               # -a3/a0

    roots = np.linalg.eigvals(C)   # shape (n, 3), complex

    omega_1 = roots[:, 0]
    omega_2 = roots[:, 1]
    omega_3 = roots[:, 2]
    return omega_1, omega_2, omega_3


# dmin = 2e-9
# alpha = np.sqrt(1+e**(-2))/dmin

k_tot = np.linspace(1e-7,10,1000000)

omega_VP,omega_conv_VP = omega.VPs_growth(k_tot)

omega_div_1, omega_div_2, omega_div_3 = omega_EVP(lamda_div, k_tot)
omega_conv_1, omega_conv_2, omega_conv_3 = omega_EVP(lamda_conv, k_tot)

omega_nVP_div, omega_nVP_conv = omega.VP_nonlocal(k_tot, l = 10e3)

dark2_color = plt.get_cmap('Set1')


fig, ((ax1,ax2)) = plt.subplots(2, 1, sharex = True, figsize = (4, 6), constrained_layout =True)

axs = [ax1,ax2]


ax1.set_title('Divergence')
ax1.plot(k_tot, omega_VP, color = 'royalblue', label = 'VP')
ax1.plot(k_tot, omega_nVP_div, color = omega.original_cmap(omega.norm(10e3)), label = 'nVP')
ax1.plot(k_tot, omega_div_1, color = 'gray', label = 'EVP-1')
ax1.plot(k_tot, omega_div_2, color = 'forestgreen', label = 'EVP-2')
ax1.plot(k_tot, omega_div_3, color = 'crimson', label = 'EVP-3')
ax1.set_yscale('symlog',linthresh = 1e-3)
ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(omega.every_two_symlog))
ax1.yaxis.set_minor_locator(SymmetricalLogLocator(
        base=10, linthresh=1e-3, subs=np.arange(1, 10)
    ))


ax2.set_title('Convergence')
ax2.plot(k_tot, omega_conv_1, color = 'gray')
ax2.plot(k_tot, omega_conv_2, color = 'forestgreen')
ax2.plot(k_tot, omega_conv_3, color = 'crimson')
ax2.plot(k_tot, omega_conv_VP, color = 'goldenrod')
ax2.plot(k_tot, omega_nVP_conv,color = omega.original_cmap(omega.norm(10e3)))

labels = ['(a)','(b)']
for i, ax in enumerate(axs):
    ax.grid(alpha = 0.6)
    ax.set_xscale('log')
    ax.grid(alpha = 0.5)
    ax.text(0.039, 0.955, labels[i], transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top', horizontalalignment='left')
    
leg = ax1.legend(
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),   # just outside ax1
    frameon=False,
    handlelength=0,
    handletextpad=0
)
for handle, text in zip(leg.legend_handles, leg.get_texts()):
    text.set_color(handle.get_color())
    
    
# plt.xscale('log')
# plt.yscale('symlog',linthresh = 1e-3)
fig.supylabel(r'Re{$\omega$} (s$^{-1}$)')
fig.supxlabel(r'$k$ (m)')
plt.savefig('EVP_omega.png')