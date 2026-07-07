"""
This is a script which plots all of the relevant growth rates. 

It includes:
-> Viscous-Plastic Rheology 
-> nonlocal Viscous-Plastic Rheology
-> Elastic-Viscous-Plastic Rheology

The code for figures is also included. 
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from growth_rates import *
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.ticker import FixedLocator
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib as mpl
from matplotlib.colors import ListedColormap


okabe_ito = ['#CC79A7','#0072B2','#009E73','#E69F00','#D55E00']

# mpl.colormaps.register(ListedColormap(okabe_ito, name='okabe_ito'))

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

# Colors for plot
# start, end = 0.4, 1
start,end = 0,0.5
original_cmap = plt.get_cmap('Blues_r')
# 2. Sample 256 colors from that specific range
colors = original_cmap(np.linspace(start, end, 5))

original_cmap = plt.get_cmap('magma')
norm = LogNorm(vmin=10e3*0.5, vmax=250e3*1.2)

plot_presentation = 0
plot_paper = 1

styles = {
    "presentation": plot_presentation,
    "paper": plot_paper,
}


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate paper- or presentation-style figures."
    )

    parser.add_argument(
        "--mode",
        choices=list(styles.keys()),
        default="paper",
        help="Styling target (default: paper).",
    )
    return parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
 
    # if args.list or args.figure is None:
    #     print("Available styles:")
    #     for name in styles:
    #         print(f"  - {name}")
    #     if args.figure is None and not args.list:
    #         parser.print_usage()
    #         return 1
    #     return 0
    print(args.mode)
    
    # k = np.linspace(1e-8,1, 10000000)
    k = np.logspace(-8, 0, 3000) 
    l_nl = np.array([10,50,100,150,200])*1e3

    #--------------------------------------------------------
    #Computing the growth rates 
    omega_vp_div, omega_vp_conv = VPs_growth(k)
    omega_nonlocal_div,_ = VP_nonlocal(k)
    omega_nVP_viscous,bnloverk = VP_nonlocal_viscous(k)
    omega_sVP_viscous,bloverk = VP_local_viscous(k)
    
    omega_EVP_div = omega_EVP(lamda_div, k)
    omega_EVP_conv = omega_EVP(lamda_conv, k)
    omega_EVP_viscous,_ = EVP_viscous(k)
    
    
    if args.mode == 'paper':
        fig, (ax1, ax2, ax3) = plt.subplots(
            3, 1,
            figsize=(4.7,8.7),
            sharex = True,
            constrained_layout = True
        )

        axs = [ax1,ax2,ax3]

        plt.rcParams.update({
                    "font.size": 10,
                    "axes.titlesize": 10,
                    "axes.labelsize": 10,
                    "legend.fontsize": 10,
                    "xtick.labelsize": 8,
                    "ytick.labelsize": 8,
                })

        ax2.plot(k, omega_vp_div, color="#999999", lw=1.5, alpha = 0.8,label='VP')
        ax1.plot(k, omega_sVP_viscous, color='#999999', lw=1.5, label='VP', alpha = 0.8)
        ax1.plot(k, omega_EVP_viscous, color = '#4d4d1a', label = 'EVP')
        ax1.axvspan(1e-8, np.sqrt(4*rho*Pstarstar*h0**2*(1+C_s*A0)/nu_max**2), alpha=0.3, color='darkgreen')
        ax1.axvspan(np.sqrt(4*rho*Pstarstar*h0**2*(1+C_s*A0)/nu_max**2), 1,alpha=0.3, color='crimson')
        print(np.sqrt(4*rho*Pstarstar*h0**2*(1+C_s*A0)/nu_max**2))
        for i, l in enumerate(l_nl):
            ax2.plot(k, VP_nonlocal(k, l)[0], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
            # ax1.plot(k, VP_nonlocal_viscous(k, l)[0], color=original_cmap(norm(l)), lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
            ax1.plot(k, VP_nonlocal_viscous(k, l)[0], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))

        linestyles = ['-', '-.', '--']

        count = 0
        for omega, ls in zip(omega_EVP_div, linestyles):
            if count == 0:
                ax2.plot(k, omega, color='#4d4d1a', linestyle=ls, label = 'EVP')
            else:
                ax2.plot(k, omega, color='#4d4d1a', linestyle=ls)
                
            count+=1
        for omega, ls in zip(omega_EVP_conv, linestyles):
            ax3.plot(k, omega, color='#4d4d1a', linestyle=ls)

        for ax in axs:
            ax.grid(True, alpha=0.5, lw=0.5)
            ax.tick_params(length=3, width=0.8)
            
        labels = ['(a)','(b)','(c)']
        for i, ax in enumerate(axs):
                ax.grid(alpha = 0.5)
                ax.text(0.039, 0.955, labels[i], transform=ax.transAxes,
                    fontsize=10,
                    verticalalignment='top', horizontalalignment='left')

        # fig.canvas.draw()  # force layout computation
        ax_pos = ax2.get_position()
        center = (ax_pos.x0 + ax_pos.x1) / 2
        fig.supxlabel(r'$k$ (m$^{-1}$)', x=center-0.015)
        # fig.supxlabel(r'$k$ (m$^{-1}$)',x = 0.465)
        ax1.set_ylabel(r'$\text{Re}\{\omega_+\}$ (s$^{-1}$)')
        ax2.set_ylabel(r'$\text{Re}\{\omega\}$ (s$^{-1}$)')
        ax3.set_ylabel(r'$\text{Re}\{\omega\}$ (s$^{-1}$)')
        # ax2.legend(fontsize=8, frameon=False)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax1.set_xscale('log')
        ax1.set_title('Viscous')
        ax2.set_title('Plastic-Divergence')
        ax3.set_title('Plastic-Convergence')
        fig.align_ylabels()
        # ax1.set_yscale('symlog', linthresh=1e-6)
        # ax1.yaxis.set_minor_locator(SymmetricalLogLocator(
        #     base=10, linthresh=1e-6, subs=np.arange(1, 10)
        # ))
        # ax1.set_ylim(-1e-2, 1e-6)
        # ax1.yaxis.set_major_locator(FixedLocator(ticks_neg))
        leg = ax2.legend(
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),   # just outside ax1
            frameon=False,
            handlelength=0,
            handletextpad=0
        )
        for handle, text in zip(leg.legend_handles, leg.get_texts()):
            text.set_color(handle.get_color())
        # ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        # ax1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two_symlog))


        # ax2.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10))
        # ax2.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two))
        fig.set_constrained_layout_pads(hspace=0.08)
        fig.savefig("dispersion_VP_nonlocal.png", bbox_inches="tight")


        plt.close(fig)
    return 0




if __name__ == "__main__":
    sys.exit(main())
    
    
# #------------------------------------------
# # Figures
# if Figure:
    

#     plt.figure()
#     plt.plot(k,bnloverk, color = 'r')
#     plt.plot(k,bloverk, color = 'b')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.savefig('bnloverk.png')


#     #-----------------------------------------
#     #Fig 2 in paper 



#     #---------------------------------------------------------------------
#     # figure for presentation

#     #----------------------------------------------------------
#     if Presentation:
#         # Viscous
#         fig1, ax1 = plt.subplots(1,1,figsize = (5,3))
#         ax1.spines['top'].set_visible(False)
#         ax1.spines['right'].set_visible(False)
#         ax1.plot(k, omega_sVP_viscous, color='goldenrod', lw=1.5, label='VP', alpha = 0.8)

#         for i, l in enumerate(l_nl):
#             ax1.plot(k, VP_nonlocal_viscous(k, l)[0], color=colors[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
            
#         # ax1.axvspan(1e-7, np.sqrt(4*rho*Pstarstar*h0**2*(1+C*A0)/nu_max**2), alpha=0.3, color='darkgreen')
#         # ax1.axvspan(np.sqrt(4*rho*Pstarstar*h0**2*(1+C*A0)/nu_max**2), 1,alpha=0.3, color='crimson')
#         ax1.set_yscale('symlog', linthresh=1e-6)
#         ax1.set_xscale('log')

#         ax1.yaxis.set_minor_locator(SymmetricalLogLocator(
#             base=10, linthresh=1e-6, subs=np.arange(1, 10)
#         ))
#         ax1.set_ylim(-1e-2, 1e-6)
#         ax1.set_xlabel(r'$k$ (m$^{-1}$)', fontsize=14)
#         ax1.annotate(r"$\text{Re}\{\omega_+\}$"+'\n'+"(s$^{-1}$)", xy=(-0.25, 0.8), xytext=(0, 3),
#                         xycoords="axes fraction", textcoords="offset points",
#                         ha="left", fontsize=14)

#         leg = ax1.legend(
#             loc='center left',
#             bbox_to_anchor=(0.73, 0.48),   # just outside ax1
#             frameon=False,
#             handlelength=0,
#             handletextpad=0,
#             fontsize = 13
#         )
#         for handle, text in zip(leg.legend_handles, leg.get_texts()):
#             text.set_color(handle.get_color())
#         ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

#         ax1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two_symlog))
#         fig1.savefig('viscous_regime_omega_plus.png')



#         #----------------------------------------------------------
#         # Plastic
#         fig2, ax2 = plt.subplots(1,1,figsize = (5,3))
#         ax2.spines['top'].set_visible(False)
#         ax2.spines['right'].set_visible(False)
#         ax2.plot(k, omega_vp_div, color='goldenrod', lw=1.5, label='VP', alpha = 0.8)

#         for i, l in enumerate(l_nl):
#             ax2.plot(k, VP_nonlocal(k, l)[1], color=colors[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
            
#         # ax1.axvspan(1e-7, np.sqrt(4*rho*Pstarstar*h0**2*(1+C*A0)/nu_max**2), alpha=0.3, color='darkgreen')
#         # ax1.axvspan(np.sqrt(4*rho*Pstarstar*h0**2*(1+C*A0)/nu_max**2), 1,alpha=0.3, color='crimson')
#         ax2.set_yscale('log')
#         ax2.set_xscale('log')

#         # ax1.yaxis.set_minor_locator(SymmetricalLogLocator(
#         #     base=10, linthresh=1e-6, subs=np.arange(1, 10)
#         # ))
#         # ax1.set_ylim(-1e-2, 1e-6)
#         ax2.set_xlabel(r'$k$ (m$^{-1}$)', fontsize=14)
#         ax2.annotate(r"$\text{Re}\{\omega_+\}$"+'\n'+"(s$^{-1}$)", xy=(-0.27, 0.8), xytext=(0, 3),
#                         xycoords="axes fraction", textcoords="offset points",
#                         ha="left", fontsize=14)

#         leg = ax2.legend(
#             loc='center left',
#             bbox_to_anchor=(0.01, 0.6),   # just outside ax1
#             frameon=False,
#             handlelength=0,
#             handletextpad=0,
#             fontsize = 13
#         )

#         for handle, text in zip(leg.legend_handles, leg.get_texts()):
#             text.set_color(handle.get_color())
#         ax2.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10))
#         ax2.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(every_two))

#         fig2.savefig('plastic_regime_omega_plus.png')

