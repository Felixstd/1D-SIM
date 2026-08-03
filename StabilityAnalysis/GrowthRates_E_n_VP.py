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
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


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

def plot_paper(k, l_nl, 
               omega_vp_div, omega_sVP_viscous, 
               omega_EVP_viscous, omega_EVP_div, 
               omega_EVP_conv):
    """
    This is a function to produce the papers figure for the 
    growth rates of the (n)VP and (n)EVP models. 

    Args:
        k (_type_): _description_
        l_nl (_type_): _description_
        omega_vp_div (_type_): _description_
        omega_sVP_viscous (_type_): _description_
        omega_EVP_viscous (_type_): _description_
        omega_EVP_div (_type_): _description_
        omega_EVP_conv (_type_): _description_
    """
    
    plt.rcParams.update({
                "font.size": 10,
                "axes.titlesize": 10,
                "axes.labelsize": 10,
                "legend.fontsize": 10,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
            })
    
    fig, (ax1, ax2) = plt.subplots(
            2, 1,
            figsize=(4.7,5.7),
            sharex = True,
            constrained_layout = True
        )

    axs = [ax1,ax2]

    #----------------------------------------------------------
    # - VP Rheology 
    #----------------------------------------------------------
    
    #-------------------------
    # - Viscous
    #-------------------------
    # ax1.plot(k, omega_sVP_viscous, color='k', lw=1.5, label='VP', alpha = 0.8)
    
    ax1.plot(k, omega_vp_div, color="k", lw=1.5, alpha = 0.8,label='VP')
    
    # ax1.axvspan(1e-8, np.sqrt(4*rho*Pstarstar*h0**2*(1+C_s*A0)/nu_max**2), alpha=0.3, color='darkgreen')
    # ax1.axvspan(np.sqrt(4*rho*Pstarstar*h0**2*(1+C_s*A0)/nu_max**2), 1,alpha=0.3, color='crimson')
    
    #----------------------------------------------------------
    # - nVP Rheology 
    #----------------------------------------------------------
    for i, l in enumerate(l_nl):
        #-------------------------
        # - Viscous
        #-------------------------
        # ax1.plot(k, VP_nonlocal_viscous(k, l)[0], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
        
        #-------------------------
        # - Plastic
        #-------------------------
        ax1.plot(k, VP_nonlocal(k, l)[0], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
        ax2.plot(k, VP_nonlocal(k, l)[1], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
    
    #----------------------------------------------------------
    # - EVP Rheology 
    #----------------------------------------------------------
    
    #-------------------------
    # - Viscous
    #-------------------------
    # for omega_EVP in omega_EVP_viscous:
    #     # print(omega_)
    #     ax1.plot(k, omega_EVP, color = 'k', linestyle = '--',label = 'EVP')
    
    #-------------------------
    # - Plastic
    #-------------------------
    for i, l in enumerate(l_nl):  
        count = 0
        for omega in omega_EVP_div:
            if count == 0:
                ax1.plot(k, omega, color='k', linestyle= '--', label = 'EVP')
            else:
                ax1.plot(k, omega, color='k', linestyle= '--')
                
            count+=1
            
        for omega in omega_EVP_conv:
            idx_pos = np.where(omega.real>0)[0]
            if len(idx_pos)>0:
                ax2.plot(k, omega, color='k' , linestyle= '--')
                
    #----------------------------------------------------------
    # - nEVP Rheology 
    #----------------------------------------------------------
    
    for i, l in enumerate(l_nl):
        
        #-------------------------
        # - Viscous
        #-------------------------
        # omega_nEVP_visc = omega_nEVP_viscous(k, l)
        # for omega in omega_nEVP_visc:
        #     ax1.plot(k, omega, linestyle = '--', color=okabe_ito[i], lw=1.2)
            
        #-------------------------
        # - Plastic
        #-------------------------
        omega_nEVP_div, omega_nEVP_conv = omega_nEVP_plastic(k, l)
        
        for omega in omega_nEVP_div:
            ax1.plot(k, omega, linestyle = '--', color=okabe_ito[i], lw=1.2)
            
        for omega in omega_nEVP_conv:
            idx_pos = np.where(omega.real>0)[0]
            if len(idx_pos)>0:
                ax2.plot(k, omega, linestyle = '--', color=okabe_ito[i], lw=1.2)
                
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
    fig.supylabel(r'$\text{Re}\{\omega\}$ (s$^{-1}$)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    # ax1.set_ylim(-1, 1e-6)
    # ax1.set_title('Viscous')
    ax1.set_title('Divergence')
    ax2.set_title('Convergence')
    # ax3.set_yscale('log')
    fig.align_ylabels()

        
    entries = [
    ('VP',                'k',  '-'),
    ('EVP',               'k',  '--'),
    (r'$l_c = 10\,\mathrm{km}$',  okabe_ito[0], None),
    (r'$l_c = 50\,\mathrm{km}$',  okabe_ito[1], None),
    (r'$l_c = 100\,\mathrm{km}$', okabe_ito[2], None),
    ]

    # Position in axes-fraction coords (tune x0, y0, dy to taste)
    x_text  = 1.05   # left edge of text, just outside axes
    x_line0 = 1.20   # start of handle line
    x_line1 = 1.3   # end of handle line
    y0, dy  = 0.68, -0.1

    for i, (label, color, ls) in enumerate(entries):
        y = y0 + i*dy
        ax1.text(x_text, y, label, color=color, ha='left', va='center',
                transform=ax1.transAxes, fontsize=10)
        if ls is not None:
            ax1.plot([x_line0, x_line1], [y, y], color=color, lw=1.5,
                    linestyle=ls, transform=ax1.transAxes, clip_on=False)
        
        
    fig.set_constrained_layout_pads(hspace=0.08)
    fig.savefig("dispersion_VP_nonlocal_a005_plastic.png", bbox_inches="tight")
    
    
        
    fig, axs = plt.subplots(
            3, 2,
            figsize=(4.7,8.7),
            sharex = True,
            constrained_layout = True
        )

    # axs = [ax1,ax2,ax3]

    #----------------------------------------------------------
    # - EVP Rheology 
    #----------------------------------------------------------
    
    #-------------------------
    # - Viscous
    #-------------------------
    for omega_EVP in [omega_EVP_viscous[-1]]:
        axs[0,0].plot(k, omega_EVP.real, color = 'k', linestyle = '--',label = 'EVP')
        axs[0,1].plot(k, omega_EVP.imag, color = 'k', linestyle = '--',label = 'EVP')
    
    #-------------------------
    # - Plastic
    #-------------------------
    for i, l in enumerate(l_nl):  
        count = 0
        for omega in omega_EVP_div:
            if count == 0:
                axs[1,0].plot(k, omega, color='k', linestyle= '--', label = 'EVP')
            else:
                axs[1,0].plot(k, omega, color='k', linestyle= '--')
                
            count+=1
            
        for omega in omega_EVP_conv:
            idx_pos = np.where(omega.real>0)[0]
            if len(idx_pos)>0:
                axs[2,0].plot(k, omega, color='k' , linestyle= '--')
                
    #----------------------------------------------------------
    # - nEVP Rheology 
    #----------------------------------------------------------
    
    for i, l in enumerate(l_nl):
        
        #-------------------------
        # - Viscous
        #-------------------------
        # omega_nEVP_visc = omega_nEVP_viscous(k, l)
        # for omega in omega_nEVP_visc:
        #     axs[0,1].plot(k, omega.imag, linestyle = '--', color=okabe_ito[i], lw=1.2)
        #     axs[0,0].plot(k, omega.real, linestyle = '--', color=okabe_ito[i], lw=1.2)
        #-------------------------
        # - Plastic
        #-------------------------
        omega_nEVP_div, omega_nEVP_conv = omega_nEVP_plastic(k, l)
        
        for omega in omega_nEVP_div:
            axs[1,0].plot(k, omega, linestyle = '--', color=okabe_ito[i], lw=1.2)
            
        for omega in omega_nEVP_conv:
            idx_pos = np.where(omega.real>0)[0]
            if len(idx_pos)>0:
                axs[2,0].plot(k, omega, linestyle = '--', color=okabe_ito[i], lw=1.2)
                
    # for ax in axs:
    #     ax.grid(True, alpha=0.5, lw=0.5)
    #     ax.tick_params(length=3, width=0.8)
        
    labels = ['(a)','(b)','(c)','(d)','(e)','(f)']
    print(axs.flatten())
    for i, ax in enumerate(axs.flatten()):
            ax.grid(alpha = 0.5)
            ax.text(0.039, 0.955, labels[i], transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top', horizontalalignment='left')
            ax.set_xscale('log')
            ax.tick_params(length=3, width=0.8)

    # fig.canvas.draw()  # force layout computation
    ax_pos = ax2.get_position()
    center = (ax_pos.x0 + ax_pos.x1) / 2
    fig.supxlabel(r'$k$ (m$^{-1}$)', x=center-0.015)
    axs[1,0].set_ylabel(r'$\text{Re}\{\omega\}$ (s$^{-1}$)')
    axs[1,0].set_yscale('log')
    axs[0,0].set_yscale('symlog', linthresh = 1e-6)
    axs[0,1].set_yscale('symlog', linthresh = 1e-6)
    # ax1.set_ylim(-1, 1e-8)
    # ax1.set_title('Viscous')
    # ax2.set_title('Plastic-Divergence')
    # ax3.set_title('Plastic-Convergence')
    # ax3.set_yscale('log')
    fig.align_ylabels()

        
    # entries = [
    # ('VP',                'k',  '-'),
    # (r'$l_c = 10\,\mathrm{km}$',  okabe_ito[0], None),
    # (r'$l_c = 50\,\mathrm{km}$',  okabe_ito[1], None),
    # (r'$l_c = 100\,\mathrm{km}$', okabe_ito[2], None),
    # ('EVP',               'k',  '--'),
    # ]

    # # Position in axes-fraction coords (tune x0, y0, dy to taste)
    # x_text  = 1.05   # left edge of text, just outside axes
    # x_line0 = 1.20   # start of handle line
    # x_line1 = 1.3   # end of handle line
    # y0, dy  = 0.68, -0.1

    # for i, (label, color, ls) in enumerate(entries):
    #     y = y0 + i*dy
    #     axs[1,0].text(x_text, y, label, color=color, ha='left', va='center',
    #             transform=ax2.transAxes, fontsize=10)
    #     if ls is not None:
    #          axs[1,0].plot([x_line0, x_line1], [y, y], color=color, lw=1.5,
    #                 linestyle=ls, transform=ax2.transAxes, clip_on=False)
        
        
    # fig.set_constrained_layout_pads(hspace=0.08)
    fig.savefig("dispersion_EVP_nonlocal_a005.png", bbox_inches="tight")

def plot_presentation(k, l_nl, 
               omega_vp_div, omega_sVP_viscous, 
               omega_EVP_viscous, omega_EVP_div, 
               omega_EVP_conv):
    
    plt.rcParams.update({
                "font.size": 12,
                "axes.titlesize": 12,
                "axes.labelsize": 12,
                "legend.fontsize": 12,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10,
            })
    
    #-----------------------------------------------
    # viscous regime 
    #-----------------------------------------------
    
    fig1 = plt.figure(1)
    
    plt.close(1)
    
    #-----------------------------------------------
    # Plastic regime 
    #-----------------------------------------------
    
    fig2 = plt.figure(2, figsize = (6,3))
    ax2 = plt.axes()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.plot(k, omega_vp_div, color="k", lw=1.5, alpha = 0.8,label='VP')
    for i, l in enumerate(l_nl):
        ax2.plot(k, VP_nonlocal(k, l)[0], color=okabe_ito[i], lw=1.2,label=r'$l_c = {:.0f}\,\mathrm{{km}}$'.format(l/1e3))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    # ax2.spines['bottom'].set_alpha(0.7)
    # ax2.spines['left'].set_alpha(0.7)
    # ax2.tick_params(axis='both', alpha = 0.7)
    leg = ax2.legend(
        loc='upper left',
        bbox_to_anchor=(0.01, 1),   # just outside ax1
        frameon=False,
        handlelength=0,
        handletextpad=0
    )
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        text.set_color(handle.get_color())
        
    ax2.set_xlabel(r'$k$ (m$^{-1}$)')
    
        
    ax2.set_ylabel('Re$\{\omega_+\}$\n(s$^{-1}$)', rotation = 0, 
                   ha = 'left')
    ax2.yaxis.set_label_coords(-0.22, 0.8) 
    fig2.savefig('growth_plastic_presentation.png')
    plt.close(2)
    
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
 
    k = np.logspace(-8, 0, 3000) 
    l_nl = np.array([10,50,100])*1e3

    #--------------------------------------------------------
    #Computing the growth rates 
    omega_vp_div, omega_vp_conv = VPs_growth(k)
    omega_nonlocal_div,_ = VP_nonlocal(k)
    omega_nVP_viscous,bnloverk = VP_nonlocal_viscous(k)
    omega_sVP_viscous,bloverk = VP_local_viscous(k)
    
    omega_EVP_div = omega_EVP_plastic(lamda_div, k)
    omega_EVP_conv = omega_EVP_plastic(lamda_conv, k)
    omega_EVP_visc = omega_EVP_viscous(k)
    
    if args.mode == 'paper':
        plot_paper(k, l_nl, 
                    omega_vp_div, omega_sVP_viscous, 
                    omega_EVP_visc, omega_EVP_div, 
                    omega_EVP_conv)
    
    elif args.mode == 'presentation':
        plot_presentation(k, l_nl, 
                    omega_vp_div, omega_sVP_viscous, 
                    omega_EVP_viscous, omega_EVP_div, 
                    omega_EVP_conv)
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

