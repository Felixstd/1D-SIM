import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as colors

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def plot_mechanical_energy(dates, datadict_energy, Parameters, figdir, expno):
    
    Plat_tot, Pfric_s_tot, Pfric_R_tot, Ph_tot, Pw_tot, Ppot_tot = datadict_energy.values()
    
    
    Nx = np.shape(Plat_tot)[1]
    X = np.arange(0, Parameters.Nx)*Parameters.dx/1e3
    dates = dates*Parameters.dt/(60*60)
    
    #---- Setting the colorbar ----#
    norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
    cmap = plt.cm.viridis
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    fig_all_Plat, ax_all_Plat = plt.subplots(figsize = (5, 4))

    fig_all_Pfric_s, ax_all_Pfric_s = plt.subplots(figsize = (5, 4))
    
    fig_all_Pfric_R, ax_all_Pfric_R = plt.subplots(figsize = (5, 4))
    # fig_all_Pfric_s = plt.figure(2, figsize = (5, 4))
    # ax_all_Pfric_s = fig_all_Pfric_s.axes()
    
    # fig_all_Pfric_R = plt.figure(3, figsize = (5, 4))
    # ax_all_Pfric_R = fig_all_Pfric_R.axes()
    
    
    for k, date in enumerate(dates):
        
        print('Plotting Energy (Time in space): ', date)
        
        Plat, Pfric_s, Pfric_R, Ph, Pw, Ppot = Plat_tot[k], \
                                                Pfric_s_tot[k], \
                                                Pfric_R_tot[k], \
                                                Ph_tot[k], \
                                                Pw_tot[k], \
                                                Ppot_tot[k]
        
        color = cmap(k/10)
        ax_all_Plat.plot(X[1:-1], Plat[:-1], label = '{}'.format(date), color = color)
        ax_all_Pfric_s.plot(X[1:-1], Pfric_s[:-1], label = '{}'.format(date), color = color)
        ax_all_Pfric_R.plot(X[1:-1], Pfric_R[:-1], label = '{}'.format(date), color = color)
    
        # ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax_all_Plat.grid()
    ax_all_Plat.set_xlabel('x (km)')
    ax_all_Plat.set_ylabel(r'$P_{lat}$ (W/m$^2$)')
    fig_all_Plat.colorbar(sm, ax = ax_all_Plat, label = 'Time (hr)')
    fig_all_Plat.savefig(figdir+'time_plot_Plat_{}.png'.format(expno))
    
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_Pfric_R.grid()
    ax_all_Pfric_R.set_xlabel(r'$x$ (km)')
    ax_all_Pfric_R.set_ylabel(r'$P_{fric}^R$ (W/m$^2$)')
    # ax_all_h.set_aspect('equal', 'datalim')
    fig_all_Pfric_R.colorbar(sm, ax = ax_all_Pfric_R, label = 'Time (hr)')
    fig_all_Pfric_R.savefig(figdir+'time_plot_PfricR_{}.png'.format(expno))
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_Pfric_s.grid()
    ax_all_Pfric_s.set_xlabel('$x$ (km)')
    ax_all_Pfric_s.set_ylabel(r'$P_{fric}^s$ (W/m$^2$)')
    sm = cm.ScalarMappable(cmap = cmap,
                           norm = norm)
    fig_all_Pfric_s.colorbar(sm, ax = ax_all_Pfric_s, 
                       label = 'Time (hr)')
    fig_all_Pfric_s.savefig(figdir+'time_plot_PfricS_{}.png'.format(expno))
    plt.close()
    
    
def plot_energy(datadict_energy_Tavg,num_plas, num_visc,number_mixed_viszeta, number_mixed_viseta, Parameters, figdir):
    
    
    Plat_Tavg, Pfric_S_Tavg, Pfric_R_Tavg,Pfric_S_visc_Tavg, Pfric_R_visc_Tavg, \
        Pfric_S_plas_Tavg, Pfric_R_plas_Tavg, P_fric_R_pZ_vE_Tavg, P_fric_R_vZ_pE_Tavg , \
            Ph_Tavg, Pw_Tavg, Ppot_Tavg = datadict_energy_Tavg.values()

    
    
    
    energy = [Plat_Tavg, Pfric_S_Tavg, Pfric_R_Tavg, Ph_Tavg, Pw_Tavg, Ppot_Tavg] 
    labels = [r'$P_{lat}$ (W/m$^2$)', r'$P_{f}^S$ (W/m$^2$)', r'$P_{f}^R$ (W/m$^2$)', r'$P_{h}$ (W/m$^2$)', r'$P_{w}$ (W/m$^2$)', r'$P_{pot}$ (W/m$^2$)']
    time = np.arange(1, Parameters.nstep)*Parameters.dt/(60*60)
    
    sum_energy = np.sum(energy, axis = 0)
    
    for i, E in enumerate(energy):
        print('Plotting energy (Time average): ', i)
        plt.figure(figsize = (5, 4))
        ax1 = plt.axes()
        
        ax1.plot(time, E[0], color = 'k')
        ax1.grid()
        ax1.set_xlabel(r'Time (hr)')
        ax1.set_ylabel(labels[i])
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.savefig(figdir+'energy_{}.png'.format(i))

    fig = plt.figure()
    # plt.plot(time, Ppot_Tavg[0], label = r'$P_{pot}$ ')
    # plt.plot(time, Pfric_R_Tavg[0], label = r'$P_{fric}^R$')
    plt.plot(time, Pfric_R_visc_Tavg[0], label = r'$P_{fric, vis}^R$ ')
    plt.plot(time, Pfric_R_plas_Tavg[0], label = r'$P_{fric, plas}^R$ ')

    # plt.plot(time, Pfric_R_plas_Tavg[0] + Pfric_R_visc_Tavg[0])
    plt.plot(time, P_fric_R_pZ_vE_Tavg[0], label = r'$P_{fric, plas \; \zeta}^R$')
    plt.plot(time, P_fric_R_vZ_pE_Tavg[0], label = r'$P_{fric, plas \; \eta}^R$')
    fig.legend(loc = 'outside upper center', bbox_to_anchor = (1.1, 0.9))
    plt.xlabel(r'Time (hr)')
    plt.ylabel(r'P (W/m$^2$)')
    plt.savefig(figdir+'energy_internal_stresses.png')
    
    fig = plt.figure(figsize = (7, 4))
    ax = plt.axes()
    plt.plot(time, Pfric_R_visc_Tavg[0], label = r'$P_{fric, vis}^R$ ')
    plt.plot(time, Pfric_R_plas_Tavg[0], label = r'$P_{fric, plas}^R$ ')
    plt.plot(time, Ppot_Tavg[0], label = r'$P_{pot}^R$ ')
    plt.plot(time, Ppot_Tavg[0] + Pfric_R_Tavg[0],  label = r'$P_{pot} + P_{fric}^R$')
    ax.set_yscale('symlog', linthresh=1e-1)
    fig.legend(loc = 'outside upper center', bbox_to_anchor = (1, 0.9))
    plt.xlabel(r'Time (hr)')
    plt.ylabel(r'P (W/m$^2$)')
    plt.savefig(figdir+'energy_internal_stresses_tot.png')
    
    fig = plt.figure()
    plt.plot(time, Pfric_R_plas_Tavg[0], label = r'$P_{fric, plas}^R$ ')
    plt.plot(time, Pfric_R_Tavg[0], label = r'$P_{fric}^R$')    
    plt.plot(time, Pfric_R_visc_Tavg[0], label = r'$P_{fric, vis}^R$ ')
    plt.plot(time, Pfric_R_visc_Tavg[0]+Pfric_R_plas_Tavg[0], linestyle = '--')
    fig.legend(loc = 'outside upper center', bbox_to_anchor = (1.1, 0.9))
    plt.xlabel(r'Time (hr)')
    plt.ylabel(r'P (W/m$^2$)')
    plt.savefig(figdir+'energy_fric_visc_plas.png')
    
    
    plt.figure()
    ax = plt.axes()
    plt.plot(time, sum_energy[0], label = r'$\sum P$ ')
    # plt.plot(time, Ph_Tavg[0], label = r'$P_{h}$')
    plt.plot(time, Ppot_Tavg[0], label = r'$P_{pot}$')
    # plt.plot(time, Pw_Tavg[0], label = r'$P_{w}$')
    plt.plot(time, np.sum(energy[0:3], axis =0)[0], label = r'$P_i$')
    plt.yscale('symlog', linthresh = 1)
    plt.legend()
    plt.xlabel(r'Time (hr)')
    plt.ylabel(r'Power' +'\n(W/m$^2$)', rotation = 0, ha = 'left')
    ax.yaxis.set_label_coords(-0.34,0.82)
    plt.savefig(figdir+'energy_tot.png')
    
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(time, num_plas, label = r'$N_{plas}$')
    # plt.plot(time, num_visc+num_plas)
    plt.plot(time, num_visc,  label = r'$N_{visc}$')
    # plt.plot(time,number_mixed_viszeta, label = r'$N_{visc \; \zeta, plas \; \eta}$')
    # plt.plot(time, number_mixed_viseta, label = r'$N_{visc \; \eta, plas \; \zeta }$')
    plt.legend(loc = 'best')
    plt.xlabel(r'Time (hr)')
    plt.ylabel('Number'+'\nof grid cells', rotation = 0, ha = 'left')
    ax.yaxis.set_label_coords(-0.37,0.78)
    # plt.grid()
    plt.savefig(figdir+'number_regime.png')
    
  