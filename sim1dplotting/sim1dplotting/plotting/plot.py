import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl
import cmocean
# import scienceplots

# plt.style.use('science')
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def normalize_stresses(sigI, sigII, p, eps = 1e-12):
    
    sigI_norm = np.full_like(sigI, np.nan)
    sigII_norm = np.full_like(sigII, np.nan)

    # Create a mask where p >= 1e-12
    # print(np.where(p >= 1e-20))
    mask = p >= 1e-12

    # Apply the mask to compute the normalized values
    sigI_norm[mask] = sigI[mask] / p[mask]
    sigII_norm[mask] = sigII[mask] / p[mask]
    
    return sigI_norm, sigII_norm

            
def histograms(ax, var, mu0, mu_infty,varname, I = False, I_0 = 0):
        

    min, max = np.min(var), np.max(var)
    bins = np.linspace(min, max, 1000)
        
        
    # fig, ax1 = plt.subplots(1, 1, figsize = (12, 6))
    
    if I: 
        var_nonmax=var[var < max] 
        max = np.max(var_nonmax)
        bins = np.linspace(min, max, 1000)
        ax.axvline(I_0, color = 'k', label = r'$I_0$')
        plt.legend()
        

    ax.hist(var.flatten(), bins, color = 'b')
    ax.set_xlim(mu0-0.05, 0.15)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel(varname)
    ax.set_ylabel('Counts')

    # plt.savefig(figdir+expno+'/'+filename+'{}.png'.format(date))
    
        
def plot_individual_time(dates, expno, data_dict, dx, dt, figdir, mu_0, mu_infty, Parameters, Dissipation = False):
    
        
    if Parameters.dissipation and Parameters.GC:
        divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, Wdissip_tot, P_tot, muI_tot = data_dict.values()
    elif Parameters.GC and Parameters.dissipation == False:
        divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, P_tot, muI_tot = data_dict.values()
    elif Parameters.GC == False and Parameters.dissipation == True:
        divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, Wdissip_tot= data_dict.values()
    else:
        divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot = data_dict.values()
    
        
    Nx = np.shape(divergence_tot)[1]
    X = np.arange(0, Nx)*dx/1e3
    
    
    fig_all_A = plt.figure(1, figsize = (5, 4))
    ax_all_A = plt.axes()
    
    fig_all_h = plt.figure(2, figsize = (5, 4))
    ax_all_h = plt.axes()
    
    fig_all_u = plt.figure(3, figsize = (5, 4))
    ax_all_u = plt.axes()
    
    fig_all_div = plt.figure(5, figsize = (5, 4))
    ax_all_div= plt.axes()
    
    # print(dates)
    
    if Parameters.GC:
        fig_all_mu = plt.figure(4)
        ax_all_mu = plt.axes()
    
    dates = dates*dt/(60*60)
    
    
    norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
    cmap = plt.cm.viridis
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        divergence, sig = divergence_tot[k], sigma_tot[k]
        sig_norm = sigma_norm_tot[k]
        h, A = h_tot[k], A_tot[k]
        u = u_tot[k]
        eta, zeta = eta_tot[k], zeta_tot[k]
        eta_max = np.max(eta)
        zeta_max = np.max(zeta)
        
        if Parameters.GC:
            muI = muI_tot[k]
            P = P_tot[k]
        
        if Dissipation:
            w_dissip = Wdissip_tot[k]
        
        
        #----- Deformation Plots -----#
        
        fig, (ax1)= plt.subplots(1, 1, sharex = True, figsize = (10, 12))
        
        
        # fig, (ax1, ax2, ax3, ax4) = plt.subplots(2, 2, sharex = True, figsize = (10,12))
        
        ax1.plot(X[1:-1], u[:-1], color = 'r')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('u (m/s)')
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_u_{}_{}.png'.format(k, expno))
        
        fig, (ax1)= plt.subplots(1, 1, sharex = True, figsize = (5, 4))
             
        ax1.plot(X[1:-1], eta[1:-1], color = 'r', label = r'$\eta$')
        ax1.plot(X[1:-1], zeta[1:-1], color = 'b', label = r'$\zeta$')
        ax1.set_xlabel('x (km)')
        ax1.set_ylabel(r'Viscosity (Ns/m)')
        # ax1.set_xlim(Nx)
        ax1.grid()
        plt.legend()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_eta_{}_{}.png'.format(k, expno))
        
        fig, (ax1)= plt.subplots(1, 1, sharex = True, figsize = (7, 3))
        
        colors_sig = []
        for i in range(1, len(eta)-1):
            eta_i = eta[i]
            zeta_i = zeta[i]
            
            if eta_i == eta_max and zeta_i == zeta_max :
                colors_sig.append('r')
            
            elif eta_i != eta_max and zeta_i == zeta_max :
                
                colors_sig.append('b')
                
            elif eta_i == eta_max and zeta_i != zeta_max :
                
                colors_sig.append('orange')
                
            else:
                colors_sig.append('g')
                
        
        ax1.scatter(divergence[1:-1], sig_norm[1:-1], c = colors_sig)
        ax1.set_xlabel(r'$\dot{\epsilon}_I$ (1/s)')
        ax1.set_ylabel(r'$\sigma_{11}/P$')
        ax1.set_xscale('symlog', linthresh = 1e-8)
        ax1.set_xlim(-1e-3, 1e-3 )
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        # ax1.set_aspect(1./ax1.get_data_ratio())        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_sigma_{}_{}.png'.format(k, expno))
        
        if Dissipation:
            fig, (ax4)= plt.subplots(1, 1, sharex = True, figsize = (10, 12))
            
            ax4.plot(X[1:-1], w_dissip[1:-1], color = 'r')
            ax4.set_xlabel('X (km)')
            ax4.set_ylabel(r'$W_\sigma$ (N(ms)$^{-1}$)')
            # ax1.set_xlim(Nx)
            ax4.grid()
            # ax1.set_aspect('equal')
            # ax4.set_aspect(1./ax1.get_data_ratio())        
            
            # plt.axis('equal')
            
            fig.savefig(figdir+'ice_wdissip_{}_{}.png'.format(k, expno))
            
        if Parameters.GC:
            fig, (ax4)= plt.subplots(1, 1, sharex = True, figsize = (5, 5))
            
            criteria_conv = (-muI-1)*P/1e9
            criteria_div = (muI-1)*P/1e9
            
            # ax4.plot(X[1:-1], criteria_conv[:-1], color = 'r')
            # ax4.plot(X[1:-1], criteria_div[:-1], color = 'g')
            # ax4.plot(X[1:-1], divergence[1:-1], color = 'k')
            print(criteria_div[:-1][divergence[1:-1] > 0])
            # ax4.scatter( divergence[1:-1][divergence[1:-1] < 0], criteria_conv[:-1][divergence[1:-1] < 0])
            ax4.scatter( divergence[1:-1][divergence[1:-1] > 0], criteria_div[:-1][divergence[1:-1] > 0])
            # ax4.plot(X[1:-1], P[:-1]/P.max(), color = 'b')
            ax4.set_xlabel('X (km)')
            ax4.set_ylabel(r'$criteria$')
            # ax4.set_xlim(150, 350)
            ax4.grid()
            # ax1.set_aspect('equal')
            # ax4.set_aspect(1./ax1.get_data_ratio())        
            
            # plt.axis('equal')
            
            fig.savefig(figdir+'ice_criteria_{}_{}.png'.format(k, expno))
        
        #--- Plot of A and h ---# 
        fig, (ax1, ax2)= plt.subplots(1, 2, sharex = True, figsize = (10, 12))
        
        ax1.plot(X[1:-1], A[1:-1], color = 'r')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('A')
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())

        ax2.plot(X[1:-1], h[1:-1], color = 'r')
        ax2.set_xlabel(r'$x$ (km)')
        ax2.set_ylabel(r'$h$ (m)')
        ax2.grid()
        # ax2.set_aspect('equal')
        ax2.set_aspect(1./ax2.get_data_ratio())
        
        color = cmap(k/10)
        
        
        ax_all_A.plot(X[1:-1], A[1:-1], label = '{}'.format(date), color = color)
        ax_all_u.plot(X[1:-1], u[:-1], label = '{}'.format(date), color = color)
        ax_all_h.plot(X[1:-1], h[1:-1], label = '{}'.format(date), color = color)
        
        ax_all_div.plot(X[1:-1], divergence[1:-1], label = '{}'.format(date), color = color)
        
        if Parameters.GC:
            ax_all_mu.plot(X[1:-1], muI[0:-1], label = '{}'.format(date), color = color)
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_A_h_{}_{}.png'.format(k, expno))
    
    # ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_A.grid()
    ax_all_A.set_xlabel('X (km)')
    ax_all_A.set_ylabel('A')
    fig_all_A.colorbar(sm, ax = ax_all_A, label = 'Time (hr)')
    fig_all_A.savefig(figdir+'time_plot_A_{}.png'.format(expno))
    
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # ax_all_h.grid()
    ax_all_h.set_xlabel(r'$x$ (km)')
    ax_all_h.set_ylabel(r'$h$ (m)')
    # ax_all_h.set_aspect('equal', 'datalim')
    fig_all_h.colorbar(sm, ax = ax_all_h, label = 'Time (hr)')
    fig_all_h.savefig(figdir+'time_plot_h_{}.png'.format(expno))
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # ax_all_u.grid()
    ax_all_u.set_xlabel('$x$ (km)')
    ax_all_u.set_ylabel('$u$ (m/s)')
    sm = cm.ScalarMappable(cmap = cmap,
                           norm = norm)
    fig_all_u.colorbar(sm, ax = ax_all_u, 
                       label = 'Time (hr)')
    fig_all_u.savefig(figdir+'time_plot_u_{}.png'.format(expno))
    plt.close()
    
    if Parameters.GC:
        ax_all_mu.grid()
        ax_all_mu.set_xlabel('$x$ (km)')
        ax_all_mu.set_ylabel(r'$\mu$')
        sm = cm.ScalarMappable(cmap = cmap,
                            norm = norm)
        fig_all_mu.colorbar(sm, ax = ax_all_mu, 
                        label = 'Time (hr)')
        fig_all_mu.savefig(figdir+'time_plot_mu_{}.png'.format(expno))
        plt.close()
    
    ax_all_div.grid()
    ax_all_div.set_xlabel('$x$ (km)')
    ax_all_div.set_ylabel(r'$\mu$')
    sm = cm.ScalarMappable(cmap = cmap,
                           norm = norm)
    fig_all_div.colorbar(sm, ax = ax_all_div, 
                       label = 'Time (hr)')
    fig_all_div.savefig(figdir+'time_plot_div_{}.png'.format(expno))
    plt.close()
        
        
def plot_variable(dx, dt, time, var, label, norm, cmap, expno, figdir, varname) : 
    
    var_shape = np.shape(var)
    # print(var_shape)
    X = np.arange(0, var_shape[1])*dx/1e3
    time_plot = time*dt/(60*60)
    # print(time_plot)
    # print(X)

    fig = plt.figure()
    ax = plt.axes()
    pc = plt.pcolormesh(time_plot, X, var.T, cmap = cmap, norm = norm)
    fig.colorbar(pc, ax = ax, label = label)
    plt.xlabel('Time (hr)')
    plt.ylabel('x (km)')
    plt.savefig(figdir+'Evolution_{}_{}.png'.format(expno, varname))
    plt.close()
        
    
def plot_velocities(max_vel, max_capping,colors, dt, tstep, expno, figdir, regime):
    
    time_plot = tstep*dt/(60*60)

    fig = plt.figure(figsize = (5, 4))
    ax = plt.axes()
    # ax2 = ax.twinx()
    
    # labels = []
    for i, maxvel in enumerate(max_vel):
        # print(i)
        
        mantissa, exponent = f"{ max_capping[i]:.2e}".split('e')
    
        # label = r"$\nu_{{max}} = {} \times 10^{{{}}}$".format(mantissa, int(exponent))
        

        pc = ax.plot(time_plot, maxvel,color = colors[i])
        ax.tick_params(axis='y',  color = colors[i],labelcolor=colors[i])
    
    # pc = plt.plot(time_plot, min_vel, color = 'b')
    fig.legend(loc='lower center', 
            labelcolor='linecolor',
            bbox_to_anchor=(0.5, 0.12))
    # plt.grid()
    ax.set_xlabel('Time (hr)')
    ax.set_ylabel('max{$u_i$} (m/s)')
    ax.grid()
    # ax.set_ylim(0, 10)
    # ax2.set_ylim(0, 0.01)
    # ax.set_yticks(np.linspace(0, 10, 5))
    # ax2.set_yticks(np.linspace(0, 0.01, 5))
    plt.savefig(figdir+'maxmin_velocities_{}.png'.format(regime))
    plt.close()
    
    
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
    
    
def plot_initial_conditions(A_tot, h_tot, u_tot, Parameters, figdir, expno):
    
    		#----- Initial Conditions ------#
    A_init = A_tot[0]
    h_init = h_tot[0]
    u_init = u_tot[0]
    X = np.arange(0, Parameters.Nx)*Parameters.dx/1e3
    plt.figure()
    ax1 = plt.axes()
    ax2 = ax1.twinx()
    # ax1.plot(X[1:-1], A_init[1:-1], color = 'k')
    # print(X)
    ax1.plot(X[1:-1], h_init[1:-1], color = 'k')
    ax1.set_ylabel(r'$h$ (m)',color = 'k' )
    ax1.set_xlabel(r'x (m)',color = 'k' )
    ax2.plot(X[1:-1], u_init[:-1], linestyle = '-', color = 'b')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax2.set_ylabel(r'$u$ (m/s)', color='b')
    ax2.tick_params(axis='y', labelcolor='b')
    plt.savefig(figdir+'init_conditions_{}.png'.format(expno))
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
    
    
    