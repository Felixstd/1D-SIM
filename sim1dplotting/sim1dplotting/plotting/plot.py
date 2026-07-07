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
    

    divergence_tot = data_dict['divergence_dates']
    h_tot          = data_dict['h_dates']
    A_tot          = data_dict['A_dates']
    u_tot          = data_dict['u_dates']
    sigma_tot      = data_dict['sigma_dates']
    sigma_norm_tot = data_dict['signorm_dates']
    zeta_tot       = data_dict['zeta_dates']
    eta_tot        = data_dict['eta_dates']
    
    if Parameters.GC:
        P_tot   = data_dict['P_dates']
        muI_tot = data_dict['muI_dates'] 
        
    if Parameters.dissipation:
        Wdissip_tot = data_dict['Wsigma_dates']
    
    if Parameters.strength:
        P_tot = data_dict['P_dates']
        
        
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
    
def plot_initial_conditions(A_tot, h_tot, u_tot, X, figdir, expno):
    
    		#----- Initial Conditions ------#
    A_init = A_tot[0]
    h_init = h_tot[0]
    u_init = u_tot[0]
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
    
def plot_figures_paper(dates, 
                       datadict_lowres, datadict_highres, 
                       dx_lr, dx_hr, dt, 
                       figdir, figname, 
                       tanh = False, 
                       highres_3 = False,
                       datadict_hr_highres = None, 
                       dx_hr2 = None):
    

    if highres_3: 
        h_tot_hr2, A_tot_hr2, u_tot_hr2 = (datadict_hr_highres[k] for k in ('h_dates', 'A_dates', 'u_dates'))
        # print(np.shape(h_tot_hr2))
        Nx_hr2 = np.shape(h_tot_hr2)[1]
        X_hr2 = np.arange(Nx_hr2) * dx_hr2 / 1e3
        Nx_hr2_u = np.shape(u_tot_hr2)[1]
        X_hr2_u = np.arange(Nx_hr2_u) * dx_hr2 / 1e3
        
    h_tot_lr, A_tot_lr, u_tot_lr = (datadict_lowres[k] for k in ('h_dates', 'A_dates', 'u_dates'))
    h_tot_hr, A_tot_hr, u_tot_hr = (datadict_highres[k] for k in ('h_dates', 'A_dates', 'u_dates'))
    
    Nx_lr = np.shape(h_tot_lr)[1]
    Nx_hr = np.shape(h_tot_hr)[1]
    Nx_lr_u = np.shape(u_tot_lr)[1]
    Nx_hr_u = np.shape(u_tot_hr)[1]
    
    X_lr = np.arange(Nx_lr) * dx_lr / 1e3
    X_hr = np.arange(Nx_hr) * dx_hr / 1e3
    X_lr_u= np.arange(Nx_lr_u) * dx_lr / 1e3
    X_hr_u = np.arange(Nx_hr_u) * dx_hr / 1e3
    
    dates = dates*dt/(60*60)
    
    norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
    # cmap = plt.cm.viridis
    cmap = plt.get_cmap('viridis')
    # sm = cm.ScalarMappable(cmap=cmap, norm=norm)


    # instead of: cmap = plt.get_cmap('viridis_r')
    base_cmap = plt.get_cmap('viridis_r')
    cmap = truncate_colormap(base_cmap, 0.1,1) 
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    k_idx = [0,2,4,6,8,10,12,14,16,18,19]
    plt.rcParams.update({
        "font.size": 12,
        # "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    
    if highres_3:
        fig, ((ax1, ax3, ax5), (ax2, ax4, ax6)) = plt.subplots(2,3,figsize = (13,7),
                                               sharex = True,
                                               constrained_layout = True)
        
        axs = [ax1,ax2,ax3,ax4,ax5,ax6]    
        
        
        for k, date in enumerate(np.append(dates[::2],dates[-1])):
        # for k, date in enumerate(dates[-2:-1]):

            k = k_idx[k]
            print('Plotting: ', date)
            
            h_lr, A_lr,u_lr = h_tot_lr[k], A_tot_lr[k], u_tot_lr[k]
            h_hr, A_hr,u_hr = h_tot_hr[k], A_tot_hr[k], u_tot_hr[k]
            h_hr2, A_hr2,u_hr2 = h_tot_hr2[k], A_tot_hr2[k], u_tot_hr2[k]
            
            color = cmap(norm(date))
            
            
            ax1.plot(X_lr[1:-1], h_lr[1:-1], label = '{}'.format(date), color = color)
            ax2.plot(X_lr_u, u_lr, label = '{}'.format(date), color = color)
            
            ax3.plot(X_hr[1:-1], h_hr[1:-1], label = '{}'.format(date), color = color)
            ax4.plot(X_hr_u, u_hr, label = '{}'.format(date), color = color)
            
            ax5.plot(X_hr2[1:-1], h_hr2[1:-1], label = '{}'.format(date), color = color)
            ax6.plot(X_hr2_u, u_hr2, label = '{}'.format(date), color = color)
            
        
        # ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # labels = ['(a)','(b)','(c)','(d)','(e)','(f)']
        labels = ['(a)','(d)','(b)','(e)','(c)','(f)']
        
        for i, ax in enumerate(axs):
            ax.grid()
            ax.text(0.04, 0.955, labels[i], transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top', horizontalalignment='left')
            
        
            
            
        fig.supxlabel('x (km)')
        ax1.set_ylabel('$A, h$ (m)')
        ax2.set_ylabel(r'$u$ (ms$^{-1}$)')
        fig.align_ylabels()
        
        if tanh:
            ax1.set_title(r'Local')
            ax3.set_title(r'Nonlocal')
            
            ax2.set_ylim(-1.5,1.5)
            ax4.set_ylim(-1.5,1.5)
            
        elif highres_3: 
            ax1.set_title(r'$\Delta x = 1$ km')
            ax3.set_title(r'$\Delta x = 500$ m')
            ax5.set_title(r'$\Delta x = 250$ m')
            
            ax2.set_ylim(-0.41,0.41 )
            ax4.set_ylim(-0.41,0.41 )
            ax6.set_ylim(-0.41,0.41 )
        
        else:
            ax1.set_title(r'$\Delta x = 1$ km')
            ax3.set_title(r'$\Delta x = 500$ m')
            
            
            ax2.set_ylim(-0.6,0.6)
            ax4.set_ylim(-0.6,0.6)
            

        fig.colorbar(sm, ax = [ax5,ax6], label = 'Time (hr)')
        fig.savefig(figdir+figname)
        
        plt.close()
        
    else:
    
        fig, ((ax1, ax3),(ax2,ax4)) = plt.subplots(2,2,figsize = (9,6),
                                                sharex = True,
                                                constrained_layout = True)
        
        axs = [ax1,ax2,ax3,ax4]    
        
        
        for k, date in enumerate(np.append(dates[::2],dates[-1])):
        # for k, date in enumerate(dates[-2:-1]):

            k = k_idx[k]
            print('Plotting: ', date)
            
            h_lr, A_lr,u_lr = h_tot_lr[k], A_tot_lr[k], u_tot_lr[k]
            h_hr, A_hr,u_hr = h_tot_hr[k], A_tot_hr[k], u_tot_hr[k]

            
            color = cmap(k/len(dates))
            
            
            ax1.plot(X_lr[1:-1], h_lr[1:-1], label = '{}'.format(date), color = color)
            ax2.plot(X_lr[2:-1], u_lr[1:-1], label = '{}'.format(date), color = color)
            
            ax3.plot(X_hr[1:-1], h_hr[1:-1], label = '{}'.format(date), color = color)
            ax4.plot(X_hr[2:-1], u_hr[1:-1], label = '{}'.format(date), color = color)
            
        
        # ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        labels = ['(a)','(b)','(c)','(d)']
        
        for i, ax in enumerate(axs):
            ax.grid()
            ax.text(0.04, 0.955, labels[i], transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top', horizontalalignment='left')
            
        
            
            
        fig.supxlabel('x (km)')
        ax1.set_ylabel('$A, h$ (m)')
        ax2.set_ylabel(r'$u$ (ms$^{-1}$)')
        fig.align_ylabels()
        
        if tanh:
            ax1.set_title(r'Local')
            ax3.set_title(r'Nonlocal')
            
            ax2.set_ylim(-1.5,1.5)
            ax4.set_ylim(-1.5,1.5)
        
        else:
            ax1.set_title(r'$\Delta x = 1$ km')
            ax3.set_title(r'$\Delta x = 500$ m')
            
            ax2.set_ylim(-0.6,0.6)
            ax4.set_ylim(-0.6,0.6)
        
        

        fig.colorbar(sm, ax = [ax3,ax4], label = 'Time (hr)')
        fig.savefig(figdir+figname)
        
        plt.close()
        
        
        
def plot_E_VP_resolution(dates, data_exps, Parameters, 
                           figdir, figname):
    """
    This is a function to plot different models at different resolution. 
    
    This figure should plot the standard VP and EVP. 
    
    The data is organised as:
        1. VP-1000m
        2. VP-500m
        3. VP-250m
        4. EVP-1000m
        5. EVP-500m
        6. EVP-250m

    Args:
        dates (_type_): _description_
        data_exps (_type_): _description_
        Parameters (_type_): _description_
        figdir (_type_): _description_
        figname (_type_): _description_
    """
    
    plt.rcParams.update({
        "font.size": 12,
        # "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    
    dates = dates*Parameters.dt/(60*60)
    
    # norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
    # instead of: cmap = plt.get_cmap('viridis_r')
    base_cmap = plt.get_cmap('viridis_r')
    cmap = truncate_colormap(base_cmap, 0.1,1) 
    norm = Normalize(vmin=dates[0], vmax=dates[-1])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    
    k_idx = [0,2,4,6,8,10,12,14,16,18,19]
    
    fig, axs2d = plt.subplots(4,3,figsize = (12,9),
                                               sharex = True,
                                               sharey='row',
                                               constrained_layout = True)
    axs = axs2d.T.ravel()
    # exp index -> (h_axis, u_axis)
    exp_axes = {1: (axs[0], axs[1]),
                2: (axs[4], axs[5]),
                3: (axs[8], axs[9]), 
                4: (axs[2], axs[3]),
                5: (axs[6], axs[7]),
                6: (axs[10], axs[11])}   # fixed: was colliding with exp 2 before

    date_positions = np.append(dates[::2], dates[-1])
    colors = cmap(np.linspace(0, 1, len(dates))) 
    
    for i, exp in enumerate(Parameters.expno):
            print(i)
            ax_h, ax_u = exp_axes[i+1]
            data = data_exps[i]
            h_tot, u_tot, x = data['h'], data['u'], data['x']

            for pos, date in enumerate(date_positions):
                k = k_idx[pos]
                h, u = h_tot[k], u_tot[k]
                ax_h.plot(x[1:-1]/1e3, h[1:-1], label=str(date), color=colors[k])
                ax_u.plot(x[1:-1]/1e3, u[:-1], color=colors[k]) 
                
    labels = ['(a)','(d)','(g)','(j)','(b)','(e)','(h)','(k)','(c)','(f)','(i)','(l)']
        
    for i, ax in enumerate(axs):
        ax.grid()
        ax.text(0.04, 0.955, labels[i], transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top', horizontalalignment='left')
            
        
            
            
    fig.supxlabel('x (km)')
    axs[0].set_ylabel('$A, h$ (m)')
    axs[1].set_ylabel(r'$u$ (ms$^{-1}$)')
    axs[2].set_ylabel('$A, h$ (m)')
    axs[3].set_ylabel(r'$u$ (ms$^{-1}$)')
    axs[0].set_title(r'$\Delta x = 1$ km')
    axs[4].set_title(r'$\Delta x = 500$ m')
    axs[8].set_title(r'$\Delta x = 250$ m')
    fig.align_ylabels()
    fig.colorbar(sm, ax = axs[8:], 
                 aspect = 40,
                 label = 'Time (hr)')
    
    fig.savefig(figdir+figname)
    
    
        
        
def make_figure_A(Parameters, data_exps, n_lines):
    fig2, ax3 = plt.subplots(1, 1, figsize=(5, 3), facecolor='none')

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_facecolor('none')

    # Collect all line data first
    all_lines = []
    for i, exp in enumerate(Parameters.expno):
        data = data_exps[i]
        A_tot = data['A']
        x = data['x']
        if int(Parameters.lscale[i]) == 0:
            all_lines.append((x, A_tot, 'goldenrod', r'VP', True))   # True = is VP
        else:
            all_lines.append((x, A_tot, colors[i-1], r'$l_c = {}$ km'.format(int(Parameters.lscale[i])), False))

    # Plot visible lines (up to n_lines)
    for i, (x, A_tot, color, label, is_vp) in enumerate(all_lines):
        if i < n_lines:
            if is_vp:
                ax3.plot(x[1:-1]/1e3, A_tot[-1][1:-1], color=color,
                        linewidth=1.6, alpha=0.7, zorder=3, label=label)
            else:
                ax3.plot(x[1:-1]/1e3, A_tot[-1][1:-1], color=color,
                        linewidth=1.2, label=label)
        else:
            # Invisible, no label — keeps axes fixed
            if is_vp:
                ax3.plot(x[1:-1]/1e3, A_tot[-1][1:-1], alpha=0, linewidth=1.6, zorder=3)
            else:
                ax3.plot(x[1:-1]/1e3, A_tot[-1][1:-1], alpha=0, linewidth=1.2)

    leg = ax3.legend(
        loc='center left',
        bbox_to_anchor=(0.77, 0.5),
        frameon=False,
        handlelength=0,
        handletextpad=0
    )
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        text.set_color(handle.get_color())

    ax3.set_xlabel('x (km)')
    ax3.set_ylabel(r'$A$')

    # fig2.savefig(f'A_end_newP_presentation_{n_lines}.png', bbox_inches='tight', transparent=True)
    plt.close()
    
def plot_sensitivity_lc(Parameters, data_exps,
                        figname,
                        number_plots = 2):
    """
    Plotting function to analyse the sensitivity of the 
    cooperative length scale on the fractures. 

    Args:
        Parameters (_type_): _description_
        data_exps (_type_): _description_
        number_plots (int, optional): _description_. Defaults to 2.
    """

    labels = ['(a)','(b)','(c)','(d)']
    
    start, end = 0.4, 1
    start,end = 0,0.5
    original_cmap = plt.get_cmap('magma')
    norm = LogNorm(vmin=np.min(Parameters.lscale[1:])*0.5, vmax=np.max(Parameters.lscale)*1.2)

    # 2. Sample 256 colors from that specific range
    colors = original_cmap(np.linspace(start, end, len(data_exps)))
    
    print(original_cmap(250))
    
    if number_plots == 2:
        fig,((ax1,ax2)) = plt.subplots(2,1, 
                                       figsize=(6,5), 
                                       sharex = True, 
                                       constrained_layout = True)
        axs = [ax1,ax2]
    else:
        fig,((ax3,ax4),(ax1,ax2)) = plt.subplots(2,2, figsize = (9,6), sharex = True,constrained_layout = True)
        axs = [ax3,ax4,ax1,ax2]
    
    plt.rcParams.update({
        "font.size": 12,
        # "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    
    
    for i, ax in enumerate(axs):
        ax.set_axisbelow(True)
        ax.grid(zorder = 0,alpha = 0.5)
        ax.text(0.03, 0.97, labels[i], transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top', horizontalalignment='left')	
        
    for i, exp in enumerate(Parameters.expno):
        data = data_exps[i]
        A_tot = data['A']
        h_tot = data['h']
        P_tot = data['base']['P_dates']
        u_tot = data['u']
        x = data['x']
        
        if int(Parameters.lscale[i]) == 0:
            ax1.plot(x[1:-1]/1e3, A_tot[-1][1:-1], color = 'goldenrod', linewidth = 1.6, alpha = 0.7,zorder = 0,label = r'VP',)
            ax2.plot(x[1:-1]/1e3, u_tot[-1][0:-1], color = 'goldenrod', linewidth = 1.6, alpha = 0.7,zorder = 0)
        # elif int(Parameters.lscale[i]) == 100:
        #     ax1.plot(x[1:-1]/1e3, A_tot[-1][1:-1], color = colors[i-1], linewidth = 1.2,zorder = 3,label = r'$l_c = {}$ km'.format(int(Parameters.lscale[i])))
        #     ax2.plot(x[1:-1]/1e3, u_tot[-1][0:-1], color = colors[i-1], linewidth = 1.2,zorder = 3)
        else:

            ax1.plot(x[1:-1]/1e3, A_tot[-1][1:-1], color = original_cmap(norm(Parameters.lscale[i])), linewidth = 1.2,zorder = 2,label = r'$l_c = {}$ km'.format(int(Parameters.lscale[i])))
            ax2.plot(x[1:-1]/1e3, u_tot[-1][0:-1], color = original_cmap(norm(Parameters.lscale[i])), linewidth = 1.2,zorder = 2)
              
        if number_plots == 4:

            if int(Parameters.lscale[i]) == 0:
                ax3.plot(x[1:-1]/1e3, A_tot[-15][1:-1], color = 'goldenrod', linewidth = 1.6, alpha = 0.7,zorder = 0)
                ax4.plot(x[1:-2]/1e3, u_tot[-15][1:-1], color = 'goldenrod', linewidth = 1.6, alpha = 0.7,zorder = 0)
            
            elif int(Parameters.lscale[i]) == 100:
                ax3.plot(x[1:-1]/1e3, A_tot[-15][1:-1], color =colors[i-1], linewidth = 1.2,zorder = 3)
                ax4.plot(x[1:-2]/1e3, u_tot[-15][1:-1], color =colors[i-1], linewidth = 1.2,zorder = 3)
            
            else:
                ax3.plot(x[1:-1]/1e3, A_tot[-15][1:-1], color = colors[i-1], linewidth = 1.2,zorder = 2)
                ax4.plot(x[1:-2]/1e3, u_tot[-15][1:-1], color = colors[i-1], linewidth = 1.2,zorder = 2)
    
    # fig.legend(bbox_to_anchor = (1.18,0.87), labelcolor = 'linecolor')
    leg = ax1.legend(
            loc='center left',
            bbox_to_anchor=(1.03, 0.5),   # just outside ax1
            frameon=False,
            handlelength=0,
            handletextpad=0
        )
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        text.set_color(handle.get_color())

    ax2.set_ylim(-0.6,0.6)
    ax1.set_title(r'$t = 53$ (hr)')
    ax1.set_ylabel(r'$A, h$ (m)',fontsize =12)
    ax.set_yscale('symlog',linthresh = 1e-5)
    
    
    if number_plots ==2:
        ax2.set_ylabel(r'$u$ (ms$^{-1}$)',fontsize =12)
        fig.supxlabel('x (km)',x = 0.44,y = -0.05)
        
    elif number_plots == 4:
        ax3.set_title(r'$t = 14$ (hr)')
        ax3.set_ylabel(r'$u$ (ms$^{-1}$)',fontsize =12)
        ax4.set_ylim(-0.4,0.4)
    
    
    # ax2.set_xlabel('x (km)')
    


    fig.align_ylabels()
    plt.savefig(figname)
    plt.close(2)

def plot_sensitivity_regul(dates,Parameters, data_exps):
    
    original_cmap = plt.get_cmap('magma')
    norm = LogNorm(vmin=10*0.5, vmax=250*1.2)
    
    dt = Parameters.dt

    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize = (9,6),constrained_layout = True,sharex =True)
			# fig,((ax1,ax2)) = plt.subplots(2,1, figsize=(4.5,5), sharex = True)

    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    # print(np.where(np.array(Parameters.lscale)>0)[0][0])
    axs = [ax1,ax2,ax3,ax4]
    labels_panel = ['(a)','(b)','(c)','(d)']
    # labels_regul = ['VP','nVP','tanh','EVP']
    # labels_regul = ['tanh','EVP']
    # colors = ['goldenrod',original_cmap(norm(Parameters.lscale[np.where(np.array(Parameters.lscale)>0)[0][0]])), 
    #           'royalblue','darkgreen']
    
    # colors = ['royalblue','darkgreen']
            
    dates = dates*dt/(60*60)
    
    norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
    cmap = plt.cm.viridis
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    for i, ax in enumerate(axs):
        ax.grid(zorder = 0,alpha = 0.5)
        ax.text(0.03, 0.97, labels_panel[i], transform=ax.transAxes,
        	fontsize=10,
        	verticalalignment='top', horizontalalignment='left')
    
    for i, exp in enumerate(Parameters.expno):
        data = data_exps[i]
        A_tot = data['A']
        x = data['x']

        k_idx = [0,2,4,6,8,10,12,14,16,18,19]

        for k, date in enumerate(np.append(dates[::2],dates[-1])):
            k = k_idx[k]
            color = cmap(k/len(dates))
            
            if Parameters.lscale[i] == 0:
                ax1.plot(x[1:-1]/1e3, A_tot[k][1:-1], alpha = 0.8, color = color)
                
            elif Parameters.lscale[i] == 1:
                ax3.plot(x[1:-1]/1e3, A_tot[k][1:-1], alpha = 0.8, color = color)
            elif Parameters.lscale[i] == 2:
                ax2.plot(x[1:-1]/1e3, A_tot[k][1:-1], alpha = 0.8, color = color)
            else:
                ax4.plot(x[1:-1]/1e3, A_tot[-k][1:-1], alpha = 0.8, color = color)

            
    
    # fig.legend(bbox_to_anchor = (1.18,0.87), labelcolor = 'linecolor')
    leg = ax1.legend(
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),   # just outside ax1
            frameon=False,
            handlelength=0,
            handletextpad=0
        )
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        text.set_color(handle.get_color())

    
    # ax1.set_xlabel('x (km)')
    fig.colorbar(sm, ax = [ax2,ax4], label = 'Time (hr)')
    fig.supxlabel('x (km)')
    ax1.set_title(r'$\Delta x = 1$ km ')
    ax2.set_title(r'$\Delta x = 500$ km ')
    fig.supylabel(r'$A, h$ (m)',fontsize =12)
    
    plt.savefig('./Figures/A_end_RP_regul.png')
    plt.close(2)

def plot_ridgebuilding(Parameters, data_exps):
    """
    Function to make the plots for the ridge building 
    experiments. A linear fit on the ridge is computed for every 
    lines with VP and nVP to compare with the theoretical profile. 

    Args:
        Parameters (_type_): _description_
        data_exps (_type_): _description_
    """
    
    original_cmap = plt.get_cmap('magma')
    norm = Normalize(vmin=5, vmax=250)

    file_fits = open("./output_fits_ridging.txt", "w")

    fig,(ax1) = plt.subplots(1,1, figsize = (8,4),constrained_layout = True,sharex =True)


    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })    
    
    for i, exp in enumerate(Parameters.expno):
        data = data_exps[i]
        h_tot = data['h']
        
        
        x = data['x']
        
        l_scale = Parameters.lscale[i]

        color = original_cmap(norm(l_scale))

        indice = np.where(h_tot[-1][1:-1] > 1.05)
        if len(indice[0]) < 1:
            continue
        else:
            indice = indice[0][-1]
            
        if indice > 498:
            continue
        else:
            fit = np.polyfit(x[1:indice], h_tot[-1][1:indice],deg=1,cov=True) 
            
        if Parameters.lscale[i] == 0:
            ax1.plot(x[1:-1]/1e3, h_tot[-1][1:-1], alpha = 0.8, color = 'royalblue', label = 'VP')
            m = (2*1.3*1.2e-3*(5)**2/(27.5e3*(np.sqrt(1+2**(-2))+1)))
            theory = -m*x[1:indice]+fit[0][1]
            ax1.plot(x[1:indice]/1e3, theory, color = 'gray', linestyle = '--', label = 'theory',zorder = 3)
            file_fits.write('Theory slope: '+str(m)+'\n')
        
        else:
            ax1.plot(x[1:-1]/1e3, h_tot[-1][1:-1], alpha = 0.8, 
                    color = color, label = r'$l_c = {}$ km'.format(int(Parameters.lscale[i])))

        err_slope = fit[1][0,0]
        
        file_fits.write(f'l_c: {Parameters.lscale[i]} slope: {fit[0][0]} err_slope: {err_slope}\n')
        
        
    leg = ax1.legend(
            loc='center left',
            bbox_to_anchor=(0.81, 0.7),   # just outside ax1
            frameon=False,
            handlelength=0,
            handletextpad=0
        )
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        text.set_color(handle.get_color())

    ax1.grid(alpha = 0.5)
    fig.supxlabel('x (km)',x = 0.52)
    fig.supylabel(r'$h$ (m)',fontsize =12)
    plt.xlim(0,2000)
    plt.savefig('./Figures/ridge_building_lscale2.png')
    plt.close(2)
    file_fits.close()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    