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
    
        
def plot_individual_time(dates, expno, data_dict, dx, figdir, mu_0, mu_infty, MuPhi = True, Dissipation = False):
    
        
    if Dissipation:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, Wdissip_tot = data_dict.values()
    else:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot = data_dict.values()
    
        
    Nx = np.shape(divergence_tot)[1]
    X = np.arange(0, Nx)*dx/1e3
    fig_all_A = plt.figure(1)
    ax_all_A = plt.axes()
    
    fig_all_h = plt.figure(2, figsize = (5, 4))
    ax_all_h = plt.axes()
    
    fig_all_u = plt.figure(3)
    ax_all_u = plt.axes()
    
    print(dates)
    
    norm = colors.Normalize(vmin=0, vmax=dates[-1]/(60*60))  # assuming k ranges from 0 to 10
    cmap = plt.cm.viridis
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        divergence = divergence_tot[k]
        h, A = h_tot[k], A_tot[k]
        u = u_tot[k]
        eta, zeta = eta_tot[k], zeta_tot[k]
        
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
        
        fig.savefig(figdir+'ice_u_{}_{}.png'.format(date, expno))
        
        fig, (ax1)= plt.subplots(1, 1, sharex = True, figsize = (10, 12))
        
        ax1.plot(X[1:-1], eta[1:-1], color = 'r')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel(r'$\eta$ (Ns/m)')
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_eta_{}_{}.png'.format(date, expno))
        
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
            
            fig.savefig(figdir+'ice_wdissip_{}_{}.png'.format(date, expno))
        
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
        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_A_h_{}_{}.png'.format(date, expno))
    
    # ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_A.grid()
    ax_all_A.set_xlabel('X (km)')
    ax_all_A.set_ylabel('A')
    fig_all_A.colorbar(sm, ax = ax_all_A)
    fig_all_A.savefig(figdir+'time_plot_A_{}.png'.format(expno))
    
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_h.grid()
    ax_all_h.set_xlabel(r'$x$ (km)')
    ax_all_h.set_ylabel(r'$h$ (m)')
    # ax_all_h.set_aspect('equal', 'datalim')
    fig_all_h.colorbar(sm, ax = ax_all_h, label = 'Time (hr)')
    fig_all_h.savefig(figdir+'time_plot_h_{}.png'.format(expno))
    plt.close()
    
    # ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_u.grid()
    ax_all_u.set_xlabel('$x$ (km)')
    ax_all_u.set_ylabel('$u$ (m/s)')
    sm = cm.ScalarMappable(cmap = cmap,
                           norm = norm)
    fig_all_u.colorbar(sm, ax = ax_all_u, 
                       label = 'Time (hr)')
    fig_all_u.savefig(figdir+'time_plot_u_{}.png'.format(expno))
    plt.close()
        
        
def plot_variable(dx, dt, time, var, label, norm, cmap, expno, figdir, varname) : 
    
    var_shape = np.shape(var)
    # print(var_shape)
    X = np.arange(0, var_shape[1])*dx/1e3
    time_plot = time*dt/(60*60)
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
    ax2 = ax.twinx()
    
    # labels = []
    for i, maxvel in enumerate(max_vel):
        print(i)
        
        mantissa, exponent = f"{ max_capping[i]:.2e}".split('e')
    
        label = r"$\nu_{{max}} = {} \times 10^{{{}}}$".format(mantissa, int(exponent))
        
        if i > 0:
            ax2.plot(time_plot, maxvel,color = colors[i], label = label)
            ax2.tick_params(axis='y', color = colors[i], labelcolor=colors[i])
            
        else:
            pc = ax.plot(time_plot, maxvel,color = colors[i], label = label)
            ax.tick_params(axis='y',  color = colors[i],labelcolor=colors[i])
    
    # pc = plt.plot(time_plot, min_vel, color = 'b')
    fig.legend(loc='lower center', 
            labelcolor='linecolor',
            bbox_to_anchor=(0.5, 0.12))
    # plt.grid()
    ax.set_xlabel('Time (hr)')
    ax.set_ylabel('max{$u_i$} (m/s)')
    ax.grid()
    ax.set_ylim(0, 10)
    ax2.set_ylim(0, 0.01)
    ax.set_yticks(np.linspace(0, 10, 5))
    ax2.set_yticks(np.linspace(0, 0.01, 5))
    plt.savefig(figdir+'maxmin_velocities_{}.png'.format(regime))
    plt.close()
    
    
    
    