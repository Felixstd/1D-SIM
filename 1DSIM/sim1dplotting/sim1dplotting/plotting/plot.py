import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean
# import scienceplots

# plt.style.use('science')

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
    
        
def plot_individual_time(dates, expno, data_dict, dx, figdir, mu_0, mu_infty, MuPhi = True):
    
        
    h_tot, A_tot, u_tot = data_dict.values()
    
        
    Nx = np.shape(h_tot)[1]
    X = np.arange(0, Nx)*dx

    fig_all_A = plt.figure(1)
    ax_all_A = plt.axes()
    
    fig_all_h = plt.figure(2)
    ax_all_h = plt.axes()
    
    for k, date in enumerate(dates):
        
        print('Plotting: ', date)
        
        # divergence = divergence_tot[k]
        h, A = h_tot[k], A_tot[k]
        u = u_tot[k]
        
        
        #----- Deformation Plots -----#
        
        fig, (ax1)= plt.subplots(1, 1, sharex = True, figsize = (10, 12))
        
        ax1.plot(X[1:-1], u[1:-1], color = 'r')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('u (m/s)')
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_u_{}_{}.png'.format(date, expno))
        
        #--- Plot of A and h ---# 
        fig, (ax1, ax2)= plt.subplots(1, 2, sharex = True, figsize = (10, 12))
        
        ax1.plot(X[1:-1], A[1:-1], color = 'r')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('A')
        # ax1.set_xlim(Nx)
        ax1.grid()
        # ax1.set_aspect('equal')
        ax1.set_aspect(1./ax1.get_data_ratio())

        ax2.plot(X[1:-1], 1-h[1:-1], color = 'r')
        ax2.set_xlabel('X (km)')
        ax2.set_ylabel('1-h (m)')
        ax2.grid()
        # ax2.set_aspect('equal')
        ax2.set_aspect(1./ax2.get_data_ratio())
        
        ax_all_A.plot(X[1:-1], A[1:-1], label = '{}'.format(date))
        ax_all_h.plot(X[1:-1], 1-h[1:-1], label = '{}'.format(date))
        
        
        # plt.axis('equal')
        
        fig.savefig(figdir+'ice_A_h_{}_{}.png'.format(date, expno))
    
    ax_all_A.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_A.grid()
    ax_all_A.set_xlabel('X (km)')
    ax_all_A.set_ylabel('A')
    fig_all_A.savefig(figdir+'time_plot_A_{}.png'.format(expno))
    
    ax_all_h.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax_all_h.grid()
    ax_all_h.set_xlabel('X (km)')
    ax_all_h.set_ylabel('1-h (m)')
    fig_all_h.savefig(figdir+'time_plot_h_{}.png'.format(expno))
        
        
def plot_variable(dx, time, var, label, norm, cmap, expno, figdir, varname) : 
    
    var_shape = np.shape(var)
    
    X = np.arange(0, var_shape[1])*dx
    
    fig = plt.figure()
    ax = plt.axes()
    pc = plt.pcolormesh(time, X, var.T, cmap = cmap, norm = norm)
    fig.colorbar(pc, ax = ax, label = label)
    plt.xlabel('Time (hr)')
    plt.ylabel('x (km)')
    plt.savefig(figdir+'Evolution_{}_{}.png'.format(expno, varname))
        
    
    
    
    