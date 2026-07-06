import numpy as np
from scipy.signal import chirp, find_peaks, peak_widths
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


def find_FractureWidth(data_exps, l_scale, 
                       figdir = './Figures/Figs_sens_lscale/',
                       rel_height = 0.5):
    """
    Function used to compute the width of the fracture 
    for different length scale.
    
    It computes the fracture width at 25,50,75 % of the 


    Args:
        data_exps (_type_): _description_
        l_scale (_type_): _description_
    """
    rel_heights = np.linspace(0.25,0.75,3)
    
    n_exps = len(data_exps)
    n_rh = len(rel_heights)

    fracture_widths = np.full((n_exps, n_rh), np.nan)
    fracture_A = np.full((n_exps, n_rh), np.nan)
    fracture_errors = np.full((n_rh-1,n_exps), np.nan)
    
    for k, data in enumerate(data_exps):
        
        A = data['A'][-1]
        A_init = data['A'][0]
    
        peaks, _ = find_peaks(-A)
        
        for i,rel_height in enumerate(rel_heights):

            results_half = peak_widths(-A, peaks, rel_height=rel_height)
            widths, width_heights, left_ips, right_ips = results_half
            idx = np.argmax(widths)
            
            fracture_widths[k,i] = widths[idx]
            fracture_A[k,i] = -width_heights[idx]
            
        
        plt.figure()
        plt.plot(A_init)
        plt.plot(A)
        plt.plot(peaks, A[peaks], "x")
        # plt.hlines(*results_half[1:], color="C2")
        plt.hlines(-width_heights[idx], left_ips[idx], right_ips[idx], color="C2")
        plt.savefig(figdir+'A_last_{}_75.png'.format(l_scale[k]))

    fracture_errors[0,:] = abs(fracture_widths[:,1]-fracture_widths[:,0])
    fracture_errors[1,:] = abs(fracture_widths[:,1]-fracture_widths[:,2])

    color_w = 'darkgreen'
    l_c_ticks = np.linspace(0,250,6)
    plt.figure(figsize = (5,3))
    ax = plt.axes()
    ax.grid(zorder = 0 ,alpha = 0.7)

    sc = ax.scatter(l_scale, fracture_widths[:,1], 
                    c=fracture_A[:,1], cmap='viridis',
                    vmin = np.min(fracture_A[:,1]),vmax=np.max(fracture_A[:,1]),
                    zorder = 3)
    
    
    
    ax.errorbar(l_scale, fracture_widths[:,1],yerr = fracture_errors,fmt='none',    
    elinewidth=1,     # Makes the vertical error bar line thin
    capsize=4,        # Adds the horizontal bar (cap) at the ends
    capthick=1      ,  # Makes the horizontal bar thin)
    ecolor = 'black',
    zorder = 2
    )
    ax.set_yscale('log')
    ax.set_ylabel('Width (km)')
    cbar = plt.colorbar(sc)
    cbar.set_label(r'$A$')
    ax.set_xlabel(r'$l_c$ (km)' )
    
    ax.axhspan(xmin = 0,xmax=250,ymin=0,ymax=200, color = 'gray', 
                     edgecolor = None,alpha = 0.2,zorder = 0)
    
    ax.set_xticks(l_c_ticks)
    
    print(fracture_widths)
    

    plt.savefig(figdir+'localization_lscale.png')
