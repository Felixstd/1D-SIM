import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')



#------- Viscous Experiment (1) -------#
datadict_viscous = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_exp_1.npy', 
                              allow_pickle=True).item()
datadict_viscous_energy = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_energy_saved_1.npy', 
                              allow_pickle=True).item()
datadict_viscous_energy_conv = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_energy_saved_7.npy', 
                              allow_pickle=True).item()





divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, dt, dx, tstep, num_regime = datadict_viscous.values()

values = list(datadict_viscous_energy.values())
Plat_Tavg, Pfric_S_Tavg, Pfric_R_Tavg, Pfric_S_visc_Tavg, Pfric_R_visc_Tavg, \
Pfric_S_plas_Tavg, Pfric_R_plas_Tavg, P_fric_R_pZ_vE_Tavg, P_fric_R_vZ_pE_Tavg, \
Ph_Tavg, Pw_Tavg, Ppot_Tavg = values

values_conv = list(datadict_viscous_energy_conv.values())
Plat_Tavg_conv, Pfric_S_Tavg_conv, Pfric_R_Tavg_conv, Pfric_S_visc_Tavg_conv, Pfric_R_visc_Tavg_conv, \
Pfric_S_plas_Tavg_conv, Pfric_R_plas_Tavg_conv, P_fric_R_pZ_vE_Tavg_conv, P_fric_R_vZ_pE_Tavg_conv, \
Ph_Tavg_conv, Pw_Tavg_conv, Ppot_Tavg_conv = values_conv
            
energy = [Plat_Tavg, Pfric_S_Tavg, Pfric_R_Tavg, Ph_Tavg, Pw_Tavg, Ppot_Tavg] 
energy_conv = [Plat_Tavg_conv, Pfric_S_Tavg_conv, Pfric_R_Tavg_conv, Ph_Tavg_conv, Pw_Tavg_conv, Ppot_Tavg_conv] 


labels = [r'$P_{lat}$ (W/m$^2$)', r'$P_{f}^S$ (W/m$^2$)', r'$P_{f}^R$ (W/m$^2$)', r'$P_{h}$ (W/m$^2$)', r'$P_{w}$ (W/m$^2$)', r'$P_{pot}$ (W/m$^2$)']
time = np.arange(0, 299999)*dt/(60*60)

sum_energy = np.sum(energy, axis = 0)
sum_energy_conv = np.sum(energy_conv, axis = 0)


fig, (ax1, ax2) = plt.subplots(1,2, figsize = (11, 4), sharex = True)

axs = [ax1, ax2]

ax1.plot(time, sum_energy[0], label = r'$\sum P$ ')
ax1.plot(time, Ph_Tavg[0], label = r'$P_{h}$')
ax1.plot(time, Ppot_Tavg[0], label = r'$P_{pot}$')
ax1.plot(time, Pw_Tavg[0], label = r'$P_{w}$')
ax1.plot(time, np.sum(energy[0:3], axis =0)[0], label = r'$P_i$')
ax1.set_title('Divergence')


ax2.plot(time, sum_energy_conv[0])
ax2.plot(time, Ph_Tavg_conv[0])
ax2.plot(time, Ppot_Tavg_conv[0])
ax2.plot(time, Pw_Tavg_conv[0])
ax2.plot(time, np.sum(energy_conv[0:3], axis =0)[0])
ax2.set_title('Convergence')


fig.legend(loc = 'outside upper center', bbox_to_anchor = (0.95, 0.9))
for ax in axs:
    
    ax.grid()
fig.supxlabel(r'Time (hr)')
fig.supylabel(r'P (W/m$^2$)', x = 0.05)
plt.savefig('Viscous_energy_tot.png')

Nx = np.shape(divergence_tot)[1]
X = np.arange(0, Nx)*dx

norm = colors.Normalize(vmin=0, vmax=tstep[-1])  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (11, 4), sharex = True)

axs = [ax1, ax2]
c = 0
for i, ax in enumerate(axs):
    
    ax.grid()
 
    divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, dt, dx, dates, num_regime = datadict_viscous.values()
    Nx = np.shape(divergence_tot)[1]
    # ax.set_aspect('auto')
    
    X = np.arange(0, Nx)*dx
    
    for k, date in enumerate(dates):
        divergence = divergence_tot[k]
        h, A = h_tot[k], A_tot[k]
        u = u_tot[k]
        eta, zeta = eta_tot[k], zeta_tot[k]
        
        color = cmap(k/10)
        
        if c > 0:
            ax.plot(X[1:-1], u[:-1], label = '{}'.format(date), color = color)
        else:
            ax.plot(X[1:-1], h[1:-1], label = '{}'.format(date), color = color)
       
    c +=1
    
fig.supxlabel(r'$x$ (km)', y = 0.01, x = 0.43)
fig.align_ylabels()
ax1.set_ylabel(r'$h$ (m)')
fig.colorbar(sm, ax=axs, label = 'Time (hr)')
ax2.set_ylabel(r'$u$ (m/s)')
fig.savefig('Viscous_numerical.png', dpi = 500)



#------ VP Experiments -------#
datadict_lowres = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_exp_2.npy', 
                              allow_pickle=True).item()
datadict_highres = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_exp_8.npy', 
                              allow_pickle=True).item()

divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, dt_low, dx_low, tstep_low, num_regime_low = datadict_lowres.values()
divergence_tot_high, h_tot_high, A_tot_high, u_tot_high, eta_tot_high, zeta_tot_high, dt_high, dx_high, tstep_high, num_regime_high = datadict_highres.values()

#---- u and h ----#
Nx_low = np.shape(divergence_tot)[1]
X_low = np.arange(0, Nx_low)*dx_low

norm = colors.Normalize(vmin=0, vmax=tstep_low[-1])  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (11, 8), sharex = True)

axs = [ax1, ax2, ax3, ax4]
c = 0
for i, ax in enumerate(axs):
    
    ax.grid()
 
    if i % 2 == 0:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, dt, dx, dates, num_regime = datadict_lowres.values()
        Nx = np.shape(divergence_tot)[1]
		
    else:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, dt, dx, dates, num_regime = datadict_highres.values()
        Nx = np.shape(divergence_tot)[1]
    
    X = np.arange(0, Nx)*dx
    
    for k, date in enumerate(dates):
        divergence = divergence_tot[k]
        h, A = h_tot[k], A_tot[k]
        u = u_tot[k]
        eta, zeta = eta_tot[k], zeta_tot[k]
        
        color = cmap(k/10)
        
        if c > 1:
            ax.plot(X[1:-1], u[:-1], label = '{}'.format(date), color = color)
        else:
            ax.plot(X[1:-1], h[1:-1], label = '{}'.format(date), color = color)
       
    c +=1
    
fig.supxlabel(r'$x$ (km)', y = 0.05, x = 0.43)
fig.align_ylabels()
ax1.set_ylabel(r'$h$ (m)')
ax1.set_title(r'$\Delta x = $ 1km')
ax2.set_title(r'$\Delta x = $ 500m')
fig.colorbar(sm, ax=axs, label = 'Time (hr)')
ax3.set_ylabel(r'$u$ (m/s)')
fig.savefig('VP_numerical.png', dpi = 500)

#---- number in viscous and plastic -----# 
tstep_tot = np.arange(0, np.shape(num_regime[0])[0])
print(tstep_tot)
plt.figure()
plt.plot(tstep_tot/(60*60)*dt, num_regime_low[0]/num_regime_low[1], linestyle = '-', color = 'r', label = r'$\Delta x = 1$ km')
plt.plot(tstep_tot/(60*60)*dt, num_regime_high[0]/num_regime_high[1], linestyle = '-', color = 'b',  label = r'$\Delta x = 500$ m')
plt.grid()
plt.legend()
plt.xlabel('Time (hr)')
plt.ylabel(r'$N_{visc}/N_{plas}$')
plt.savefig('number_regime_vicous.png')



