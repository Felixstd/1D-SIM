import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import matplotlib.cm as pltcm

from sim1dplotting.utils.TimeUtilities import TimeUtility
from sim1dplotting.data.readnamelist import namelist
from sim1dplotting.data import read_data
from sim1dplotting.plotting import plot
from sim1dplotting.analysis import analysis


import configparser
import cmocean as cm
import warnings
import os


plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')


config_exp = configparser.ConfigParser()
config_exp.read('./namelistAnalysis')




Parameters = namelist(configuration_exp  = config_exp['Experiment'], 
                      configuration_rheo = config_exp['Rheology'], 
                      configuration_fig  = config_exp['Figures'], 
                      configuration_time = config_exp['Time'])

warnings.filterwarnings("ignore")

# For 1 to 13
#--- Dates of the simulation ---#
#---   They are in seconds   ---#

b_exp = []
A_exp = []
u_exp = []
zeta_exp = []
eta_exp = []
maxvel_exp = []

for exp in Parameters.expno: 
	print('Reading Experiment', exp)
	Dates_Config = TimeUtility(configuration_time = config_exp['Time'], 
                           configuration_fig  = config_exp['Figures'])
	tstep, dates_time = (TimeUtility.read_time(Dates_Config, expno = exp))
	# tstep = [int(1)]
	# dates_time = [0, 1]
	# tstep = tstep[:2]
	# print(tstep)

	if not os.path.isdir(Parameters.figdir+str(exp)):
		os.mkdir(Parameters.figdir+str(exp))
	
	if Parameters.read_all: 
		#---------- READING DATA ----------#
  
		datadict = read_data.read_data(exp, 
									int(Parameters.dt), 
									int(Parameters.dx/1e3), 
									Parameters.solv, 
									Parameters.imex, 
									Parameters.adv, 
									tstep, 
									Parameters.outputdir, MuPhi = Parameters.muphi, Dissipation = Parameters.dissipation)
		if Parameters.dx < 1:
			if Parameters.Nx == 1002:
				Parameters.dx= 1
		if Parameters.dissipation:
			divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, Wdissip_tot = datadict.values()
		else:
			divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot = datadict.values()
   
		#----- Initial Conditions ------#
		A_init = A_tot[0]
		h_init = h_tot[0]
		u_init = u_tot[0]
		X = np.arange(0, Parameters.Nx)*Parameters.dx/1e3
		plt.figure()
		ax1 = plt.axes()
		ax2 = ax1.twinx()
		# ax1.plot(X[1:-1], A_init[1:-1], color = 'k')
		ax1.plot(X[1:-1], h_init[1:-1], color = 'k')
		ax1.set_ylabel(r'$h$ (m)',color = 'k' )
		ax1.set_xlabel(r'x (m)',color = 'k' )
		ax2.plot(X[1:-1], u_init[:-1], linestyle = '--', color = 'b')
		ax2.set_ylabel(r'$u$ (m/s)', color='b')
		ax2.tick_params(axis='y', labelcolor='b')
		plt.savefig('init_conditions.png')

	if Parameters.savevar:
		datadict_saved = datadict.copy()
		datadict_saved['timestep'] = int(Parameters.dt)
		datadict_saved['resolution'] = int(Parameters.dx/1e3)
		datadict_saved['dates'] = dates_time
		np.save('SavedExperiments/datadict_exp_{}.npy'.format(exp), datadict_saved)
# X = np.linspace(0, nx*dx, nx)

# 
	if Parameters.plotfields:
   
		plot.plot_variable(Parameters.dx, 
							Parameters.dt,
							tstep, 
							np.array(A_tot), 
							r'$A$', 
							colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), 
							cm.cm.ice, 
							exp,
							Parameters.figdir+str(exp)+'/', 
							'A')
		plot.plot_variable(Parameters.dx, 
							Parameters.dt,
							tstep, 
							np.array(h_tot), 
							r'$h$ (m)', 
							colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), 
							cm.cm.ice, 
							exp,
							Parameters.figdir+str(exp)+'/', 
							'h')


#---------- Plotting ----------#
	if Parameters.plotsingle:
		plot.plot_individual_time(tstep, 
                            		exp, 
                              		datadict, 
                                	Parameters.dx,
									Parameters.figdir+str(exp)+'/', 
         							0, 
                					0,
                     				MuPhi = Parameters.muphi,
                         			Dissipation = Parameters.dissipation)
  
   
	if Parameters.maxvelocities:
     
	
		file = "/aos/home/fstdenis/1D-SIM/output_post_files/min_max_vel_{}.out".format(exp)
		max_velocities, min_velocities, tstep = read_data.read_maxvelocities(file)
		print(len(max_velocities))
		maxvel_exp.append(max_velocities)
   
		
 
if Parameters.maxvelocities:
	plot.plot_velocities(maxvel_exp, Parameters.maxcapping, [mpl.colormaps['Set1'].colors[2],mpl.colormaps['Set1'].colors[0]] , Parameters.dt, tstep, exp, './', 'viscous')




#------ Viscous Experiments -------#
datadict_lowcapping = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_exp_58.npy', 
                              allow_pickle=True).item()
datadict_highcapping = np.load('/aos/home/fstdenis/1D-SIM/sim1dplotting/SavedExperiments/datadict_exp_59.npy', 
                              allow_pickle=True).item()

divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, Wdissip_tot, _, _, _ = datadict_lowcapping.values()
dx = datadict_lowcapping['resolution']
dates = datadict_lowcapping['dates']
print(dates)
Nx = np.shape(divergence_tot)[1]
X = np.arange(0, Nx)*dx

norm = colors.Normalize(vmin=0, vmax=dates[-1])  # assuming k ranges from 0 to 10
cmap = plt.cm.viridis
sm = pltcm.ScalarMappable(cmap=cmap, norm=norm)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (11, 8), sharex = True)

axs = [ax1, ax2, ax3, ax4]
c = 0
for i, ax in enumerate(axs):
    
    ax.grid()
 
    if i % 2 == 0:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, Wdissip_tot, dt, dx, dates = datadict_lowcapping.values()
    else:
        divergence_tot, h_tot, A_tot, u_tot, eta_tot, zeta_tot, Wdissip_tot, dt, dx, dates= datadict_highcapping.values()
        
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
ax1.set_title(r'$\nu_{max} = 1.25 \times 10^{8}$ Nsm$^{{-1}}$')
ax2.set_title(r'$\nu_{max} = 1.25 \times 10^{12}$ Nsm$^{{-1}}$')
fig.colorbar(sm, ax=axs, label = 'Time (hr)')
ax3.set_ylabel(r'$u$ (m/s)')
fig.savefig('viscous_numerical.png', dpi = 500)

