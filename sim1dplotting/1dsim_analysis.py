import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import matplotlib.cm as pltcm

from sim1dplotting.utils.TimeUtilities import TimeUtility
from sim1dplotting.data.readnamelist import namelist
from sim1dplotting.data import read_data
from sim1dplotting.data import load_data
from sim1dplotting.plotting import plot
from sim1dplotting.analysis import analysis


import configparser
import cmocean as cm
import warnings
import os

colors_list1 = mpl.color_sequences['Dark2']

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

#--------------------------------------------------------------
# Reading the namelist
#--------------------------------------------------------------
config_exp = configparser.ConfigParser()
config_exp.read('./namelistAnalysis')

Parameters = namelist(configuration_exp  = config_exp['Experiment'], 
                      configuration_rheo = config_exp['Rheology'], 
                      configuration_fig  = config_exp['Figures'], 
                      configuration_time = config_exp['Time'])

Dates_Config = TimeUtility(configuration_time = config_exp['Time'], 
						   configuration_fig  = config_exp['Figures'])

warnings.filterwarnings("ignore")

# For 1 to 13
#--- Dates of the simulation ---#
#---   They are in seconds   ---#

b_exp = []
A_exps = []
x_exps = []
u_exp = []
zeta_exp = []
eta_exp = []
maxvel_exp = []

#--------------------------------------------------------------
# Looping through the experiments
#--------------------------------------------------------------

for ind, exp in enumerate(Parameters.expno): 
	print('Reading Experiment', exp)
 
 
#--------------------------------------------------------------
# Reading the time for exp
#--------------------------------------------------------------
	print('Reading time')
 
	tstep, dates_time = (TimeUtility.read_time(Dates_Config, expno = exp))
 
	if not os.path.isdir(Parameters.figdir+str(exp)):
		os.mkdir(Parameters.figdir+str(exp))
	
	if Parameters.read_all: 
		 
		data = load_data.load_experiment_data(Parameters, exp, ind, tstep)


		if Parameters.savevar:
			datadict_saved = data['base'].copy()
			datadict_saved['timestep'] = float(Parameters.dt)
			datadict_saved['resolution'] = float(Parameters.dx[ind]/1e3)
			datadict_saved['dates'] = dates_time
			# datadict_saved['Number_Regimes'] = [number_viscous, number_plastic, number_mixed_viszeta, number_mixed_viseta]
			np.save('SavedExperiments/datadict_exp_{}.npy'.format(exp), datadict_saved)
			datadict_energy_saved = data['energy_avg'].copy()
			np.save('SavedExperiments/datadict_energy_saved_{}.npy'.format(exp), datadict_energy_saved)

	if Parameters.maxvelocities:
		file = "/aos/home/fstdenis/1D-SIM/output_post_files/min_max_vel_{}.out".format(exp)
		max_velocities, min_velocities, tstep = read_data.read_maxvelocities(file)
		print(len(max_velocities))
		maxvel_exp.append(max_velocities)
		
		datadict_saved = {'max': max_velocities, 'min': min_velocities, 'tstep': tstep}
		np.save('SavedExperiments/datadict_maxvel_{}.npy'.format(exp), datadict_saved)
# 
	if Parameters.plotfields:
		
		plot.plot_variable(Parameters.dx[ind], 
							Parameters.dt,
							tstep, 
							np.array(data["A"]), 
							r'$A$', 
							colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), 
							cm.cm.ice, 
							exp,
							Parameters.figdir+str(exp)+'/', 
							'A')
		plot.plot_variable(Parameters.dx[ind], 
							Parameters.dt,
							tstep, 
							np.array(data["h"]), 
							r'$h$ (m)', 
							colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), 
							cm.cm.ice, 
							exp,
							Parameters.figdir+str(exp)+'/', 
							'h')


#---------- Plotting ----------#
	if Parameters.plotsingle:
     
		plot.plot_initial_conditions(data["A"], data["h"], data["u"], data["x"]/1e3, Parameters.figdir+str(exp)+'/', exp)
    
		plot.plot_individual_time(tstep, 
                            		exp, 
                              		data["base"], 
                                	Parameters.dx[ind],
                                 	Parameters.dt,
									Parameters.figdir+str(exp)+'/', 
         							0, 
                					0,
                     				Parameters,
                         			Dissipation = Parameters.dissipation)

		if Parameters.mechanical_energy:
			plot.plot_mechanical_energy(
                            tstep, 
                            data["energy"], 
                            Parameters, 
                            Parameters.figdir+str(exp)+'/', 
                            exp
                            )
			# plot.plot_energy(datadict_energy_Tavg,number_plastic, number_viscous, number_mixed_viszeta, number_mixed_viseta, Parameters, Parameters.figdir+str(exp)+'/')
			

   

if Parameters.maxvelocities:
	plot.plot_velocities(maxvel_exp, Parameters.maxcapping, [mpl.colormaps['Set1'].colors[2],mpl.colormaps['Set1'].colors[0]] , Parameters.dt, tstep, exp, './', 'viscous')


if len(Parameters.expno) > 1:

	#------------------------------
	# Theoretical profile 
	#------------------------------
 
	A0 = 1
	C = 20
	ell = 2
	ell_2 = np.sqrt(1+ell**(-2))
	rho = 900
	Pstar = 2.75e4

 
	omega_sqrt = np.sqrt(Pstar*np.exp(-C*(1-A0))*(1+C*A0)*(ell_2-1)/rho)
	time = np.linspace(0, dates_time[1]*3600, 1000)
 
    
	plt.figure(1)
	print(dates_time)
	for i, exp in enumerate(Parameters.expno):
     
		A_tot = A_exps[i]
		x = x_exps[i]
		idx_mid = np.where((x < 275e3) & (x > 225e3))[0]
		A_min = np.zeros(len(A_tot))
  
		omega = omega_sqrt*(1/Parameters.dx[i])
		A_theo = A0*np.exp(-omega*time)
		
	
		for j,A in enumerate(A_tot):
			A_mid = A[idx_mid]
			A_min[j] = np.min(A_mid)
		print(A_min)
		plt.plot(dates_time[1:], A_min, color = colors_list1[i], label = r'$\Delta x = {}$ m'.format(int(Parameters.dx[i])))
		# plt.plot(time/3600, A_theo, color = colors_list1[i], linestyle = '--')
	plt.legend()
	plt.xlabel('Time (hr)')
	plt.ylabel(r'min($A$)')
	plt.savefig('min_A_exps.png')
	plt.close(1)
