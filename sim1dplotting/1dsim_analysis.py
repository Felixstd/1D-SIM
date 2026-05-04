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
	Dates_Config = TimeUtility(configuration_time = config_exp['Time'], 
							  configuration_fig  = config_exp['Figures'])
	tstep, dates_time = (TimeUtility.read_time(Dates_Config, expno = exp))
	print(tstep)
	tstep = tstep[1:]
	if not os.path.isdir(Parameters.figdir+str(exp)):
		os.mkdir(Parameters.figdir+str(exp))
	
	if Parameters.read_all: 
		#---------- READING DATA ----------#

		if Parameters.mechanical_energy:
			datadict, datadict_energy = read_data.read_data(
				exp, 
				int(Parameters.dt_read), 
				int(Parameters.dx[ind]/1e3), 
				Parameters.solv, 
				Parameters.imex, 
				Parameters.adv, 
				tstep, 
				Parameters.outputdir, 
				Dissipation=Parameters.dissipation, 
				Energy=Parameters.mechanical_energy, 
				GC = Parameters.GC
			)

			datadict_energy_Tavg = read_data.read_avg_energy(
       			exp, 
          		int(Parameters.dt_read), 
            	int(Parameters.dx[ind]/1e3), 
             	Parameters.solv, 
              	Parameters.imex, 
               	Parameters.adv, 
                Parameters.nstep, 
                Parameters.outputdir)
   
			number_viscous, number_plastic = read_data.read_visc_plas_file('outputs_VP/num_visc_plas_{}.out'.format(exp))
   
			number_mixed_viszeta, number_mixed_viseta = read_data.read_visc_plas_file('outputs_VP/num_mixed_visc_plas_{}.out'.format(exp))
		else:
			datadict = read_data.read_data(
				exp, 
				int(Parameters.dt_read), 
				int(Parameters.dx[ind]/1e3), 
				Parameters.solv, 
				Parameters.imex, 
				Parameters.adv, 
				tstep, 
				Parameters.outputdir, 
				Dissipation=Parameters.dissipation
			)

		if Parameters.dx[ind] < 1:
			if Parameters.Nx[ind] == 627:
				Parameters.dx[ind] = 800
			if Parameters.Nx[ind]== 802:
				Parameters.dx[ind] = 625
			if Parameters.Nx[ind] == 202:
				Parameters.dx[ind]= 2500
    
			if Parameters.Nx[ind] == 1002:
				Parameters.dx[ind]= 500
			if Parameters.Nx[ind] == 2002:
				Parameters.dx[ind] = 250
			if Parameters.Nx[ind] == 4002:
				Parameters.dx[ind] = 125
			if Parameters.Nx[ind] == 5002:
				Parameters.dx[ind] = 100
			if Parameters.Nx[ind] == 8002:
				Parameters.dx[ind] = 62.5
    
		x_exp = np.linspace(0, Parameters.Nx[ind]*Parameters.dx[ind], int(Parameters.Nx[ind]))	
    
		if Parameters.dissipation and Parameters.GC:
			divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, Wdissip_tot, P_tot, muI_tot = datadict.values()
		elif Parameters.GC and Parameters.dissipation == False:
			divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, P_tot, muI_tot = datadict.values()
		elif Parameters.GC == False and Parameters.dissipation == True:
			divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot, Wdissip_tot= datadict.values()
		else:
			divergence_tot, h_tot, A_tot, u_tot, sigma_tot, sigma_norm_tot, eta_tot, zeta_tot = datadict.values()
		
		plot.plot_initial_conditions(A_tot, h_tot, u_tot, x_exp/1e3, Parameters.figdir+str(exp)+'/', exp)
  
		A_exps.append(A_tot)
		x_exps.append(x_exp)


		if Parameters.savevar:
			datadict_saved = datadict.copy()
			datadict_saved['timestep'] = float(Parameters.dt)
			datadict_saved['resolution'] = float(Parameters.dx[ind]/1e3)
			datadict_saved['dates'] = dates_time
			datadict_saved['Number_Regimes'] = [number_viscous, number_plastic, number_mixed_viszeta, number_mixed_viseta]
			np.save('SavedExperiments/datadict_exp_{}.npy'.format(exp), datadict_saved)
			datadict_energy_saved = datadict_energy_Tavg.copy()
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
							np.array(A_tot), 
							r'$A$', 
							colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), 
							cm.cm.ice, 
							exp,
							Parameters.figdir+str(exp)+'/', 
							'A')
		plot.plot_variable(Parameters.dx[ind], 
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
                            datadict_energy, 
                            Parameters, 
                            Parameters.figdir+str(exp)+'/', 
                            exp
                            )
			plot.plot_energy(datadict_energy_Tavg,number_plastic, number_viscous, number_mixed_viszeta, number_mixed_viseta, Parameters, Parameters.figdir+str(exp)+'/')
			

   

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
