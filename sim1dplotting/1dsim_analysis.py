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
zeta_exp = []
eta_exp = []
maxvel_exp = []
data_exps = []

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
	# print(tstep)
	if not os.path.isdir(Parameters.figdir+str(exp)):
		os.mkdir(Parameters.figdir+str(exp))

#--------------------------------------------------------------
# Reading data for the experiment
#--------------------------------------------------------------
	if Parameters.read_all: 
		 
		data = load_data.load_experiment_data(Parameters, Parameters.solv[ind],exp, ind, tstep)
		data_exps.append(data)

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


plt.close()

# plot.plot_E_VP_resolution(tstep, data_exps, Parameters, 
#                            './', './VP_EVP_res.png')
# plot.plot_laplacian_p(data_exps, Parameters)
plot.plot_models_resolution(tstep, data_exps, Parameters, 
                           './', './VP_EVP_res_u.png',var = 'u')
    
analysis_regul  = False
analysis_plot = False
ridge_building = False

if (len(Parameters.expno) >= 1):
    
    
    
	# plot.plot_figures_paper(tstep, 
	# 						data_exps[0]["base"], data_exps[1]["base"], 
	# 						Parameters.dx[0], Parameters.dx[1], Parameters.dt, 
	# 						'./', 'nVP_resolution_omega1.5.png', 
    #    						highres_3=True, 
    #          				datadict_hr_highres=data_exps[2]["base"], 
    #              			dx_hr2=Parameters.dx[2]	)
		
	if (Parameters.helmP):  
      
		plot.plot_sensitivity_lc(Parameters,data_exps,
						'Figures/Figs_sens_lscale/sensitivity_lc.png', number_plots=2)
  
	if analysis_regul:

		plot.plot_sensitivity_regul(tstep,Parameters,data_exps)
  
	if ridge_building: 
		plot.plot_ridgebuilding(Parameters,data_exps)
