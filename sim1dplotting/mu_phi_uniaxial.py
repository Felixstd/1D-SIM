import numpy as np
from sim1dplotting.data import read_data
from sim1dplotting.plotting import plot
from sim1dplotting.analysis import analysis
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean as cm
import warnings

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--expno', dest='expno', type=str, help='Experiment number')
args = parser.parse_args()
expno = str(args.expno)

warnings.filterwarnings("ignore")

# For 1 to 13
#--- Dates of the simulation ---#
#---   They are in seconds   ---#

outputdir = "/aos/home/fstdenis/1D-SIM/output/"
figdir = '/aos/home/fstdenis/1D-SIM/Experiments/'+expno+'/'

dates = []
with open('../output/time_run.{}'.format(expno)) as f:
        lines = f.readlines()
        for line in lines:
                dates.append(int(line))

dates = np.asarray(dates[:-1])
print('Time steps', dates.T)


dt = 1
time = dates*dt/(60*60)

dx = 10
solv = 1
imex  = 0
adv = 1


Nx = 100

N_transect = 101
muphi = 0


#---------- READING DATA ----------#
datadict = read_data.read_data(expno, dt, dx, solv, imex, adv, dates, outputdir, MuPhi = muphi)
divergence_tot, h_tot, A_tot, u_tot = datadict.values()

nx = np.shape(A_tot[0])[0]
X = np.linspace(0, nx*dx, nx)

# 
plot.plot_variable(dx, time, np.array(A_tot), r'$A$', colors.SymLogNorm(vmin=0, vmax=1, linthresh=0.1), cm.cm.ice, expno,figdir, 'A')

#---------- Analysing Wind Forcing----------#
# analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

#---------- Analysing Invariants ----------#
# analysis.invariants(datadict, N_transect, time, figdir+expno+'/', muphi)

#---------- Plotting ----------#
plot.plot_individual_time(dates, expno, datadict, dx, figdir, 0, 0, MuPhi = muphi)

#----- Initial Conditions ------#
A_init = A_tot[0]
h_init = h_tot[0]
u_init = u_tot[0]

plt.figure()
ax1 = plt.axes()
ax2 = ax1.twinx()
ax1.plot(X, A_init, color = 'b')
ax1.set_ylabel(r'$A, \; h$ (m)',color = 'k' )
ax2.plot(X[1:], u_init, linestyle = '--', color = 'k')
ax2.set_ylabel(r'$u$ (m/s)', color='b')
ax2.tick_params(axis='y', labelcolor='b')
plt.savefig('init_conditions.png')

