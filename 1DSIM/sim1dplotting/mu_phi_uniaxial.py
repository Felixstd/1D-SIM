import numpy as np
from sim1dplotting.data import read_data
from sim1dplotting.plotting import plot
from sim1dplotting.analysis import analysis
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean as cm
import warnings

warnings.filterwarnings("ignore")


expno = '01'
outputdir = "/aos/home/fstdenis/1D-SIM/1DSIM/output/"
figdir = '/aos/home/fstdenis/1D-SIM/1DSIM/Experiments/'+expno+'/'
dates = []
with open('../output/time_run.{}'.format(expno)) as f:
        lines = f.readlines()
        for line in lines:
                dates.append(int(line))

dates = np.asarray(dates[:-1])

dt = 10
time = dates*dt/(60*60)

dx = 10
solv = 1
imex  = 0
adv = 1


Nx = 100

N_transect = 101
muphi = 0


# mu_0 = 0.1

# mu_infty = 0.9

#---------- READING DATA ----------#
datadict = read_data.read_data(expno, dt, dx, solv, imex, adv, dates, outputdir, MuPhi = muphi)
h_tot, A_tot, u_tot = datadict.values()




plot.plot_variable(dx, time, np.array(A_tot), r'$A$ (m)', colors.Normalize(vmin=0, vmax=1e-7), cm.cm.ice, expno,figdir, 'h')

#---------- Analysing Wind Forcing----------#
# analysis.wind_forcing(datadict, N_transect, dy, Ny, time, figdir+expno+'/', muphi)

#---------- Analysing Invariants ----------#
# analysis.invariants(datadict, N_transect, time, figdir+expno+'/', muphi)

#---------- Plotting ----------#
plot.plot_individual_time(dates, expno, datadict, dx, figdir, 0, 0, MuPhi = muphi)


