import numpy as np
from sim1dplotting.data import read_data
from sim1dplotting.plotting import plot
from sim1dplotting.analysis import analysis
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cmocean as cm
import warnings

warnings.filterwarnings("ignore")

# For 1 to 13
#--- Dates of the simulation ---#
#---   They are in seconds   ---#
 

dates = np.array([720,
        1440,
        2160,
        2880,
        3600,
        4320,
        5040,
        5760,
        6480])

dates = np.array([  0,
           7,
          14,
          21,
          28,
          35,
          42,
          49,
          56,
          63])

dates = np.array([864,
        1728,
        2592,
        3456,
        4320,
        5184,
        6048,
        6912, 
        7776])

# dates = np.array([ 0,
#        86400,
#       172800,
#       259200,
#       345600,
#       432000,
#       518400,
#       604800,
#       691200,
#       777600,
#       864000])

# dates = np.array([17280,
#        34560,
#        51840,
#        69120,
#        86400,
#       103680,
#       120960,
#       138240,
#       155520,
#       172800])

# dates = np.array([0,
#         2592,
#         5184,
#         7776,
#        10368,
#        12960,
#        15552,
#        18144,
#        20736,
#        23328])
# dates = np.array([0,
#         1728,
#         3456,
#         5184,
#         6912,
#         8640,
#        10368,
#        12096,
#        13824,
#        15552])

# dates = np.array([0,
#          720,
#         1440,
#         2160,
#         2880,
#         3600,
#         4320,
#         5040,
#         5760,
#         6480])

dates = np.array([
        8640,
       17280,
       25920,
       34560,
       43200,
       51840,
       60480,
       69120,
       77760])


#        25920
#        51840
#        77760
#       103680
#       129600
#       155520
#       181440
#       207360
#       233280
#       259200
      

# dates = np.array([1728,
#         3456,
#         5184,
#         6912,
#         8640,
#        10368,
#        12096,
#        13824,
#        15552])

        #    0
        #  864
        # 1728
        # 2592
        # 3456
        # 4320
        # 5184
        # 6048
        # 6912
        # 7776
        # 8640



expno = '01'
outputdir = "/aos/home/fstdenis/1D-SIM/1DSIM/output/"
figdir = '/aos/home/fstdenis/1D-SIM/1DSIM/Experiments/'+expno+'/'
dates = []
with open('../output/time_run.{}'.format(expno)) as f:
        lines = f.readlines()
        for line in lines:
                dates.append(int(line))

dates = np.asarray(dates[:-1])
print(dates)
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


