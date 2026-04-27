import numpy as np
import matplotlib.pyplot as plt

import scienceplots
import warnings
warnings.filterwarnings('ignore')

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def find_roots(a, b, c, d):
    
    coefficients = np.array([a, b, c, d])

    roots = np.roots(coefficients)
    
    return roots

rhoice = 900
h0 = 1
