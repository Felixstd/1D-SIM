import numpy as np


dt = 15
N= 100
e = 2 
h0 = 1
A0 = 1
C_s = 20
Pstar = 27.5e3
Pstarstar = 27.5e3*np.exp(-C_s*(1-A0))
rhoi = 900
T = 0.36*dt
alpha = 1/((1+e**(-2))*T)
alpha_e = (1+e**(-2))
lamda_div =(np.sqrt(1+e**(-2))-1)/2
lamda_conv = (-np.sqrt(1+e**(-2))-1)/2
rho = 900
nu_max = 1e12
P0 = Pstarstar*h0