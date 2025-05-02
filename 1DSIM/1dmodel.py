import numpy as np 
import matplotlib.pyplot as plt

# Solving for the 1d advection equation using the FCT method


delta_x = 1e-2
L = 3
N = int(2*L/delta_x)
delta_t = 0.01
T = 10
Nt = int(T/delta_t)
x = np.linspace(-L, L, N)
print(x)

u = np.zeros(N)
A = np.zeros((Nt, N))

#initial condition
A[0, abs(x)<1] = (1-x[abs(x)<1]**2)**(1/20)
u[0] = 0
u[abs(x)<1] = x[abs(x)<1]

u[-1] = 0 
u[0] = 0 
    
u_half_pos = (u[1:] + u[:-1]) / 2

u_half_pos = np.zeros(N-1)
for i in range(0, N-1):
    u_half_pos[i] = (u[i+1] + u[i])/2

u_half_neg = np.zeros(N-1)
for i in range(1, N):
    u_half_neg[i-1] = (u[i] + u[i-1])/2
    



#half step


    







