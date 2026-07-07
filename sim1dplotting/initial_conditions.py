import numpy as np
import matplotlib.pyplot as plt

nx = 500
Deltax = 1000
sharpness = 0.0000001**(-0.1)
u = np.zeros(nx)
u2 = np.zeros(nx)
h = np.zeros(nx)
A = np.zeros(nx)
A_2 = np.zeros(nx)
A_3 = np.zeros(nx)

sharpness2 = 1 / 20000  


# sharpness2 = 1 
# window = 0.5d0 * (tanh(sharpness * (ramp + L/2)) - tanh(sharpness * (ramp - L/2)))

for i in range(1, nx):
    
    i_center1 = nx/2+1
    i_center2 = nx/2        
    half_width = (2*nx/3 - nx/3)/2 

    val = (np.tanh((i - nx/3) / sharpness) 
         - np.tanh((i - (2*nx/3+1)) / sharpness))
    A[i] = min(max(val, 0.0) / 2, 1.0)
    
    A_2[i] = (np.tanh((i - (i_center1 - half_width) )/ sharpness)
         - np.tanh((i - (i_center1 + half_width)) / sharpness))/2
    
    A_3[i] = (np.tanh((i - (i_center2 - half_width) )/ sharpness)
         - np.tanh((i - (i_center2 + half_width)) / sharpness))/2
    
    # if A[i] < 1e-4:
    #     A[i] = 0
    if A[i] > 0.0:
        h[i] = A[i]

for i in range(2, nx):    
    L = nx* Deltax
    x = (i+1) * Deltax
    ramp = (x - L/2)
    window1 = 1 / (np.exp((abs(-ramp))/10000))
    window = 0.5 * (np.tanh(sharpness2 * (ramp + 80e3)) - np.tanh(sharpness2 * (ramp - 80e3)))
    u[i] = (ramp/1e5)*window
    # u2[i] = ramp/0.65e6*window1
    
    plateau = 40e3             
    window = 1/2 * (np.tanh((ramp + 80e3)*sharpness2) - np.tanh((ramp - 80e3)*sharpness2))
    u2[i] = np.sign(ramp)*np.maximum(np.abs(ramp) - plateau, 0) / 1e5 * window
    
    # ramp = (x - L/2) / 1e2
    # window = 1 / (1 + np.exp((abs(ramp) - 800)/100))
    # u2[i] = ramp*window/1000000
    


half_width = (nx - 1) /6
i_lo = nx/2 - half_width+1
i_hi = nx/2 + half_width    

sharpness = 1e-1    

# for i in range(2, nx):
    # A_2[i] = (np.tanh(sharpness * (i - i_lo)) - np.tanh(sharpness * (i - i_hi))) / 2
    # if A_2[i] < 1e-4:
    #     A_2[i] = 0
    # A_2[i] = min(max(A_2[i], 0.0), 1.0)
    
print(len(np.where(A_2> 0)[0]),len(np.where(A> 0)[0]))
# print(A_2)




print(u[249],u[250],u[251])
print(len(np.where(A[:251] > 0)[0]),len(np.where(A[251:] > 0)[0]))
plt.figure(figsize = (7,5))
ax1 = plt.axes()
# ax = plt.twinx()
ax1.plot(u+1, color  ='b', label = 'tanh')
# plt.axvline(250)
# ax.plot(u2, label = 'exp')
ax1.plot(A,color = 'r', label = 'def', alpha = 0.5)
# ax1.plot(A_2,color = 'g',/ label = 'nx+1', alpha = 0.5)
# ax1.plot(A_3,color = 'b', label = 'nx', alpha = 0.7)
plt.xlabel('x (km)')
plt.legend()
ax1.legend()
ax1.set_ylabel('A,h')
# ax.set_ylabel('u')
plt.savefig('initcond.png',dpi = 500)


arr1 = np.array([
    0.0000000000000000,
    0.19987156275717063,
    0.25444947738507540,
    0.30791237081489131,
    0.35137383318556337,
    0.37592365070271572,
    0.37592365070271572,
    0.35137383318556337,
    0.30791237081489131,
    0.25444947738507540,
    0.19987156275717063,
    0.0000000000000000
])
x = np.arange(0.5,10.5, 1)
arr2 = np.array([
    0.0000000000000000,
   -8.9107621629081786e-06,
   -7.7078734446758588e-06,
   -5.2802682880658960e-06,
   -1.8835744242790017e-06,
    1.8835744242790017e-06,
    5.2802682880658960e-06,
    7.7078734446758588e-06,
    8.9107621629081786e-06,
    0.0000000000000000,
    0.0000000000000000
])

plt.figure()
plt.plot(x,arr1[1:-1], color ='r')
plt.plot(arr2*1e6,color ='b')
plt.savefig('testinitcond.png')


arr_u = np.array([
    0.0000000000000000,
   -9.0741128951386061e-03,
   -8.4632825359629120e-03,
   -6.6403677026784912e-03,
   -3.6721244361146966e-03,
    0.0000000000000000,
    3.6721244361146966e-03,
    6.6403677026784912e-03,
    8.4632825359629120e-03,
    9.0741128951386061e-03,
    0.0000000000000000
])
x_a = np.arange(1.5,11.5, 1)
x_u = np.arange(1,12, 1)

arr_a = np.array([
    0.0000000000000000,
    0.19987156275717063,
    0.25444947738507540,
    0.30791237081489131,
    0.35137383318556337,
    0.37592365070271572,
    0.37592365070271572,
    0.35137383318556337,
    0.30791237081489131,
    0.25444947738507540,
    0.19987156275717063,
    0.0000000000000000
])
plt.figure()
plt.plot(x_a,arr_a[1:-1], color ='r')
plt.plot(x_u,arr_u+0.35,color ='b')
plt.savefig('testinitcond.png')