import numpy as np
import matplotlib.pyplot as plt

Nx = 400
x = np.arange(0, Nx+1,1)
hinit = np.exp(-((x-Nx/2)**2/2500)**5)
idx = np.where(hinit>0.2)[0]
uinit = x/2e3*(np.tanh((x - idx[0])/10) - np.tanh((x - idx[-1])/10))

left = -x /1e3 * (np.tanh((x - idx[0])/1e1) - np.tanh((x - idx[int(len(idx)/2)])/1e1))
right = x /1e3 * (np.tanh((x - idx[int(len(idx)/2)])/1e1) - np.tanh((x - idx[-1])/1e1))
uinit = left + right

uinit = np.cos(np.pi*x/(2*Nx))
uinit[hinit < 0.2] = 0


x = np.linspace(-3, 3, 100)
u = np.zeros_like(x)

u[abs(x) < 1]= x[abs(x) < 1]
u[abs(x) > 1] = 0

plt.figure()
plt.plot(uinit)
plt.plot(hinit)
plt.savefig('test.png')

import numpy as np
import matplotlib.pyplot as plt

# Grid
Nx = 200
x = np.arange(Nx+1)

# hinit is a sharp step
hinit = np.exp(-((x - Nx/2)**2 / 5000)**5)

# Raw u profile
u_raw = np.linspace(-2.5, 2.5, Nx+1)

# Hard cutoff
u_hard = u_raw.copy()
u_hard[hinit < 0.2] = 0

# 1) Tanh smoothing
width = 0.06
mask_tanh = 0.1 * (1 + np.tanh((hinit - 0.2) / width))
u_tanh = u_raw * mask_tanh
u_tanh[hinit < 1e-2] = 0

# 2) Smoothstep
def smoothstep(x, edge0, edge1):
    t = np.clip((x - edge0) / (edge1 - edge0), 0, 1)
    return t * t * (3 - 2 * t)

mask_smooth = smoothstep(hinit, 0.15, 0.25)
u_smooth = u_raw * mask_smooth

# 3) Gaussian decay
sigma = 0.05
mask_gauss = np.exp(-((hinit - 0.2)**2) / (2 * sigma**2))
u_gauss = u_raw * mask_gauss


print(hinit, u_tanh)
# Plot
plt.figure(figsize=(8,5))
# plt.plot(x, u_hard, label='Hard cutoff', color='k', lw=2)
plt.plot(x, u_tanh, label='Tanh smooth', lw=2)
# plt.plot(x, u_smooth, label='Smoothstep', lw=2)
# plt.plot(x, u_gauss, label='Gaussian', lw=2)
plt.xlabel('x')
plt.ylabel('uinit')
plt.legend()
plt.grid(True)
plt.title('Smoothing uinit based on hinit step profile')
plt.show()


plt.savefig('test.png')