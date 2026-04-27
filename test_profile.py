# import numpy as np
# import matplotlib.pyplot as plt

# # ! Parameters
# # real(8), parameter :: L = 1d0   ! period in h
# # real(8) :: phi, w

# # phi = modulo(h(i), L) / L   ! normalized phase in [0,1)
# # w   = sin(pi*phi)           ! smooth bump: 0→1→0

# # utp(i) = ((i - nx/2d0)/1d5) * w

# #             if (nx .eq. 500) then
# #                A(i) = min(max(tanh((real(i)-real(nx)/3)/0.0000001**(-0.1)) &
# #                         - tanh((real(i)-2*real(nx)/3)/0.0000001**(-0.1)), 0d0)/2, 1d0)
# #                if (A(i) .gt. 0d0) then  
# #                   h(i) = A(i)
                  
                  
# # h = np.zeros(500)
# x = np.arange(0, 500)

# h = (np.tanh((x-500/3)/0.0000001**(-0.1)) - np.tanh((x-2*500/3)/0.0000001**(-0.1)))/2
# # phi = np.mod(h, 1)
# # w = np.sin(np.pi*phi)

# # # t = h(i) / H                     ! t in (0,1)
# # w = np.exp( -1 / (h*(1 - h)) )
# # tri = 2*(np.abs(phi-1/2))
# # # w = 1-np.tanh((tri-1/2)/0.01)
# # print(w)
# # u = (x-250)*w


# # a, b = 500/3, 2*500/3   # interval
# # eps = 0.3          # smoothing width

# # # smooth step using tanh
# # def smooth_step(x, edge, eps, sign=1):
# #     return 0.5*(1 + np.tanh(sign*(x - edge)/eps))

# # # window: turns on at a, off at b
# # w = smooth_step(x, a, eps, +1) * smooth_step(x, b, eps, -1)

# # ramp = x

# # # smooth window [a,b]
# # def smooth_cut(x, a, b, eps):
# #     left  = 0.5*(1 + np.tanh((x - a)/eps))   # turns on at a
# #     right = 0.5*(1 - np.tanh((x - b)/eps))   # turns off at b
# #     return left * right

# # w = smooth_cut(x, a, b, eps)
# # y = ramp * w   # smoothed ramp



# hmin, hmax = 0.00000000001, 1.0   # active window in h
# eps = 0.2               # smoothing width

# # linear ramp in i
# ramp = (x - 500/2) / 1e5
# # ramp = (i - nx/2) / 1e5

# hmin, hmax = 0.00000001, 1.0  # start and end of smoothing
# eps = 0.3              # smoothing width

# # smooth transitions: 0 below hmin, 0 above hmax
# w = np.ones_like(h)
# w = w * 0.5*(1 + np.tanh((h - hmin)/eps))  # smooth turn-on at hmin
# w = w * 0.5*(1 - np.tanh((h - hmax)/eps))  # smooth turn-off at hmax

# # final profile
# utp = ramp * w
# # final profile
# # utp = ramp * w

# i = np.arange(0, 500)
# nx = len(i)

# # linear ramp passing through 0
# ramp = (i - nx/2) / 1e5

# # define half-sine transition region
# left_edge = 300 # start of smooth left transition
# right_edge = 300 # start of smooth right transition

# w = np.ones_like(i, dtype = np.float32)
# # define symmetric window
# fade = int(500)  # number of points for smooth edges on each side
# idx_left = np.arange(fade)
# # w[idx_left] = np.sin((np.pi/2) * idx_left / fade)

# # right half-sine
# idx_right = np.arange(fade)
# # w[-fade + idx_right] = np.sin((np.pi/2) * (fade - idx_right) / fade)

# # additional tanh smoothing on both edges
# eps = 12.0
# x = np.arange(nx)

# # left tanh on top
# # w[:fade] *= 0.5 * (1 + np.tanh((x[:fade] - fade/2)/eps))

# # right tanh on top
# # w[-fade:] *= 0.5 * (1 + np.tanh((-(x[-fade:] - (nx - fade/2)))/eps))

# linear_fraction = 1/100        # central linear region
# fade_fraction = (1 - linear_fraction)/2  # fraction for each edge
# fade = int(fade_fraction * nx)      # points for each edge
# linear_start = fade
# linear_end = nx - fade

# # create window
# w = np.ones(nx)

# # left transition (0 -> 1)
# x_left = np.arange(fade)
# w[:fade] = 0.5 * (1 - np.cos(np.pi * x_left / fade))

# # right transition (1 -> 0)
# x_right = np.arange(fade)
# w[-fade:] = 0.5 * (1 - np.cos(np.pi * (fade - x_right) / fade))

# # central linear region stays 1
# w[linear_start:linear_end] = 1.0

# y = np.linspace(-3, 3, 500)

# nx = 800
# i = np.arange(0, nx)
# ramp = (i - nx/2) / 1e2
# eps = 0.1  # controls sharpness of the transition

# # sigmoid taper: goes from 1 in center to 0 near edges
# window = 1 / (1 + np.exp((np.abs(ramp) - 0.8)/eps))

# h = (np.tanh((i-nx/3)/0.0000001**(-0.1)) - np.tanh((i-2*nx/3)/0.0000001**(-0.1)))/2
# # apply window
# y_smooth = ramp * window/1000

# # # apply window
# # utp = ramp * w
# # # final profile
# # utp = ramp * w
# # # apply window to ramp
# # utp = ramp * w
# # final profile
# # utp = ramp * w*100


import numpy as np
import matplotlib.pyplot as plt

# Parameters
W = 10  # decay width, adjust for matching
x = np.linspace(0, 1000, 1000)

# Double tanh function centered at 500
f_adj = np.tanh((x - (500 - 1000/6))/W) - np.tanh((x - (500 + 1000/6))/W)

# y1, y2, y3 adjusted
y1 = -(x - 500)
y2 = 1 / (1 + np.exp((np.abs(y1) - 1000/6)/W))
y3 = y1 * y2 / 100

# Plot
plt.figure(figsize=(10,6))
plt.plot(x, f_adj, label='Double tanh')
plt.plot(x, y2, label='y2 (sigmoid)')
plt.plot(x, y3, label='y3 = y1*y2/100')
plt.xlabel('x')
plt.ylabel('Function value')
plt.title('Adjusted functions centered at 500 with same decay')
plt.legend()
plt.grid(True)
plt.savefig('testprofile.png')

# fig = plt.figure()
# ax = plt.axes()
# ax2 = ax.twinx()
# ax.plot(y_smooth)
# ax.axvline(250, color = 'k')
# # plt.plot(w)
# ax2.plot(h)
# plt.savefig('testprofile.png')