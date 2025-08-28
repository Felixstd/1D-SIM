#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')
def omega_viscous_constant(k, zeta_max):

    a = rhoice*h0
    
    b = zeta_max*k**2
    
    c = Pstar*k**2*h0+0j
    
    omega_plus = 1/(2*a)*(-b + np.sqrt(b**2 - 4*a*c))
    omega_minus = 1/(2*a)*(-b - np.sqrt(b**2 - 4*c*a))
    
    max_omega = h0*Pstar/zeta_max

    return omega_plus, omega_minus, max_omega

def omega_plastic_VP(k, u0, h0, e):
    
    gamma_div = np.sqrt(1+e**(-2))/2 - 1/2
    gamma_conv = np.sqrt(1+e**(-2))/2 - 1/2
    
    a = rhoice*u0
    
    b = -1J*rhoice*h0*u0*k
    
    c_div = -h0*k**2*Pstarstar*gamma_div
    
    omega_min = (-b-np.sqrt(b**2-4*a*c_div))/(2*a)
    
    return omega_min

def omega_viscous_dhdx_dAdx(k, dh0dx, dA0dx, du0dx, u0, zeta_max):

    a = rhoice*h0
    
    b = -(zeta_max*k**2 + rhoice*h0*du0dx + 1j*u0*rhoice*h0*k)
    
    c = zeta_max*du0dx*k**2+Pstarstar*(h0*k**2 - 2*C*dh0dx*dA0dx + C*h0*A0*k**2 + C**2*h0*(dA0dx)**2) + \
        1j*(u0*zeta_max*k**3 - Pstarstar*(2*C*h0*dA0dx*k +k*dh0dx+ C*A0*dh0dx*k + C**2*A0*h0*k*dA0dx))
    
    omega_plus = 1/(2*a)*(-b + np.sqrt(b**2 - 4*a*c))
    omega_minus = 1/(2*a)*(-b - np.sqrt(b**2 - 4*a*c))
    
    
    omega_3 = du0dx + 1j*u0*k
    
    max = h0/zeta_max*(Pstarstar + A0*C*Pstarstar - 2*rhoice*u0**2) + du0dx
    
    return omega_plus, omega_minus, omega_3, max

def omega_viscous_dhdx(k, dh0dx, du0dx, u0, zeta_max):
    a = rhoice*h0
    
    b = -(zeta_max*k**2 + rhoice*h0*du0dx + 1j*u0*rhoice*h0*k)
    
    c = zeta_max*du0dx*k**2-Pstarstar*(h0*k**2) + \
        1j*(u0*zeta_max*k**3 + Pstarstar* k*dh0dx)
    
    omega_plus = 1/(2*a)*(-b + np.sqrt(b**2 - 4*a*c))
    omega_minus = 1/(2*a)*(-b - np.sqrt(b**2 - 4*a*c))
    
    omega_max = -h0/zeta_max*(Pstarstar) + du0dx
    
    return omega_plus, omega_minus, omega_max
    

def frequency_mub(k, mub, dA0dx, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    
    Pstarstar = Pstar*np.exp(-C*(1-A0))
    mub_conv = -mub-1
    mub_div = mub-1
    lambdastar_div = Pstarstar*mub_div*C/(rhoice*h0)
    lambdastar_conv = Pstarstar*mub_conv*C/(rhoice*h0)
    
    omega_plus_div = np.sqrt(lambdastar_div*(dA0dx**2*C - A0*k**2) + 1j*k*dA0dx*(1+C*A0))
    omega_minus_div = -np.sqrt(lambdastar_div*(dA0dx**2*C - A0*k**2) + 1j*k*dA0dx*(1+C*A0))
    
    omega_plus_conv= np.sqrt(lambdastar_conv*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    omega_minus_conv = -np.sqrt(lambdastar_conv*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    
    return omega_plus_div, omega_minus_div, omega_plus_conv, omega_minus_conv

zeta_max_tot = np.array([1e8, 1e9, 1e10, 1e11, 1e12])
e = 0
rhoice = 900
h0, A0, C= 1, 1, 20


Pstar = 27.5e3
Pstarstar = Pstar*np.exp(-C*(1-A0))

if e < 1:
    zeta_max_e_tot = zeta_max_tot
else:
    zeta_max_e_tot = zeta_max_tot*(1+e**(-2))

dh0dx_arr = np.power(10.0, np.arange(-8, -1, 1))
dh0dx_arr = np.concatenate([np.array([0.0]), dh0dx_arr])
dA0dx_arr = dh0dx_arr
print(len(dA0dx_arr))

k = np.linspace(1e-7, 1e0, int(1e7))

zeta_max_0 = 1e12
# omega_plus_0, omega_minus_0, argument0 = omega_viscous_dhdx(k, 0, 1e12)

#%%
fig = plt.figure(figsize = (5, 4))
for i, zeta_max_e in enumerate(zeta_max_e_tot):
    zeta_max = zeta_max_tot[i]
    # omega_plus_0, omega_minus_0, argument0 = omega_viscous_dhdx(k, 0, zeta_max_e)
    mantissa, exponent = f"{zeta_max:.1e}".split('e')
    
    label = r"$\zeta_{{max}} = {} \times 10^{{{}}}$".format(mantissa, int(exponent))
    # plt.plot(k, omega_minus_0, label = label)
    
fig.legend(loc='upper center', 
            labelcolor='linecolor',  
            bbox_to_anchor=(1.06, 0.9))
# plt.axis('equal')
plt.grid()
plt.xscale('log')
plt.ylabel(r'$Re\{\omega_-\}$ (s$^{-1}$)')
plt.xlabel(r'$k$ (1/m)')
plt.title(r'$\partial h_0 / \partial x = 0$, $e = {}$'.format(e))
plt.savefig('omega_vp_e{}_dhdx0.png'.format(e))

# fig , (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (6, 4))
#%%

dx = 1e4

dh0dx = 1/dx
dA0dx = 1/dx
u0 = 0.1
du0dx =u0/dx


# fig = plt.figure(figsize = (8, 4))
# ax = plt.axes()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (7, 5), sharex = True)  
fig2, (ax3) = plt.subplots(1, 1, figsize = (8, 4), sharex = True)  
zeta_max_tot = np.array([1.25e8, 1.25e9, 1.25e12])
colors = ['darkred', 'forestgreen', 'royalblue']

colors = mpl.colormaps['Set1'].colors
# zeta_max_tot = np.array([1e8])
# colors = ['k']
labels = []
for i, zeta_max in enumerate(zeta_max_tot):
    
    # omega_plus_dhdx_dadx, omega_minus_dhdx_dadx, omega_dudx, max_omega =  omega_viscous_dhdx_dAdx(k, dh0dx, 0, 0, 0, zeta_max)
    # omega_plus_dhdx, omega_minus_dhdx, max_omega =  omega_viscous_dhdx(k, dh0dx, du0dx, u0, zeta_max)
    # omega_plus_dhdx_negdu0dx, omega_minus_dhdx_negdu0dx, max_omega_negdu0dx =  omega_viscous_dhdx(k, dh0dx, -du0dx, u0, zeta_max)
    omega_plus_constant, omega_minus_constant, max_omega = omega_viscous_constant(k, zeta_max)
    
    mantissa, exponent = f"{zeta_max:.2e}".split('e')
    label = r"$\nu_{{max}} = {} \times 10^{{{}}}$ Nsm$^{{-1}}$".format(mantissa, int(exponent))
    # p = ax.plot(k, omega_minus_dhdx, color = colors[i], label = label)
    p = ax1.plot(k, omega_plus_constant, color = colors[i], label = label)
    p = ax2.plot(k, np.imag(omega_plus_constant), color = colors[i])
    
    P = ax3.plot(k, omega_plus_constant, color = colors[i], label = label, linewidth = 2)
    # p = ax3.plot(k, omega_plus_constant, color = colors[i])
    
    # p = ax.plot(k, omega_minus_dhdx_negdu0dx, color = colors[i], linestyle = '--', label = label)
    # ax.plot(k, omega_minus_dhdx, color ='b', linestyle = '--', label = label)
    labels.append(mlines.Line2D([], [], color = p[0].get_color(), marker='None', linestyle='None',
                            label=label))
    # ax.axhline(max_omega)
    
for ax in [ax1, ax2]:
    
    ax.grid()
    ax.set_xscale('log')
    
    
ax1.set_ylabel(r'$\text{Re}\{\omega_+\}$ (s$^{-1}$)', fontsize = 15)
ax2.set_ylabel(r'$\text{Im}\{\omega_+\}$ (s$^{-1}$)', fontsize = 15)
ax3.set_ylabel(r'$\text{Re}\{\omega_+\}$ (s$^{-1}$)', fontsize = 15)
# ax3.set_ylabel(r'$\text{Im}\{\omega_+\}$ (s$^{-1}$)', fontsize = 15)
fig.supxlabel(r'$k$ (m$^{-1}$)', fontsize = 15)
fig.align_ylabels()
fig.legend(loc='upper center', 
                handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.07, 0.9), fontsize = 15)

ax3.set_xscale('log')
fig2.supxlabel(r'$k$ (m$^{-1}$)', fontsize = 15)
fig2.legend(loc='lower right', 
                handles=labels, labelcolor='linecolor',  bbox_to_anchor=(0.9, 0.12), fontsize = 15)
# ax.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}, dA_0/dx = {} \; \text{{m}}^{{-1}}$'.format(du0dx, dh0dx, dA0dx))
# ax.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}$'.format(du0dx, dh0dx))
fig.savefig('omega_viscous_constant_plus.png')
fig2.savefig('omega_viscous_constant_plus_real.png')

omega_minus_zeta_min_max = []

#%%
#------ Plastic ---------#

omega_min = omega_plastic_VP(k, u0, h0, 2)

fig = plt.figure(figsize = (8, 4))
ax = plt.axes()

ax.plot(k, omega_min, color = 'k')

ax.grid()
ax.set_xscale('log')
ax.set_ylabel(r'$\text{Re}\{\omega_-\}$ (s$^{-1}$)', fontsize = 15)
ax.set_xlabel(r'$k$ (m$^{-1}$)', fontsize = 15)
# ax.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}, dA_0/dx = {} \; \text{{m}}^{{-1}}$'.format(du0dx, dh0dx, dA0dx))
# ax.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}$'.format(du0dx, dh0dx))
fig.savefig('omega_plastic_VP.png')

# for i, zeta_max_e in enumerate(zeta_max_e_tot):
#     zeta_max = zeta_max_tot[i]
#     labels = []
#     print('Plotting zeta_max: ', zeta_max/1e8)
    
    
#     fig , (ax1, ax2) = plt.subplots(2, 1, figsize = (6, 4))
    
#     fig2, (ax3, ax4) = plt.subplots(2, 1, figsize = (6, 4))
    
#     fig3, (ax5, ax6) = plt.subplots(2, 1, figsize = (6, 4))

#     for c, dh0dx in enumerate(dh0dx_arr): 
#         dA0dx = dA0dx_arr[c]
        
#         omega_plus_dhdx, omega_minus_dhdx, argument = omega_viscous_dhdx(k, 0, zeta_max_e)
#         omega_plus_dhdx_dadx, omega_minus_dhdx_dadx =  omega_viscous_dhdx_dAdx(k, dA0dx, dh0dx, zeta_max_e)
#         omega_plus_dadx, omega_minus_dadx =  omega_viscous_dhdx_dAdx(k, 0, dA0dx, zeta_max_e)
        
#         omega_plus_div_mub, omega_minus_div_mub, omega_plus_conv_mub, omega_minus_conv_mub = frequency_mub(k, 10, dA0dx)
#         if ((i < 1) or (i > 3)) and (c > 6 ) : 
#             print('here')
#             omega_minus_zeta_min_max.append(omega_minus_div_mub)
            
#         # break
        
#         if c < 1:  
#     #         p = ax1.plot(k, omega_plus_dhdx, color = 'k', label = 'dhdx = {}'.format(dh0dx))
#     #         ax2.plot(k, omega_minus_dhdx, color = 'k')
            
#             p = ax3.plot(k, omega_plus_dhdx_dadx, color = 'k', label = 'dhdx = dAdx = {}'.format(dh0dx))
#             ax4.plot(k, omega_minus_dhdx_dadx, color = 'k')
            
#     #         p = ax5.plot(k, omega_plus_dadx, color = 'k', label = 'dhdx = dAdx = {}'.format(dh0dx))
#     #         ax6.plot(k, omega_minus_dadx, color = 'k')
#     #         # ax3.plot(k, argument)
#         else:
#     #         p = ax1.plot(k, omega_plus_dhdx, label = 'dhdx = {}'.format(dh0dx))
#     #         ax2.plot(k, omega_minus_dhdx)
            
#             p = ax3.plot(k, omega_plus_dhdx_dadx, label = 'dhdx = dAdx = {}'.format(dh0dx))
#             ax4.plot(k, omega_minus_dhdx_dadx)
            
#     #         p = ax5.plot(k, omega_plus_dadx, label = 'dhdx = dAdx = {}'.format(dh0dx))
#     #         ax6.plot(k, omega_minus_dadx)
#     #         # ax3.plot(k, argument)
            
#         mantissa, exponent = f"{dh0dx:.1e}".split('e')
#         label = r'$\partial h_0/ \partial x = {} \times 10^{{{}}}$'.format(mantissa, int(exponent))
#         labels.append(mlines.Line2D([], [], color = p[0].get_color(), marker='None', linestyle='None',
#                             label=label))

#     # # ax1.plt.plot(k, omega_minus_0, color = 'k', label = )

#     for ax in [ax3, ax4]:
#         ax.grid()
#         ax.set_xscale('log')
        
#     # ax1.set_ylabel(r'$Re\{\omega_+\}$ (s$^{-1}$)')
#     # ax2.set_ylabel(r'$Re\{\omega_-\}$ (s$^{-1}$)')
#     # ax3.grid()
#     # ax3.set_xscale('log')
#     ax3.set_ylabel(r'$Re\{\omega_+\}$ (s$^{-1}$)')
#     # ax4.set_ylabel(r'$Re\{\omega_-\}$ (s$^{-1}$)')
    
#     # ax5.set_ylabel(r'$Re\{\omega_+\}$ (s$^{-1}$)')
#     # ax6.set_ylabel(r'$Re\{\omega_-\}$ (s$^{-1}$)')
#     # ax3.set_ylabel(r'$atan(B/A)$')
#     mantissa, exponent = f"{zeta_max:.1e}".split('e')
#     title = r"$\zeta_{{max}} = {} \times 10^{{{}}}$, $e = {}$".format(mantissa, int(exponent), e)
#     # fig.suptitle(title)
#     # fig.align_ylabels()
    
#     fig2.suptitle(title)
#     fig2.align_ylabels()
    
#     # fig3.suptitle(title)
#     # fig3.align_ylabels()
    
    

#     # fig.legend(loc='upper center', 
#     #             handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.02, 0.9))
#     fig2.legend(loc='upper center', 
#                 handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.02, 0.9))

#     # fig3.legend(loc='upper center', 
#     #             handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.02, 0.9))
    
#     # fig.supxlabel(r'$k$ (m$^{-1}$)')
#     fig2.supxlabel(r'$k$ (m$^{-1}$)')
#     # fig3.supxlabel(r'$k$ (m$^{-1}$)')
    
#     # fig.savefig('omega_vp_zeta{}_e{}.png'.format(exponent, e))
#     fig2.savefig('omega_vp_zeta{}_e{}_dAhdx.png'.format(exponent, e))
#     # fig3.savefig('omega_vp_zeta{}_e{}_dAdx.png'.format(exponent, e))
#     plt.close()
#     plt.clf()
    
# # fig, (ax1, ax2) = plt.subplots(2, 1)

# # for ax in [ax1, ax2]:
# #     ax.grid()
# #     ax.set_xlabel('log')
    
# fig = plt.figure(figsize = (10, 5))
# ax1 = plt.axes()

# labels = [mlines.Line2D([], [], color = 'k', marker='None', linestyle='None',
#                             label=r'$\nu_{max} = 10^8$ (Nsm$^{-1})$'), 
#           mlines.Line2D([], [], color = 'r', marker='None', linestyle='None',
#                             label=r'$\nu_{max} = 10^{12}$ (Nsm$^{-1})$')]

# ax1.plot(k, omega_minus_zeta_min_max[0], color = 'k', label=r'$\nu_{max} = 10^8$ (N s m$^{-1})$')
# ax1.plot(k, omega_minus_zeta_min_max[1], color = 'r', label=r'$\nu_{max} = 10^{12}$ (N s m$^{-1})$')
# ax1.fill_between(k[omega_minus_zeta_min_max[0]< 0], omega_minus_zeta_min_max[0][omega_minus_zeta_min_max[0]< 0], 0, color = 'purple')
# plt.xscale('log')
# fig.legend(handles = labels, loc='upper left', labelcolor='linecolor',
#           bbox_to_anchor = (0.03, 0.8), fontsize = 15)
# # plt.grid()
# plt.yticks(fontsize = 15)
# plt.xticks(fontsize = 15)
# ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
# fig.supxlabel('$k$ (1/m)', fontsize = 15)
# fig.supylabel(r'$\text{Re}\{\omega_-\} (1/s)$', fontsize = 15)
# plt.savefig('viscous_testmub.png')
