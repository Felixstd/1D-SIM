import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import scienceplots
import warnings
warnings.filterwarnings('ignore')

plt.style.use('science')

def friction(I, mu0, muinf, I0 = 1e-3):

    return mu0 + (muinf-mu0)/(I0/(I+1e-12) + 1)

def inertial(dudx, P, h, rho = 900, dmean = 1e3):
    
    return dmean*dudx*np.sqrt(rho*h/P)

def strength(h, A, C = 20, Pstar = 27.5e3):
    
    return Pstar*h*np.exp(-C*(1-A))

def dilantancy(mu, phi = 20*np.pi/180):
    
    return (mu-np.tan(phi))/(1+mu*np.tan(phi))

def sigma11(P, mu, mub, dudx, dilat):
    
    eta = 1e12*np.tanh(mu*P/(abs(dudx+1e-12)*1e12))
    zeta = 2e8*P*np.tanh(mub/(abs(dudx+1e-12)*2e8))
    
    # sig11 = P*(np.sign(dudx)*(mu/2 + mub)-1)
    sig11 = P*(dudx*(eta/2+zeta)-1)
    sig11_dilat = P*(np.sign(dudx)*(mu/2 + mub)-1) - (mub - mu/2)*dilat*P
    return sig11, sig11_dilat


def frequency(k, mu0, muinf, mub, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    
    delta_mu = muinf - mu0
    lamda = Pstar*np.exp(-C*(1-A0))
    Gamma = h0 * np.sqrt(lamda*rhoice)*(dmean*delta_mu/(2*I0))
    gamma_plus = (mu0/2 + mub) -1
    gamma_minus = -(mu0/2 + mub) -1
    
    #Divergence
    omega_plus_d = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1)+0j))/(2*rhoice*h0)
    omega_minus_d = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1)+0j))/(2*rhoice*h0)

    #convergence
    omega_plus_c = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)+0j))/(2*rhoice*h0)
    omega_minus_c = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)+0j))/(2*rhoice*h0)
    return [omega_plus_d, omega_minus_d], [omega_plus_c, omega_minus_c], \
        [-h0*lamda*gamma_plus*(C*A0+1)/Gamma, Gamma*k**2/(2*rhoice*h0)]
        
def frequency_constant_state(k, mu0, muinf, mub, u0, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):

    delta_mu = muinf - mu0
    lamda = Pstar*np.exp(-C*(1-A0))
    Gamma = h0 * np.sqrt(lamda*rhoice)*(dmean*delta_mu/(2*I0))
    gamma_plus = (mu0/2 + mub) -1
    gamma_minus = -2*(mu0/2 + mub) -1
    
    omega_plus_c = 1/(2*rhoice*h0)*(Gamma*k**2 + np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
        4*rhoice*h0*k**2*(gamma_minus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma)))
    omega_minus_c = 1/(2*rhoice*h0)*(Gamma*k**2 - np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
        4*rhoice*h0*k**2*(gamma_minus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma)))
    
    omega_plus_d = 1/(2*rhoice*h0)*(Gamma*k**2 + np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
    4*rhoice*h0*k**2*(gamma_plus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma)))
    omega_minus_d = 1/(2*rhoice*h0)*(Gamma*k**2 - np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
    4*rhoice*h0*k**2*(gamma_plus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma)))
    
    omega_max = (Gamma*k**2)*1/(2*rhoice*h0)
    racine_c = (1/(2*rhoice*h0))*np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
        4*rhoice*h0*k**2*(gamma_minus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma))
    racine_d = (1/(2*rhoice*h0))*np.sqrt((Gamma*k**2 + 1j*rhoice*u0*h0*k)**2 + \
        4*rhoice*h0*k**2*(gamma_plus*lamda*h0*(1+C*A0) - 1j*k*u0*Gamma))
    
    return [omega_plus_d, omega_minus_d], [omega_plus_c, omega_minus_c], omega_max, racine_c, racine_d

k = np.linspace(1e-6, 1, 1000000)
mu0 = .1
muinf = 0.8
u0 = 0.1
# mub = [0.01]
mub = [0.1, 0.5, 1, 10]
dudx = np.linspace(-1/100, 1/100, 100)
A = np.ones_like(dudx)
h = np.ones_like(dudx)

ice_strength = strength(h, A)
Inertial_num = inertial(dudx, ice_strength, h)
friction_coef = friction(Inertial_num, mu0, muinf)
tanpsi = dilantancy(friction_coef)
stress1, s1_d = sigma11(ice_strength, friction_coef, 0.5, dudx, tanpsi)
stress2, s2_d = sigma11(ice_strength, friction_coef, 10, dudx, tanpsi)
stress3, s3_d = sigma11(ice_strength, friction_coef, 0.1, dudx, tanpsi)




colors = ['k', 'g', 'r', 'b']
fig, (ax1, ax2)= plt.subplots(1, 2, sharex = True, figsize = (9, 4))
fig1, (ax3, ax4)= plt.subplots(1, 2, sharex = True, figsize = (9, 4))
labels = []
for i, b in enumerate(mub):
    omega_divergence, omega_convergence, max_values = frequency(k,  mu0, muinf, b)
    omega_divergence_u0, omega_convergence_u0, omega_max, \
        racine_div, racine_conv = frequency_constant_state(k,  mu0, muinf, b, u0)
    
    
    ax1.plot(k, omega_divergence[1], color = colors[i], linestyle = '--')
    # ax3.plot(k, omega_divergence_u0[1], color = colors[i], linestyle = '--', alpha = 0.5)
    # ax3.plot(k, omega_max, color = 'orange', linestyle = '--', alpha = 0.5)
    # ax4.plot(k, omega_max, color = 'orange', linestyle = '--', alpha = 0.5)
    # ax3.plot(k, omega_max - racine_div, color = colors[i], linestyle = '-', alpha = 0.5)
    # ax3.plot(k, omega_max - racine_conv, color = colors[i], linestyle = '-', alpha = 0.5)
    ax3.plot(k, omega_convergence_u0[1], color = colors[i], linestyle = '-')
    ax3.plot(k, omega_divergence_u0[1], color = colors[i], linestyle = '--')
    # ax3.set_ylim(-0.1, 0.1)
    if i == 3:
        ax1.axhline(max_values[0])
    ax1.plot(k, omega_convergence[1], color = colors[i],linestyle = '-')
    ax2.plot(k, omega_divergence[0], color = colors[i], linestyle = '--')
    ax2.plot(k, omega_convergence[0], color = colors[i],linestyle = '-')
    
    # ax3.plot(k, omega_convergence_u0[1], color = colors[i],linestyle = '-', alpha = 0.5)
    ax4.plot(k, omega_convergence_u0[0], color = colors[i], linestyle = '--', alpha = 0.5)
    ax4.plot(k, omega_divergence_u0[0], color = colors[i],linestyle = '-', alpha = 0.5)
    labels.append(mlines.Line2D([], [], color=colors[i], marker='None', linestyle='None',
                        label=r'$\mu_b = {}$'.format(b)))

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    ax.set_xscale('log')

# ax1.set_ylim(-0.01, 0)
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='--',
                        label='Divergence')] + labels
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Convergence')] + labels

fig.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig1.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig.supxlabel(r'$k$ (1/m)')
fig1.supxlabel(r'$k$ (1/m)')

ax1.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)')
ax2.set_ylabel(r'$\text{Re}\{\omega_+\}$ (1/s)')
ax3.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)')
ax4.set_ylabel(r'$\text{Re}\{\omega_+\}$ (1/s)')
fig.savefig('omega.png')
fig1.savefig('omega_u0.png')


omega_divergence_1,omega_convergence_1, _, _, _ = frequency_constant_state(k,  mu0, muinf, 0.1, u0)
omega_divergence_2,omega_convergence_2, _, _, _= frequency_constant_state(k,  mu0, muinf, 1, u0)

diff = abs(omega_divergence_2[0] - omega_divergence_1[0])/omega_divergence_1[0]

plt.figure()
plt.plot(k, diff)
plt.xscale('log')
plt.xlabel(r'$k$ (1/m)')
plt.ylabel(r'error (big - small)')
plt.savefig('diff.png')



plt.figure(figsize = (8, 5))
ax0 = plt.axes()

plt.plot(dudx, stress3, color = 'g', label = r'$\mu_b = 0.5$')
# plt.plot(dudx, s3_d, color = 'g', label = r'$\mu_b = 0.5$')
# plt.plot(dudx, stress1, color = 'k', label = r'$\mu_b = 1$')
# plt.plot(dudx, stress2, color = 'r', label = r'$\mu_b = 10$')
plt.grid()
# ax0.set_xscale('symlog', linthresh = 1e-7)
plt.xlabel(r'$\frac{\partial u}{\partial x}$ (m/s)')
plt.ylabel(r'$\sigma_{11}$ (N/m)')
plt.legend()
plt.savefig('stress.png')