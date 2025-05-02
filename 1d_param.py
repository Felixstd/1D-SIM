import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import scienceplots

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
    
    sig11 = P*(np.sign(dudx)*(mu/2 + mub)-1)
    sig11_dilat = P*(np.sign(dudx)*(mu/2 + mub)-1) - (mub - mu/2)*dilat*P
    return sig11, sig11_dilat


def frequency(k, mu0, muinf, mub, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    
    delta_mu = muinf - mu0
    lamda = Pstar*np.exp(-C*(1-A0))
    Gamma = h0 * np.sqrt(lamda*rhoice)*(dmean*delta_mu/(2*I0))
    gamma_plus = (mu0/2 + mub) -1
    gamma_minus = -(mu0/2 + mub) -1
    
    #Divergence
    omega_plus_d = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1)))/(2*rhoice*h0)
    omega_minus_d = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1)))/(2*rhoice*h0)
    # print((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1))
    #convergence
    omega_plus_c = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)))/(2*rhoice*h0)
    omega_minus_c = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)))/(2*rhoice*h0)
    
    
    
    return [omega_plus_d, omega_minus_d], [omega_plus_c, omega_minus_c], \
        [-4*rhoice*h0**2*lamda*gamma_plus*(C*A0+1)/Gamma, -4*rhoice*h0**2*lamda*gamma_minus*(C*A0+1)/Gamma]


k = np.linspace(0, np.pi/5e3, 100000)
mu0 = .1
muinf = 0.8
mub = [0.01, 0.5, 1, 2]
dudx = np.linspace(-1/1000, 1/1000, 1000)
A = np.ones_like(dudx)
h = np.ones_like(dudx)

ice_strength = strength(h, A)
Inertial_num = inertial(dudx, ice_strength, h)
friction_coef = friction(Inertial_num, mu0, muinf)
tanpsi = dilantancy(friction_coef)
stress1, s1_d = sigma11(ice_strength, friction_coef, 1, dudx, tanpsi)
stress2, s2_d = sigma11(ice_strength, friction_coef, 10, dudx, tanpsi)
stress3, s3_d = sigma11(ice_strength, friction_coef, 0.5, dudx, tanpsi)




colors = ['k', 'g', 'r', 'b']
fig, (ax1, ax2)= plt.subplots(1, 2, sharex = True, figsize = (9, 4))
labels = []
for i, b in enumerate(mub):
    omega_divergence, omega_convergence, max_values = frequency(k,  mu0, muinf, b)
    ax1.plot(k, omega_divergence[1], color = colors[i], linestyle = '--')
    if i == 3:
        ax1.axhline(max_values[0])
    ax1.plot(k, omega_convergence[1], color = colors[i],linestyle = '-')
    ax2.plot(k, omega_divergence[0], color = colors[i], linestyle = '--')
    ax2.plot(k, omega_convergence[0], color = colors[i],linestyle = '-')
    labels.append(mlines.Line2D([], [], color=colors[i], marker='None', linestyle='None',
                        label=r'$\mu_b = {}$'.format(b)))

ax1.grid()

ax2.grid()
plt.xscale('log')
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='--',
                        label='Divergence')] + labels
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Convergence')] + labels

fig.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig.supxlabel(r'$k$ (1/m$^2$)')
ax1.set_ylabel(r'$\omega_-$ (1/s)')
ax2.set_ylabel(r'$\omega_+$ (1/s)')
plt.savefig('omega.png')


plt.figure(figsize = (8, 5))
ax0 = plt.axes()

# plt.plot(dudx, stress3, color = 'g', label = r'$\mu_b = 0.5$')
plt.plot(dudx, s3_d, color = 'g', label = r'$\mu_b = 0.5$')
# plt.plot(dudx, stress1, color = 'k', label = r'$\mu_b = 1$')
# plt.plot(dudx, stress2, color = 'r', label = r'$\mu_b = 10$')
plt.grid()
ax0.set_xscale('symlog', linthresh = 1e-6)
plt.xlabel(r'$\frac{\partial u}{\partial x}$ (m/s)')
plt.ylabel(r'$\sigma_{11}$ (N/m)')
plt.legend()
plt.savefig('stress.png')