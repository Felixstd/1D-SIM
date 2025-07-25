import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import scienceplots
from scipy.optimize import root_scalar
import warnings
warnings.filterwarnings('ignore')

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

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
    """
    This is for the plastic state with the base state velocity being 0 and the base state 
    thickness and concentration are constants. 

    Args:
        k (_type_): _description_
        mu0 (_type_): _description_
        muinf (_type_): _description_
        mub (_type_): _description_
        rhoice (int, optional): _description_. Defaults to 900.
        h0 (int, optional): _description_. Defaults to 1.
        A0 (int, optional): _description_. Defaults to 1.
        C (int, optional): _description_. Defaults to 20.
        Pstar (_type_, optional): _description_. Defaults to 27.5e3.
        dmean (_type_, optional): _description_. Defaults to 1e3.
        I0 (_type_, optional): _description_. Defaults to 1e-3.

    Returns:
        _type_: _description_
    """
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
        
def frequency_velocity(k, mu0, muinf, mub, du0dx, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    
    """
    This is for the plastic regime with du0dx = C but dh0dx and da0dx = 0

    Returns:
        _type_: _description_
    """
    
    delta_mu = muinf - mu0
    lamda = Pstar*np.exp(-C*(1-A0))
    Ibar = dmean*du0dx*np.sqrt(rhoice/lamda)
    Gamma = h0 * np.sqrt(lamda*rhoice)*(dmean*delta_mu*I0/(2*(I0+Ibar)**2))
    
    
    mu_0 = mu0 + (muinf-mu0)/(I0/Ibar + 1)
    
    
    gamma_plus = (mu_0/2 + mub) -1
    gamma_minus = -(mu_0/2 + mub) -1
    Gamma_star_plus = Gamma * C * du0dx/2
    Gamma_star_minus = Gamma * C * du0dx/2
    
    
    
    #Divergence
    omega_plus_d = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1) + 4*rhoice*h0**2*k**2*A0 * Gamma_star_plus +0j))/(2*rhoice*h0)
    omega_minus_d = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_plus*(C*A0+1)+4*rhoice*h0**2*k**2*A0 * Gamma_star_plus+0j))/(2*rhoice*h0)

    #convergence
    omega_plus_c = (Gamma*k**2 + np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)+ 4*rhoice*h0*k**2*A0 * Gamma_star_minus+0j))/(2*rhoice*h0)
    omega_minus_c = (Gamma*k**2 - np.sqrt((Gamma*k**2)**2 + 4*rhoice*h0**2*k**2*lamda*gamma_minus*(C*A0+1)+4*rhoice*h0*k**2*A0 * Gamma_star_minus+0j))/(2*rhoice*h0)
    
    
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

def frequency_VP(k, A0 = 1, C = 20, rhoice = 900, e =2):
    
    omega_conv_plus = np.sqrt(k**2*(1+C*A0)/rhoice * (-np.sqrt(1+e**(-2))-1)/2) + 0j
    
    omega_div_plus = np.sqrt(k**2*(1+C*A0)/rhoice * (np.sqrt(1+e**(-2))-1)/2) + 0j
    omega_div_minus = -np.sqrt(k**2*(1+C*A0)/rhoice * (np.sqrt(1+e**(-2))-1)/2) + 0j
    
    return omega_conv_plus, omega_div_plus, omega_div_minus

def frequency_mub(k, mub, dA0dx, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    
    Pstarstar = Pstar*np.exp(-C*(1-A0))
    mub_conv = -mub-1
    mub_div = mub-1
    lambdastar_div = Pstarstar*mub_div*C/(rhoice*h0)
    lambdastar_conv = Pstarstar*mub_conv*C/(rhoice*h0)
    
    omega_plus_div = np.sqrt(lambdastar_div*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    omega_minus_div = -np.sqrt(lambdastar_div*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    
    omega_plus_conv= np.sqrt(lambdastar_conv*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    omega_minus_conv = -np.sqrt(lambdastar_conv*(dA0dx**2*C - A0*k**2 + 1j*k*dA0dx*(1+C*A0)))
    
    return omega_plus_div, omega_minus_div, omega_plus_conv, omega_minus_conv


def frequency_dhdx_dudx(k, mub, dh0dx, du0dx, sign_u, rhoice = 900, h0 = 1, A0 = 1, C = 20, Pstar = 27.5e3, dmean = 1e3, I0 = 1e-3):
    """
    This the dispersion relation for the plastic regime considering dh0dx and du0dx as constants and A = A0 (a constant)
    The dispersion relation has the form 
    
    A*w^2 + Bw + C = 0
    Args:
        k (_type_): _description_
        mub (_type_): _description_
        dA0dx (_type_): _description_
        rhoice (int, optional): _description_. Defaults to 900.
        h0 (int, optional): _description_. Defaults to 1.
        A0 (int, optional): _description_. Defaults to 1.
        C (int, optional): _description_. Defaults to 20.
        Pstar (_type_, optional): _description_. Defaults to 27.5e3.
        dmean (_type_, optional): _description_. Defaults to 1e3.
        I0 (_type_, optional): _description_. Defaults to 1e-3.
    """
    
    
    #------ Defining Constants for B and C -------#
    
    Pstarstar = Pstar*np.exp(-C*(1-A0))
    deltamu = muinf - mu0
    Ibar = dmean*np.sqrt(rhoice/Pstarstar)*abs(du0dx)
    print(mub)
    
    I0_prime = deltamu*I0/(I0 + Ibar)**2
    
    mustar_divergence = mub + (mu0 + deltamu/(I0/Ibar +1))/2 - 1
    mustar_convergence = -(mub + (mu0 + deltamu/(I0/Ibar +1))/2) - 1
    
    alpha = Pstarstar*I0_prime*dmean*np.sqrt(rhoice/Pstarstar)*dh0dx
    Gamma = dmean/2*np.sqrt(rhoice*Pstarstar)*I0_prime*h0
    
    
    
    #----- Defining A -------#
    A = rhoice*h0
    
    #------- Defining B -------#
    B_divergence = -rhoice*h0*du0dx - Gamma*k**2 +alpha*sign_u*k - 1j*h0*u0*rhoice*k
    B_convergence = -rhoice*h0*du0dx - Gamma*k**2 - alpha*sign_u*k - 1j*h0**2*rhoice*k
    #------- Defining C -----#
    C_divergence = alpha*sign_u*k*du0dx + Pstarstar*mustar_divergence*h0*k**2 + Gamma*k**2*du0dx + \
                    1j*(-alpha*sign_u*k*h0 - Pstarstar*mustar_divergence*dh0dx*k**2 + Gamma*k**3*u0)
                    
    C_convergence = +alpha*sign_u*k*du0dx + Pstarstar*mustar_convergence*h0*k**2 + Gamma*k**2*du0dx + \
                    1j*(alpha*sign_u*k*h0 - Pstarstar*mustar_convergence*dh0dx*k**2 + Gamma*k**2*du0dx)
                    
    omega_plus_divergence = (-B_divergence + np.sqrt(B_divergence**2 - 4*(A)*(C_divergence)))/(2*A)
    omega_minus_divergence = (-B_divergence - np.sqrt(B_divergence**2 - 4*(A)*(C_divergence)))/(2*A) 
    
    omega_plus_convergence = (-B_convergence + np.sqrt(B_convergence**2 - 4*(A)*(C_convergence)))/(2*A)
    omega_minus_convergence = (-B_convergence - np.sqrt(B_convergence**2 - 4*(A)*(C_convergence)))/(2*A) 
    
    lim_big_k = Pstarstar*mustar_divergence/Gamma - h0*rhoice*u0**2/Gamma + du0dx
    
    
    
    return [omega_plus_divergence, omega_minus_divergence], [omega_plus_convergence, omega_minus_convergence], lim_big_k
    
    
    
    
    
k = np.linspace(1e-6, 1e0, 1000000)
mu0 = 0.2
muinf = 0.8
u0 = 0.01
du0dx = 0.01/1e2
dh0dx = 1/1e2
# mub = [0.01]
# mub = [0.5, 10]
mub = [0.1, 10]
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




colors = ['k', 'r', 'b']
fig, (ax1, ax2)= plt.subplots(1, 2, sharex = True, figsize = (9, 4))
fig1, (ax3, ax4)= plt.subplots(1, 2, sharex = True, figsize = (9, 4))
fig1= plt.figure(figsize = (7, 4))
ax3 = plt.axes()
labels = []
for i, b in enumerate(mub):
    omega_divergence, omega_convergence, max_values = frequency(k,  mu0, muinf, b)
    omega_divergence_u0, omega_convergence_u0, omega_max, \
        racine_div, racine_conv = frequency_constant_state(k,  mu0, muinf, b, u0)
    
    omega_divergence_u0, omega_convergence_u0, omega_max = frequency_velocity(k,  mu0, muinf, b, 1e8, Pstar=27.5e03/2)
    
    omega_divergence_du0dx_dh0dx, omega_convergence_du0dx_dh0dx, max_lim_k_div = frequency_dhdx_dudx(k, b, dh0dx, du0dx, 1)
    # omega_plus_div_mub, omega_minus_div_mub, omega_plus_conv_mub, omega_minus_conv_mub = frequency_mub(k, mub, dA0dx)

    ax1.plot(k, omega_divergence[1], color = colors[i], linestyle = '--')
    # ax3.plot(k, omega_divergence_u0[1], color = colors[i], linestyle = '--', alpha = 0.5)
    # ax3.plot(k, omega_max, color = 'orange', linestyle = '--', alpha = 0.5)
    # ax4.plot(k, omega_max, color = 'orange', linestyle = '--', alpha = 0.5)
    # ax3.plot(k, omega_max - racine_div, color = colors[i], linestyle = '-', alpha = 0.5)
    # ax3.plot(k, omega_max - racine_conv, color = colors[i], linestyle = '-', alpha = 0.5)
    ax3.plot(k, omega_divergence_du0dx_dh0dx[1], color = colors[i], linestyle = '-')
    ax3.plot(k, omega_divergence_du0dx_dh0dx[1], color = colors[i], linestyle = '--')
    ax3.axhline(max_lim_k_div, zorder = 0, linestyle = '--', color = colors[i])
    # ax3.set_ylim(-0.1, 0.1)
    if i == 3:
        ax1.axhline(max_values[0])
    ax1.plot(k, omega_convergence[1], color = colors[i],linestyle = '-')
    ax2.plot(k, omega_divergence[0], color = colors[i], linestyle = '--')
    ax2.plot(k, omega_convergence[0], color = colors[i],linestyle = '--')
    
    # ax3.plot(k, omega_convergence_u0[1], color = colors[i],linestyle = '-', alpha = 0.5)
    # ax4.plot(k, omega_convergence_u0[0], color = colors[i], linestyle = '--', alpha = 0.5)
    # ax4.plot(k, omega_divergence_u0[0], color = colors[i],linestyle = '-', alpha = 0.5)
    labels.append(mlines.Line2D([], [], color=colors[i], marker='None', linestyle='None',
                        label=r'$\mu_b = {}$'.format(b)))

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid()
    ax.set_xscale('log')
omega_vp_conv, omega_vp_div_plus, omega_vp_div_minus = frequency_VP(k)

# ax3.plot(k, omega_vp_conv, color = 'g')
# ax3.plot(k, omega_vp_div_plus, color = 'darkviolet', linestyle = '--')
# ax3.plot(k, omega_vp_div_minus, color = 'darkviolet', linestyle = '--')
# ax1.set_ylim(-0.01, 0)
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='--',
                        label='Divergence')] + labels
labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Convergence')] + labels #+ [mlines.Line2D([], [], color='darkviolet', marker='None', linestyle='None',
                        # label='VP')]

fig.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1, 0.9))
fig1.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9), fontsize = 15)
fig.supxlabel(r'$k$ (1/m)', fontsize = 15)
fig1.supxlabel(r'$k$ (1/m)', fontsize = 15)
# ax3.set_ylim(-0.01, 0.05)
ax1.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)')
ax2.set_ylabel(r'$\text{Re}\{\omega_+\}$ (1/s)')
ax3.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)', fontsize = 15)
ax3.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}$'.format(du0dx, dh0dx))
# ax4.set_ylabel(r'$\text{Re}\{\omega_+\}$ (1/s)')
fig.savefig('omega.png')
fig1.savefig('omega_du0dx_dh0dx.png')


omega_divergence_1,omega_convergence_1, _, _, _ = frequency_constant_state(k,  mu0, muinf, 0.1, u0)
omega_divergence_2,omega_convergence_2, _, _, _= frequency_constant_state(k,  mu0, muinf, 1, u0)

diff = abs(omega_divergence_2[0] - omega_divergence_1[0])/omega_divergence_1[0]

plt.figure()
plt.plot(k, diff)
plt.xscale('log')
plt.xlabel(r'$k$ (1/m)')
plt.ylabel(r'error (big - small)')
plt.savefig('diff.png')


print(np.where(np.diff(np.real(omega_divergence_du0dx_dh0dx[1]))> 0))
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

plt.figure()
plt.plot(k[1:], np.diff(np.real(omega_divergence_du0dx_dh0dx[1])))
plt.xscale('log')
plt.savefig('test2.png')