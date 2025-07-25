import sympy as sp
from scipy.optimize import root_scalar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def constants_dh0dx(mub):
    
    Pstarstar = Pstar*np.exp(-C*(1-A0_c))
    deltamu = muinf - mu0
    Ibar = dmean*np.sqrt(rhoice/Pstarstar)*abs(u0x_c)
    deltamu_star = deltamu*I0/(I0 + Ibar)**2
    
    mubar = mu0 + deltamu/(I0/Ibar + 1)

    alpha_c = (Pstarstar * deltamu_star * dmean * h0x_c)/(2)*np.sqrt(rhoice/Pstarstar)
    gamma_star_c  = Pstarstar*(mub + mubar/2 - 1)
    Gamma_c  = dmean*h0_c*np.sqrt(Pstarstar*rhoice)/2*deltamu_star  
    
    return alpha_c, gamma_star_c, Gamma_c


def constants_dA0dx(mub):
    
    Pstarstar = Pstar*np.exp(-C*(1-A0_c))
    Pstarstar_c = Pstarstar*C*h0_c
    deltamu = muinf - mu0
    Ibar = dmean*np.sqrt(rhoice/Pstarstar)*abs(u0x_c)
    deltamu_star = deltamu*I0/(I0 + Ibar)**2
    mubar = mu0 + deltamu/(I0/Ibar + 1)
    
    Gamma_c  = dmean*h0_c*np.sqrt(Pstarstar*rhoice)/2*deltamu_star  
    gamma_div = (mub + mubar/2) -1
    beta_div = Pstarstar_c*gamma_div - (Gamma_c*C/2)*abs(u0x_c)
    lamda_div = gamma_div*Pstarstar_c*C*A0x_c - Gamma_c*C**2/4*abs(u0x_c)*A0x_c
    alpha_div = dmean*np.sqrt(rhoice/Pstarstar)*((Gamma_c*C)/(Ibar+I0)**2*abs(u0x_c)*A0x_c + Pstarstar_c * deltamu_star/2*A0x_c)
    Gamma_star = C*Gamma_c*A0x_c/2
    
    gamma_conv = -(mub + mubar/2) -1
    beta_conv = Pstarstar_c*gamma_conv - (Gamma_c*C/2)*abs(u0x_c)
    alpha_conv = dmean*np.sqrt(rhoice/Pstarstar)*((Gamma_c*C)/(Ibar+I0)**2*abs(u0x_c)*A0x_c - Pstarstar_c * deltamu_star/2*A0x_c)
    lamda_conv = gamma_conv*Pstarstar_c*C*A0x_c - Gamma_c*C**2/4*abs(u0x_c)*A0x_c
    return beta_div, lamda_div, alpha_div, Gamma_star , Gamma_c, gamma_conv, beta_conv, alpha_conv, lamda_conv
    
    
def dispersion_relation_divergence_dh0dx(alpha_c, gamma_star_c, Gamma_c):
    """
    This the dispersion relation for divergence for 
    h = h_0(x) + h'
    A = A_0
    u = u_0(x) + u'

    Args:
        alpha_c (_type_): _description_
        gamma_star_c (_type_): _description_
        Gamma_c (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    k = sp.symbols('k', real = True)
    rhoi, h0, h0x, u0, u0x, alpha, gamma_star, Gamma = \
        sp.symbols('rhoi, h0,h0x, u0, u0x, alpha, gamma_star, Gamma', 
        real = True)

    params = {
        rhoi: rhoice,
        h0: h0_c,
        u0x: u0x_c,
        h0x: h0x_c,
        alpha: alpha_c,
        gamma_star: gamma_star_c,
        Gamma: Gamma_c,
        u0: u0_c
    }

    A = rhoi*h0
    B = (-rhoi * h0 * u0x - Gamma * k**2 + k*alpha - sp.I*(u0 * h0 *rhoi * k))
    C = gamma_star * h0 * k**2 - alpha*u0x*k + Gamma*u0x*k**2 + sp.I*(u0*Gamma*k**3 - gamma_star*h0x*k - u0*alpha*k**2)

    omega = (-B - sp.sqrt(B**2 - 4*A*C))/(2*A)

    real_omega = sp.re(omega)
    
    w = sp.lambdify(k, real_omega.subs(params), modules=['numpy'])

    max_value_real_omega =  gamma_star_c * h0_c / Gamma_c + u0x_c
    
    return real_omega, w, max_value_real_omega, params
    

def dispersion_dA0dx(beta_c, lamda_c, alpha_c, Gamma_c, gamma_star, C_c = 20):
    
    """
    This the dispersion relation for divergence for 
    h = h_0
    A = A_0(x)+h'
    u = u_0(x) + u'

    Args:
        alpha_c (_type_): _description_
        gamma_star_c (_type_): _description_
        Gamma_c (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    k = sp.symbols('k', real = True)
    rhoi, h0, A0, A0x, u0, u0x, alpha_star, lambda_star, Gamma, beta = \
        sp.symbols('rhoi, h0,A0, A0x, u0, u0x, alpha_star, lambda_star, Gamma, beta', 
        real = True)

    params = {
        rhoi: rhoice,
        h0: h0_c,
        u0x: u0x_c,
        A0x: A0x_c,
        alpha_star: alpha_c - gamma_star,
        lambda_star: lamda_c - alpha_c*C_c/2*abs(u0x_c),
        Gamma: Gamma_c,
        beta: beta_c,
        u0: u0_c,
        A0:A0_c
    }

    A = rhoi*h0
    B = (alpha_star*k - Gamma*k**2 - rhoi*h0*u0x - sp.I*(k*u0*rhoi*h0))
    C = (beta*A0*k**2 + Gamma*k**2*u0x - lambda_star *A0x - alpha_star*k*u0x + sp.I*(Gamma*k**3*u0 - A0*lambda_star*k-beta*k*A0x - alpha_star*k**2*u0))

    omega = (-B - sp.sqrt(B**2 - 4*A*C))/(2*A)

    real_omega = sp.re(omega)
    
    w = sp.lambdify(k, real_omega.subs(params), modules=['numpy'])

    max_value_real_omega =  A0_c*beta_c/Gamma_c + u0x_c
    
    return real_omega, w, max_value_real_omega, params

#----- Defining constants -------#
dx = 1e4
rhoice = 900
h0_c = 1
A0_c = 1
C = 20
Pstar = 27.5e3
dmean = 1e3
I0 = 1e-3
u0_c = 0.1
u0x_c = u0_c/dx
# u0x_c = 0
h0x_c = h0_c/dx
# h0x_c = 0
A0_x = 1
A0x_c = A0_c/dx
mu0 = 0.2
muinf = 0.8
mub = 10

mub_tot = [0.5,10]

k = sp.symbols('k', real = True)
k_vals = np.linspace(1e-6, 1e3, int(1e7))

colors = ['b', 'r']
fig1 = plt.figure(figsize = (9, 4))
ax1 = plt.axes()
labels = []

fig2 = plt.figure(figsize = (9, 4))
ax2 = plt.axes()

fig3 = plt.figure(figsize = (9, 4))
ax3 = plt.axes()

axes = [ax1, ax2, ax3]

for i, mub in enumerate(mub_tot):
    print('mub: ', mub)
    alpha_c, gamma_star_c, Gamma_c = constants_dh0dx(mub)
    real_omega_div, w_div, max_value_div, params_div = dispersion_relation_divergence_dh0dx(alpha_c, gamma_star_c, Gamma_c )
    
    
    beta_div_da0dx, lamda_div_da0dx, alpha_div_da0dx, Gamma_star_da0dx , Gamma_c_da0dx, \
        gamma_conv_da0dx, beta_conv_da0dx, alpha_conv_da0dx, lamda_conv_da0dx  = constants_dA0dx(mub)
        
    print(gamma_conv_da0dx, beta_conv_da0dx, alpha_conv_da0dx, lamda_conv_da0dx)

    real_omega_div_da0dx, w_div_da0dx, max_value_div_da0dx, params_div_da0dx = \
        dispersion_dA0dx(beta_div_da0dx, lamda_div_da0dx, alpha_div_da0dx, Gamma_star_da0dx , Gamma_c_da0dx)
        
    real_omega_conv_da0dx, w_conv_da0dx, max_value_conv_da0dx, params_conv_da0dx = \
        dispersion_dA0dx(beta_conv_da0dx, lamda_conv_da0dx, alpha_conv_da0dx, Gamma_star_da0dx , Gamma_c_da0dx)
    
    
    omega_k = w_div(k_vals)
    omega_k_da0dx = w_div_da0dx(k_vals)
    
    omega_k_da0dx_conv = w_conv_da0dx(k_vals)
    
    drealw_dk = sp.diff(real_omega_div, k)
    domega_dk = sp.lambdify(k, drealw_dk.subs(params_div), modules=['numpy'])
    
    try:
        sol = root_scalar(domega_dk, method='bisect', bracket=[1e-3, k_vals[146529]])
        ax1.axhline(w_div(sol.root), color='g', linestyle='--')
    except:
        print('No Root')
    
    dwdk = domega_dk(k_vals)
    
    
    ax1.plot(k_vals, omega_k, color = colors[i])
    ax1.axhline(max_value_div, color='k', linestyle='--')
    
    ax3.plot(k_vals, omega_k_da0dx, color = colors[i])
    ax3.axhline(max_value_conv_da0dx, color='k', linestyle='--')
    
    
    
    
    ax2.plot(k_vals, dwdk, color = colors[i])
    labels.append(mlines.Line2D([], [], color=colors[i], marker='None', linestyle='None',
                        label=r'$\mu_b = {}$'.format(mub)))


print(np.where(dwdk == np.max(dwdk)))

# # Use root finding to solve f(k) = 0
# sol = root_scalar(f, method='bisect')
# print(sol.root)
# real_w_func = sp.lambdify(sol.root, real_omega.subs(params), modules=['numpy'])
# print(real_w_func)

# k_vals = np.linspace(1e-6, 1e0, int(1e9))
# omega_k = w_div(k_vals)




#------ Plotting --------#



labels = [mlines.Line2D([], [], color='k', marker='None', linestyle='-',
                        label='Divergence')] + labels

fig1.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9), fontsize = 15)
ax1.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)', fontsize = 15)
ax1.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}$'.format(u0x_c, h0x_c))

fig3.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9), fontsize = 15)
ax3.set_ylabel(r'$\text{Re}\{\omega_-\}$ (1/s)', fontsize = 15)
ax3.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dA_0/dx = {} \; \text{{m}}^{{-1}}$'.format(u0x_c, A0x_c))

fig2.legend(loc='upper center', 
               handles=labels, labelcolor='linecolor',  bbox_to_anchor=(1.05, 0.9), fontsize = 15)
ax2.set_ylabel(r'$d\text{Re}\{\omega_-\}/dk$ (ms)$^{-1}$', fontsize = 15)
ax2.set_title(r'$du_0/dx = {} \; \text{{s}}^{{-1}}, dh_0/dx = {}$'.format(u0x_c, h0x_c))

for ax in axes:
    ax.set_xlabel("$k$ (1/m)")
    ax.set_xscale("log")
    ax.grid()


fig1.savefig('omega_dh0dx_du0dx.png')
fig2.savefig('domegadk_dh0dx_du0dx.png')
fig3.savefig('omega_dA0dx_du0dx.png')

