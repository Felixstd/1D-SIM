import numpy as np
import matplotlib.pyplot as plt
from seaiceparameters import *

plt.style.use('/aos/home/fstdenis/1D-SIM/sim1dplotting/science.mplstyle')

def VPs_growth(k):
    """
    Function to plot the plastic VP
    growth rates. In both convergence and divergence

    Args:
        k (array): wavenumbers

    Returns:
        omega_div,omega_conv: growth rates in div and conv. 
    """
    
    omega_div = np.sqrt((Pstarstar*(1+C_s*A0)/rho)*lamda_div*k**2)
    
    omega_conv = np.sqrt((Pstarstar*(1+C_s*A0)/rho)*lamda_conv*k**2).astype(complex)
    
    return omega_div,omega_conv

def VP_nonlocal(k, l = 40e3):
    """
    Function to plot the plastic nVP
    growth rates. In both convergence and divergence

    Args:
        k (array): wavenumbers

    Returns:
        omega_div,omega_conv: growth rates in div and conv. 
    """
    l0 = l*A0   
    omega_div = k*np.sqrt(lamda_div/(rho*h0)*((2*P0+(C_s*A0-1)*Pstarstar*h0))*(1/(l0**2*k**2+1)))
    omega_conv = k*np.sqrt(lamda_conv/(rho*h0)*((2*P0+(C_s*A0-1)*Pstarstar*h0))*(1/(l0**2*k**2+1)).astype(complex))
    return omega_div,omega_conv


def VP_nonlocal_viscous(k, l=40e3):
    """
    Function to plot the viscous nVP
    growth rates. It only outputs the positive growth rate

    Args:
        k (array): wavenumbers

    Returns:
        omega: positive root growth rate 
    """
    alpha = 2/(l**2)*(P0-Pstarstar)
    beta = 4*rho*h0*((Pstarstar*h0/l**2*(1+C_s*A0))+alpha*A0)/nu_max**2
    arg = 1 - beta/(k**2*(1/l**2+k**2))
    omega = nu_max*k**2/(2*rho*h0)*(-1 + np.sqrt(arg.astype(complex)))
    return omega, beta/(k**2*(1/l**2+k**2))

def VP_local_viscous(k):
    """
    Function to plot the viscous VP
    growth rates. It only outputs the positive growth rate

    Args:
        k (array): wavenumbers

    Returns:
        omega: positive root growth rate 
    """
    arg = 1 - 4*rho*Pstarstar*h0**2*(1+C_s*A0)/(nu_max*k)**2
    omega = nu_max*k**2/(2*rho*h0)*(-1 + np.sqrt(arg.astype(complex)))
    
    x = 4 * rho * Pstarstar * h0**2 * (1 + C_s * A0) / (nu_max * k)**2
    arg = 1 - x

    # Stable rewrite of nu_max*k**2/(2*rho*h0) * (-1 + sqrt(arg))
    # using -1+sqrt(1-x) = -x/(1+sqrt(1-x)); the k**2 cancels analytically.
    B = 2 * Pstarstar * h0 * (1 + C_s * A0) / nu_max
    omega = -B / (1 + np.sqrt(arg.astype(complex)))

    return omega, arg

    # return omega, 4*rho*Pstarstar*h0**2*(1+C_s*A0)/(nu_max*k)**2


def omega_EVP(lamda, k_tot):
    k_tot = np.asarray(k_tot, dtype=np.float64)
    n = len(k_tot)

    gamma = lamda * Pstarstar * alpha * (1 - C_s * A0) / rhoi
    d = -gamma * k_tot**2          # constant term for each k, coeffs = [1, alpha, 0, d]

    # build a stack of companion matrices, shape (n, 3, 3)
    C = np.zeros((n, 3, 3), dtype=np.float64)
    C[:, 1, 0] = 1.0
    C[:, 2, 1] = 1.0
    C[:, 0, 0] = -alpha           # -a1/a0
    C[:, 0, 1] = 0.0              # -a2/a0
    C[:, 0, 2] = -d               # -a3/a0

    roots = np.linalg.eigvals(C)   # shape (n, 3), complex

    omega_1 = roots[:, 0]
    omega_2 = roots[:, 1]
    omega_3 = roots[:, 2]
    return omega_1, omega_2, omega_3

def EVP_viscous(k):
    """
    Function to plot the viscous VP
    growth rates. It only outputs the positive growth rate

    Args:
        k (array): wavenumbers

    Returns:
        omega: positive root growth rate 
    """
    arg = 1 - 4*rho*Pstarstar*h0**2*(1+C_s*A0)*(alpha*T-1)/(nu_max*k)**2
    omega = nu_max*k**2/(2*rho*h0*(alpha*T-1))*(-1 + np.sqrt(arg.astype(complex)))

    x = 4 * rho * Pstarstar * h0**2 * (1 + C_s * A0) / (nu_max * k)**2
    arg = 1 - x

    # Stable rewrite of nu_max*k**2/(2*rho*h0) * (-1 + sqrt(arg))
    # using -1+sqrt(1-x) = -x/(1+sqrt(1-x)); the k**2 cancels analytically.
    B = 2 * Pstarstar * h0 * (1 + C_s * A0) / nu_max
    omega = -B / (1 + np.sqrt(arg.astype(complex)))

    return omega, arg
    