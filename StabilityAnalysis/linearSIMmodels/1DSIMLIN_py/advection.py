import numpy as np

from global_constants import *
from numba import jit, prange 

@jit(nopython=True)
def calc_flux(ui, uip1, Tim1, Ti, Tip1):
    """
    Function used to compute the flux coming in and out 
    of a grid cell. 
    
    It's combined with the upwind advection method where the 
    x-derivative in 
    
    h_t + (uh)_x = 0
    
    depend on the sign of the velocity. It's either a forward 
    or backward derivative. 
    

    Args:
        ui (_type_): _description_
        uip1 (_type_): _description_
        Tim1 (_type_): _description_
        Ti (_type_): _description_
        Tip1 (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    if (ui > 0) :                                                                        
        flux1 = ui*Tim1
    else:
        flux1 = ui*Ti


    if (uip1 < 0):                                                                             
        flux2 = uip1*Ti
    else:
        flux2 = uip1*Tip1

    flux = flux2 - flux1
    
    return flux

@jit(nopython=True)
def upwind_scheme(u, h, N, concentration = False) : 
    """
    
    This is a function used to solve the continuity equation using 
    the upwind scheme. It advects forward in time using the 
    fluxes computed with the last function. 
    
    hnp1 = hnm1 - dt/dx*Flux

    Args:
        u (_type_): _description_
        h (_type_): _description_
        N (_type_): _description_
        concentration (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    
    hout = np.zeros_like(h)
    for i in prange(1, N-1):
        
        flux = calc_flux(u[i], u[i+1], h[i-1], h[i], h[i+1])
        
        hout[i] = h[i] - DtoverDx*flux
        hout[i] = max(hout[i], 0)
        if concentration:
            hout[i] = min(hout[i], 1)

    hout[0] = 0
    hout[-1] = 0
    
    return hout