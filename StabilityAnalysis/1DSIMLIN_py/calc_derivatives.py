import numpy as np

from numba import jit

@jit(nopython=True)
def first_derivative(x, deltax, N) : 
    """
    Function to evaluate the first derivative of an array x. 
    
    This is the formula for a second order accuracy. 

    Args:
        x (_type_): _description_
        deltax (_type_): _description_
        N (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    derivative = np.zeros_like(x)
    
    derivative[1:N-1] = (x[2:] - x[:-2])/(2*deltax)
    
    derivative[0] = (-3*x[0] + 4*x[1] - x[2])/(2*deltax)
    derivative[-1] = (3*x[-1] - 4*x[-2] + x[-3])/(2*deltax)
    
    return derivative

@jit(nopython=True)
def second_derivative(x, deltax, N):
    """
    Function to evaluate the second derivative of an array x. 
    
    This is the formula for a second order accuracy. 

    Args:
        x (_type_): _description_
        deltax (_type_): _description_
        N (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    sec_derivative = np.zeros_like(x)
    
    sec_derivative[1:N] = (x[2:] - 2*x[1:N] + x[:-2])/(deltax**2)
    sec_derivative[0] = (-x[3] + 4*x[2] - 5*x[1] + 2*x[0])/deltax**2
    sec_derivative[-1] = (-x[-4] + 4*x[-3] - 5*x[-2] + 2*x[-1])/deltax**2
    
    
    return sec_derivative