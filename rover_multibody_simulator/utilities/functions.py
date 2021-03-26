# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:36:38 2021

@author: Matteo
"""

def step5(x, x_0, y_0, x_max, y_max):
    """
    This function approximate the Heaviside step function

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    x_0 : TYPE
        DESCRIPTION.
    y_0 : TYPE
        DESCRIPTION.
    x_max : TYPE
        DESCRIPTION.
    y_max : TYPE
        DESCRIPTION.

    Returns
    -------
    c : TYPE
        DESCRIPTION.

    """
    
    if x < x_0:
        c = y_0
        
    elif x > x_max:
        c = y_max
        
        
    else:
        a = y_max - y_0;
        delta = (x-x_0)/(x_max-x_0);
        c = y_0 + a*delta**3*(10-15*delta+6*delta**2)
    
    

    return c