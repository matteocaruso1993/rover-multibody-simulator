# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:36:38 2021

@author: Matteo
"""
import numpy as np

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



def wrap2pi(angles):
    
    
    new_angles = (angles + np.pi) % (2*np.pi) -np.pi
    
    return new_angles





def findContinguousContactRegions(input_array):
    """
    Parameters
    ----------
    input_array : A numpy array of size (n,)
        This is an array of booleans.

    Returns
    -------
    groups : List of list. Each list represent a set of indices of at least 
            2 contiguous True elements.
        

    """
    groups = list()
    tmp = list()
    for i in range(input_array.shape[0]):
        if input_array[i]:
            tmp.append(i)
            
        else:
            if tmp and len(tmp)>=2:
                groups.append(tmp)
            
            tmp = list()
    
    if len(tmp) >= 2:
        groups.append(tmp)
    
    return groups
        