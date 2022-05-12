# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:36:38 2021

@author: Matteo
"""
import numpy as np
from numba import jit
import sobol
import os
import pandas as pd

@jit(nopython=True)
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



def generateWheelsPoints(R,theta):
    
    points = np.zeros((theta.shape[0],3))
    
    z = R*np.cos(theta)
    y = R*np.sin(theta)
    points[:,-1] = z
    points[:,1] = y
    return points


def clip(v, min_v, max_v):
    if v > max_v:
        v = max_v
    elif v < min_v:
        v = min_v
    else:
        pass
    return v


def generateDOE(in_dict, n_points = 100, method='random'):
    cols = len(list(in_dict.keys()))
    
    out = np.zeros((n_points,cols), dtype=float)
    
    for n, key in enumerate(in_dict.keys()):
        out[:,n] = in_dict[key]['default']
    
    mask = [in_dict[key]['variable'] for key in in_dict.keys()]
    col_idx = np.where(mask)[0]
    #print(col_idx)
    dims = col_idx.shape[0]
    
    if method == 'sobol':
        out_1 = sobol.sample(dims, n_points)
    elif method == 'random':
        out_1 = np.random.rand(n_points,dims)
    else:
        raise ValueError('Invalid Method')
        
    
    out[:,col_idx] = out_1
        
        
    counter=0
    for key in in_dict.keys():
        if in_dict[key]['variable']:
            in_dict[key]['values'] = in_dict[key]['min'] + out[:,counter]*(in_dict[key]['max'] - in_dict[key]['min'])
        else:
            in_dict[key]['values'] = out[:,counter]
        counter+=1
        
    return in_dict


def saveVars(filename, var_list):
    with open(filename, 'a') as f:
        string_to_write = ', '.join(map(str,var_list))
        f.write('\n' + string_to_write)
        
    f.close()
    
    
def saveSimulationData(sim, destination_folder):
    soluzione = sim.ode_sol
    out = np.vstack((soluzione['t'],soluzione['y']))
    header = ['t']
    for c in sim.gen_coord:
        header.append(c.name)
            
    for c in sim.gen_speeds:
        header.append(c.name)
            
    my_string = ','.join(map(str,header))
        
    num = 0
    while True:
        filename = os.path.join(destination_folder,'sim_py_v' + str(num) + '_transpose_con_header.csv')
        if os.path.exists(filename):
            num += 1
        else:
            break
                    
        
    np.savetxt(filename, out.T, delimiter=',', header=my_string)
    
    return os.path.split(filename)


def saveConfig(sim, destination_filename):
    with open(destination_filename,'w') as f:
        sim.config.write(f)
        
    f.close()
    
    
    
    

def loadExperimentalData(path_to_folder, sim_type, pin_id):
    if sim_type == 'dead':
        prefix = 'DD_'
    elif sim_type == 'skew':
        prefix = 'SD_'
    else:
        raise ValueError('Invalid sim type')
        
    valid_pin_id = [1,2,3,4]
    if pin_id not in valid_pin_id:
        raise ValueError('Invalid pin ID')
        
    filename_x = os.path.join(path_to_folder, prefix + str(pin_id) + '_X.csv')
    filename_y = os.path.join(path_to_folder, prefix + str(pin_id) + '_Y.csv')     
    data_x = pd.read_csv(filename_x, skiprows=[0,1])
    data_y = pd.read_csv(filename_y, skiprows=[0,1])
    
    out = {'Pin ' +str(pin_id): {'x': data_x, 'y': data_y}}
    
    return out


def loadAllExperimentData(path_to_folder, sim_type):
    pin_id = [1,2,3,4]
    data_full = dict()
    
    for ID in pin_id:
        data = loadExperimentalData(path_to_folder=path_to_folder, sim_type=sim_type, pin_id=ID)
        data_full = {**data_full, **data}
    
    return data_full




def computeEnergyLike(in_array):
    return np.sum(in_array**2)




    
    

    
    
