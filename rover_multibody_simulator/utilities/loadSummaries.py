# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:34:17 2021

@author: Matteo
"""

import os
from rover_multibody_simulator.utilities import postprocess
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt


def loadSummaries(path_to_folder_simulations, idx_list):
    #load_data
    data_loaded = dict()
    doe_to_load = idx_list

    for idx in doe_to_load:
        path_to_folder = os.path.join(path_to_folder_simulations, str(idx))
        path_to_summary = os.path.join(path_to_folder, 'Summary_con_RMS.csv')
        
        loaded_summ = postprocess.loadSummaryData(path_to_summary)
        data_loaded[str(idx)] = loaded_summ
        
    return data_loaded


def createInterpSurface(dict_of_loaded_summaries,var1, var2, fitness, n_points = 100, mode='grid'):
    #usiamo griddata
    var1_np = np.empty((0,))
    var2_np = np.empty((0,))
    fitness_np = np.empty((0,))
    
    for key in dict_of_loaded_summaries:
        var1_np = np.hstack((var1_np, dict_of_loaded_summaries[key][var1].to_numpy()))
        var2_np = np.hstack((var2_np, dict_of_loaded_summaries[key][var2].to_numpy()))
        fitness_np = np.hstack((fitness_np, dict_of_loaded_summaries[key][fitness].to_numpy()))
    
        
    x_min = var1_np.min()
    x_max = var1_np.max()
    y_min = var2_np.min()
    y_max = var2_np.max()
    
    
    X,Y = np.meshgrid(np.linspace(x_min, x_max, n_points), np.linspace(y_min, y_max, n_points))
    
    if mode=='grid':
        interp = interpolate.griddata(np.hstack((var1_np.reshape(var1_np.shape[0],-1), var2_np.reshape(var2_np.shape[0],-1))),\
                                  fitness_np, (X,Y), method='cubic')
    elif mode=='interpolant':
        print(var1_np.shape)
        interp = interpolate.SmoothBivariateSpline(np.array([var1_np]), np.array([var2_np.flatten()]), np.array([fitness_np.flatten()]))
        
    return interp, X, Y
    
    
    
    
    
    

if __name__ == '__main__':
    idx_to_load = [11,14,15,17,18,19]
    sim_path = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\DesignOfExperiments'
    test = loadSummaries(sim_path, idx_to_load)
    
    interp, X, Y = createInterpSurface(test, ' f*g', ' c_fr1', 'MSE', 2000, mode='interpolant')
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x_new = np.linspace(0.1,1,1000)
    y_new = np.linspace(100,2000,1000)
    X_new,Y_new = np.meshgrid(x_new, y_new)
    #ax.plot_surface(X,Y,interp(X[0,:],Y[:,0]))
    ax.plot_surface(X_new,Y_new,interp(x_new,y_new))
    ax.set_xlabel('f*g')
    ax.set_ylabel('c')
    ax.set_zlabel('MSE')
    