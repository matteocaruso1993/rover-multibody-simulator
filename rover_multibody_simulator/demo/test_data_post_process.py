# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:59:08 2021

@author: Matteo
"""
import os
import matplotlib.pyplot as plt
import pandas as pd
from rover_multibody_simulator.utilities import functions, postprocess
from tqdm import tqdm
from copy import deepcopy

show_exp = False
sim_type = 'dead'
plot_results = True
show_var_coverage = True
generate_new_csv = False
correct_offset = True
debug_offset = False

if __name__ == '__main__':
    path_to_csv_summary = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\DesignOfExperiments\6\Summary.csv'
    foldername = os.path.dirname(path_to_csv_summary)
    path_to_experimental_data = os.path.dirname(os.path.dirname(foldername))
    
    exp_data = functions.loadAllExperimentData(path_to_experimental_data, sim_type=sim_type)
        
    #load simulation summary
    sim_summary = postprocess.create_temporary_copy_csv(path_to_csv_summary)
        
    #Preprocess Experimental data
    exp_data = postprocess.preprocessExperimentalData(exp_data)
    if show_exp:
        plt.plot(exp_data['Pin 1']['y_proc_valid'].Time, exp_data['Pin 1']['y_proc_valid'].Amplitude)
       
    simulation_data = dict()
        
    for sim_id in tqdm(range(len(sim_summary)), desc='Loading simulation data...'):
        file_to_load = os.path.join(foldername, sim_summary['#filename'][sim_id])
        #print(file_to_load)
        sim_data_tmp = postprocess.loadSummaryData(file_to_load)
        simulation_data[str(sim_id)] = sim_data_tmp
       
    if plot_results:
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(exp_data['Pin 1']['y_proc_valid'].Time, exp_data['Pin 1']['y_proc_valid'].Amplitude)
        IDs = list(simulation_data.keys())
        for id_s in IDs:
            ax.plot(simulation_data[id_s]['# t'], simulation_data[id_s]['z'], label=id_s)
            ax.grid()
        #ax.legend()
        
    #Compute error
    rms_list = list()
    mse_list = list()
    R2_list = list()
        
    for id_sim in tqdm(simulation_data.keys(), desc='Generating new summary...'):
        y_off = postprocess.computeVerticalOffset(simulation_data[id_sim], exp_data, 'Pin 1', 'y')
        exp_tmp = deepcopy(exp_data['Pin 1']['y_proc_valid'][['Time','Amplitude']].to_numpy())
        if correct_offset:
            exp_tmp[:,1] += y_off
            
        if debug_offset:
            exp_tmp_debug = deepcopy(exp_data)
            exp_tmp_debug['Pin 1']['y_proc_valid']['Amplitude'] += y_off
            print(postprocess.computeVerticalOffset(simulation_data[id_sim], exp_tmp_debug, 'Pin 1', 'y'))
            
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(simulation_data[id_sim][['# t','z']].to_numpy(), exp_tmp, debug=False, scale=1000, mean = True)
        rms_list.append(rms)
        mse_list.append(mse)
        R2_list.append(r2)
        
    if generate_new_csv:
        sim_summary['RMS'] = rms_list
        sim_summary['MSE'] = mse_list
        sim_summary['R2'] = R2_list
        sim_summary.to_csv(os.path.join(foldername, 'Summary_con_RMS.csv'))