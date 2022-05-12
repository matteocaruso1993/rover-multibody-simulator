# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 13:10:52 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
from rover_multibody_simulator.utilities import postprocess
import re
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm

sim = RoverSimulator()
s = 'data/wrappers/simple_model'
sys.path.insert(0,os.path.join(os.getcwd(),*s.split('/')))
sim.loadLambdaFunctions(model_name='simple_model', wrapper_folder='data/wrappers/simple_model')
sim.loadGroundPoints()
sim.loadFrictionModel()


doe_folder_path = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\DesignOfExperiments'
folders_to_open = [name for name in os.listdir(doe_folder_path) if os.path.isdir(os.path.join(doe_folder_path,name))]

pin_names = ['P_{SA_r}', 'P_{SA_l}']
pin_idx = list()
for pin_name in pin_names:
    for n, pt in enumerate(sim.points):
        if pt.name == pin_name:
            pin_idx.append(n)
            break

pins = dict(zip(pin_names, pin_idx))
#pins['P_{SA_r}']['name'] = 'z_pin_r'


for folder in tqdm(folders_to_open, desc='Preprocessing data'):
    files = os.listdir(os.path.join(doe_folder_path,folder))
    valid_filename = list()
    for file in files:
        if not 'sim_py_v' in file:
            continue
        else:
            valid_filename.append(file)
    
    
    #Start to open files
    new_folder_name = os.path.join(doe_folder_path,folder,'WithPinsHeight')
    if not os.path.exists(new_folder_name):
        os.makedirs(new_folder_name)
    
    
    
    
    
    for filename in valid_filename:
        data = postprocess.loadSummaryData(os.path.join(doe_folder_path,folder,filename))
        num = re.findall(r'\d+',filename)[0]
        
        config_file = os.path.join(doe_folder_path, folder, 'configs', 'config_v' + num +'.ini')
        if os.path.exists(config_file):
            sim.config_path = config_file
        
        try:
            sim.loadConfig(output=False)
            sim._createParamMapping()
                                 
            par = sim._map_parameters['vars-to-sub']
            z_pin_r = list()
            z_pin_l = list()
            
        
            for i in range(data.shape[0]):
                var = data.iloc[i, 1:16]
                for key in pins.keys():
                    pt_c = sim.lambdified_points[pins[key]](*np.hstack((var, par)))
                    if key == 'P_{SA_r}':
                        z_pin_r.append(pt_c[0,2])
                    else:
                        z_pin_l.append(pt_c[0,2])
                    
            
            new_data = np.hstack((data.to_numpy(), 
                                  np.array(z_pin_r).reshape(-1,1), np.array(z_pin_l).reshape(-1,1)))
            
            original_names = list(data.keys()) + ['z_pin_r','z_pin_l']
            data1 = pd.DataFrame(new_data, columns=original_names)
            
            data1.to_csv(os.path.join(new_folder_name, filename))
        except:
            print('Skipping DOE:\t'+folder+', simulation:\t'+num)
            
    
    #break
    
    
    
    
