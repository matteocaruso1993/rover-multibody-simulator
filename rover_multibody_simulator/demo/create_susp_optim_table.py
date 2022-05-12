# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:41:59 2022

@author: Matteo
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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

save_table = True



root = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\DesignOfExperimentsSuspensions\k6'

files_list = os.listdir(root)
folders_list = [f for f in files_list if os.path.isdir(os.path.join(root,f))]

ID_springs = ['fr1','fr2','br1','br2','fl1','fl2','bl1','bl2']
vars_to_extract = ['z','chi','psi']
header = ID_springs + vars_to_extract


M_id = [' M0_' + item for item in ID_springs]

for folder in folders_list:
    new_root = os.path.join(root,folder)
    path_to_csv_summary = os.path.join(new_root,'Summary.csv')
    
    print('Loading Summary from:\t %s\n'%path_to_csv_summary)
    sim_summary = postprocess.create_temporary_copy_csv(path_to_csv_summary)
    
    
    for sim in tqdm(range(sim_summary.shape[0]),desc='loading simulation files...'):
        sim_to_load = os.path.join(new_root, sim_summary['#filename'][sim])
        sim_data_tmp = postprocess.loadSummaryData(sim_to_load)
        springs_values = sim_summary[M_id].iloc[sim,:]
        
        if not vars_to_extract:
            vars_to_extract = sim_data_tmp.keys().to_list()
            vars_to_extract.remove(' #t')
        
        arr = np.zeros((len(vars_to_extract),))
        for n,key in enumerate(vars_to_extract):
            arr[n] = sim_data_tmp[key].iloc[-1]
            
        pd_to_append = pd.DataFrame(np.expand_dims(np.hstack((springs_values.to_numpy(), arr)),axis=0), columns=[s.strip(' ') for s in M_id] + vars_to_extract)
        
        if 'final_pd' not in locals():
            final_pd = deepcopy(pd_to_append)
        else:
            final_pd = final_pd.append(pd_to_append,ignore_index=True)
            


if save_table:
    final_pd.to_csv(os.path.join(root, 'data_table.csv'))
                
        

