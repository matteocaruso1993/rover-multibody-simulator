# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 21:54:58 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
from rover_multibody_simulator.utilities import functions
import sys
import os
import numpy as np
from tqdm import tqdm
import datetime
import time
import re


sim_step_time = .01
sim_time = 1

debug = False
skip_confirm = True
animate = True

ID_springs = ['fr1','fr2','br1','br2','fl1','fl2','bl1','bl2']

M_id = ['M0_' + item for item in ID_springs]
C_id = ['c_' + item for item in ID_springs]


#foldername = r'C:\Users\Matteo\Desktop\RoverMatlab\Experiment\DesignOfExperiments'
foldername = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\DesignOfExperiments'
subfolder_id = 0

vars_to_save = ['#filename','config_name','date','sim_time','computation_time','rover_mass','drop_height', 'ground_damping', 'ground_stiffness', 'exp', 'f*g',
                'k_spring']

vars_to_save += M_id + C_id
    


if __name__ == '__main__':  
    while True:
        if os.path.exists(os.path.join(foldername, str(subfolder_id))):
            subfolder_id += 1
        else:
            os.makedirs(os.path.join(foldername, str(subfolder_id)))
            os.makedirs(os.path.join(foldername, str(subfolder_id), 'configs'))
            break
        
        
    
    sim = RoverSimulator()
    s = 'data/wrappers/simple_model'
    sys.path.insert(0,os.path.join(os.getcwd(),*s.split('/')))
    sim.loadLambdaFunctions(model_name='simple_model', wrapper_folder='data/wrappers/simple_model')
    sim.loadGroundPoints()
    sim.loadFrictionModel()
    
    coeff = {'g1': 0,
             'g2': 0,
             'g3': 0,
             'g4': 0,
             'g5': 0,
             'g6': 0}
    
    
    sim.friction.setCoeffiecients(*list(coeff.values()))
    fric_prop = sim.friction.getFrictionProperties()
    
    vars_to_save += list(coeff.keys()) + list(fric_prop.keys())
    
    with open(os.path.join(foldername,str(subfolder_id),'Summary.csv'), 'w') as f:
        string_to_write = ', '.join(map(str,vars_to_save))
        f.write(string_to_write)
    f.close()
        
    
    
    gen_coord = len(sim.gen_coord)*[0]
    
    variables_to_change = {'damping': {'section_name':'Ground Properties', 'min':10, 'max':100, 'default':sim.config.getfloat('Ground Properties','damping'),'variable':True},
                           'height': {'section_name':None, 'min':0.2762, 'max':0.30, 'default':0.2832, 'variable': False},
                           'c_springs': {'section_name': None, 'min':.01, 'max':1, 'default': sim.config.getfloat('Springs Definition','c_fr1'), 'variable': True},
                           'c_scale': {'section_name': None, 'min':.2, 'max':1, 'default': 1, 'variable': False},
                           'stiffness': {'section_name': 'Ground Properties','min':900000, 'max':1100000, 'default':sim.config.getfloat('Ground Properties','stiffness'), 'variable':False},
                           'exp': {'section_name': 'Ground Properties', 'min':1, 'max':2.2, 'default':sim.config.getfloat('Ground Properties','exp'), 'variable':True},
                           'f_x_E':{'section_name': 'Springs Definition','min':100, 'max':2000, 'default':sim.config.getfloat('Springs Definition','E_fr1')*sim.config.getfloat('Springs Definition','f_fr1'), 'variable':True}}
    
    variables_to_change = functions.generateDOE(variables_to_change,600)
    
    
    for i in tqdm(range(variables_to_change['damping']['values'].shape[0]), desc='Solving'):
        gen_coord[2] = variables_to_change['height']['values'][i]
        gen_coord[4] = 0.0202
        gen_coord[5] = 0.0042
        gen_coord[6] = -0.0189
        sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
        dg = variables_to_change['damping']['values'][i]
        kg = variables_to_change['stiffness']['values'][i]
        e = variables_to_change['exp']['values'][i]
        
        E = np.sqrt(variables_to_change['f_x_E']['values'][i])
        f = np.sqrt(variables_to_change['f_x_E']['values'][i])
        
        
        sim.config.set('Ground Properties','damping',str(dg))
        sim.config.set('Ground Properties','damping',str(dg))
        sim.config.set('Ground Properties','exp',str(e))
        
        
        
        scalec = variables_to_change['c_scale']['values'][i]
        c_spring_new = variables_to_change['c_springs']['values'][i]
        for ID in ID_springs:
            #Change Spring damping
            #for n in [1,2]:
            
            if ID[-1] == '1':
                sim.config.set('Springs Definition', 'c_'+ ID, str(c_spring_new))
            elif ID[-1]=='2':
                sim.config.set('Springs Definition', 'c_'+ ID, str(c_spring_new*scalec))
                
            sim.config.set('Springs Definition', 'E_'+ ID, str(E))
            sim.config.set('Springs Definition', 'f_'+ ID, str(f))
                
        
        
        sim._createParamMapping()
        
        time_start = time.time()
        sim.simulate(sim_step_time, sim_time, out=False)
        time_end = time.time()
        
        filename = functions.saveSimulationData(sim, os.path.join(foldername, str(subfolder_id)))
        
        time_date = datetime.datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
        
        #vars_to_save = ['filename','date','rover_mass','drop_height', 'ground_damping', 'ground_stiffness', 'exp', 'f*g',
        #        'k_spring']
        f_times_g = sim.config.getfloat('Springs Definition','E_fr1')*sim.config.getfloat('Springs Definition','f_fr1')
        #Save config file
        sim_num = re.findall(r'\d+', filename[1])
        config_name = os.path.join(foldername,str(subfolder_id), 'configs', 'config_v'+ sim_num[0] +'.ini')
        functions.saveConfig(sim, config_name)
        var_data = [filename[1], os.path.split(config_name)[1],time_date, sim_time,time_end-time_start, sim.config.get('Model Description','rover_mass'), gen_coord[2], dg, kg, e, f_times_g, sim.config.getfloat('Springs Definition','k_fr1')]
        
        M_to_save = [sim.config.get('Springs Definition',item) for item in M_id]
        c_to_save = [sim.config.get('Springs Definition',item) for item in C_id]
        
        var_data += M_to_save + c_to_save + list(coeff.values()) + list(fric_prop.values())
        
        
        functions.saveVars(os.path.join(foldername,str(subfolder_id),'Summary.csv'), var_data)
        
        
        
        sim.reset()
        
    
    