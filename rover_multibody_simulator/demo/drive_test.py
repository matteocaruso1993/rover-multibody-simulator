# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 15:29:59 2021

@author: Matteo Caruso
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
from rover_multibody_simulator.utilities import functions
import argparse
import sys
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt

sim_step_time = .001
sim_time = 7#☺.2 #2
#sim_time = 4 #1ms
#sim_time = 3 #2ms
drop_height = .5

debug = False
skip_confirm = True
animate = False
print_simulation_output = True
save_config = True



if __name__ == '__main__':
    sim = RoverSimulator()
    s = 'data/wrappers/full_model'
    sys.path.insert(0,os.path.join(os.getcwd(),*s.split('/')))
    sim.loadLambdaFunctions(model_name='full-parametric-model', wrapper_folder='data/wrappers/full_model')
    #sim.loadConfig()
    sim.loadGroundPoints()
    sim.loadFrictionModel()
    
    sim.config.set('Ground Interpolator','kx','1')
    sim.config.set('Ground Interpolator','ky','1')
    
    
    sim.changeTerrain('ramp20_v1')
    
    #Initialize wheel controllers
    #sim.wheel_controller.setSteerTarget([0,0,0,0])
    #sim.wheel_controller.setDriveTarget([-4,-4,-4,-4])
    #sim.wheel_controller.setDriveTarget(4*[-2*11.76470588235294])
    steer_idx = sim.wheel_controller.getWheelIndexes()['steer']
    drive_idx = sim.wheel_controller.getWheelIndexes()['drive']
    
    
    sim.loadInitialConditions('test_stiffer','full', subs_vel=True)
    
    gen = sim.current_gen_coord
    spe = sim.current_gen_speed
    gen[1] = 0 #0.3m/s
    #gen[1] -= 0.5 #1m/s
    #gen[1] -= 1 #2m/s 
    #gen[2] += 0.01
    sim.setInitialConditions(gen, len(gen)*[1e-3])
    
    #sim.setConstantSpeed(steer_idx + drive_idx, 4*[0] + 4*[-2*11.76470588235294]) #2 m/s
    sim.setConstantSpeed(steer_idx + drive_idx, 4*[0] + 4*[-4]) #0.3 m/s
    #sim.setConstantSpeed(steer_idx + drive_idx, 4*[0] + 4*[-1*11.76470588235294]) #1 m/s
    
    
    sim.config.set('Ground Properties','damping','40')
    sim.config.set('Ground Properties','damping','100')
    sim.config.set('Ground Properties','exp','1.5')
    sim.config.set('Ground Properties','stiffness','1e8')
    
    
    
    coeff = {'g1': 1.8097017139445475,
             'g2': 222.22459836588484,
             'g3': 185.03061961596708,
             'g4': 0.21,
             'g5': 204.91424590571341,
             'g6': 0} #ms 0.3 md 0.2
    
    
    ID_springs = ['fr','br','fl','bl']
    stiffness_scale = 1
    preload_scale = .9
    f_times_e = 10000
    scalec = 1#.9      #100
    c_spring_new = .1*scalec
    scale_second_module = 1
    for ID in ID_springs:
        #Change Spring damping
        e = np.sqrt(f_times_e)
        f = np.sqrt(f_times_e)
        
        for n in [1,2]:
            if n == 1:
                sim.config.set('Springs Definition', 'c_'+ ID+str(n), str(c_spring_new))
            elif n==2:
                sim.config.set('Springs Definition', 'c_'+ ID+str(n), str(c_spring_new*scale_second_module))
            sim.config.set('Springs Definition', 'k_' + ID+str(n), str(stiffness_scale*sim.config.getfloat('Springs Definition', 'k_fr1')))
            sim.config.set('Springs Definition', 'M0_' + ID+str(n), str(preload_scale*sim.config.getfloat('Springs Definition', 'M0_' + ID+str(n))))
            sim.config.set('Springs Definition', 'E_' + ID+str(n), str(e))
            sim.config.set('Springs Definition', 'f_' + ID+str(n), str(f))
            
            
            print('Current spring damping is: \t {:10.6f}'.format(sim.config.getfloat('Springs Definition', 'c_'+ ID+str(n))))
            print('Current spring stiffness is: \t {:10.6f}'.format(sim.config.getfloat('Springs Definition', 'k_'+ ID+str(n))))
            print('Current spring preload is: \t {:10.6f}'.format(sim.config.getfloat('Springs Definition', 'M0_'+ ID+str(n))))
    
    if not debug:
        sim.friction.setCoeffiecients(*list(coeff.values()))
        print(sim.friction.getCoefficients())
        print('Current floor damping is: \t {:10.6f}'.format(sim.config.getfloat('Ground Properties','damping')))
        if not skip_confirm:    
            input('Press Enter to start the simulation')
        sim._createParamMapping()
        sim.simulate(sim_step_time, sim_time, out = print_simulation_output)
        
        
        
        soluzione = sim.ode_sol
        out = np.vstack((soluzione['t'],soluzione['y']))
        header = ['t']
        for c in sim.gen_coord:
            header.append(c.name)
            
        for c in sim.gen_speeds:
            header.append(c.name)
            
        my_string = ','.join(map(str,header))
        
        num = 8
        while True:
            filename = 'sim_py_v' + str(num) + '_transpose_con_header.csv'
            if os.path.exists(filename):
                num += 1
            else:
                break
                    
        
        np.savetxt(filename, out.T, delimiter=',', header=my_string)
        
        cwd = os.getcwd();
        save_dir = r'C:\Users\Matteo\Desktop\RoverMatlab\Experiment'
        
        if save_config:
            functions.saveConfig(sim, 'config.ini')
        
        # shutil.copyfile(os.path.join(cwd,filename), os.path.join(save_dir, filename))
        
        # with open(os.path.join(save_dir,'README.txt'), 'a') as f:
        #     string_to_save1 = 'v'+str(num)+" "+ 'k_ground=%10.6f, c_ground=%10.6f, c_springs=%10.6f,'%(sim.config.getfloat('Ground Properties','stiffness'),sim.config.getfloat('Ground Properties','damping'), c_spring_new) + " "
        #     string_to_save2 = ''
        #     for key in sim.friction.getFrictionProperties().keys():
        #         string_to_save2 += key + "=" + str(sim.friction.getFrictionProperties()[key]) + ", "
                
        #     f.write('\n' + string_to_save1 + string_to_save2[:-2])
            
        # f.close()
            
        
        
        if animate:
            sim.animateMotion()