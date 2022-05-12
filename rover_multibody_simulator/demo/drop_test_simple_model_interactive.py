# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 21:47:48 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
import argparse
import sys
import os
import numpy as np
import shutil

sim_step_time = .01
sim_time = 1
drop_height = .5

debug = True
skip_confirm = True
animate = False
print_simulation_output = True



if __name__ == '__main__':
    sim = RoverSimulator()
    s = 'data/wrappers/simple_model'
    sys.path.insert(0,os.path.join(os.getcwd(),*s.split('/')))
    sim.loadLambdaFunctions(model_name='simple_model', wrapper_folder='data/wrappers/simple_model')
    #sim.loadConfig()
    sim.loadGroundPoints()
    sim.loadFrictionModel()
    gen_coord = len(sim.gen_coord)*[0]
    sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
    current_wheel_height = sim.getLowestWheelPoint()[1]
    gen_coord[2] = drop_height - current_wheel_height
    #sim.config.set('Simulator', 'integrator','Euler-Explicit')
    
    #sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
    #sim.simulate(sim_step_time, sim_time)
    h_offset = -0.0028492199999999995 #-0.00960784#7/1000
    #gen_coord[2] = 0.2762 + h_offset
    #gen_coord[4] = np.deg2rad(2.26) #0.0202*20
    #gen_coord[5] = 0.0042
    #gen_coord[6] = -0.0189
    gen_coord[2] = 0.28104992000000006
    gen_coord[4] = 0.03944444109507184#0.0202
    gen_coord[5] = 0.0042
    gen_coord[6] =  0.0189
    sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
    sim.config.set('Ground Properties','damping','40')
    sim.config.set('Ground Properties','damping','100')
    #sim.config.set('Ground Properties','stiffness','1000000')
    """
    
    coeff = {'g1': 1.8097017139445475,
             'g2': 222.22459836588484,
             'g3': 185.03061961596708,
             'g4': 0.21,
             'g5': 204.91424590571341,
             'g6': 0} #ms 0.3 md 0.2"""
    
    coeff = {'g1': 0,
             'g2': 0,
             'g3': 0,
             'g4': 0,
             'g5': 0,
             'g6': 0}#0"""
    
    """
    coeff = {'g1': 6.2854693753104245,
             'g2': 230.8206645931289,
             'g3': 213.26760031726008,
             'g4': 0.3993,
             'g5': 199.6578038398171,
             'g6': 0} #ms= .4 md= .3"""
    
    ID_springs = ['fr','br','fl','bl']
    stiffness_scale = 1
    preload_scale = 1
    f_times_e = 10000
    scalec = 1#.9      #100
    c_spring_new = .01*scalec
    scale_second_module = .5
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
        print('The rover centre of mass current height is: \t {:10.6f} [m]'.format(gen_coord[2]))
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
        
        shutil.copyfile(os.path.join(cwd,filename), os.path.join(save_dir, filename))
        
        with open(os.path.join(save_dir,'README.txt'), 'a') as f:
            string_to_save1 = 'v'+str(num)+" "+ 'k_ground=%10.6f, c_ground=%10.6f, c_springs=%10.6f,'%(sim.config.getfloat('Ground Properties','stiffness'),sim.config.getfloat('Ground Properties','damping'), c_spring_new) + " "
            string_to_save2 = ''
            for key in sim.friction.getFrictionProperties().keys():
                string_to_save2 += key + "=" + str(sim.friction.getFrictionProperties()[key]) + ", "
                
            f.write('\n' + string_to_save1 + string_to_save2[:-2])
            
        f.close()
            
        
        
        if animate:
            sim.animateMotion()