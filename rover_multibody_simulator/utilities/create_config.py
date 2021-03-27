# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 22:24:08 2021

@author: Matteo
"""

import configparser
import math
import os

config = configparser.ConfigParser()

config.add_section('Model Description')

config['Model Description']['rover_mass'] = '1'
config['Model Description']['rover_inertia_xx'] = '0.0075'
config['Model Description']['rover_inertia_yy'] = '0.0075'
config['Model Description']['rover_inertia_zz'] = '0.0075'

config['Model Description']['arm_offset'] = '0.3'
config['Model Description']['arm_cm_distance'] = '0.02'
config['Model Description']['link_lenght'] = '0.2'
config['Model Description']['hinge1_len'] = str(0.2*math.cos(math.radians(20)))
config['Model Description']['hinge1_height'] = str(0.2*math.sin(math.radians(20)))
config['Model Description']['swing_arm_mass'] = '0.1'
config['Model Description']['swing_arm_inertia_xx'] = '0.0005'
config['Model Description']['swing_arm_inertia_yy'] = '0.0005'
config['Model Description']['swing_arm_inertia_zz'] = '0.0005'
config['Model Description']['swing_arm_inertia_xy'] = '0'
config['Model Description']['swing_arm_inertia_xz'] = '0'
config['Model Description']['swing_arm_inertia_yz'] = '0'
config['Model Description']['link_cm_dist'] = '0.1'
config['Model Description']['link_mass'] = '0.02'
config['Model Description']['link_inertia_xx'] = '0.00025'
config['Model Description']['link_inertia_yy'] = '0.00025'
config['Model Description']['link_inertia_zz'] = '0.00025'
config['Model Description']['wheel_offset'] = '0.1'
config['Model Description']['wheel_mass'] = '0.05'
config['Model Description']['wheel_inertia_xx'] = '0.0075'
config['Model Description']['wheel_inertia_yy'] = '0.0075'
config['Model Description']['wheel_inertia_zz'] = '0.0075'
config['Model Description']['wheel_radius'] = '.1'

config['Model Description']['simplified'] = 'True'
config['Model Description']['debugging'] = 'False'
config['Model Description']['action-reaction'] = 'True'
config['Model Description']['initial symbols substitution'] = 'True'



config.add_section('Springs Definition')
config['Springs Definition']['k_fr1'] = '10.5'
config['Springs Definition']['M0_fr1'] = '1'
config['Springs Definition']['E_fr1'] = '1.5'
config['Springs Definition']['f_fr1'] = '133'
config['Springs Definition']['c_fr1'] = '5e-3'

config['Springs Definition']['k_fr2'] = '10.5'
config['Springs Definition']['M0_fr2'] = '1'
config['Springs Definition']['E_fr2'] = '1.5'
config['Springs Definition']['f_fr2'] = '133'
config['Springs Definition']['c_fr2'] = '5e-3'

config['Springs Definition']['k_br1'] = '10.5'
config['Springs Definition']['M0_br1'] = '1'
config['Springs Definition']['E_br1'] = '1.5'
config['Springs Definition']['f_br1'] = '133'
config['Springs Definition']['c_br1'] = '5e-3'

config['Springs Definition']['k_br2'] = '10.5'
config['Springs Definition']['M0_br2'] = '1'
config['Springs Definition']['E_br2'] = '1.5'
config['Springs Definition']['f_br2'] = '133'
config['Springs Definition']['c_br2'] = '5e-3'

config['Springs Definition']['k_fl1'] = '10.5'
config['Springs Definition']['M0_fl1'] = '1'
config['Springs Definition']['E_fl1'] = '1.5'
config['Springs Definition']['f_fl1'] = '133'
config['Springs Definition']['c_fl1'] = '5e-3'

config['Springs Definition']['k_fl2'] = '10.5'
config['Springs Definition']['M0_fl2'] = '1'
config['Springs Definition']['E_fl2'] = '1.5'
config['Springs Definition']['f_fl2'] = '133'
config['Springs Definition']['c_fl2'] = '5e-3'

config['Springs Definition']['k_bl1'] = '10.5'
config['Springs Definition']['M0_bl1'] = '1'
config['Springs Definition']['E_bl1'] = '1.5'
config['Springs Definition']['f_bl1'] = '133'
config['Springs Definition']['c_bl1'] = '5e-3'

config['Springs Definition']['k_bl2'] = '10.5'
config['Springs Definition']['M0_bl2'] = '1'
config['Springs Definition']['E_bl2'] = '1.5'
config['Springs Definition']['f_bl2'] = '133'
config['Springs Definition']['c_bl2'] = '5e-3'




config.add_section('Ground Properties')
config['Ground Properties']['stiffness'] = '10000'
config['Ground Properties']['damping'] = '10'
config['Ground Properties']['exp'] = '1.5'
config['Ground Properties']['max depth'] = '1e-4'



config.add_section('Plotting')
config['Plotting']['center mass'] = 'True'
config['Plotting']['live animation'] = 'True'
config['Plotting']['save animation'] = 'True'
config['Plotting']['animation name'] = 'auto'
config['Plotting']['plot contact point'] = 'True' #For debugging porpouse
config['Plotting']['animation fps'] = '30'





config.add_section('Simulator')
config['Simulator']['simulation_time'] = '5'
config['Simulator']['eom_method'] = 'Kane'
config['Simulator']['integrator'] = 'ipv - RK45'#'odeint'    #[odeint, ipv - RK45, ipv - Adams, 'Euler-Explicit']
config['Simulator']['step size'] = '1e-5'

#Used only for odeint
config['Simulator']['min step size'] = 'auto'
config['Simulator']['max stiff order'] = 'auto'
config['Simulator']['max non-stiff order'] = 'auto'
config['Simulator']['lambdifier'] = 'lambdify' #options = ['lambdify', theano, autowrap - f2py, autowrap - cython]
config['Simulator']['autowrap save files'] = 'True'
config['Simulator']['autowrap save dir'] = 'data/wrappers'


file_path = os.path.dirname(__file__)
parent_folder = os.path.abspath(os.path.join(file_path,'..'))
   
if not os.path.isdir(os.path.join(parent_folder,'data','config')):
    os.mkdir(os.path.join(parent_folder,'data','config'))
    
    
with open(os.path.join(parent_folder, 'data','config','config.ini'), 'w') as configfile:
    config.write(configfile)



