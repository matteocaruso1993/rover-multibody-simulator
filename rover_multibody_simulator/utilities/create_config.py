# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 22:24:08 2021

@author: Matteo
"""

import configparser
import math
import os
import json


#model = 'simple_model'
model = 'full-parametric-model'
copy_to_model = True

config = configparser.ConfigParser()

config.add_section('Model Description')


kg_mm2kg_m = 1e-6
leg_scale_mass = 1.0304335182393098



#ROVER BODY PROPERTIES
config['Model Description']['rover_mass'] = '1.267' #'0.92'
config['Model Description']['rover_inertia_xx'] = str(56118.91/10*kg_mm2kg_m)
config['Model Description']['rover_inertia_yy'] = str(4372.79*kg_mm2kg_m)
config['Model Description']['rover_inertia_zz'] = str(8383.09*kg_mm2kg_m)
config['Model Description']['rover_inertia_xy'] = '0' # TO CHECK!
config['Model Description']['rover_inertia_xz'] = '0' # TO CHECK!
config['Model Description']['rover_inertia_yz'] = '0' # TO CHECK!







#BOGIE BODY PROPERTIES
config['Model Description']['arm_offset_x'] = str(0.1328161 + 0.49128/1000)#TO CHECK!
config['Model Description']['arm_offset_y'] = '0.0' #TO CHECK! # 0 per semplicità
config['Model Description']['arm_offset_z'] = '0.0' #TO CHECK! #20mm forse, 15.5080mm
config['Model Description']['arm_cm_distance'] = '-0.004087'
config['Model Description']['link_lenght'] = str('0.11081344461055073')#str(0.1101876+0.002)#'0.1091876'   #TO CHECK!
config['Model Description']['link_angle'] = str(math.radians(13.818282155900022))    #14.44
config['Model Description']['hinge1_len'] = str(float(config['Model Description']['link_lenght'])\
                                                *math.cos(float(config['Model Description']['link_angle']))) #TO CHECK!
config['Model Description']['hinge1_height'] = str(float(config['Model Description']['link_lenght'])\
                                                *math.sin(float(config['Model Description']['link_angle']))) #TO CHECK!
config['Model Description']['swing_arm_mass'] = str(leg_scale_mass*0.168615)
config['Model Description']['swing_arm_inertia_xx'] = str(92.454*kg_mm2kg_m)
config['Model Description']['swing_arm_inertia_yy'] = str(740.29055*kg_mm2kg_m)
config['Model Description']['swing_arm_inertia_zz'] = str(690.80801*kg_mm2kg_m)
config['Model Description']['swing_arm_inertia_xy'] = str(0.20764*kg_mm2kg_m)
config['Model Description']['swing_arm_inertia_xz'] = str(0)
config['Model Description']['swing_arm_inertia_yz'] = str(1.10413*kg_mm2kg_m)


#CENTRAL LINK PROPERTIES
config['Model Description']['central_link_cm_dist_x'] = '0.0'
config['Model Description']['central_link_cm_dist_y'] = '0.050808'  #'0.048808'
config['Model Description']['central_link_cm_dist_z'] = '0.009641'
config['Model Description']['central_link_lenght'] = str(.1)#str(0.09139+0.004)
config['Model Description']['central_link_mass'] = str(leg_scale_mass*0.0766)
config['Model Description']['central_link_inertia_xx'] = str(39.7567*kg_mm2kg_m)
config['Model Description']['central_link_inertia_yy'] = str(83.8321*kg_mm2kg_m)
config['Model Description']['central_link_inertia_zz'] = str(58.4776*kg_mm2kg_m)
config['Model Description']['central_link_inertia_xy'] = str(2.49E-02*kg_mm2kg_m)
config['Model Description']['central_link_inertia_xz'] = str(-3.729*kg_mm2kg_m)
config['Model Description']['central_link_inertia_yz'] = str(-3.11E-02*kg_mm2kg_m)


#TERMINAL LINK PROPERTIES
config['Model Description']['terminal_link_cm_dist_x'] = '0.06471'
config['Model Description']['terminal_link_cm_dist_y'] = '0.1413'
config['Model Description']['terminal_link_cm_dist_z'] = '-0.01137'
config['Model Description']['terminal_link_lenght'] = str(0.16739+0.002)  #'0.16639'
config['Model Description']['terminal_link_mass'] = str(leg_scale_mass*0.38157)
config['Model Description']['terminal_link_inertia_xx'] = str(39.7567*kg_mm2kg_m)
config['Model Description']['terminal_link_inertia_yy'] = str(83.8321*kg_mm2kg_m)
config['Model Description']['terminal_link_inertia_zz'] = str(58.4776*kg_mm2kg_m)
config['Model Description']['terminal_link_inertia_xy'] = str(2.49E-02*kg_mm2kg_m)
config['Model Description']['terminal_link_inertia_xz'] = str(-3.729*kg_mm2kg_m)
config['Model Description']['terminal_link_inertia_yz'] = str(-3.11E-02*kg_mm2kg_m)


#config['Model Description']['link_cm_dist'] = '0.1'
#config['Model Description']['link_mass'] = '0.02'
#config['Model Description']['link_inertia_xx'] = '0.00025'
#config['Model Description']['link_inertia_yy'] = '0.00025'
#config['Model Description']['link_inertia_zz'] = '0.00025'



#WHEEL BODY PROPERTIES
config['Model Description']['wheel_offset_x'] = str(0.087 + 2/1000)#'#'0.10577' #!!!!!!!!!!!
config['Model Description']['wheel_offset_y'] = str(0.166394 +1.4/1000) #!!!!!!!!!!!!!!!
config['Model Description']['wheel_offset_z'] = str(-0.02713 - 5/1000-0.6/1000) #!!!!!!!!!!!!!!!!!!
config['Model Description']['wheel_mass'] = str(leg_scale_mass*0.2242)
config['Model Description']['wheel_inertia_xx'] = str(1108.7749*kg_mm2kg_m)
config['Model Description']['wheel_inertia_yy'] = str(675.2718*kg_mm2kg_m)
config['Model Description']['wheel_inertia_zz'] = str(675.2718*kg_mm2kg_m)
config['Model Description']['wheel_radius'] = '.085'#str(.085+0.005791804035500003) #'.0943709' #'.085'

config['Model Description']['simplified'] = 'False'
config['Model Description']['debugging'] = 'False'
config['Model Description']['action-reaction'] = 'True'
config['Model Description']['initial symbols substitution'] = 'False'
config['Model Description']['wheel-torque-control-debug'] = 'False'


#TORSIONAL SPRINGS PROPERTIES
config.add_section('Springs Definition')
config['Springs Definition']['k_fr1'] = '9.1875' #10.5
config['Springs Definition']['M0_fr1'] = '2.1384'# '1'
config['Springs Definition']['E_fr1'] = '1.5'
config['Springs Definition']['f_fr1'] = '133'
config['Springs Definition']['c_fr1'] = '5e-3'
config['Springs Definition']['M0_fr1_stdv'] = '0.0434'

config['Springs Definition']['k_fr2'] = '9.1875'
config['Springs Definition']['M0_fr2'] = '1.4842'
config['Springs Definition']['E_fr2'] = '1.5'
config['Springs Definition']['f_fr2'] = '133'
config['Springs Definition']['c_fr2'] = '5e-3'
config['Springs Definition']['M0_fr2_stdv'] = '0.0505'

config['Springs Definition']['k_br1'] = '9.1875'
config['Springs Definition']['M0_br1'] = '2.2756'
config['Springs Definition']['E_br1'] = '1.5'
config['Springs Definition']['f_br1'] = '133'
config['Springs Definition']['c_br1'] = '5e-3'
config['Springs Definition']['M0_br1_stdv'] = '0.0838'

config['Springs Definition']['k_br2'] = '9.1875'
config['Springs Definition']['M0_br2'] = '1.2415'
config['Springs Definition']['E_br2'] = '1.5'
config['Springs Definition']['f_br2'] = '133'
config['Springs Definition']['c_br2'] = '5e-3'
config['Springs Definition']['M0_br2_stdv'] = '0.0357'

config['Springs Definition']['k_fl1'] = '9.1875'
config['Springs Definition']['M0_fl1'] = '2.329'
config['Springs Definition']['E_fl1'] = '1.5'
config['Springs Definition']['f_fl1'] = '133'
config['Springs Definition']['c_fl1'] = '5e-3'
config['Springs Definition']['M0_fl1_stdv'] = '0.0985'

config['Springs Definition']['k_fl2'] = '9.1875'
config['Springs Definition']['M0_fl2'] = '1.3572'
config['Springs Definition']['E_fl2'] = '1.5'
config['Springs Definition']['f_fl2'] = '133'
config['Springs Definition']['c_fl2'] = '5e-3'
config['Springs Definition']['M0_fl2_stdv'] = '0.0341'


config['Springs Definition']['k_bl1'] = '9.1875'
config['Springs Definition']['M0_bl1'] = '1.9174'
config['Springs Definition']['E_bl1'] = '1.5'
config['Springs Definition']['f_bl1'] = '133'
config['Springs Definition']['c_bl1'] = '5e-3'
config['Springs Definition']['M0_bl1_stdv'] = '0.0771'


config['Springs Definition']['k_bl2'] = '9.1875'
config['Springs Definition']['M0_bl2'] = '1.6447'
config['Springs Definition']['E_bl2'] = '1.5'
config['Springs Definition']['f_bl2'] = '133'
config['Springs Definition']['c_bl2'] = '5e-3'
config['Springs Definition']['M0_bl2_stdv'] = '0.0589'



config.add_section('Ground Properties')
config['Ground Properties']['stiffness'] = '1e6' #10000
config['Ground Properties']['damping'] = '1e2'
config['Ground Properties']['exp'] = '1.5'
config['Ground Properties']['max depth'] = '1e-4'


config.add_section('Ground Interpolator')
config['Ground Interpolator']['kx'] = '3'
config['Ground Interpolator']['ky'] = '3'
config['Ground Interpolator']['s'] = '0'



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
config['Simulator']['integrator'] = 'ipv - RK23'#ipv - RK23' #'ipv - Radau'#'odeint'    #[odeint, ipv - RK45, ipv - Adams, 'Euler-Explicit']
config['Simulator']['step size'] = '1e-5'
config['Simulator']['rtol'] = '1e-3'
config['Simulator']['atol'] = '1e-3'

#Used only for odeint
config['Simulator']['min step size'] = 'auto'
config['Simulator']['max stiff order'] = 'auto'
config['Simulator']['max non-stiff order'] = 'auto'
config['Simulator']['lambdifier'] = 'lambdify' #options = ['lambdify', theano, autowrap - f2py, autowrap - cython]
config['Simulator']['autowrap save files'] = 'True'
config['Simulator']['autowrap save dir'] = 'data/wrappers/full_model'

config.add_section('Friction Properties')
config['Friction Properties']['load_model'] = 'True'
config['Friction Properties']['g1'] = '2.0430'
config['Friction Properties']['g2'] = '239.1462'
config['Friction Properties']['g3'] = '120.0534'
config['Friction Properties']['g4'] = '0.5'
config['Friction Properties']['g5'] = '84.4111'
config['Friction Properties']['g6'] = '0'

config.add_section('Controller')
config['Controller']['frequency'] = '1000'


config.add_section('Wheel-Steer-Controller')
config['Wheel-Steer-Controller']['K1'] = '5' #0.3
config['Wheel-Steer-Controller']['K2'] = '5'
config['Wheel-Steer-Controller']['K3'] = '5'
config['Wheel-Steer-Controller']['K4'] = '5'
config['Wheel-Steer-Controller']['D1'] = '.2'
config['Wheel-Steer-Controller']['D2'] = '.2'
config['Wheel-Steer-Controller']['D3'] = '.2'
config['Wheel-Steer-Controller']['D4'] = '.2'
config['Wheel-Steer-Controller']['I1'] = '.2'
config['Wheel-Steer-Controller']['I2'] = '.2'
config['Wheel-Steer-Controller']['I3'] = '.2'
config['Wheel-Steer-Controller']['I4'] = '.2' #0.05




config.add_section('Wheel-Drive-Controller')
config['Wheel-Drive-Controller']['K1'] = '7'
config['Wheel-Drive-Controller']['K2'] = '7'
config['Wheel-Drive-Controller']['K3'] = '7'
config['Wheel-Drive-Controller']['K4'] = '7'
config['Wheel-Drive-Controller']['D1'] = '0.0000'
config['Wheel-Drive-Controller']['D2'] = '0.0000'
config['Wheel-Drive-Controller']['D3'] = '0.0000'
config['Wheel-Drive-Controller']['D4'] = '0.0000'
config['Wheel-Drive-Controller']['I1'] = '0'
config['Wheel-Drive-Controller']['I2'] = '0'
config['Wheel-Drive-Controller']['I3'] = '0'
config['Wheel-Drive-Controller']['I4'] = '0'



file_path = os.path.dirname(__file__)
parent_folder = os.path.abspath(os.path.join(file_path,'..'))
   
if not os.path.isdir(os.path.join(parent_folder,'data','config')):
    os.makedirs(os.path.join(parent_folder,'data','config'))
    
    
with open(os.path.join(parent_folder, 'data','config','config.ini'), 'w') as configfile:
    config.write(configfile)
    
    
if copy_to_model:
    with open(os.path.join(parent_folder, 'data','model',model,'config.ini'), 'w') as configfile:
        config.write(configfile)



