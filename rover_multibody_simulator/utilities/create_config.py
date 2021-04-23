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


kg_mm2kg_m = 1e-6



#ROVER BODY PROPERTIES
config['Model Description']['rover_mass'] = '0.92'
config['Model Description']['rover_inertia_xx'] = str(56118.91*kg_mm2kg_m)
config['Model Description']['rover_inertia_yy'] = str(4372.79*kg_mm2kg_m)
config['Model Description']['rover_inertia_zz'] = str(8383.09*kg_mm2kg_m)
config['Model Description']['rover_inertia_xy'] = '0' # TO CHECK!
config['Model Description']['rover_inertia_xz'] = '0' # TO CHECK!
config['Model Description']['rover_inertia_yz'] = '0' # TO CHECK!







#BOGIE BODY PROPERTIES
config['Model Description']['arm_offset_x'] = '0.1328161' #TO CHECK!
config['Model Description']['arm_offset_y'] = '0.0' #TO CHECK! # 0 per semplicità
config['Model Description']['arm_offset_z'] = '0.0' #TO CHECK! #20mm forse, 15.5080mm
config['Model Description']['arm_cm_distance'] = '-0.004087'
config['Model Description']['link_lenght'] = str('0.11081344461055073')#str(0.1101876+0.002)#'0.1091876'   #TO CHECK!
config['Model Description']['link_angle'] = str(math.radians(13.818282155900022))    #14.44
config['Model Description']['hinge1_len'] = str(float(config['Model Description']['link_lenght'])\
                                                *math.cos(float(config['Model Description']['link_angle']))) #TO CHECK!
config['Model Description']['hinge1_height'] = str(float(config['Model Description']['link_lenght'])\
                                                *math.sin(float(config['Model Description']['link_angle']))) #TO CHECK!
config['Model Description']['swing_arm_mass'] = '0.168615'
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
config['Model Description']['central_link_mass'] = '0.0766'
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
config['Model Description']['terminal_link_mass'] = '0.38157'
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
config['Model Description']['wheel_offset_x'] = '0.10577'
config['Model Description']['wheel_offset_y'] = '0.16639'
config['Model Description']['wheel_offset_z'] = '-0.02713'
config['Model Description']['wheel_mass'] = '0.2242'
config['Model Description']['wheel_inertia_xx'] = str(1108.7749*kg_mm2kg_m)
config['Model Description']['wheel_inertia_yy'] = str(675.2718*kg_mm2kg_m)
config['Model Description']['wheel_inertia_zz'] = str(675.2718*kg_mm2kg_m)
config['Model Description']['wheel_radius'] = str(.085+0.005791804035500003) #'.0943709' #'.085'

config['Model Description']['simplified'] = 'False'
config['Model Description']['debugging'] = 'False'
config['Model Description']['action-reaction'] = 'True'
config['Model Description']['initial symbols substitution'] = 'False'
config['Model Description']['wheel-torque-control-debug'] = 'False'


#TORSIONAL SPRINGS PROPERTIES
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
config['Ground Properties']['stiffness'] = '1e6' #10000
config['Ground Properties']['damping'] = '1e2'
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
config['Simulator']['integrator'] = 'ipv - RK23'#ipv - RK23' #'ipv - Radau'#'odeint'    #[odeint, ipv - RK45, ipv - Adams, 'Euler-Explicit']
config['Simulator']['step size'] = '1e-5'

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
config['Controller']['frequency'] = '100'


config.add_section('Wheel-Steer-Controller')
config['Wheel-Steer-Controller']['K1'] = '.3'
config['Wheel-Steer-Controller']['K2'] = '.3'
config['Wheel-Steer-Controller']['K3'] = '.3'
config['Wheel-Steer-Controller']['K4'] = '.3'
config['Wheel-Steer-Controller']['D1'] = '.05'
config['Wheel-Steer-Controller']['D2'] = '.05'
config['Wheel-Steer-Controller']['D3'] = '.05'
config['Wheel-Steer-Controller']['D4'] = '.05'
config['Wheel-Steer-Controller']['I1'] = '.05'
config['Wheel-Steer-Controller']['I2'] = '.05'
config['Wheel-Steer-Controller']['I3'] = '.05'
config['Wheel-Steer-Controller']['I4'] = '.05'




config.add_section('Wheel-Drive-Controller')
config['Wheel-Drive-Controller']['K1'] = '.1'
config['Wheel-Drive-Controller']['K2'] = '.1'
config['Wheel-Drive-Controller']['K3'] = '.1'
config['Wheel-Drive-Controller']['K4'] = '.1'
config['Wheel-Drive-Controller']['D1'] = '0'
config['Wheel-Drive-Controller']['D2'] = '0'
config['Wheel-Drive-Controller']['D3'] = '0'
config['Wheel-Drive-Controller']['D4'] = '0'
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



