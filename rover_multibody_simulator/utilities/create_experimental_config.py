# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 23:19:06 2021

@author: Matteo
"""

import configparser
import os
import numpy as np


config = configparser.ConfigParser()
config.add_section('General')
config['General']['scale'] = '10.07383'
config.add_section('Pin 1 Y')
vertical_offset = -7/1000
config['Pin 1 Y']['n_initial_data_to_skip'] = '5'
config['Pin 1 Y']['n_final_data_to_skip'] = '5'
config['Pin 1 Y']['y_offset'] = str(0.2832 + vertical_offset)
config['Pin 1 Y']['time_offset'] = '0.77'


file_path = os.path.dirname(__file__)
parent_folder = os.path.abspath(os.path.join(file_path,'..'))
   
if not os.path.isdir(os.path.join(parent_folder,'data','config')):
    os.makedirs(os.path.join(parent_folder,'data','config'))
    
    
with open(os.path.join(parent_folder, 'data','config','config_experimental.ini'), 'w') as configfile:
    config.write(configfile)