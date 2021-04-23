# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 22:40:19 2021

@author: Matteo
"""
import os
from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n', help='Model name', required=True, dest='model_name', action='store')
parser.add_argument('-t', help='Model type', required=True, dest='model_type', action='store', choices=['full','simple'])
parser.add_argument('-s', help='Initial Substitution', required=True, dest='initial_sub', action='store',choices = ['True', 'False'])

args = parser.parse_args()


def main(model_name=''):
    if not model_name:
        raise ValueError('Model name must be specified')
    sim = RoverSimulator()
    tmp_str = sim.config.get('Simulator', 'autowrap save dir')
    tmp_dvd = tmp_str.split('/')
    tmp_dvd[-1] = model_name
    sim.config.set('Simulator', 'autowrap save dir', os.path.join(*tmp_dvd))
    if args.model_type == 'full':
        sim.config.set('Model Description', 'simplified', 'False')
        sim.config.set('Model Description', 'debugging', 'False')
    else:
        sim.config.set('Model Description', 'simplified', 'True')
        sim.config.set('Model Description', 'debugging', 'False')
        
    sim.config.set('Model Description', 'initial symbols substitution', args.initial_sub)
    
    sim.initialize()
    sim.lambdifyAll('autowrap - cython')
    sim.formEquationsOfMotion('autowrap - cython')
    sim.saveLambdaFunctions(model_name=model_name)
    
    
    
    
    
    
    
if __name__ == '__main__':
    main(args.model_name)