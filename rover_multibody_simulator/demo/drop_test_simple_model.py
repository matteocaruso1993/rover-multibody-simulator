# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 18:05:24 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator
import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('-d', type=float, help='Drop height of the closest wheel to the ground', required=True, dest='drop_height', action='store')
parser.add_argument('-t', type=float, help='Simulation time', required=True, dest='sim_time', action='store')
parser.add_argument('-st', type=float, help = 'Simulation step time', required=False, dest='sim_step_time', action='store', default=.01)
parser.add_argument('-v', type=bool, help = 'Show video', required=False, default=False, dest='show_video',action='store')

args = parser.parse_args()

def main():
    sim = RoverSimulator()
    s = 'data/wrappers/simple_model'
    sys.path.insert(0,os.path.join(os.getcwd(),*s.split('/')))
    sim.loadLambdaFunctions(model_name='simple_model', wrapper_folder='data/wrappers/simple_model')
    sim.loadGroundPoints()
    sim.loadFrictionModel()
    gen_coord = len(sim.gen_coord)*[0]
    sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
    current_wheel_height = sim.getLowestWheelPoint()[1]
    gen_coord[2] = args.drop_height - current_wheel_height
    print('The rover centre of mass current height is: \t {:10.6f} [m]'.format(gen_coord[2]))
    input('Press Enter to start the simulation')
    sim.setInitialConditions(gen_coord, len(gen_coord)*[0])
    sim.simulate(args.sim_step_time, args.sim_time)
    
    if args.show_video:
        sim.animateMotion(override_save='yes')
    

if __name__ == '__main__':
    #print(sys.path)
    main()

