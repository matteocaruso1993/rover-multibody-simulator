# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 11:08:26 2021

@author: Matteo
"""

from rover_multibody_simulator import four_ws_rover_dynamic_simulator as sim
import json


if __name__ == '__main__':
    simulator = sim.RoverSimulator()
    simulator.initialize()
    simulator.lambdifyAll()
    simulator.formEquationsOfMotion('autowrap - f2py')
    gen_coord = len(simulator.gen_coord)*[0]
    #gen_coord[:3] = 3*[0]
    #gen_coord[2] = 1
    #gen_coord[0] = 0
    
    simulator.loadGroundPoints()
    simulator.changeTerrain('sine_terrain')
    
    with open('initial_conditions_flat_surface.json','r') as f:
        pos = json.load(f)
    
    gen_coord = pos['generalized coordinates']
    
    simulator.setInitialConditions(gen_coord,len(simulator.gen_coord)*[1e-8])
    simulator.loadFrictionModel()
    simulator.setDrivingTorque(*4*[1])
    simulator.simulate(0.02,2)
    # in_coord = simulator.getCurrentState()
    # gen_coord1 = in_coord[:,0]
    # simulator.setDrivingTorque(*4*[5])
    # simulator.setInitialConditions(gen_coord1, len(gen_coord1)*[1e-6])
    # simulator.reset()
    # simulator.simulate(0.02,3)
    #simulator.animateMotion()

