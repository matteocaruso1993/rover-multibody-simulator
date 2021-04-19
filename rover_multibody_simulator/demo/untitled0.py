# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 22:34:11 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator

if __name__ == '__main__':
    sim = RoverSimulator()
    sim.initialize()
    sim.lambdifyAll()
    sim.formEquationsOfMotion('autowrap - f2py')
    gen_coord = len(sim.gen_coord)*[0]
    gen_coord[2] = 0.7027821962519036
    sim.setInitialConditions(gen_coord,len(sim.gen_coord)*[0])
    sim.loadGroundPoints()
    sim.simulate(.02, 1)