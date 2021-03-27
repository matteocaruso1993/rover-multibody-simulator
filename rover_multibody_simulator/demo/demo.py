# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 11:08:26 2021

@author: Matteo
"""

from rover_multibody_simulator import four_ws_rover_dynamic_simulator as sim


if __name__ == '__main__':
    simulator = sim.RoverSimulator()
    simulator.initialize()
    simulator.lambdifyAll()
    simulator.formEquationsOfMotion('autowrap - f2py')
    gen_coord = 15*[0]
    gen_coord[:3] = 3*[0]
    gen_coord[2] = 1
    gen_coord[0] = 0
    simulator.createGound(0.0,-0.2,0)
    simulator.setInitialConditions(gen_coord,15*[0])
    simulator.simulate(0.02,1)
    simulator.animateMotion()