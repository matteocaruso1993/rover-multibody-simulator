# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 16:09:06 2021

@author: Matteo
"""

from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator

if __name__ == '__main__':
    sim = RoverSimulator()
    sim.config.set('Model Description', 'debugging', 'true')
    sim.initialize()
    sim.lambdifyAll()
    sim.formEquationsOfMotion('autowrap - f2py')
    gen_coord = len(sim.gen_coord)*[0]
    #gen_coord[:3] = 3*[0]
    #gen_coord[2] = 1
    #gen_coord[0] = 0
    sim.setInitialConditions(gen_coord,len(sim.gen_coord)*[0])
    sim.loadGroundPoints()
    sim.simulate(.1, 2)
    sim.animateMotion()