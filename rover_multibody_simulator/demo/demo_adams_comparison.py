<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 22:34:11 2021

@author: Matteo
"""

import numpy as np
from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator

if __name__ == '__main__':
    sim = RoverSimulator()
    sim.initialize()
    sim.lambdifyAll('autowrap - cython')
    #sim.formEquationsOfMotion('autowrap - cython')
    #sim.loadLambdaFunctions(model_name='full-no-parametric')
    #gen_coord = len(sim.gen_coord)*[0]
    #gen_coord[2] = 0.7027821962519036# -0.4
    #gen_coord[2] = 0.7027821962519036-0.39037973397488457
    #sim.wheel_controller.setSteerTarget(4*[np.pi/2])
    #sim.setInitialConditions(gen_coord,len(sim.gen_coord)*[0])
    #sim.loadGroundPoints()
    #sim.loadFrictionModel()
    #sim.simulate(.02, 1)
    #sim.animateMotion(override_save='no-save')
=======
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 22:34:11 2021

@author: Matteo
"""

import numpy as np
from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator

if __name__ == '__main__':
    sim = RoverSimulator()
    sim.initialize()
    sim.lambdifyAll()
    sim.formEquationsOfMotion('autowrap - f2py')
    gen_coord = len(sim.gen_coord)*[0]
    #gen_coord[2] = 0.7027821962519036 -0.4
    sim.wheel_controller.setSteerTarget(4*[np.pi/2])
    sim.setInitialConditions(gen_coord,len(sim.gen_coord)*[0])
    sim.loadGroundPoints()
    sim.loadFrictionModel()
    #sim.simulate(.02, 1)
>>>>>>> b6d4c17251816852f4d7e9bf65d38a15bd77234e
