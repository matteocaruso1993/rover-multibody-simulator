# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
from rover_multibody_simulator.four_ws_rover_dynamic_simulator_full_par import RoverSimulator

def main():
    sim = RoverSimulator()
    sim.initialize()
    sim.lambdifyAll('autowrap - cython')
    sim.formEquationsOfMotion('autowrap - f2py')
    sim.saveLambdaFunctions('full-parametric-model')
    



if __name__ == '__main__':
    if os.getcwd() != os.path.abspath(os.path.dirname(__file__)):
        os.chdir(os.path.abspath(os.path.dirname(__file__)))
    main()