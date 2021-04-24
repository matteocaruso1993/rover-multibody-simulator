# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:53:43 2021

@author: Matteo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, dual_annealing

from geneticalgorithm import geneticalgorithm as ga


class Friction:
    def __init__(self):
        self._coefficients = dict()
        self._friction_properties = dict()
        
        self.reset()
        
        
    def reset(self):
        self._coefficients['g1'] = 0
        self._coefficients['g2'] = 0
        self._coefficients['g3'] = 0
        self._coefficients['g4'] = 0
        self._coefficients['g5'] = 0
        self._coefficients['g6'] = 0
        
    def setCoeffiecients(self, g1, g2, g3, g4, g5, g6):
        self._coefficients['g1'] = g1
        self._coefficients['g2'] = g2
        self._coefficients['g3'] = g3
        self._coefficients['g4'] = g4
        self._coefficients['g5'] = g5
        self._coefficients['g6'] = g6
        
        self._computeFrictionProperties()
        
        
        
    def getCoefficients(self, in_type = 'dict'):
        if in_type == 'dict':
            return self._coefficients
        elif in_type == 'list':
            return list(self._coefficients.values())
        else:
            raise ValueError('Invalid input argument')
        
    def computeFrictionCoeffiecient(self, in_vel):
        #mu = g1*(tanh(g2*in_vel) - tanh(g3*in_vel)) + g4*tanh(g5*in_vel) + g6*in_vel;
        c = list(self._coefficients.values())
        
        mu = c[0]*(np.tanh(c[1]*in_vel) - np.tanh(c[2]*in_vel)) + c[3]*np.tanh(c[4]*in_vel) + c[5]*in_vel
        return mu
        
    def showFrictionCoefficient(self, vel_domain = [-0.2, 0.2], samples = 100):
        v = np.linspace(vel_domain[0], vel_domain[1], samples)
        mu = self.computeFrictionCoeffiecient(v)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(v,mu)
        ax.set_xlabel('v [m/s]')
        ax.set_ylabel('\mu ')
        plt.show()
        
    def getVeltoFrictionMapping(self, vel_in):
        if not isinstance(vel_in, np.ndarray):
            raise ValueError('Invalid input argument! Expected a np array')
        mu = self.computeFrictionCoeffiecient(vel_in)
        return np.array((vel_in, mu)).T
    
    
    def getFrictionProperties(self):
        return self._friction_properties
    
    def setFrictionProperties(self, mu_s, mu_d, v_s, v_d, show_result = False):
        
        self._friction_properties['mu static'] = mu_s
        self._friction_properties['v static'] = v_s
        
        self._friction_properties['mu dynamic'] = mu_d
        self._friction_properties['v dynamic'] = v_d
        
        bounds = [(0,1e3),(0,1e3), (0,1e3), (0,1e3)]
        
        vel_test = np.linspace(0,v_d,1000)
        
        def optimFunction(x):
            fun = np.zeros((4,), dtype=float)
            
            fun[0] = x[0]*(np.tanh(x[1]*v_s)-np.tanh(x[2]*v_s))+ mu_d*np.tanh(x[3]*v_s)-mu_s
            fun[1] = x[0]*(np.tanh(x[1]*v_d)-np.tanh(x[2]*v_d))+ mu_d*np.tanh(x[3]*v_d)-mu_d
            vn = v_d + 0.01
            fun[2] = x[0]*(np.tanh(x[1]*vn)-np.tanh(x[2]*vn))+ mu_d*np.tanh(x[3]*vn)-mu_d
            vn += 0.01
            fun[3] = x[0]*(np.tanh(x[1]*vn)-np.tanh(x[2]*vn))+ mu_d*np.tanh(x[3]*vn)-mu_d
            
            err = np.sum(fun**2)
            
            if x[1] < x[2]:
                err += 100000
    
            if np.any(x < 0):
                err += 100000
            
            g = [x[0], x[1],x[2], mu_d, x[3], 0]
            
            self.setCoeffiecients(*g)
            mu = self.computeFrictionCoeffiecient(vel_test)
            pos_mu = mu[mu> mu_s]
            err += np.sum(100*pos_mu)
                
            return err
        
        
        #model = differential_evolution(optimFunction, bounds, tol=1e-6, maxiter=10000, popsize=100)
        model = dual_annealing(optimFunction, bounds, maxiter=10000, x0 = np.ones(4,) )
        algorithm_param = {'max_num_iteration': 1000,\
                   'population_size':100,\
                   'mutation_probability':0.1,\
                   'elit_ratio': 0.01,\
                   'crossover_probability': 0.5,\
                   'parents_portion': 0.3,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv':None}
        
        #model = ga(function = optimFunction, dimension = 4,variable_type = 'real', variable_boundaries=np.array(bounds), algorithm_parameters = algorithm_param)
        #model.run()
        
        g = [model.x[0], model.x[1],model.x[2], mu_d, model.x[3], 0]
        self.setCoeffiecients(*g)
        
        if show_result:
            self.showFrictionCoefficient()
        
        return model
                
        
        
        
        
        
    def _computeFrictionProperties(self):
        v = np.linspace(0,0.2,1000)
        mapp = self.getVeltoFrictionMapping(v)
        
        mu_s_idx = np.argmax(mapp[:,1])
        
        self._friction_properties['mu static'] = mapp[mu_s_idx,1]
        self._friction_properties['v static'] = mapp [mu_s_idx,0]
        
        new_arr = mapp[mu_s_idx:,:]
        
        mu_d_idx = np.argmin(new_arr[:,1])
        
        self._friction_properties['mu dynamic'] = np.mean(new_arr[mu_d_idx:,1])
        self._friction_properties['v dynamic'] = new_arr[mu_d_idx,0]
        
        
        

if __name__ == '__main__':
    p = [2.0430, 239.1462, 120.0534, 0.5000, 84.4111, 0]
    f = Friction()
    r = f.setFrictionProperties(0.8, 0.5, 0.005, 0.02, show_result=True)
    #f.setCoeffiecients(*p)
    #prop = f.getFrictionProperties()
