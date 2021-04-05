# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:53:43 2021

@author: Matteo
"""

import numpy as np
import matplotlib.pyplot as plt


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
    f.setCoeffiecients(*p)
    prop = f.getFrictionProperties()
        