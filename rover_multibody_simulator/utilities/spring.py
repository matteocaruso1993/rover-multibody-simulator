# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:31:38 2021

@author: Matteo
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('font')


class Spring:
    def __init__ (self, k, T0, r):
        self.k = k
        self.T0 = T0
        self.r = r
        self._approximation_factors = dict()
        
        
    def setApproximationFactors(self, E, f):
        self._approximation_factors['E'] = E
        self._approximation_factors['f'] = f
        self._approximation_factors['scaling_factor'] = self.r 
        
    
    def computeResponse(self, in_angle):
        M = -(2*self.T0*self.r/np.pi)*np.arctan(self._approximation_factors['E']*self._approximation_factors['f']*in_angle) - self.r**2*self.k*in_angle
        return M
        
    
    def showResponse(self, domain = [-np.pi/4, np.pi/4], step = np.deg2rad(.0001)):
        angle_domain = np.hstack((np.arange(domain[0], domain[1], step), domain[1]))
        
        M = self.computeResponse(angle_domain)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'$\beta [rad]$', fontsize=14)
        ax.tick_params(labelsize=14)
        ax.set_ylabel('M [Nm]', fontsize=14)
        ax.set_xlim(-0.1,0.1)
        ax.set_ylim(-1.5,1.5)
        ax.plot(angle_domain, -np.sign(angle_domain)*self.T0*self.r - self.r**2*self.k*angle_domain, label='Theoric')
        
        
        test = [100,500,1000,10000]
        for val in test:
            self.setApproximationFactors(np.sqrt(val), np.sqrt(val))
            M = self.computeResponse(angle_domain)
            ax.plot(angle_domain, M, '--',label='Approximation ($f$ = %d)'%val)
        
        #ax.plot(angle_domain, M, label='Approximation (E = %d)'%val)
        ax.legend()
        fig.tight_layout()
        fig.savefig('spring_response.png', dpi=600)
        
        
        


if __name__ == '__main__':
    s = Spring(17800, 40, 0.025)
    s.setApproximationFactors(100, 100)
    s.showResponse()
        
        
    