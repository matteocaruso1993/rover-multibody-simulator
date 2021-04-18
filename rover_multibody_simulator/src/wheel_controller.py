# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:31:09 2021

@author: Matteo
"""
import numpy as np

class RoverWheelControl:
    def __init__(self, n_wheels = 4):
        self._n_wheels = n_wheels
        self._steer_gain = None
        self._drive_gain = None
        
        self._steer_target = None #Steer angle
        self._drive_target = None #Wheel velocity
        
        self.__initialize()
        
    
    def __initialize(self):
        self._steer_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        
    def setSteerGains(self, in_array):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            self._steer_gain[i,i] = in_array[i]
            
    def getSteerGains(self):
        return self._steer_gain
            
        
    def setDriveGains(self, in_array):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            self._drive_gain[i,i] = in_array[i]
            
            
    def getDriveGains(self):
        return self._drive_gain
        
        
    def setSteerTarget(self, in_array):
        self._steer_target = in_array
        
    def setDriveTarget(self, in_array):
        self._drive_target = in_array
        
    
    
    def reset(self):
        self.__initialize()
        
        
    def computeTorqueCorrection(self, in_vect):
        correction_steer = np.zeros((self._n_wheels,1))
        correction_drive = np.zeros((self._n_wheels,1))
        
        in_steer = in_vect[:self._n_wheels]
        in_drive = in_vect[self._n_wheels]
        
        if self._steer_target is not None:
            #Compute correction
            err_steer = in_steer - self._steer_target
            correction_steer = -np.dot(self._steer_gain, err_steer)
            
        if self._drive_target is not None:
            #Compute correction
            err_drive = in_drive - self._drive_target
            correction_drive = -np.dot(self._drive_gain, err_drive)
            
        
        out = np.hstack((correction_steer.T, correction_drive.T))
        return out
    
    
    

if __name__ == '__main__':
    control = RoverWheelControl()
    wheel_current_speed = np.zeros((4,))
    wheel_current_steer = np.zeros((4,))
    steer_gains = [1,1,1,1]
    drive_gains = [1,1,1,1]
    
    control.setSteerGains(steer_gains)
    control.setDriveGains(drive_gains)
    
    control.setSteerTarget(4*[np.pi/2])
    control.setDriveTarget(4*[2])
    
    
    sim_time = 10
    time_step = .001
    
    time_array = np.arange(0,sim_time, time_step)
    history = list()
    corr_history = list()
    curr_state = np.hstack((wheel_current_steer, wheel_current_speed))
    history.append(curr_state.tolist())
    for t in time_array:
        corr = control.computeTorqueCorrection(curr_state)
        #print(curr_state)
        curr_state += corr*time_step
        
        history.append(curr_state.tolist())
    
    