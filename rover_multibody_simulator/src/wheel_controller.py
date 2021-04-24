<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:31:09 2021

@author: Matteo
"""
import numpy as np
from ..utilities.functions import clip

class RoverWheelControl:
    def __init__(self, n_wheels = 4, controller_frequency = 10):
        self._n_wheels = n_wheels
        self._steer_gain = None
        self._steer_d_gain = None
        self._steer_i_gain = None
        
        
        self._drive_gain = None
        self._drive_d_gain = None
        self._drive_i_gain = None
        
        self._steer_target = None #Steer angle
        self._drive_target = None #Wheel velocity
        self.__wheel_index = None #Wheels index
        
        
        
        #Add controller working frequency
        self._sampling_period = 1/controller_frequency
        self._current_control_time = 0
        self._last_control_action = None
        self._last_error = None
        self._integrator_term = None
        self._initial_step = True
        
        self._saturation_block = None
        
        self.__initialize()
        
    
    def __initialize(self):
        
        self._steer_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._steer_d_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_d_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._steer_i_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_i_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._last_control_action = np.zeros((2*self._n_wheels))
        self._last_error = np.zeros((2*self._n_wheels))
        self._integrator_term = np.zeros((2*self._n_wheels))
        self._current_control_time = 0
        
        self._saturation_block = dict()
        self._saturation_block['steer'] = [-float('inf'), float('inf')]
        self._saturation_block['drive'] = [-float('inf'), float('inf')]
        
        
        
    def setSteerGains(self, in_array, gain_type):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            if gain_type == 'P':
                self._steer_gain[i,i] = in_array[i]
            elif gain_type == 'D':
                self._steer_d_gain[i,i] = in_array[i]
            elif gain_type == 'I':
                self._steer_i_gain[i,i] = in_array[i]
            else:
                raise ValueError('Unexpected input argument for gain_type. Must be either: P, D, I')
            
    def getSteerGains(self):
        return {'P': self._steer_gain, 'D': self._steer_d_gain, 'I': self._steer_i_gain}
            
        
    def setDriveGains(self, in_array, gain_type):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            
            if gain_type == 'P':
                self._drive_gain[i,i] = in_array[i]
            
            elif gain_type == 'D':
                self._drive_d_gain[i,i] = in_array[i]
            elif gain_type == 'I':
                self._drive_i_gain[i,i] = in_array[i]
            else:
                raise ValueError('Unexpected input argument for gain_type. Must be either: P, D, I')
                
    def setSaturation(self, min_v, max_v, type_v):
        
        if min_v > max_v:
            raise ValueError('Minimum value must be less than the maximum one')
        
        if type_v != 'steer' and type_v != 'drive':
            raise ValueError('Input type argument must be either: steer or drive')
            
            
        self._saturation_block[type_v] = [min_v, max_v]
        
    def getSaturation(self):
        return self._saturation
        
            
            
            
    def getDriveGains(self):
        return {'P': self._drive_gain, 'D': self._drive_d_gain, 'I': self._drive_i_gain}
        
        
    def setSteerTarget(self, in_array):
        self._steer_target = in_array
        
    def setDriveTarget(self, in_array):
        self._drive_target = in_array
        
    def setWheelIndexes(self, steer_idx, drive_idx):
        self.__wheel_index = dict()
        self.__wheel_index['steer'] = steer_idx
        self.__wheel_index['drive'] = drive_idx
        
    def getWheelIndexes(self):
        return self.__wheel_index
    
    
    
    def reset(self):
        self.__initialize()
        
        
        
    def computeTorqueCorrection(self, in_vect):
        correction_steer = np.zeros((self._n_wheels,))
        correction_drive = np.zeros((self._n_wheels,))
        err_der_steer = np.zeros((self._n_wheels,))
        err_der_drive = np.zeros((self._n_wheels,))
        
        in_steer = in_vect[:self._n_wheels]
        in_drive = in_vect[self._n_wheels:]
        
        if self._steer_target is not None:
            #Compute correction
            err_steer = in_steer - self._steer_target
            
            #Proportional
            correction_steer = -np.dot(self._steer_gain, err_steer)
            
            #Derivative
            
            err_der_steer = (-self._last_error[:self._n_wheels] + err_steer)/self._sampling_period
            
            correction_steer -= np.dot(self._steer_d_gain, err_der_steer)
            
            #Integral term
            self._integrator_term[:self._n_wheels] += err_steer*self._sampling_period
            correction_steer -= np.dot(self._steer_i_gain, self._integrator_term[:self._n_wheels])
            self._last_error[:self._n_wheels] = err_steer
            
            for i, val in enumerate(correction_steer):
                correction_steer[i] = clip(val)
                
            
            
        if self._drive_target is not None:
            #Compute correction
            err_drive = in_drive - self._drive_target
            
            
            #Proportional
            correction_drive = -np.dot(self._drive_gain, err_drive)
            
            #Derivative
            err_der_drive = (-self._last_error[self._n_wheels:] + err_drive)/self._sampling_period
            correction_drive -= np.dot(self._drive_d_gain, err_der_drive)
            
            #Integral term
            self._integrator_term[self._n_wheels:] += err_drive*self._sampling_period
            correction_drive -= np.dot(self._drive_i_gain, self._integrator_term[self._n_wheels:])
            self._last_error[self._n_wheels:] = err_drive
            
            for i, val in enumerate(correction_drive):
                correction_drive[i] = clip(val)
        
        out = np.hstack((correction_steer, correction_drive))
        
        
        
        #Update controller time
        self._current_control_time += self._sampling_period
        
        #Update control action
        self._last_control_action = out
        
        
        
        return out
    
    
        
    
    
    

if __name__ == '__main__':
    control = RoverWheelControl(controller_frequency=5)
    wheel_current_speed = np.zeros((4,))
    wheel_current_steer = np.zeros((4,))
    steer_gains = [1,1,1,1]
    drive_gains = [1,1,1,1]
    
    control.setSteerGains(steer_gains,'P')
    control.setSteerGains(steer_gains,'D')
    control.setSteerGains(steer_gains,'I')
    
    control.setDriveGains(drive_gains, 'P')
    control.setDriveGains(drive_gains, 'D')
    control.setDriveGains(drive_gains, 'I')
    
    control.setSteerTarget(4*[np.pi/2])
    control.setDriveTarget(4*[2])
    
    
    sim_time = 10
    time_step = .001
    
    time_array = np.arange(0,sim_time, time_step)
    history = list()
    corr_history = list()
    curr_state = np.hstack((wheel_current_steer, wheel_current_speed))
    history.append(curr_state.tolist())
    corr = np.zeros((8,))
    for t in time_array:
        if control._current_control_time + control._sampling_period >= t:
            corr = control.computeTorqueCorrection(curr_state)
        #print(curr_state)
        curr_state += corr*time_step
        
        history.append(curr_state.tolist())
    
=======
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 13:31:09 2021

@author: Matteo
"""
import numpy as np
from ..utilities.functions import clip

class RoverWheelControl:
    def __init__(self, n_wheels = 4, controller_frequency = 10):
        self._n_wheels = n_wheels
        self._steer_gain = None
        self._steer_d_gain = None
        self._steer_i_gain = None
        
        
        self._drive_gain = None
        self._drive_d_gain = None
        self._drive_i_gain = None
        
        self._steer_target = None #Steer angle
        self._drive_target = None #Wheel velocity
        self.__wheel_index = None #Wheels index
        
        
        
        #Add controller working frequency
        self._sampling_period = 1/controller_frequency
        self._current_control_time = 0
        self._last_control_action = None
        self._last_error = None
        self._integrator_term = None
        self._initial_step = True
        
        self._saturation_block = None
        
        self.__initialize()
        
    
    def __initialize(self):
        
        self._steer_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._steer_d_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_d_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._steer_i_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        self._drive_i_gain = np.zeros((self._n_wheels, self._n_wheels), dtype='float64')
        
        self._last_control_action = np.zeros((2*self._n_wheels))
        self._last_error = np.zeros((2*self._n_wheels))
        self._integrator_term = np.zeros((2*self._n_wheels))
        self._current_control_time = 0
        
        self._saturation_block = dict()
        self._saturation_block['steer'] = [-float('inf'), float('inf')]
        self._saturation_block['drive'] = [-float('inf'), float('inf')]
        
        
        
    def setSteerGains(self, in_array, gain_type):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            if gain_type == 'P':
                self._steer_gain[i,i] = in_array[i]
            elif gain_type == 'D':
                self._steer_d_gain[i,i] = in_array[i]
            elif gain_type == 'I':
                self._steer_i_gain[i,i] = in_array[i]
            else:
                raise ValueError('Unexpected input argument for gain_type. Must be either: P, D, I')
            
    def getSteerGains(self):
        return {'P': self._steer_gain, 'D': self._steer_d_gain, 'I': self._steer_i_gain}
            
        
    def setDriveGains(self, in_array, gain_type):
        if not isinstance(in_array, list) and not isinstance(in_array, np.ndarray):
            raise ValueError('Unexpected input argument')
            
        n = len(in_array)
        if n != self._n_wheels:
            raise ValueError('Input array must have size equal to the numbers of wheels')
            
        for i in range(n):
            
            if gain_type == 'P':
                self._drive_gain[i,i] = in_array[i]
            
            elif gain_type == 'D':
                self._drive_d_gain[i,i] = in_array[i]
            elif gain_type == 'I':
                self._drive_i_gain[i,i] = in_array[i]
            else:
                raise ValueError('Unexpected input argument for gain_type. Must be either: P, D, I')
                
    def setSaturation(self, min_v, max_v, type_v):
        
        if min_v > max_v:
            raise ValueError('Minimum value must be less than the maximum one')
        
        if type_v != 'steer' and type_v != 'drive':
            raise ValueError('Input type argument must be either: steer or drive')
            
            
        self._saturation_block[type_v] = [min_v, max_v]
        
    def getSaturation(self):
        return self._saturation
        
            
            
            
    def getDriveGains(self):
        return {'P': self._drive_gain, 'D': self._drive_d_gain, 'I': self._drive_i_gain}
        
        
    def setSteerTarget(self, in_array):
        self._steer_target = in_array
        
    def setDriveTarget(self, in_array):
        self._drive_target = in_array
        
    def setWheelIndexes(self, steer_idx, drive_idx):
        self.__wheel_index = dict()
        self.__wheel_index['steer'] = steer_idx
        self.__wheel_index['drive'] = drive_idx
        
    def getWheelIndexes(self):
        return self.__wheel_index
    
    
    
    def reset(self):
        self.__initialize()
        
        
        
    def computeTorqueCorrection(self, in_vect):
        correction_steer = np.zeros((self._n_wheels,))
        correction_drive = np.zeros((self._n_wheels,))
        err_der_steer = np.zeros((self._n_wheels,))
        err_der_drive = np.zeros((self._n_wheels,))
        
        in_steer = in_vect[:self._n_wheels]
        in_drive = in_vect[self._n_wheels:]
        
        if self._steer_target is not None:
            #Compute correction
            err_steer = in_steer - self._steer_target
            
            #Proportional
            correction_steer = -np.dot(self._steer_gain, err_steer)
            
            #Derivative
            
            err_der_steer = (-self._last_error[:self._n_wheels] + err_steer)/self._sampling_period
            
            correction_steer -= np.dot(self._steer_d_gain, err_der_steer)
            
            #Integral term
            self._integrator_term[:self._n_wheels] += err_steer*self._sampling_period
            correction_steer -= np.dot(self._steer_i_gain, self._integrator_term[:self._n_wheels])
            self._last_error[:self._n_wheels] = err_steer
            
            for i, val in enumerate(correction_steer):
                correction_steer[i] = clip(val)
                
            
            
        if self._drive_target is not None:
            #Compute correction
            err_drive = in_drive - self._drive_target
            
            
            #Proportional
            correction_drive = -np.dot(self._drive_gain, err_drive)
            
            #Derivative
            err_der_drive = (-self._last_error[self._n_wheels:] + err_drive)/self._sampling_period
            correction_drive -= np.dot(self._drive_d_gain, err_der_drive)
            
            #Integral term
            self._integrator_term[self._n_wheels:] += err_drive*self._sampling_period
            correction_drive -= np.dot(self._drive_i_gain, self._integrator_term[self._n_wheels:])
            self._last_error[self._n_wheels:] = err_drive
            
            for i, val in enumerate(correction_drive):
                correction_drive[i] = clip(val)
        
        out = np.hstack((correction_steer, correction_drive))
        
        
        
        #Update controller time
        self._current_control_time += self._sampling_period
        
        #Update control action
        self._last_control_action = out
        
        
        
        return out
    
    
        
    
    
    

if __name__ == '__main__':
    control = RoverWheelControl(controller_frequency=5)
    wheel_current_speed = np.zeros((4,))
    wheel_current_steer = np.zeros((4,))
    steer_gains = [1,1,1,1]
    drive_gains = [1,1,1,1]
    
    control.setSteerGains(steer_gains,'P')
    control.setSteerGains(steer_gains,'D')
    control.setSteerGains(steer_gains,'I')
    
    control.setDriveGains(drive_gains, 'P')
    control.setDriveGains(drive_gains, 'D')
    control.setDriveGains(drive_gains, 'I')
    
    control.setSteerTarget(4*[np.pi/2])
    control.setDriveTarget(4*[2])
    
    
    sim_time = 10
    time_step = .001
    
    time_array = np.arange(0,sim_time, time_step)
    history = list()
    corr_history = list()
    curr_state = np.hstack((wheel_current_steer, wheel_current_speed))
    history.append(curr_state.tolist())
    corr = np.zeros((8,))
    for t in time_array:
        if control._current_control_time + control._sampling_period >= t:
            corr = control.computeTorqueCorrection(curr_state)
        #print(curr_state)
        curr_state += corr*time_step
        
        history.append(curr_state.tolist())
    
>>>>>>> b6d4c17251816852f4d7e9bf65d38a15bd77234e
    