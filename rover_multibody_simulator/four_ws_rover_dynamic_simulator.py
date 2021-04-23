# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 22:45:22 2021

@author: Matteo Matteo Caruso
@email: matteo.caruso@phd.units.it

"""



from __future__ import print_function, division
from sympy import symbols, simplify, lambdify, atan, diff, Dummy, Matrix
from sympy.utilities.codegen import codegen
from sympy.utilities.autowrap import autowrap
from sympy.physics.mechanics import dynamicsymbols, find_dynamicsymbols, ReferenceFrame, Point, inertia, RigidBody, Lagrangian, LagrangesMethod, mprint, KanesMethod, msubs
from sympy.physics.vector import init_vprinting
from sympy.printing.theanocode import theano_function
from scipy.integrate import odeint, solve_ivp
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, tri
import datetime
import numpy as np
import time
from tqdm import tqdm
import configparser
import os
import dill
import json

from .utilities.functions import step5, generateWheelsPoints, findContinguousContactRegions
from .src.friction import Friction





init_vprinting(use_latex=None,pretty_print=False)


class RoverSimulator:
    def __init__(self,print_latex=False, config_path = None):
        
        
        self.frames = list()
        self.lambdified_frames = list()
        self.lambdified_frames_ang_vel = list()
        self.lambdified_frames_ang_acc = list()
        
        
        self.points = list()
        
        self.lambdified_points = list()
        self.lambdified_points_vel = list()
        self.lambdified_points_acc = list()
        
        self.lambdified_partial_velocities = dict()
        
        self.contact_parametric = dict()
        self.friction = None
              
        
        self.bodies = list()
        self.torques = list()
        self.__driving_torques = list()
        self.__current_driving_torque = None
        
        
        if config_path is None:
            self.config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','config','config.ini')
            #print(self.config_path)
        else:
            self.config_path = config_path
        self.method = 'Lagrange'
        
        
        self.gen_coord = list()
        self.gen_speeds = list()
        
        self.gravity_forces = list()
        
        self.q_initial_condition = None
        self.q_d_initial_condition = None
        
        self.current_gen_coord = None
        self.current_gen_speed = None
        
        self.lambdified = False
        
        self.eq_const = list() #For Kane Only
        self.__config_found = False
        self.config = None
        
        self.kane_method = dict()
        self.kane_method['lambda func'] = dict()
        
        self.lagrange_method = dict()
        self.lagrange_method['lambda func'] = dict()
        
        self.ode_sol = None
        
        
        self.__ground = dict()
        self.__ground['coeffs'] = dict()
        
        
        self.__fig = None
        self.__ax = None
        
        
        self.__gen_coord_history = list(list())
        self.__gen_speed_history = list(list())
        self.__time_history = list()
        
        self.__reaction_torque = list()
        
        self.__map = dict()
        
        self.__wheel_points = None
        
        self.ground_new = None
        self.__ground_interpolator = None
        
        
        
        
        
        if print_latex:
            init_vprinting(use_latex='mathjax',pretty_print=False)
        else:
            init_vprinting(use_latex=None,pretty_print=False)
        
        
        self.loadConfig()
        
    
    def setMethod(self,new_method):
        if not new_method == 'Lagrange' and not new_method == 'Kane':
            raise ValueError('Invalid Input Argument')
        
        self.method = new_method
                
    
            
    
    def loadConfig(self):
        
        config = configparser.ConfigParser()
        if not config.read(self.config_path):
            self.__config_found = False
            print('Couldn''t locate config file at the specified path! Can not proceed with pre substitution')
        else:
            self.__config_found = True
            self.config = config
            print('Config file successfully loaded')
        

    def initialize(self):
        
        """ This method initialize the class roverSimulator trying to use constants given
        the configuration file. If configuration can not be found constants will be treated as symbols.
        
        
        """
        print('Initializing')
        #Start with defining the frames needed for the rover
        inertial_frame = ReferenceFrame('I') #Inertial reference frame
        rover_pitch_frame = ReferenceFrame('R_{pitch}')
        rover_intermediate_frame = ReferenceFrame('R_{int}') #Intermediate reference frame
        rover_frame = ReferenceFrame('R') #Rover final frame
        
        
        #%% Constant Symbols definition
        rover_mass = symbols('m_R')
        d_h = symbols('d_h')
        off_swing_dst, l_hinge_1, h_hinge_1 = symbols('d_{SA_{CM}}, l_{h_1},h_{h_1}')
        mass_swing_arm = symbols('m_{SA}')
        swing_arm_i_xx, swing_arm_i_yy, swing_arm_i_zz, swing_arm_i_xy, swing_arm_i_xz,swing_arm_i_yz = symbols('I_{SA_{xx}}, I_{SA_{yy}}, I_{SA_{zz}}, I_{SA_{xy}}, I_{SA_{xz}}, I_{SA_{yz}}')
        d_CM_link, l_link = symbols('d_{L_{CM}}, l_{L}')
        d_CM_link = d_CM_link.subs({d_CM_link: l_link/2}) #Assumed centre mass in middle of the link
        mass_link, link_inertia_xx,link_inertia_yy,link_inertia_zz = symbols('m_L, I_{L_{xx}}, I_{L_{yy}}, I_{L_{zz}}')
        wheel_offset = symbols('d_W')
        wheel_mass, wheel_inertia_xx,wheel_inertia_yy, wheel_inertia_zz = symbols('m_W, I_{W_{xx}}, I_{W_{yy}}, I_{W_{zz}}')
        rover_inertia_x,rover_inertia_y,rover_inertia_z = symbols('I_{R_{xx}}, I_{R_{yy}}, I_{R_{zz}}')
        wheel_radius = symbols('w_R')
        
        E_fr1, E_br1, E_fl1, E_bl1, E_fr2, E_br2, E_fl2, E_bl2 = symbols('E_{fr1}, E_{br1}, E_{fl1}, E_{bl1}, E_{fr2}, E_{br2}, E_{fl2}, E_{bl2}')
        k_fr1, k_br1, k_fl1, k_bl1, k_fr2, k_br2, k_fl2, k_bl2 = symbols('k_{fr1}, k_{br1}, k_{fl1}, k_{bl1}, k_{fr2}, k_{br2}, k_{fl2}, k_{bl2}')
        M0_fr1, M0_br1, M0_fl1, M0_bl1, M0_fr2, M0_br2, M0_fl2, M0_bl2 = symbols('M0_{fr1}, M0_{br1}, M0_{fl1}, M0_{bl1}, M0_{fr2}, M0_{br2}, M0_{fl2}, M0_{bl2}')
        f_fr1, f_br1, f_fl1, f_bl1, f_fr2, f_br2, f_fl2, f_bl2 = symbols('f_{fr1}, f_{br1}, f_{fl1}, f_{bl1}, f_{fr2}, f_{br2}, f_{fl2}, f_{bl2}')
        c_fr1, c_br1, c_fl1, c_bl1, c_fr2, c_br2, c_fl2, c_bl2 = symbols('c_{fr1}, c_{br1}, c_{fl1}, c_{bl1}, c_{fr2}, c_{br2}, c_{fl2}, c_{bl2}')
        
        def springTorqueExpression(sim_var, k, M, E, f):
            expression = -((M*2/np.pi)*atan(sim_var*E*f) + k*sim_var)
            #expression = -(k*sim_var)
            return expression
        
        
        
        
        
        if self.__config_found:
            
            #Load rover body properties
            
            rover_mass = self.config.getfloat('Model Description','rover_mass')
            rover_inertia_x = self.config.getfloat('Model Description','rover_inertia_xx')
            rover_inertia_y = self.config.getfloat('Model Description','rover_inertia_yy')
            rover_inertia_z = self.config.getfloat('Model Description','rover_inertia_zz')
            rover_inertia_xy =self.config.getfloat('Model Description', 'rover_inertia_xy')
            rover_inertia_xz =self.config.getfloat('Model Description', 'rover_inertia_xz')
            rover_inertia_yz =self.config.getfloat('Model Description', 'rover_inertia_yz')
            
            # Load bogie parameters
            d_h = self.config.getfloat('Model Description','arm_offset')
            off_swing_dst = self.config.getfloat('Model Description','arm_cm_distance')
            l_link = self.config.getfloat('Model Description','link_lenght')
            l_hinge_1 = self.config.getfloat('Model Description','hinge1_len')
            h_hinge_1 = self.config.getfloat('Model Description','hinge1_height')
            mass_swing_arm = self.config.getfloat('Model Description','swing_arm_mass')
            swing_arm_i_xx = self.config.getfloat('Model Description','swing_arm_inertia_xx')
            swing_arm_i_yy = self.config.getfloat('Model Description','swing_arm_inertia_yy')
            swing_arm_i_zz = self.config.getfloat('Model Description','swing_arm_inertia_zz')
            swing_arm_i_xy = self.config.getfloat('Model Description','swing_arm_inertia_xy')
            swing_arm_i_xz = self.config.getfloat('Model Description','swing_arm_inertia_xz')
            swing_arm_i_yz = self.config.getfloat('Model Description','swing_arm_inertia_yz')
            
            
            # Load central link parameters
            d_CM_link_x = self.config.getfloat('Model Description','link_cm_dist')
            mass_link = self.config.getfloat('Model Description','link_mass')
            link_inertia_xx = self.config.getfloat('Model Description','link_inertia_xx')
            link_inertia_yy = self.config.getfloat('Model Description','link_inertia_yy')
            link_inertia_zz = self.config.getfloat('Model Description','link_inertia_zz')
            wheel_offset = self.config.getfloat('Model Description','wheel_offset')
            wheel_mass = self.config.getfloat('Model Description','wheel_mass')
            wheel_inertia_xx = self.config.getfloat('Model Description','wheel_inertia_xx')
            wheel_inertia_yy = self.config.getfloat('Model Description','wheel_inertia_yy')
            wheel_inertia_zz = self.config.getfloat('Model Description','wheel_inertia_zz')
            wheel_radius = self.config.getfloat('Model Description','wheel_radius')
            
            E_fr1 = self.config.getfloat('Springs Definition','E_fr1')
            E_br1 = self.config.getfloat('Springs Definition','E_br1')
            E_fl1 = self.config.getfloat('Springs Definition','E_fl1')
            E_bl1 = self.config.getfloat('Springs Definition','E_bl1')
            
            E_fr2 = self.config.getfloat('Springs Definition','E_fr2')
            E_br2 = self.config.getfloat('Springs Definition','E_br2')
            E_fl2 = self.config.getfloat('Springs Definition','E_fl2')
            E_bl2 = self.config.getfloat('Springs Definition','E_bl2')
            
            k_fr1 = self.config.getfloat('Springs Definition','k_fr1')
            k_br1 = self.config.getfloat('Springs Definition','k_br1')
            k_fl1 = self.config.getfloat('Springs Definition','k_fl1')
            k_bl1 = self.config.getfloat('Springs Definition','k_bl1')
            
            k_fr2 = self.config.getfloat('Springs Definition','k_fr2')
            k_br2 = self.config.getfloat('Springs Definition','k_br2')
            k_fl2 = self.config.getfloat('Springs Definition','k_fl2')
            k_bl2 = self.config.getfloat('Springs Definition','k_bl2')
            
            f_fr1 = self.config.getfloat('Springs Definition','f_fr1')
            f_br1 = self.config.getfloat('Springs Definition','f_br1')
            f_fl1 = self.config.getfloat('Springs Definition','f_fl1')
            f_bl1 = self.config.getfloat('Springs Definition','f_bl1')
            
            f_fr2 = self.config.getfloat('Springs Definition','f_fr2')
            f_br2 = self.config.getfloat('Springs Definition','f_br2')
            f_fl2 = self.config.getfloat('Springs Definition','f_fl2')
            f_bl2 = self.config.getfloat('Springs Definition','f_bl2')
            
            M0_fr1 = self.config.getfloat('Springs Definition','M0_fr1')
            M0_br1 = self.config.getfloat('Springs Definition','M0_br1')
            M0_fl1 = self.config.getfloat('Springs Definition','M0_fl1')
            M0_bl1 = self.config.getfloat('Springs Definition','M0_bl1')
            
            M0_fr2 = self.config.getfloat('Springs Definition','M0_fr2')
            M0_br2 = self.config.getfloat('Springs Definition','M0_br2')
            M0_fl2 = self.config.getfloat('Springs Definition','M0_fl2')
            M0_bl2 = self.config.getfloat('Springs Definition','M0_bl2')
            
            c_fr1 = self.config.getfloat('Springs Definition','c_fr1')
            c_br1 = self.config.getfloat('Springs Definition','c_br1')
            c_fl1 = self.config.getfloat('Springs Definition','c_fl1')
            c_bl1 = self.config.getfloat('Springs Definition','c_bl1')
            
            c_fr2 = self.config.getfloat('Springs Definition','c_fr2')
            c_br2 = self.config.getfloat('Springs Definition','c_br2')
            c_fl2 = self.config.getfloat('Springs Definition','c_fl2')
            c_bl2 = self.config.getfloat('Springs Definition','c_bl2')
            
            

            self.setMethod(self.config['Simulator']['eom_method'])
            
#%% Dynamics symbols definition
        theta, psi = dynamicsymbols('theta, psi') #First two rotation of the rover
        chi = dynamicsymbols('chi')
        x,y,z = dynamicsymbols('x,y,z')
        phi = dynamicsymbols('phi')
        alpha_1 = dynamicsymbols('alpha1')
        alpha_2 = dynamicsymbols('alpha2')
        delta_1 = dynamicsymbols('delta1')
        omega_1 = dynamicsymbols('omega1')
        beta_1 = dynamicsymbols('beta1')
        beta_2 = dynamicsymbols('beta2')
        delta_2 = dynamicsymbols('delta2')
        omega_2 = dynamicsymbols('omega2')
        gamma_1 = dynamicsymbols('gamma1')
        gamma_2 = dynamicsymbols('gamma2')
        delta_3 = dynamicsymbols('delta3')
        omega_3 = dynamicsymbols('omega3')
        eps_1 = dynamicsymbols('epsilon1')
        eps_2 = dynamicsymbols('epsilon2')
        delta_4 = dynamicsymbols('delta4')
        omega_4 = dynamicsymbols('omega4')
        
        if self.method == 'Lagrange':
            
            x_d,y_d,z_d = dynamicsymbols('x,y,z',1)
            theta_d, psi_d = dynamicsymbols('theta, psi',1)
            chi_d = dynamicsymbols('chi',1)
            phi_d = dynamicsymbols('phi',1)
            alpha_1_d = dynamicsymbols('alpha1',1)
            alpha_2_d = dynamicsymbols('alpha2',1)
            delta_1_d = dynamicsymbols('delta1',1)
            omega_1_d = dynamicsymbols('omega1',1)
            beta_1_d = dynamicsymbols('beta1',1)
            beta_2_d = dynamicsymbols('beta2',1)
            delta_2_d = dynamicsymbols('delta2',1)
            omega_2_d = dynamicsymbols('omega2',1)
            gamma_1_d = dynamicsymbols('gamma1',1)
            gamma_2_d = dynamicsymbols('gamma2',1)
            delta_3_d = dynamicsymbols('delta3',1)
            omega_3_d = dynamicsymbols('omega3',1)
            eps_1_d = dynamicsymbols('epsilon1',1)
            eps_2_d = dynamicsymbols('epsilon2',1)
            delta_4_d = dynamicsymbols('delta4',1)
            omega_4_d = dynamicsymbols('omega4',1)
        else:

            x_d,y_d,z_d = dynamicsymbols('u1, u2, u3')
            theta_d, chi_d, psi_d = dynamicsymbols('u4, u5, u6')
            phi_d = dynamicsymbols('u7')
            alpha_1_d = dynamicsymbols('u8')
            alpha_2_d = dynamicsymbols('u9')
            delta_1_d = dynamicsymbols('u10')
            omega_1_d = dynamicsymbols('u11')
            beta_1_d = dynamicsymbols('u12')
            beta_2_d = dynamicsymbols('u13')
            delta_2_d = dynamicsymbols('u14')
            omega_2_d = dynamicsymbols('u15')
            gamma_1_d = dynamicsymbols('u16')
            gamma_2_d = dynamicsymbols('u17')
            delta_3_d = dynamicsymbols('u18')
            omega_3_d = dynamicsymbols('u19')
            eps_1_d = dynamicsymbols('u20')
            eps_2_d = dynamicsymbols('u21')
            delta_4_d = dynamicsymbols('u22')
            omega_4_d = dynamicsymbols('u23')
        
        self.gen_coord = [x,y,z,theta, chi,psi,phi,alpha_1,alpha_2, delta_1, omega_1,beta_1,beta_2, delta_2, omega_2,
                  gamma_1,gamma_2, delta_3, omega_3,eps_1,eps_2, delta_4, omega_4]   
        
        self.gen_speeds = [x_d,y_d,z_d,theta_d, chi_d,psi_d,phi_d,alpha_1_d,alpha_2_d, delta_1_d, omega_1_d,beta_1_d,beta_2_d, delta_2_d, omega_2_d,
                  gamma_1_d,gamma_2_d, delta_3_d, omega_3_d,eps_1_d,eps_2_d, delta_4_d, omega_4_d]
        
        if self.config.getboolean('Model Description','simplified'):
            
            delta_1 = delta_2 = delta_3 = delta_4 = omega_1 =omega_2 =omega_3 = omega_4 = omega_4 = 0
            delta_1_d = delta_2_d = delta_3_d = delta_4_d = omega_1_d =omega_2_d =omega_3_d = omega_4_d = omega_4_d = 0
            
            self.gen_coord = [x,y,z,theta,chi,psi,phi,alpha_1,alpha_2,beta_1,beta_2,
                  gamma_1,gamma_2,eps_1,eps_2]
            self.gen_speeds = [x_d,y_d,z_d,theta_d,chi_d, psi_d,phi_d,alpha_1_d,alpha_2_d, beta_1_d,beta_2_d,
                  gamma_1_d,gamma_2_d, eps_1_d,eps_2_d]
            
        if self.config.getboolean('Model Description','debugging'):
            z = 1
            z_d = 0
            delta_1 = delta_2 = delta_3 = delta_4 = omega_1 =omega_2 =omega_3 = omega_4 = omega_4 = 0
            delta_1_d = delta_2_d = delta_3_d = delta_4_d = omega_1_d =omega_2_d =omega_3_d = omega_4_d = omega_4_d = 0
            
            self.gen_coord = [x,y,theta,chi, psi,phi,alpha_1,alpha_2,beta_1,beta_2,
                  gamma_1,gamma_2,eps_1,eps_2]
            self.gen_speeds = [x_d,y_d,theta_d, chi_d,psi_d,phi_d,alpha_1_d,alpha_2_d, beta_1_d,beta_2_d,
                  gamma_1_d,gamma_2_d, eps_1_d,eps_2_d]
                  
        for i in range(len(self.gen_speeds)):
            self.gen_coord[i].__class__.__module__ = '__main__'
            self.gen_speeds[i].__class__.__module__ = '__main__'
            
            self.eq_const.append(diff(self.gen_coord[i]) - self.gen_speeds[i])
            


        
#%% Rover Body
        
        rover_intermediate_frame.orient(inertial_frame, 'Axis', (theta, inertial_frame.z)) #Yaw
        rover_pitch_frame.orient(rover_intermediate_frame, 'Axis', (chi, rover_intermediate_frame.x))
        rover_frame.orient(rover_pitch_frame, 'Axis',(psi,rover_pitch_frame.y))
        
        
        rover_intermediate_frame.set_ang_vel(inertial_frame,theta_d*inertial_frame.z)
        rover_pitch_frame.set_ang_vel(rover_intermediate_frame, chi_d*rover_intermediate_frame.x)
        rover_frame.set_ang_vel(rover_pitch_frame, psi_d*rover_pitch_frame.y)
        
        O = Point('O')
        O.set_vel(inertial_frame,0) #Assign null speed to origin
        O.set_acc(inertial_frame, 0)

        rover_centre_mass = Point('P_{R_{CM}}')
        rover_centre_mass.set_pos(O,x*inertial_frame.x + y*inertial_frame.y + z*inertial_frame.z)
        

        rover_centre_mass.set_vel(inertial_frame, x_d*inertial_frame.x + y_d*inertial_frame.y + z_d*inertial_frame.z)
        rover_centre_mass.set_acc(inertial_frame, diff(x_d)*inertial_frame.x + diff(y_d)*inertial_frame.y + diff(z_d)*inertial_frame.z)
        

        rover_inertia_dyadic = inertia(rover_frame,rover_inertia_x,rover_inertia_y,rover_inertia_z)

        rover_inertia = (rover_inertia_dyadic, rover_centre_mass)

        rover_body = RigidBody('Rover Body',rover_centre_mass, rover_frame, rover_mass, rover_inertia)


        g = symbols('g')
        g = float(g.subs({g:9.81}))


        rover_grav_force_vector = -rover_mass*g*inertial_frame.z
        self.gravity_forces.append((rover_centre_mass, rover_grav_force_vector))

        second_case_example = -rover_grav_force_vector.dot(rover_centre_mass.pos_from(O))
        rover_body.potential_energy = second_case_example

#%% Right Swing Arm
        right_swing_arm_hinge_pt = Point('P_{SA_r}')
        right_swing_arm_hinge_pt.set_pos(rover_centre_mass, d_h*rover_frame.x)

        #Set linear velocity of this point (Useful to use Poisson rule and check the result with previous one)
        right_swing_arm_hinge_pt.v2pt_theory(rover_centre_mass,inertial_frame,rover_frame)
        right_swing_arm_hinge_pt.a2pt_theory(rover_centre_mass,inertial_frame,rover_frame)
        


        right_swing_arm_frame = ReferenceFrame('R_{SA_{rf}}')
        right_swing_arm_frame.orient(rover_frame, 'Axis',(phi,rover_frame.x))

        
        right_swing_arm_frame.set_ang_vel(rover_frame,phi_d*rover_frame.x)

        right_swing_centre_mass = Point('SA_{r,CM}')
        right_swing_centre_mass.set_pos(right_swing_arm_hinge_pt, -off_swing_dst*right_swing_arm_frame.z)

        #Compute right swing arm centre mass speed by assigning manual speed and using Poisson rule
        right_swing_centre_mass.v2pt_theory(right_swing_arm_hinge_pt, inertial_frame, right_swing_arm_frame)
        right_swing_centre_mass.a2pt_theory(right_swing_arm_hinge_pt, inertial_frame, right_swing_arm_frame)

        #Define now the first hinge location for the front leg
        right_front_hinge1_pt = Point('P_{rfh1}')
        right_front_hinge1_pt.set_pos(right_swing_arm_hinge_pt, l_hinge_1*right_swing_arm_frame.y -h_hinge_1*right_swing_arm_frame.z)


        #Compute hinge liner speed using poisson rule
        right_front_hinge1_pt.v2pt_theory(right_swing_arm_hinge_pt,inertial_frame,right_swing_arm_frame)
        right_front_hinge1_pt.a2pt_theory(right_swing_arm_hinge_pt,inertial_frame,right_swing_arm_frame)

        #Define now the first hinge location for the back leg
        right_back_hinge1_pt = Point('P_{rbh1}')
        right_back_hinge1_pt.set_pos(right_swing_arm_hinge_pt, -l_hinge_1*right_swing_arm_frame.y -h_hinge_1*right_swing_arm_frame.z)
        
        right_back_hinge1_pt.v2pt_theory(right_swing_arm_hinge_pt,inertial_frame,right_swing_arm_frame)
        right_back_hinge1_pt.a2pt_theory(right_swing_arm_hinge_pt,inertial_frame,right_swing_arm_frame)
        
        
        swing_arm_inertia_dyadic = inertia(right_swing_arm_frame,swing_arm_i_xx,swing_arm_i_yy,swing_arm_i_zz,swing_arm_i_xy,                                  swing_arm_i_yz, swing_arm_i_xz)


        right_swing_arm_inertia = (swing_arm_inertia_dyadic, right_swing_centre_mass)


        right_swing_arm_body = RigidBody('Right Swing Arm',right_swing_centre_mass,                                 right_swing_arm_frame, mass_swing_arm, right_swing_arm_inertia)



        swing_arm_grav_force_vector = -mass_swing_arm*g*inertial_frame.z
        
        self.gravity_forces.append((right_swing_centre_mass, swing_arm_grav_force_vector))

        potential_field_right_swing_arm = -swing_arm_grav_force_vector.dot(right_swing_centre_mass.pos_from(O))
        right_swing_arm_body.potential_energy = potential_field_right_swing_arm


#%% First intermediate link front right

        right_front_first_int_frame = ReferenceFrame('R_{FR_{i_1}}')

        right_front_first_int_frame.orient(right_swing_arm_frame,'Axis',(-atan(h_hinge_1/l_hinge_1)+alpha_1, right_swing_arm_frame.x))

        right_front_first_int_frame.set_ang_vel(right_swing_arm_frame, alpha_1_d*right_swing_arm_frame.x)
        
        #Create points and assign locations
        first_int_front_right_centre_mass = Point('P_{rf_{I1_cm}}')
        first_int_front_right_centre_mass.set_pos(right_front_hinge1_pt, d_CM_link*right_front_first_int_frame.y)
        right_front_hinge2_pt = Point('P_{rfh2}')
        right_front_hinge2_pt.set_pos(right_front_hinge1_pt, l_link*right_front_first_int_frame.y)


        first_int_front_right_centre_mass.v2pt_theory(right_front_hinge1_pt, inertial_frame, right_front_first_int_frame)
        first_int_front_right_centre_mass.a2pt_theory(right_front_hinge1_pt, inertial_frame, right_front_first_int_frame)


        right_front_hinge2_pt.v2pt_theory(right_front_hinge1_pt, inertial_frame, right_front_first_int_frame)
        right_front_hinge2_pt.a2pt_theory(right_front_hinge1_pt, inertial_frame, right_front_first_int_frame)


        link_inertia_dyadic_first_front_right = inertia(right_front_first_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)

        link_inertia_first_front_right = (link_inertia_dyadic_first_front_right, first_int_front_right_centre_mass)


        first_front_right_int_link_body = RigidBody('Front Right Intermediate Link 1',first_int_front_right_centre_mass,                                 right_front_first_int_frame, mass_link, link_inertia_first_front_right)

        link_grav_force_vector = -mass_link*g*inertial_frame.z

        self.gravity_forces.append((first_int_front_right_centre_mass, link_grav_force_vector))
        first_front_right_int_link_body.potential_energy = -link_grav_force_vector.dot(first_int_front_right_centre_mass.pos_from(O))
        
        
        torque_expr = springTorqueExpression(alpha_1, k_fr1, M0_fr1, E_fr1, f_fr1)
        torque_vector = (torque_expr -c_fr1*alpha_1_d)*right_front_first_int_frame.x
        torque = (right_front_first_int_frame, torque_vector)
        torque1 = (right_swing_arm_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)

#%% Second Inermediate link front right + Wheel
        
        
        right_front_second_int_frame = ReferenceFrame('R_{FR_{i_2}}')
        right_front_second_int_vert_frame = ReferenceFrame('R_{FR_{vert}}')
        right_front_second_int_frame.orient(right_front_first_int_frame,'Axis',(alpha_2, right_front_first_int_frame.x))
        right_front_second_int_vert_frame.orient(right_front_second_int_frame,'Axis',(atan(h_hinge_1/l_hinge_1),right_front_second_int_frame.x))

        
            
        right_front_second_int_frame.set_ang_vel(right_front_first_int_frame,alpha_2_d*right_front_first_int_frame.x)
        
        right_front_second_int_vert_frame.set_ang_vel(right_front_second_int_frame,0) #Rigid with previous frame, but just rotated

        second_int_front_right_centre_mass = Point('P_{rf_{I2_{cm}}')
        front_right_wheel_centre_location = Point('P_{wfr}')

        second_int_front_right_centre_mass.set_pos(right_front_hinge2_pt, d_CM_link*right_front_second_int_frame.y)
        front_right_wheel_centre_location.set_pos(right_front_hinge2_pt, l_link*right_front_second_int_frame.y + wheel_offset*right_front_second_int_frame.x)


        second_int_front_right_centre_mass.v2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)
        second_int_front_right_centre_mass.a2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)
        
        front_right_wheel_centre_location.v2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)
        front_right_wheel_centre_location.a2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)


        # Define a help reference frame that rotates about the z axis of the rotated one
        
        
        front_right_wheel_steer_frame = ReferenceFrame('FRS_f')
        front_right_wheel_frame = ReferenceFrame('FRW_f')
        
        front_right_wheel_steer_frame.orient(right_front_second_int_vert_frame,'Axis',(delta_1,right_front_second_int_vert_frame.z))
        front_right_wheel_frame.orient(front_right_wheel_steer_frame,'Axis',(omega_1,front_right_wheel_steer_frame.x))
        
        front_right_wheel_steer_frame.set_ang_vel(right_front_second_int_vert_frame,delta_1_d*right_front_second_int_vert_frame.z)
        front_right_wheel_frame.set_ang_vel(front_right_wheel_steer_frame, omega_1_d*front_right_wheel_steer_frame.x)
        
        
        link_inertia_dyadic_second_front_right = inertia(right_front_second_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        front_right_wheel_inertia_dyadic = inertia(front_right_wheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)



        link_inertia_second_front_right = (link_inertia_dyadic_second_front_right, second_int_front_right_centre_mass)
        front_right_wheel_inertia = (front_right_wheel_inertia_dyadic, front_right_wheel_centre_location)



        link_second_front_right_body = RigidBody('Right Front Second Intermediate Link',second_int_front_right_centre_mass,                                 right_front_second_int_frame, mass_link, link_inertia_second_front_right)
        link_second_front_right_body.potential_energy = -link_grav_force_vector.dot(second_int_front_right_centre_mass.pos_from(O))
        self.gravity_forces.append((second_int_front_right_centre_mass, link_grav_force_vector))
        
        torque_expr = springTorqueExpression(alpha_2, k_fr2, M0_fr2, E_fr2, f_fr2)
        torque_vector = (torque_expr -c_fr2*alpha_2_d)*right_front_second_int_frame.x
        torque = (right_front_second_int_frame, torque_vector)
        torque1 = (right_front_first_int_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
        
        front_right_wheel_body = RigidBody('Front Right Wheel',front_right_wheel_centre_location,                                 front_right_wheel_frame, wheel_mass, front_right_wheel_inertia)

        wheel_grav_force_vector = -wheel_mass*g*inertial_frame.z
        self.gravity_forces.append((front_right_wheel_centre_location, wheel_grav_force_vector))
        front_right_wheel_body.potential_energy = -wheel_grav_force_vector.dot(front_right_wheel_centre_location.pos_from(O))


#%% First Intermediate Link Back Right
        right_back_first_int_frame = ReferenceFrame('R_{BR_{i_1}}')
        right_back_first_int_frame.orient(right_swing_arm_frame,'Axis',(atan(h_hinge_1/l_hinge_1)+beta_1, right_swing_arm_frame.x))
            
        right_back_first_int_frame.set_ang_vel(right_swing_arm_frame,beta_1_d*right_swing_arm_frame.x)

        #Create points and assign locations
        first_int_back_right_centre_mass = Point('P_{rbI1_{cm}}')
        first_int_back_right_centre_mass.set_pos(right_back_hinge1_pt, -d_CM_link*right_back_first_int_frame.y)
        right_back_hinge2_pt = Point('P_{rbh2}')
        right_back_hinge2_pt.set_pos(right_back_hinge1_pt, -l_link*right_back_first_int_frame.y)
        
        
        first_int_back_right_centre_mass.v2pt_theory(right_back_hinge1_pt, inertial_frame, right_back_first_int_frame)
        first_int_back_right_centre_mass.a2pt_theory(right_back_hinge1_pt, inertial_frame, right_back_first_int_frame)
        
        right_back_hinge2_pt.v2pt_theory(right_back_hinge1_pt, inertial_frame, right_back_first_int_frame)
        right_back_hinge2_pt.a2pt_theory(right_back_hinge1_pt, inertial_frame, right_back_first_int_frame)


        link_inertia_dyadic_first_back_right = inertia(right_back_first_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)


        link_inertia_first_back_right = (link_inertia_dyadic_first_back_right, first_int_back_right_centre_mass)
        
        link_first_back_right_body = RigidBody('First intermediate right back link body',first_int_back_right_centre_mass,                                 right_back_first_int_frame, mass_link, link_inertia_first_back_right)
        
        link_first_back_right_body.potential_energy = -link_grav_force_vector.dot(first_int_back_right_centre_mass.pos_from(O))
        self.gravity_forces.append((first_int_back_right_centre_mass, link_grav_force_vector))
        
        torque_expr = springTorqueExpression(beta_1, k_br1, M0_br1, E_br1, f_br1)
        torque_vector = (torque_expr -c_br1*beta_1_d)*right_back_first_int_frame.x
        torque = (right_back_first_int_frame, torque_vector)
        torque1 = (right_swing_arm_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
        
#%% Second Intermediate Link Back Right + Wheel

        
        right_back_second_int_frame = ReferenceFrame('R_{BR_{I_2}}')
        right_back_second_int_vert_frame = ReferenceFrame('R_{vert_{rb}}')
        right_back_second_int_frame.orient(right_back_first_int_frame,'Axis',(beta_2, right_back_first_int_frame.x))
        right_back_second_int_vert_frame.orient(right_back_second_int_frame,'Axis',(-atan(h_hinge_1/l_hinge_1),right_back_second_int_frame.x))
        
        right_back_second_int_frame.set_ang_vel(right_back_first_int_frame,beta_2_d*right_back_first_int_frame.x)

        right_back_second_int_vert_frame.set_ang_vel(right_back_second_int_frame,0) #Rigid with previous frame, but just rotated
        
        second_int_back_right_centre_mass = Point('P_{rbI2_{cm}}')
        back_right_wheel_centre_location = Point('P_{wbr}')
        
        second_int_back_right_centre_mass.set_pos(right_back_hinge2_pt, -d_CM_link*right_back_second_int_frame.y)
        back_right_wheel_centre_location.set_pos(right_back_hinge2_pt, -l_link*right_back_second_int_frame.y + wheel_offset*right_back_second_int_frame.x)


        # Set linear velocities of these two points
        second_int_back_right_centre_mass.v2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)
        second_int_back_right_centre_mass.a2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)
        
        back_right_wheel_centre_location.v2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)
        back_right_wheel_centre_location.a2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)


        # Define a help reference frame that rotates about the z axis of the rotated one
        
        back_right_wheel_steer_frame = ReferenceFrame('BRS_f')
        back_right_wheel_steer_frame.orient(right_back_second_int_vert_frame,'Axis',(delta_2,right_back_second_int_vert_frame.z))
        back_right_wheel_steer_frame.set_ang_vel(right_back_second_int_vert_frame,delta_2_d*right_back_second_int_vert_frame.z)

        # Define now the wheel reference frame
        
        back_right_wheel_frame = ReferenceFrame('BRW_f')
        back_right_wheel_frame.orient(back_right_wheel_steer_frame,'Axis',(omega_2,back_right_wheel_steer_frame.x))
        back_right_wheel_frame.set_ang_vel(back_right_wheel_steer_frame, omega_2_d*back_right_wheel_steer_frame.x)



        link_inertia_dyadic_second_back_right = inertia(right_back_second_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        back_right_wheel_inertia_dyadic = inertia(back_right_wheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)
        
        
        link_inertia_second_back_right = (link_inertia_dyadic_second_back_right, second_int_back_right_centre_mass)
        back_right_wheel_inertia = (back_right_wheel_inertia_dyadic, back_right_wheel_centre_location)


        link_second_back_right_body = RigidBody('Right Back Second Intermediate Link',second_int_back_right_centre_mass,                                 right_back_second_int_frame, mass_link, link_inertia_second_back_right)
        link_second_back_right_body.potential_energy = -link_grav_force_vector.dot(second_int_back_right_centre_mass.pos_from(O))
        self.gravity_forces.append((second_int_back_right_centre_mass, link_grav_force_vector))

        torque_expr = springTorqueExpression(beta_2, k_br2, M0_br2, E_br2, f_br2)
        torque_vector = (torque_expr -c_br2*beta_2_d)*right_back_second_int_frame.x
        torque = (right_back_second_int_frame, torque_vector)
        torque1 = (right_back_first_int_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        

        back_right_wheel_body = RigidBody('Back Right Wheel',back_right_wheel_centre_location,                                 back_right_wheel_frame, wheel_mass, back_right_wheel_inertia)

        back_right_wheel_body.potential_energy = -wheel_grav_force_vector.dot(back_right_wheel_centre_location.pos_from(O))
        self.gravity_forces.append((back_right_wheel_centre_location, wheel_grav_force_vector))


#%% Left Swing arm


        left_swing_arm_hinge_pt = Point('P_{SA_l}')
        left_swing_arm_hinge_pt.set_pos(rover_centre_mass, -d_h*rover_frame.x)


        left_swing_arm_hinge_pt.v2pt_theory(rover_centre_mass,inertial_frame,rover_frame)
        left_swing_arm_hinge_pt.a2pt_theory(rover_centre_mass,inertial_frame,rover_frame)


        # Since the two swing arms motion are coupled, such that the rotate of same quantity but in opposite direction it is not needed to define a new generalized coordinate.

        left_swing_arm_frame = ReferenceFrame('R_{SA_l}')
        left_swing_arm_frame.orient(rover_frame, 'Axis',(-phi,rover_frame.x))

    
        #Assign manual angular velocity
        left_swing_arm_frame.set_ang_vel(rover_frame,-phi_d*rover_frame.x)
        

        left_swing_centre_mass = Point('SA_{l_{CM}}')
        left_swing_centre_mass.set_pos(left_swing_arm_hinge_pt, -off_swing_dst*left_swing_arm_frame.z)
        
        #Set center mass linear velocity
        left_swing_centre_mass.v2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)
        left_swing_centre_mass.a2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)
        
        
        #Define now the first hinge location for the front leg
        left_front_hinge1_pt = Point('P_{h_l_f1}')
        left_front_hinge1_pt.set_pos(left_swing_arm_hinge_pt, l_hinge_1*left_swing_arm_frame.y -h_hinge_1*left_swing_arm_frame.z)
        
        #Compute hinge liner speed using poisson rule
        left_front_hinge1_pt.v2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)
        left_front_hinge1_pt.a2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)


        #Define now the first hinge location for the back leg
        left_back_hinge1_pt = Point('P_{h_l_b1}')
        left_back_hinge1_pt.set_pos(left_swing_arm_hinge_pt, -l_hinge_1*left_swing_arm_frame.y -h_hinge_1*left_swing_arm_frame.z)
        
        left_back_hinge1_pt.v2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)
        left_back_hinge1_pt.a2pt_theory(left_swing_arm_hinge_pt,inertial_frame,left_swing_arm_frame)



        swing_arm_inertia_dyadic_left = inertia(left_swing_arm_frame,swing_arm_i_xx,swing_arm_i_yy,swing_arm_i_zz,swing_arm_i_xy,                                  swing_arm_i_yz, swing_arm_i_xz)
        
        swing_arm_left_inertia = (swing_arm_inertia_dyadic_left, left_swing_centre_mass)
        
        swing_arm_left_body = RigidBody('Left Swing Arm Body',left_swing_centre_mass,                                 left_swing_arm_frame, mass_swing_arm, swing_arm_left_inertia)
        swing_arm_left_body.potential_energy = -swing_arm_grav_force_vector.dot(left_swing_centre_mass.pos_from(O))

        self.gravity_forces.append((left_swing_centre_mass, swing_arm_grav_force_vector))


        #%% First intermediate link - Left Front Leg
        
        left_front_first_int_frame = ReferenceFrame('R_{FL_{i_1}}')
        #atan(h_hinge_1/l_hinge_1)
        left_front_first_int_frame.orient(left_swing_arm_frame,'Axis',(-atan(h_hinge_1/l_hinge_1)+gamma_1, left_swing_arm_frame.x))
        left_front_first_int_frame.set_ang_vel(left_swing_arm_frame, gamma_1_d*left_swing_arm_frame.x)
        
        
        #Create points and assign locations
        first_int_front_left_centre_mass = Point('P_{lfI1_{cm}}')
        first_int_front_left_centre_mass.set_pos(left_front_hinge1_pt, d_CM_link*left_front_first_int_frame.y)
        left_front_hinge2_pt = Point('P_{h_lf_2}')
        left_front_hinge2_pt.set_pos(left_front_hinge1_pt, l_link*left_front_first_int_frame.y)
        
        first_int_front_left_centre_mass.v2pt_theory(left_front_hinge1_pt, inertial_frame, left_front_first_int_frame)
        first_int_front_left_centre_mass.a2pt_theory(left_front_hinge1_pt, inertial_frame, left_front_first_int_frame)
        
        left_front_hinge2_pt.v2pt_theory(left_front_hinge1_pt, inertial_frame, left_front_first_int_frame)
        left_front_hinge2_pt.a2pt_theory(left_front_hinge1_pt, inertial_frame, left_front_first_int_frame)
        
        
        link_inertia_dyadic_first_front_left = inertia(left_front_first_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        
        
        link_inertia_first_front_left = (link_inertia_dyadic_first_front_left, first_int_front_left_centre_mass)
        
        link_first_front_left_body = RigidBody('First intermediate left front link body',first_int_front_left_centre_mass,                                 left_front_first_int_frame, mass_link, link_inertia_first_front_left)
        
        link_first_front_left_body.potential_energy = -link_grav_force_vector.dot(first_int_front_left_centre_mass.pos_from(O))
        self.gravity_forces.append((first_int_front_left_centre_mass, link_grav_force_vector))

        torque_expr = springTorqueExpression(gamma_1, k_fl1, M0_fl1, E_fl1, f_fl1)
        torque_vector = (torque_expr -c_fl1*gamma_1_d)*left_front_first_int_frame.x
        torque = (left_front_first_int_frame, torque_vector)
        torque1 = (left_swing_arm_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
#%% Second Intermediate Link - Left Front Leg

        left_front_second_int_frame = ReferenceFrame('R_{FL{i_2}}')
        left_front_second_int_vert_frame = ReferenceFrame('R_{Vert_{FL}}')
        left_front_second_int_frame.orient(left_front_first_int_frame,'Axis',(gamma_2, left_front_first_int_frame.x))
        left_front_second_int_vert_frame.orient(left_front_second_int_frame,'Axis',(atan(h_hinge_1/l_hinge_1),left_front_second_int_frame.x))
        
        
        
        left_front_second_int_frame.set_ang_vel(left_front_first_int_frame,gamma_2_d*left_front_first_int_frame.x)
        
        left_front_second_int_vert_frame.set_ang_vel(left_front_second_int_frame,0) #Rigid with previous frame, but just rotated
        
        # Remember: wheel_offset = symbols('d_W')
        second_int_front_left_centre_mass = Point('P_{lfI2_{cm}}')
        front_left_wheel_centre_location = Point('P_{wfl}')
        
        second_int_front_left_centre_mass.set_pos(left_front_hinge2_pt, d_CM_link*left_front_second_int_frame.y)
        front_left_wheel_centre_location.set_pos(left_front_hinge2_pt, l_link*left_front_second_int_frame.y - wheel_offset*left_front_second_int_frame.x)
        
        
        # Set linear velocities of these two points
        second_int_front_left_centre_mass.v2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        second_int_front_left_centre_mass.a2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        
        front_left_wheel_centre_location.v2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        front_left_wheel_centre_location.a2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        
        
        # Define a help reference frame that rotates about the z axis of the rotated one
        
        front_left_wheel_steer_frame = ReferenceFrame('FLS_f')
        front_left_wheel_steer_frame.orient(left_front_second_int_vert_frame,'Axis',(delta_3,left_front_second_int_vert_frame.z))
        front_left_wheel_steer_frame.set_ang_vel(left_front_second_int_vert_frame,delta_3_d*left_front_second_int_vert_frame.z)
        
        
        # Define now the wheel reference frame
        
        
        front_left_wheel_frame = ReferenceFrame('FLW_f')
        front_left_wheel_frame.orient(front_left_wheel_steer_frame,'Axis',(omega_3,front_left_wheel_steer_frame.x))
        front_left_wheel_frame.set_ang_vel(front_left_wheel_steer_frame, omega_3_d*front_left_wheel_steer_frame.x)
        
        
        link_inertia_dyadic_second_front_left = inertia(left_front_second_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        front_left_wheel_inertia_dyadic = inertia(front_left_wheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)
        
        
        link_inertia_second_front_left = (link_inertia_dyadic_second_front_left, second_int_front_left_centre_mass)
        front_left_wheel_inertia = (front_left_wheel_inertia_dyadic, front_left_wheel_centre_location)
        
        
        link_second_front_left_body = RigidBody('Left Front Second Intermediate Link',second_int_front_left_centre_mass,                                 left_front_second_int_frame, mass_link, link_inertia_second_front_left)
        link_second_front_left_body.potential_energy = -link_grav_force_vector.dot(second_int_front_left_centre_mass.pos_from(O))
        self.gravity_forces.append((second_int_front_left_centre_mass, link_grav_force_vector))
        
        torque_expr = springTorqueExpression(gamma_2, k_fl2, M0_fl2, E_fl2, f_fl2)
        torque_vector = (torque_expr -c_fl2*gamma_2_d)*left_front_second_int_frame.x
        torque = (left_front_second_int_frame, torque_vector)
        torque1 = (left_front_first_int_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
        
        
        front_left_wheel_body = RigidBody('Front Left Wheel',front_left_wheel_centre_location,                                 front_left_wheel_frame, wheel_mass, front_left_wheel_inertia)
        
        front_left_wheel_body.potential_energy = -wheel_grav_force_vector.dot(front_left_wheel_centre_location.pos_from(O))
        self.gravity_forces.append((front_left_wheel_centre_location, wheel_grav_force_vector))


#%% First intermediate link - Left Back Leg
        
        left_back_first_int_frame = ReferenceFrame('R_{BL{I_1}}')
        #atan(h_hinge_1/l_hinge_1)
        left_back_first_int_frame.orient(left_swing_arm_frame,'Axis',(atan(h_hinge_1/l_hinge_1)+eps_1, left_swing_arm_frame.x))
        left_back_first_int_frame.set_ang_vel(left_swing_arm_frame,eps_1_d*left_swing_arm_frame.x)
        
        #Create points and assign locations
        first_int_back_left_centre_mass = Point('P_{lbI1_{cm}}')
        first_int_back_left_centre_mass.set_pos(left_back_hinge1_pt, -d_CM_link*left_back_first_int_frame.y)
        left_back_hinge2_pt = Point('P_{lbh_2}')
        left_back_hinge2_pt.set_pos(left_back_hinge1_pt, -l_link*left_back_first_int_frame.y)
        
        
        first_int_back_left_centre_mass.v2pt_theory(left_back_hinge1_pt, inertial_frame, left_back_first_int_frame)
        first_int_back_left_centre_mass.a2pt_theory(left_back_hinge1_pt, inertial_frame, left_back_first_int_frame)
        
        left_back_hinge2_pt.v2pt_theory(left_back_hinge1_pt, inertial_frame, left_back_first_int_frame)
        left_back_hinge2_pt.a2pt_theory(left_back_hinge1_pt, inertial_frame, left_back_first_int_frame)
        
        
        link_inertia_dyadic_first_back_left = inertia(left_back_first_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        
        link_inertia_first_back_left = (link_inertia_dyadic_first_back_left, first_int_back_left_centre_mass)
        
        link_first_back_left_body = RigidBody('First intermediate left back link body',first_int_back_left_centre_mass,                                 left_back_first_int_frame, mass_link, link_inertia_first_back_left)
        
        link_first_back_left_body.potential_energy = -link_grav_force_vector.dot(first_int_back_left_centre_mass.pos_from(O))
        self.gravity_forces.append((first_int_back_left_centre_mass, link_grav_force_vector))
        
        
        torque_expr = springTorqueExpression(eps_1, k_bl1, M0_bl1, E_bl1, f_bl1)
        torque_vector = (torque_expr - c_bl1*eps_1_d)*left_back_first_int_frame.x
        torque = (left_back_first_int_frame, torque_vector)
        torque1 = (left_swing_arm_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
        
#%% End Link - Left Back Leg
        
        
        left_back_second_int_frame = ReferenceFrame('R_{BL_{I_2}}')
        left_back_second_int_vert_frame = ReferenceFrame('R_{vert_{lb}}')
        left_back_second_int_frame.orient(left_back_first_int_frame,'Axis',(eps_2, left_back_first_int_frame.x))
        left_back_second_int_vert_frame.orient(left_back_second_int_frame,'Axis',(-atan(h_hinge_1/l_hinge_1),left_back_second_int_frame.x))
        
        
        
        left_back_second_int_frame.set_ang_vel(left_back_first_int_frame,eps_2_d*left_back_first_int_frame.x)
        
        left_back_second_int_vert_frame.set_ang_vel(left_back_second_int_frame,0) #Rigid with previous frame, but just rotated
        
        second_int_back_left_centre_mass = Point('P_{lbI2_{cm}}')
        back_left_wheel_centre_location = Point('P_{wbl}')
        
        second_int_back_left_centre_mass.set_pos(left_back_hinge2_pt, -d_CM_link*left_back_second_int_frame.y)
        back_left_wheel_centre_location.set_pos(left_back_hinge2_pt, -l_link*left_back_second_int_frame.y - wheel_offset*left_back_second_int_frame.x)
        
        
        # Set linear velocities of these two points
        second_int_back_left_centre_mass.v2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        second_int_back_left_centre_mass.a2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        
        back_left_wheel_centre_location.v2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        back_left_wheel_centre_location.a2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        
        
        # Define a help reference frame that rotates about the z axis of the rotated one
        back_left_wheel_steer_frame = ReferenceFrame('BLS_f')
        back_left_wheel_steer_frame.orient(left_back_second_int_vert_frame,'Axis',(delta_4,left_back_second_int_vert_frame.z))
        back_left_wheel_steer_frame.set_ang_vel(left_back_second_int_vert_frame,delta_4_d*left_back_second_int_vert_frame.z)
        
        # Define now the wheel reference frame
        
        
        back_left_wheel_frame = ReferenceFrame('BLW_f')
        back_left_wheel_frame.orient(back_left_wheel_steer_frame,'Axis',(omega_4,back_left_wheel_steer_frame.x))
        back_left_wheel_frame.set_ang_vel(back_left_wheel_steer_frame, omega_4_d*back_left_wheel_steer_frame.x)
        
        
        link_inertia_dyadic_second_back_left = inertia(left_back_second_int_frame,link_inertia_xx,link_inertia_yy,link_inertia_zz)
        back_left_wheel_inertia_dyadic = inertia(back_left_wheel_frame, wheel_inertia_xx, wheel_inertia_yy, wheel_inertia_zz)
        
        
        link_inertia_second_back_left = (link_inertia_dyadic_second_back_left, second_int_back_left_centre_mass)
        back_left_wheel_inertia = (back_left_wheel_inertia_dyadic, back_left_wheel_centre_location)
        
        
        link_second_back_left_body = RigidBody('Left Back Second Intermediate Link',second_int_back_left_centre_mass,                                 left_back_second_int_frame, mass_link, link_inertia_second_back_left)
        link_second_back_left_body.potential_energy = -link_grav_force_vector.dot(second_int_back_left_centre_mass.pos_from(O))
        self.gravity_forces.append((second_int_back_left_centre_mass, link_grav_force_vector))
        
        torque_expr = springTorqueExpression(eps_2, k_bl2, M0_bl2, E_bl2, f_bl2)
        torque_vector = (torque_expr - c_bl2*eps_2_d)*left_back_second_int_frame.x
        torque = (left_back_second_int_frame, torque_vector)
        torque1 = (left_back_first_int_frame, -torque_vector)
        self.torques.append(torque)
        self.__reaction_torque.append(torque1)
        
        
        back_left_wheel_body = RigidBody('Back Left Wheel',back_left_wheel_centre_location,                                 back_left_wheel_frame, wheel_mass, back_left_wheel_inertia)
        
        back_left_wheel_body.potential_energy = -wheel_grav_force_vector.dot(back_left_wheel_centre_location.pos_from(O))
        self.gravity_forces.append((back_left_wheel_centre_location, wheel_grav_force_vector))

#%%  Rigid Bodied Build up and External Loads definition
# In this section we will first combine all the defined rigid bodies in a list, and then we will define the external loads/torques acting on the model.
        
        # Insert 4 Points which are the four legs end points 
        end_right_front = right_front_hinge2_pt.locatenew('P_{end_leg_fr}', l_link*right_front_second_int_frame.y)
        end_left_front = left_front_hinge2_pt.locatenew('P_{end_leg_fl}', l_link*left_front_second_int_frame.y)
        end_right_back = right_back_hinge2_pt.locatenew('P_{end_leg_br}', -l_link*right_back_second_int_frame.y)
        end_left_back = left_back_hinge2_pt.locatenew('P_{end_leg_bl}', -l_link*left_back_second_int_frame.y)
        
        end_right_front.v2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)
        end_right_front.a2pt_theory(right_front_hinge2_pt, inertial_frame, right_front_second_int_frame)
        
        end_left_front.v2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        end_left_front.a2pt_theory(left_front_hinge2_pt, inertial_frame, left_front_second_int_frame)
        
        end_right_back.v2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)
        end_right_back.a2pt_theory(right_back_hinge2_pt, inertial_frame, right_back_second_int_frame)
        
        end_left_back.v2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        end_left_back.a2pt_theory(left_back_hinge2_pt, inertial_frame, left_back_second_int_frame)
        
        
        
        # Insert 4 Points which represent the contact wheels point
        contact_right_front = front_right_wheel_centre_location.locatenew('P_{contact_right_front}', -wheel_radius*rover_frame.z)
        contact_right_back = back_right_wheel_centre_location.locatenew('P_{contact_right_back}', -wheel_radius*rover_frame.z)
        contact_left_front = front_left_wheel_centre_location.locatenew('P_{contact_left_front}', -wheel_radius*rover_frame.z)
        contact_left_back = back_left_wheel_centre_location.locatenew('P_{contact_left_back}', -wheel_radius*rover_frame.z)
        
        contact_right_front.v2pt_theory(front_right_wheel_centre_location, inertial_frame, front_right_wheel_frame)
        contact_right_front.a2pt_theory(front_right_wheel_centre_location, inertial_frame, front_right_wheel_frame)
        
        contact_right_back.v2pt_theory(back_right_wheel_centre_location, inertial_frame, back_right_wheel_frame)
        contact_right_back.a2pt_theory(back_right_wheel_centre_location, inertial_frame, back_right_wheel_frame)
        
        contact_left_front.v2pt_theory(front_left_wheel_centre_location, inertial_frame, front_left_wheel_frame)
        contact_left_front.a2pt_theory(front_left_wheel_centre_location, inertial_frame, front_left_wheel_frame)
        
        contact_left_back.v2pt_theory(back_left_wheel_centre_location, inertial_frame, back_left_wheel_frame)
        contact_left_back.a2pt_theory(back_left_wheel_centre_location, inertial_frame, back_left_wheel_frame)
        
        
        
        # Insert 4 Parametric points which lies inside the surface of the wheel. Tpo do this we need to
        # define other four reference frames which are the frames of the contact point
        
        
        #We define now 4 parametric points
        a1,a2,a3,a4 = symbols('a_1, a_2, a_3, a_4')
        rho1,rho2,rho3,rho4 = symbols('rho1, rho2, rho3, rho4')
        
        
        
        #FRAMES:
        frame_par_contact_front_right = ReferenceFrame('R_{par_cont_fr}')
        frame_par_contact_front_right.orient(front_right_wheel_steer_frame,'Axis', (a1,front_right_wheel_steer_frame.x))
        
        frame_par_contact_back_right = ReferenceFrame('R_{par_cont_br}')
        frame_par_contact_back_right.orient(back_right_wheel_steer_frame,'Axis', (a2, back_right_wheel_steer_frame.x))
        
        frame_par_contact_front_left = ReferenceFrame('R_{par_cont_fl}')
        frame_par_contact_front_left.orient(front_left_wheel_steer_frame,'Axis', (a3, front_left_wheel_steer_frame.x))
        
        frame_par_contact_back_left = ReferenceFrame('R_{par_cont_bl}')
        frame_par_contact_back_left.orient(back_left_wheel_steer_frame,'Axis', (a4, back_left_wheel_steer_frame.x))
        
        
        
        #CONTACT POINTS:
        contact_right_front_parametric = front_right_wheel_centre_location.locatenew('P_{contact_right_front_parametric}',\
                                                                                     -rho1*frame_par_contact_front_right.z)
        
        contact_right_back_parametric = back_right_wheel_centre_location.locatenew('P_{contact_right_back_parametric}',\
                                                                                     -rho2*frame_par_contact_back_right.z)
        
        contact_left_front_parametric = front_left_wheel_centre_location.locatenew('P_{contact_left_front_parametric}',\
                                                                                     -rho3*frame_par_contact_front_left.z)
            
        contact_left_back_parametric = back_left_wheel_centre_location.locatenew('P_{contact_left_back_parametric}',\
                                                                                   -rho4*frame_par_contact_back_left.z)
        
        
        contact_right_front_parametric.v2pt_theory(front_right_wheel_centre_location, inertial_frame, front_right_wheel_frame)
        contact_right_front_parametric.a2pt_theory(front_right_wheel_centre_location, inertial_frame, front_right_wheel_frame)
        
        contact_right_back_parametric.v2pt_theory(back_right_wheel_centre_location, inertial_frame, back_right_wheel_frame)
        contact_right_back_parametric.a2pt_theory(back_right_wheel_centre_location, inertial_frame, back_right_wheel_frame)
        
        contact_left_front_parametric.v2pt_theory(front_left_wheel_centre_location, inertial_frame, front_left_wheel_frame)
        contact_left_front_parametric.a2pt_theory(front_left_wheel_centre_location, inertial_frame, front_left_wheel_frame)
        
        contact_left_back_parametric.v2pt_theory(back_left_wheel_centre_location, inertial_frame, back_left_wheel_frame)
        contact_left_back_parametric.a2pt_theory(back_left_wheel_centre_location, inertial_frame, back_left_wheel_frame)
            
        

        #Define now the point velocities and accelerations
        
        
        

        self.bodies = [rover_body, right_swing_arm_body, 
                  first_front_right_int_link_body, 
                  link_second_front_right_body,
                  front_right_wheel_body,
                  link_first_back_right_body,
                  link_second_back_right_body,
                  back_right_wheel_body,
                  swing_arm_left_body,
                  link_first_front_left_body,
                  link_second_front_left_body,
                  front_left_wheel_body,
                  link_first_back_left_body,
                  link_second_back_left_body,
                  back_left_wheel_body]

        self.frames = [inertial_frame,
                       rover_intermediate_frame,
                       rover_pitch_frame,
                       rover_frame,
                       right_swing_arm_frame,
                       right_front_first_int_frame,
                       right_front_second_int_frame,
                       right_front_second_int_vert_frame,
                       front_right_wheel_steer_frame,
                       front_right_wheel_frame,
                       right_back_first_int_frame,
                       right_back_second_int_frame,
                       right_back_second_int_vert_frame,
                       back_right_wheel_steer_frame,
                       back_right_wheel_frame,
                       left_swing_arm_frame,
                       left_front_first_int_frame,
                       left_front_second_int_frame,
                       left_front_second_int_vert_frame,
                       front_left_wheel_steer_frame,
                       front_left_wheel_frame,
                       left_back_first_int_frame,
                       left_back_second_int_frame,
                       left_back_second_int_vert_frame,
                       back_left_wheel_steer_frame,
                       back_left_wheel_frame,
                       ]
        
        
        self.points = [O,
                       rover_centre_mass,
                       right_swing_arm_hinge_pt,
                       right_swing_centre_mass,
                       right_front_hinge1_pt,
                       right_back_hinge1_pt,
                       first_int_front_right_centre_mass,
                       right_front_hinge2_pt,
                       second_int_front_right_centre_mass,
                       front_right_wheel_centre_location,
                       first_int_back_right_centre_mass,
                       right_back_hinge2_pt,
                       second_int_back_right_centre_mass,
                       back_right_wheel_centre_location,
                       left_swing_arm_hinge_pt,
                       left_swing_centre_mass,
                       left_front_hinge1_pt,
                       left_back_hinge1_pt,
                       first_int_front_left_centre_mass,
                       left_front_hinge2_pt,
                       second_int_front_left_centre_mass,
                       front_left_wheel_centre_location,
                       first_int_back_left_centre_mass,
                       left_back_hinge2_pt,
                       second_int_back_left_centre_mass,
                       back_left_wheel_centre_location,
                       end_right_front,
                       end_left_front,
                       end_right_back,
                       end_left_back,
                       contact_right_front,
                       contact_right_back,
                       contact_left_front,
                       contact_left_back,                       
                       ]
        
        
        self.contact_parametric['points'] = [contact_right_front_parametric,
                                             contact_right_back_parametric,
                                             contact_left_front_parametric,
                                             contact_left_back_parametric]
        
        self.contact_parametric['frames'] = [frame_par_contact_front_right,
                                             frame_par_contact_back_right,
                                             frame_par_contact_front_left,
                                             frame_par_contact_back_left]
        
        
        self.contact_parametric['variables'] = list(zip([rho1, rho2, rho3, rho4], [a1, a2, a3, a4]))
        
        
        self.__createWheel()
        ground_pts = self.loadGroundPoints()
        ground_pts1 = np.array((ground_pts['x'],ground_pts['y'], ground_pts['z'])).T
        self.__ground_new = ground_pts1
        
        X = ground_pts1[:,0].reshape((ground_pts['rows'], ground_pts['columns']))
        Y = ground_pts1[:,1].reshape((ground_pts['rows'], ground_pts['columns']))
        Z = ground_pts1[:,2].reshape((ground_pts['rows'], ground_pts['columns']))
        
        #self.__ground_interpolator = interpolate.SmoothBivariateSpline(ground_pts[:,0],
        #                                                               ground_pts[:,1],ground_pts[:,2])
        
        self.__ground_interpolator = interpolate.RectBivariateSpline(Y[:,0],X[0,:], Z)
        
        steer_torques = list(symbols('TS_:4')) #Steer torque
        drive_torques = list(symbols('TD_:4')) #Drive torque
        
        self.__driving_torques = steer_torques + drive_torques
        self.__current_driving_torque = np.zeros((len(steer_torques)+len(drive_torques),))
        
        
        steer_frames = ['FRS_f','BRS_f', 'FLS_f', 'BLS_f']
        wheel_frames = ['FRW_f','BRW_f', 'FLW_f', 'BLW_f']
        link_frames = ['R_{FR_{i_2}}', 'R_{BR_{I_2}}', 'R_{FL{i_2}}', 'R_{BL_{I_2}}']
        
        t = list()
        
        for i in range(len(steer_torques)):
            for j, frame in enumerate(self.frames):
                if frame.name == steer_frames[i]:
                    break
            
            for p, frame in enumerate(self.frames):
                if frame.name == wheel_frames[i]:
                    break
            for o, frame in enumerate(self.frames):
                if frame.name == link_frames[i]:
                    break
                
                
            #Add link frames if action-reaction
            torque_vec = steer_torques[i]*self.frames[j].z + drive_torques[i]*self.frames[p].x
            torque_tup = (self.frames[p], torque_vec)
            torque_tup1 = (self.frames[o], -torque_vec)
            t.append(torque_tup)
            #t.append(torque_tup1)
            
        self.torques += t
            
                    
        

#%% #### Torques Build Up
#TO DO


    def lambdifyAll(self, method = None):
        """
        This method attempts on creating lamda functions for:
            - Frames orientation (input arguments q) 
            - Frames angular speed (input arguments q q_d) In this order!
            - Frames angular accelerations (input arguments q, q_d, q_dd) in this order!
            
            - Points location (input arguments q)
            - Points linear velocities (input arguments q, q_d) in this order
            - Points linear accelerations (input arguments q, q_d, q_dd) in this order!
            
        Parameters
        ----------
            method : string, optional
            Indication of which method to be used. If no arguments are passed then
            the method is chosen to be the one indicated in the config file. Default is None.
            Valid methods are in the following list [lambify, autowrap - f2py, autowrap - cython, theano]
        """
        
        if method is not None:
            if method != 'autowrap - f2py' and method != 'autowrap - cython' and method != 'lambdify' and method != 'theano':
                raise ValueError('Invalid input argument')
        # print('Attempting to Lambdify the most functions as possible....'+'\n')
        # sleep(0.1)
        else:
            method = self.config.get('Simulator','lambdifier')
            if method != 'autowrap - f2py' and method != 'autowrap - cython' and method != 'lambdify' and method != 'theano':
                raise ValueError('Invalid input argument')
                
        save_dir = None
        if self.config.getboolean('Simulator','autowrap save files'):
            save_dir = self.config.get('Simulator','autowrap save dir')

        self.__createVariablesMapping()
        
        cwd = os.getcwd()
        
        if not os.path.isdir(os.path.join(cwd,'data', 'wrappers')):
            os.makedirs(os.path.join(cwd,'data', 'wrappers'))
        
        
        
        
        acc_list = [diff(vel) for vel in self.gen_speeds]
        
        
        for i in tqdm(range(len(self.frames)), desc='Lambdifying reference frames tranformation matrices'):
            if method == 'lambdify':
                f = lambdify(self.gen_coord,self.frames[i].dcm(self.frames[0]), "numpy")
            elif method == 'theano':
                f = theano_function(self.gen_coord, [self.frames[i].dcm(self.frames[0])], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                f = autowrap(msubs(self.frames[i].dcm(self.frames[0]), self.__map['simplyfied']['mapping']['gen-coords'])\
                             ,tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'], backend = 'f2py')
            else:
                f = autowrap(msubs(self.frames[i].dcm(self.frames[0]), self.__map['dummyfied']['mapping']['gen-coords'])\
                             ,tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords'], backend = 'cython')
            self.lambdified_frames.append(f)
            
            
        tmp_coords = self.gen_coord + self.gen_speeds
        for i in tqdm(range(len(self.frames)), desc='Lambdifying reference frames angular velocities'):
            if method == 'lambdify':
                f = lambdify(tmp_coords, self.frames[i].ang_vel_in(self.frames[0]).to_matrix(self.frames[0]), "numpy")
            elif method == 'theano':
                f = theano_function(tmp_coords, [self.frames[i].ang_vel_in(self.frames[0]).to_matrix(self.frames[0])], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                expr = self.frames[i].ang_vel_in(self.frames[0]).to_matrix(self.frames[0])
                f = autowrap(msubs(expr, self.__map['simplyfied']['mapping']['gen-merged']), \
                             tempdir=save_dir, args=self.__map['simplyfied']['subs merged'], backend = 'f2py')
            else:
                expr = self.frames[i].ang_vel_in(self.frames[0]).to_matrix(self.frames[0])
                f = autowrap(msubs(expr, self.__map['dummyfied']['mapping']['gen-merged']), \
                             tempdir=save_dir, args=self.__map['dummyfied']['subs merged'], backend = 'cython')
    
            self.lambdified_frames_ang_vel.append(f)
        
        tmp_coords = self.gen_coord + self.gen_speeds + acc_list
        # for i in tqdm(range(len(self.frames)), desc='Lambdifying reference frames angular accelerations'):
        #     f = lambdify(tmp_coords, self.frames[i].ang_acc_in(self.frames[0]).to_matrix(self.frames[0]))
        #     self.lambdified_frames_ang_acc.append(f)
            
        
            
        
            
        
        for i in tqdm(range(len(self.points)), desc='Lambdifying points positions'):
            
            if method == 'lambdify':
                f = lambdify(self.gen_coord, self.points[i].pos_from(self.points[0]).to_matrix(self.frames[0]).T, "numpy")
            elif method == 'theano':
                f = theano_function(self.gen_coord, [self.points[i].pos_from(self.points[0]).to_matrix(self.frames[0]).T], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                expr = self.points[i].pos_from(self.points[0]).to_matrix(self.frames[0]).T
                f = autowrap(msubs(expr, self.__map['simplyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'], backend = 'f2py')
            else:
                expr = self.points[i].pos_from(self.points[0]).to_matrix(self.frames[0]).T
                f = autowrap(msubs(expr, self.__map['dummyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords'], backend = 'cython')
            self.lambdified_points.append(f)
            
        
        tmp_coords = self.gen_coord + self.gen_speeds
        for i in tqdm(range(len(self.points)), desc='Lambdifying points linear velocities'):
            expr = self.points[i].vel(self.frames[0]).to_matrix(self.frames[0]).T
            if method == 'lambdify':
                f = lambdify(tmp_coords, expr, "numpy")
            elif method == 'theano':
                f = theano_function(tmp_coords, [expr], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                f = autowrap(msubs(expr, self.__map['simplyfied']['mapping']['gen-merged']), \
                             tempdir=save_dir, args=self.__map['simplyfied']['subs merged'], backend = 'f2py')
            else:
                f = autowrap(msubs(expr, self.__map['dummyfied']['mapping']['gen-merged']), \
                             tempdir=save_dir, args=self.__map['dummyfied']['subs merged'], backend = 'cython')
            
            self.lambdified_points_vel.append(f)
            
        tmp_coords = self.gen_coord + self.gen_speeds + acc_list
        # for i in tqdm(range(len(self.points)), desc='Lambdifying points linear accelerations:'):
            
        #     f = lambdify(tmp_coords, self.points[i].acc(self.frames[0]).to_matrix(self.frames[0]).T)
        #     self.lambdified_points.append(f)
            
            
            
            
        contact_names = ['Front-Right-Wheel', 'Back-Right-Wheel', 'Front-Left-Wheel', 'Back-Left-Wheel']
        point_id = ['P_{contact_right_front}', 'P_{contact_right_back}','P_{contact_left_front}', 'P_{contact_left_back}']
        func_list = list()
        func_list_pos = list()
        func_list_vel = list()
        
        for k, name in enumerate(point_id):
            
            for i,point in enumerate(self.points):
                if point.name == name:
                    break
            
            # i-index of the point
            tmp_matrix = Matrix()
            
            
            tmp_partial = self.points[i].partial_velocity(self.frames[0], *self.gen_speeds)
            tmp_pos = self.lambdified_points[i]
            tmp_vel = self.lambdified_points_vel[i]
            
            for j in tqdm(range(len(tmp_partial)),desc='Lambdifying function for partial velocities of the contact point for the '+contact_names[k]):
                tmp_matrix = tmp_matrix.col_join(tmp_partial[j].to_matrix(self.frames[0]).T)
            
            #mprint(find_dynamicsymbols(tmp_matrix,reference_frame=self.frames[0]))
            
            
            if method == 'lambdify':
                f_n = lambdify(self.gen_coord, tmp_matrix, "numpy")
            elif method == 'theano':
                f_n = theano_function(self.gen_coord, [tmp_matrix], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                f_n = autowrap(msubs(tmp_matrix, self.__map['simplyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'], backend = 'f2py')
            else:
                f_n = autowrap(msubs(tmp_matrix, self.__map['dummyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords'], backend = 'cython')
            func_list.append(f_n)
            func_list_pos.append(tmp_pos)
            func_list_vel.append(tmp_vel)
            
        self.lambdified_partial_velocities = list(zip(contact_names, point_id, func_list_pos, func_list_vel, func_list))
        
        
    
        point_id = list()
        pos_list = list()
        vel_list = list()
        func_list = list()
        vel_projected_list = list()
        
        
        for n in range(len(self.contact_parametric['points'])):
            point_id.append(self.contact_parametric['points'][n].name)
            
            
            tmp_matrix = Matrix()
            
            
            tmp_partial = self.contact_parametric['points'][n].partial_velocity(self.frames[0], *self.gen_speeds)
            
            for j in range(len(tmp_partial)):
                tmp_matrix = tmp_matrix.col_join(tmp_partial[j].express(
                    self.contact_parametric['frames'][n]).to_matrix(self.contact_parametric['frames'][n]).T)
            
            
            
            expr_pos = self.contact_parametric['points'][n].pos_from(self.points[0]).to_matrix(self.frames[0]).T #Point position expression
            expr_vel = self.contact_parametric['points'][n].vel(self.frames[0]).express(self.contact_parametric['frames'][n]).to_matrix(self.frames[0]).T #Velocity expression
            part_vel_expr = tmp_matrix
            vel_tmp = self.contact_parametric['points'][n].vel(self.frames[0]).express(self.contact_parametric['frames'][n]) #velocity of point expressed in its frame
            vel_projected = vel_tmp - vel_tmp.args[0][0][2]*self.contact_parametric['frames'][n].z #In its frame
            
            
            expr_proj_vel = vel_projected.to_matrix(self.contact_parametric['frames'][n])
            
            
            symb = list(self.contact_parametric['variables'][n])
            
            if method == 'lambdify':
                f_p = lambdify(self.gen_coord + symb, expr_pos, "numpy")
                f_part = lambdify(self.gen_coord + symb, part_vel_expr, "numpy")
                tmp_coord = self.gen_coord + self.gen_speeds
                
                f_v = lambdify(tmp_coord + symb, expr_vel, "numpy")
                f_v_proj = lambdify(tmp_coord + symb, expr_proj_vel, "numpy")
            elif method == 'theano':
                f_p = theano_function(self.gen_coord + symb, [expr_pos], on_unused_input='ignore')
                f_part = theano_function(self.gen_coord+ symb, [part_vel_expr], on_unused_input='ignore')
                tmp_coord = self.gen_coord + self.gen_speeds
                
                f_v = theano_function(tmp_coord + symb, [expr_vel], on_unused_input='ignore')
                f_v_proj = theano_function(tmp_coord + symb, [expr_proj_vel], on_unused_input='ignore')
            elif method == 'autowrap - f2py':
                f_p = autowrap(msubs(expr_pos, self.__map['simplyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'] +symb, backend = 'f2py')
                    
                f_part = autowrap(msubs(part_vel_expr, self.__map['simplyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'] +symb, backend = 'f2py')
                    
                f_v = autowrap(msubs(expr_vel, self.__map['simplyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-merged'] +symb, backend = 'f2py')
                    
                f_v_proj = autowrap(msubs(expr_proj_vel, self.__map['simplyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-merged'] +symb, backend = 'f2py')
                    
            else:
                f_p = autowrap(msubs(expr_pos, self.__map['dummyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords']+symb, backend = 'cython')
                    
                f_part = autowrap(msubs(part_vel_expr, self.__map['dummyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords']+symb, backend = 'cython')
                    
                f_v = autowrap(msubs(expr_vel, self.__map['dummyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-merged']+symb, backend = 'cython')
                    
                f_v_proj = autowrap(msubs(expr_proj_vel, self.__map['dummyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-merged']+symb, backend = 'cython')
            
            
            
            
            pos_list.append(f_p)
            func_list.append(f_part)
            vel_list.append(f_v)
            vel_projected_list.append(f_v_proj)
            
        self.contact_parametric['lambda'] = list(zip(point_id,pos_list,vel_list, func_list, vel_projected_list))
            
        
    
    
        
        self.lambdified = True
    
    def formEquationsOfMotion(self, method = None):
        """
        This method is used to form the equations of motion of the rover either by
        using Euler-Lagrange method or Kane's method.

        Returns
        -------
        None.

        """
        
        if method is not None:
            if method != 'autowrap - f2py' and method != 'autowrap - cython' and method != 'lambdify' and method != 'theano':
                raise ValueError('Invalid input argument')
        # print('Attempting to Lambdify the most functions as possible....'+'\n')
        # sleep(0.1)
        else:
            method = self.config.get('Simulation','lambdifier')
            if method != 'autowrap - f2py' and method != 'autowrap - cython' and method != 'lambdify' and method != 'theano':
                raise ValueError('Invalid input argument')
                
        save_dir = None
        if self.config.getboolean('Simulator','autowrap save files'):
            save_dir = self.config.get('Simulator','autowrap save dir')
        
        
        
        external_forces = self.gravity_forces + self.torques
        if self.config.getboolean('Model Description','action-reaction'):
            external_forces += self.__reaction_torque
        
        
        
        if self.method == 'Lagrange':
            L = 0
            for body in self.bodies:
                body.potential_energy = 0
                
                L += Lagrangian(self.frames[0],body)
                
            m = LagrangesMethod(L, self.gen_coord, forcelist = external_forces, frame=self.frames[0])
            print('Generating Equations of motion with Euler-Lagrange''s method''. This may take a while. Please wait....')
            
            
            m.form_lagranges_equations()
            print('Done. Equations of motions successfully computed')
            
            self.lagrange_method['method'] = m
            print('Attempting to lambdifying the Mass Matrix and the Forcing vector. Please wait...')
            
            self.lagrange_method['lambda func']['mass matrix'] = lambdify(self.gen_coord,\
                                                                          m.mass_matrix)
                
            self.lagrange_method['lambda func']['mass matrix full'] = lambdify(self.gen_coord,\
                                                                          m.mass_matrix_full)
                
            self.lagrange_method['lambda func']['forcing vector'] = lambdify(self.gen_coord + self.gen_speeds, \
                                                                             m.forcing)
                
            self.lagrange_method['lambda func']['forcing vector full'] = lambdify(self.gen_coord + self.gen_speeds, \
                                                                             m.forcing_full)
                
            print('Done.')
            
            
            
        else:
            bodies = self.bodies
            for i in range(len(bodies)):
                bodies[i].potential_energy = 0
                
                
            #self.torques = []
            
            
            self.kane_method['method'] = KanesMethod(self.frames[0], self.gen_coord, self.gen_speeds, self.eq_const)
            print('Generating Equations of motion with Kane''s method''. This may take a while. Please wait....')
            
            
            
            self.kane_method['method'].kanes_equations(bodies, loads = external_forces)
            print('Done. Equations of motions successfully computed')
            
            print('Attempting to lambdifying the Mass Matrix and the Forcing vector. Please wait...')
            #self.kane_method['lambda func']['mass matrix'] = lambdify(self.gen_coord, \
                                                                     #self.kane_method['method'].mass_matrix)
                                                                     
                                                                     
            expr = self.kane_method['method'].mass_matrix_full
            tmp_coord = self.gen_coord
            
            if method == 'lambdify':
                print('lambda')
                self.kane_method['lambda func']['mass matrix full'] = lambdify(tmp_coord, \
                                                                      expr, "numpy")
            elif method == 'theano':
                print('theano')
                self.kane_method['lambda func']['mass matrix full'] = theano_function(tmp_coord, [expr])
            elif method == 'autowrap - f2py':
                print('fortran')
                self.kane_method['lambda func']['mass matrix full'] = autowrap(msubs(expr, self.__map['simplyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs gen-coords'], backend = 'f2py')
            else:
                print('cython')
                self.kane_method['lambda func']['mass matrix full'] = autowrap(msubs(expr, self.__map['dummyfied']['mapping']['gen-coords']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs gen-coords'], backend = 'cython')
                        
                
            #self.kane_method['lambda func']['forcing vector'] = lambdify(self.gen_coord + self.gen_speeds, \
                                                                      #self.kane_method['method'].forcing)
            
            expr = self.kane_method['method'].forcing_full
            tmp_coord = self.gen_coord + self.gen_speeds
            
            if method == 'lambdify':
                print('lambda')
                self.kane_method['lambda func']['forcing vector full'] = lambdify(tmp_coord + self.__driving_torques, \
                                                                      expr, "numpy")
            elif method == 'theano':
                print('theano')
                self.kane_method['lambda func']['forcing vector full'] = theano_function(tmp_coord + self.__driving_torques, [expr])
            elif method == 'autowrap - f2py':
                print('fortran')
                self.kane_method['lambda func']['forcing vector full'] = autowrap(msubs(expr, self.__map['simplyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['simplyfied']['subs merged'] + self.__driving_torques, backend = 'f2py')
            else:
                print('cython')
                self.kane_method['lambda func']['forcing vector full'] = autowrap(msubs(expr, self.__map['dummyfied']['mapping']['gen-merged']),\
                             tempdir=save_dir, args=self.__map['dummyfied']['subs merged'] + self.__driving_torques, backend = 'cython')
                
            print('Done.')
                
                
            
            
                
            
            
            
    def setInitialConditions(self, q, q_d):
        """
        This method assign the initial conditions of the system

        Parameters
        ----------
        q : list of values for the generalized coordinates
        q_d : list of values for the generalized speeds

        Returns
        -------
        None.

        """
        
        self.q_initial_condition = q
        self.q_d_initial_condition = q_d
        
        self.current_gen_coord = np.array(q, dtype='float64')
        self.current_gen_speed = np.array(q_d, dtype='float64')
        
                
    def showRoverConfiguration(self, q = None, mode = 'plot'):
        
        if q is None:
            q = self.q_initial_condition
            
        ani_intialized = False
            
        leg_right = ['P_{wbr}',
                     'P_{end_leg_br}',
                     'P_{rbh2}',
                     'P_{rbh1}',
                     'P_{SA_r}',
                     'P_{rfh1}',
                     'P_{rfh2}',
                     'P_{end_leg_fr}',
                     'P_{wfr}',]
        
        
        leg_left = ['P_{wbl}',
                     'P_{end_leg_bl}',
                     'P_{lbh_2}',
                     'P_{h_l_b1}',
                     'P_{SA_l}',
                     'P_{h_l_f1}',
                     'P_{h_lf_2}',
                     'P_{end_leg_fl}',
                     'P_{wfl}',]
        
        
        if mode == 'plot':
            self.__fig = plt.figure()
        
        if mode == 'animation' and self.__fig is None:
            self.__fig = plt.figure()
            ani_intialized = True
            
            


        leg_right_to_plot = np.empty((0,3))
        leg_left_to_plot = np.empty((0,3))
        
        for pt_name in leg_right:
            for i in range(len(self.points)):
                if self.points[i].name == pt_name:
                    break
            
            leg_right_to_plot = np.vstack((leg_right_to_plot, self.lambdified_points[i](*q)))
        
        
        for pt_name in leg_left:
            for i in range(len(self.points)):
                if self.points[i].name == pt_name:
                    break
            
            leg_left_to_plot = np.vstack((leg_left_to_plot, self.lambdified_points[i](*q)))
            
        
        
        rod = np.empty((0,3))
        rod_names = ['P_{SA_l}','P_{SA_r}']
        for pt_name in rod_names:
            for i in range(len(self.points)):
                if self.points[i].name == pt_name:
                    break
            
            rod = np.vstack((rod, self.lambdified_points[i](*q)))
        
        
        
        
        self.__ax = self.__fig.add_subplot(111,projection='3d')
        self.__ax.set_xlabel('x [m]')
        self.__ax.set_ylabel('y [m]')
        self.__ax.set_zlabel('z [m]')
        self.__ax.set_ylim(-2,2)
        self.__ax.set_xlim(-2,2)
        self.__ax.set_zlim(-2,2)
        self.__ax.set_title('Rover initial configuration')
        self.__ax.view_init(0, 180)
        self.__ax.dist = 7
        self.__ax.plot(leg_right_to_plot[:,0],leg_right_to_plot[:,1],leg_right_to_plot[:,2],'-or')
        self.__ax.plot(leg_left_to_plot[:,0],leg_left_to_plot[:,1],leg_left_to_plot[:,2],'-or')
        self.__ax.plot(rod[:,0], rod[:,1], rod[:,2],'-g')
        
        theta_points = np.linspace(0,2*np.pi)

        wheel_to_plot = np.zeros((theta_points.shape[0],3), dtype='float64')
        
        wheel_radius_test = self.config.getfloat('Model Description', 'wheel_radius')
        
        wheel_to_plot[:,0] = 0
        wheel_to_plot[:,1] = wheel_radius_test*np.cos(theta_points)
        wheel_to_plot[:,2] = wheel_radius_test*np.sin(theta_points)
        
        
        frame_names = ['FRS_f','BRS_f', 'FLS_f', 'BLS_f']
        wheel_centre_pt_names = ['P_{wfr}', 'P_{wbr}', 'P_{wfl}', 'P_{wbl}']
        
        
        iterable_var = list(zip(frame_names, wheel_centre_pt_names))
        
        for i in range(len(iterable_var)):
            for n, frame in enumerate(self.frames):
                if frame.name == iterable_var[i][0]:
                    break
            
            rot_matrix = self.lambdified_frames[n](*q)
            RW_rotated = np.dot(rot_matrix.T, wheel_to_plot.T)
            
            for n, point in enumerate(self.points):
                if point.name == iterable_var[i][1]:
                    break
            
            translation = self.lambdified_points[n](*q)
            
            RW_to_plot = RW_rotated.T + translation
            self.__ax.plot(RW_to_plot[:,0],RW_to_plot[:,1],RW_to_plot[:,2],'-b')
            
                

        
        
        
        if self.__ground['coeffs']:
            rover_pos = self.lambdified_points[1](*q)
            X,Y = np.meshgrid([rover_pos[0,0]-2, rover_pos[0,0]+2],[rover_pos[0,1]-2, rover_pos[0,1]+2])
            Z = self.__ground['coeffs']['a']*X + self.__ground['coeffs']['b']*Y + self.__ground['coeffs']['c']
            self.__ax.plot_surface(X,Y,Z)
            
            
        if self.config.getboolean('Plotting', 'center mass'):
            pts_to_plot = np.empty((0,3))
            
            centre_names = ['P_{R_{CM}}',
                            'SA_{r,CM}',
                            'P_{rf_{I1_cm}}',
                            'P_{rf_{I2_{cm}}',
                            'P_{rbI1_{cm}}',
                            'P_{rbI2_{cm}}',
                            'SA_{l_{CM}}',
                            'P_{lfI1_{cm}}',
                            'P_{lfI2_{cm}}',
                            'P_{lbI1_{cm}}',
                            'P_{lbI2_{cm}}',
                            ]
            
            
            for n, point in enumerate(self.points):
                if point.name in centre_names:
                    pts_to_plot = np.vstack((pts_to_plot,self.lambdified_points[n](*q)))
            
            
            self.__ax.scatter(pts_to_plot[:,0], pts_to_plot[:,1], pts_to_plot[:,2], c='k')
            
            if ani_intialized:
                return True
        
        
        
        if mode == 'plot':
            plt.show()
            
            
        
        
        
        
    def saveSimulator(self, path_to_file = None): #Non funziona
        
        
        #Need to clear lambda data bacause can't be saved
        obj_copy = self.__copy__()
        obj_copy.lambdified_frames = list()
        obj_copy.lambdified_frames_ang_vel = list()
        obj_copy.lambdified_frames_ang_acc = list()
        
        
        obj_copy.lambdified_points = list()
        obj_copy.lambdified_points_vel = list()
        obj_copy.lambdified_points_acc = list()
        
        obj_copy.lambdified_partial_velocities = dict()
        obj_copy.kane_method['lambda func'] = dict()
        obj_copy.lagrange_method['lambda func'] = dict()
        
        return obj_copy
        
        if path_to_file is None:
            file_list = os.listdir("data/models/")
            filename = "model"
            
            i=0
            while True:
                if i==0:
                    incr = ''
                else:
                    incr = '_' + str(i)
                    
                new_filename = filename + incr +'.pickle'
                
                if new_filename in file_list:
                    i+=1
                else:
                    break
            
            with open("data/models/" + new_filename, "wb") as dill_file:
                dill.dump(obj_copy, dill_file)
                
        else:
            with open(path_to_file, "wb") as f:
                dill.dump(obj_copy, f)
            
            
            
    def createGound(self, a, b, c):
        """
        This method initialize the ground below the rover
        as an inclined plane
        
        z = a*x + b*y + c

        Parameters
        ----------
        a : constant
        b : constant
        c : constamt

        Returns
        -------
        None.

        """
        
        self.__ground['coeffs']['a'] = a
        self.__ground['coeffs']['b'] = b
        self.__ground['coeffs']['c'] = c
        
        
        
    def simulate(self, time_step = None, sim_time = None):
        """
        This method is used to simulate the dynamics of the rover. It uses the
        Euler first order integration scheme

        Parameters
        ----------
        time_step : float, optional
            Integration step time. The default is None.
        sim_time : float, optional
            Duration of the simulation. The default is None.

        Returns
        -------
        None.

        """
        
        if time_step is None:
            time_step = self.config.getfloat('Simulator','step size')
        
        if sim_time is None:
            sim_time = self.config.getfloat('Simulator','step size')
        
        
        
        time_sim = np.arange(0, sim_time+time_step, time_step)
        
        #We will use odeint from scipy
        integrator = self.config.get('Simulator', 'integrator')
        
        
        def right_hand_side(t, x):
            
            print(t)
            
            #Unpack input
            self.current_gen_coord = x[:len(self.gen_coord)]
            self.current_gen_speed = x[len(self.gen_coord):]
            
            self.__gen_coord_history.append(self.current_gen_coord.tolist())
            self.__gen_speed_history.append(self.current_gen_speed.tolist())
            self.__time_history.append(t)
            
            if self.method == 'Kane':
                Mass_full = self.kane_method['lambda func']['mass matrix full'](*self.current_gen_coord)
                forcing_full = self.kane_method['lambda func']['forcing vector full'](*np.hstack((self.current_gen_coord, self.current_gen_speed, self.__current_driving_torque)))
                
                _, F = self.__checkWheelContact2()
                
                forcing_full += np.vstack((np.zeros((len(self.gen_coord),1)), F))
                
                sol = np.linalg.inv(Mass_full).dot(forcing_full) 
                

                
                dy_dt = sol.T.squeeze()
                
            
            return dy_dt
        
        
        print('Solving ODE')
        print('DIOCANE')
        if integrator == 'odeint':
            
            stiff_order = self.config.get('Simulator','max stiff order')
            if stiff_order == 'auto':
                stiff_order = 0
            else:
                stiff_order = int(float(stiff_order))
                
                
            non_stiff_order = self.config.get('Simulator','max non-stiff order')
            
            if non_stiff_order == 'auto':
                non_stiff_order = 0
            else:
                non_stiff_order = int(float(non_stiff_order))
                
            min_time_step = self.config.get('Simulator','min step size')
            
            if min_time_step == 'auto':
                min_time_step = 0.0
            else:
                min_time_step = float(min_time_step)
                    
            
            
            self.ode_sol = odeint(right_hand_side,np.hstack((self.current_gen_coord, self.current_gen_speed)),\
                                  time_sim, tfirst = True, full_output = 1)#, mxords=stiff_order, mxordn=non_stiff_order, hmin=min_time_step)
            
        elif integrator == 'Euler-Explicit':
            for t in tqdm(time_sim, desc='Solving ODEs with Euler-Explicit Scheme'):
                dy_dt = right_hand_side(t,np.hstack((self.current_gen_coord, self.current_gen_speed)))
                
                #Integrate
                y_t = dy_dt*time_step
                
                #Unpack solution
                self.current_gen_coord = y_t[:len(self.gen_coord)]
                self.current_gen_speed = y_t[len(self.gen_coord):]
                
                
            
            
        else:
            self.ode_sol = solve_ivp(right_hand_side, [time_sim[0], time_sim[-1]],\
                                     np.hstack((self.current_gen_coord, self.current_gen_speed)), method = integrator.split('ipv - ')[1])
        
            
            
            
    
    def __checkWheelContact(self):
        mask = np.zeros((len(self.lambdified_partial_velocities),), dtype = bool)
        #deltas = np.zeros((len(self.lambdified_partial_velocities),))
        
        k = self.config.getfloat('Ground Properties', 'stiffness')
        c_max = self.config.getfloat('Ground Properties', 'damping')
        d = self.config.getfloat('Ground Properties', 'exp')
        p_max = self.config.getfloat('Ground Properties','max depth')
        
        
        F = np.zeros((len(self.gen_speeds),1))
        
        
        for n, tup in enumerate(self.lambdified_partial_velocities):
            contact_pt_pos = tup[2](*self.current_gen_coord)
            z_lim = self.__ground['coeffs']['a']*contact_pt_pos[0,0] + \
                self.__ground['coeffs']['b']*contact_pt_pos[0,1] + self.__ground['coeffs']['c']
            
            if contact_pt_pos[0,2] <= z_lim:
                #print('contact!')
                #There is contact
                mask[n] = True
                delta = np.abs(z_lim - contact_pt_pos[0,2])
                
                wheel_contact_vel = tup[3](*np.hstack((self.current_gen_coord, self.current_gen_speed))) #In inertial frame
                
                c = step5(delta, 0, 0, p_max, c_max)
                
                F_d_wheel = -c*wheel_contact_vel
                F_k_wheel = np.array([[0,0,k*delta**d]])
                
                partial_velocities = tup[4](*self.current_gen_coord)
                
                F += np.dot(partial_velocities, F_d_wheel.T) + np.dot(partial_velocities, F_k_wheel.T)
                
        # if np.any(mask):
        #     index_list = np.where(mask)
        
        
        return mask, F
    
    
    
    def __checkWheelContact2(self):
        mask = None
        k = self.config.getfloat('Ground Properties', 'stiffness')
        c_max = self.config.getfloat('Ground Properties', 'damping')
        d = self.config.getfloat('Ground Properties', 'exp')
        p_max = self.config.getfloat('Ground Properties','max depth')
        wheel_radius = self.config.getfloat('Model Description','wheel_radius')
        
        F = np.zeros((len(self.gen_speeds),1))
        
        
        frame_names = ['FRS_f','BRS_f', 'FLS_f', 'BLS_f']
        wheel_centre_pt_names = ['P_{wfr}', 'P_{wbr}', 'P_{wfl}', 'P_{wbl}']
            
            
        iterable_var = list(zip(frame_names, wheel_centre_pt_names))
        wheel_pts = self.getWheelPoints()
        
        for i in range(len(iterable_var)):
                for n, frame in enumerate(self.frames):
                    if frame.name == iterable_var[i][0]:
                        break
                
                rot_matrix = self.lambdified_frames[n](*self.current_gen_coord)
                RW_rotated = np.dot(rot_matrix.T, wheel_pts.T)
                
                for n, point in enumerate(self.points):
                    if point.name == iterable_var[i][1]:
                        break
                
                translation = self.lambdified_points[n](*self.current_gen_coord)
                
                RW_to_plot = RW_rotated.T + translation
                
                #z_interp = interpolate.griddata((self.__ground_new[:,0],self.__ground_new[:,1]), self.__ground_new[:,2], RW_to_plot[:,:2],method='cubic')
                
                z_interp = self.__ground_interpolator(RW_to_plot[:,1],RW_to_plot[:,0], grid=False)
                
                
                if np.any(RW_to_plot[:,-1] < z_interp):
                    #There is contact!
                    
                    #Detect groups of contact
                    mask = RW_to_plot[:,-1] < z_interp
                    
                    
                    regions = findContinguousContactRegions(mask)
                    
                    for region in regions:
                        x_cm = np.mean(RW_to_plot[region,0])
                        y_cm = np.mean(RW_to_plot[region,1])
                        z_cm = np.mean(RW_to_plot[region,2])
                        
                        v = np.array([x_cm,y_cm,z_cm])
                        v_t = v -translation
                        v_rot = rot_matrix.dot(v_t.T)
                        
                        rho = np.sqrt(v_rot[1]**2 + v_rot[2]**2)[0]
                        a = np.arctan(v_rot[1]/np.abs(v_rot[2]))[0]
                        
                        #print(rho)
                        #print(a)
                        
                        delta = wheel_radius - rho
                        
                        c = step5(delta, 0, 0, p_max, c_max)
                        
                        F_k_wheel = np.array([[0,0,k*delta**d]])
                        
                        sym = np.array([rho,a])
                        
                        partial = self.contact_parametric['lambda'][i][3](*np.hstack((self.current_gen_coord, sym)))
                        vel = self.contact_parametric['lambda'][i][2](*np.hstack((self.current_gen_coord, self.current_gen_speed, sym)))
                        F_d_wheel = -c*vel
                        
                        F_d_wheel[0,:2] = 0 #Check
                        
                        F += np.dot(partial, F_d_wheel.T) + np.dot(partial, F_k_wheel.T)
                        
                        #Add Friction component
                        
                        v_f = vel[0,:2]
                        
                        v_norm = np.sqrt(v_f[0]**2 + v_f[1]**2) #Vector norm
                        
                        mu = self.friction.computeFrictionCoeffiecient(v_norm)
                        #mu/=5
                        #mu=0.1
                        #print(mu)
                        #print(mu)
                        F_friction = -mu*(F_k_wheel[0,-1])# + F_d_wheel[0,-1]) #Modulus
                        #F_friction = -mu*10
                        F_fr_vec = np.hstack((F_friction*v_f/v_norm, [0]))
                        #print(F_fr_vec)
                        
                        F += np.dot(partial, np.array([F_fr_vec]).T)
                        
                        
                    
                
        return mask, F        
                
            
             
             
                
                
    def reset(self, reloadConfig = False):
        """
        Method used to reset the simulation

        Parameters
        ----------
        reloadConfig : bool, optional
            Parameter used to indcate the method wheter to reload the config file, since it may be changed.
            The default is False.

        Returns
        -------
        None.

        """
        #First thing to reset is the current gen coords and gen speeds
        
        self.current_gen_coord = np.array(self.q_initial_condition, dtype='float64')
        self.current_gen_speed = np.array(self.q_d_initial_condition, dtype='float64')
        
        self.__time_history = list()
        self.__gen_coord_history = list(list())
        self.__gen_speed_history = list(list())
        
        self.ode_sol = None
        
        if reloadConfig:
            self.loadConfig()
            
            
    
    def __createVariablesMapping(self):
        """
        This method initilize the mapping process which will map the dynamicsymbols used 
        in the ODEs with more simple symbols. Two mappings are created, and just one is 
        used depending on the code gen mode selected.

        Returns
        -------
        None.

        """
        n = str(len(self.gen_coord))
        
        dummy_var_gen_coord = [Dummy() for i in self.gen_coord]
        dummy_var_gen_speed = [Dummy() for i in self.gen_speeds]
        
        simple_var_gen_coord = list(symbols('o_sd:' + n))
        simple_var_gen_speed = list(symbols('p_sd:' + n))
        
        dummy_dict = dict()
        dummy_dict['original gen-coords'] = self.gen_coord
        dummy_dict['original gen-speeds'] = self.gen_speeds
        dummy_dict['subs gen-coords'] = dummy_var_gen_coord
        dummy_dict['subs gen-speeds'] = dummy_var_gen_speed  
        dummy_dict['original merged'] = self.gen_coord + self.gen_speeds
        dummy_dict['subs merged'] = dummy_var_gen_coord + dummy_var_gen_speed
        
        dummy_dict['mapping'] = dict()
        
        dummy_dict['mapping']['gen-coords'] = dict(zip(self.gen_coord, dummy_var_gen_coord))
        dummy_dict['mapping']['gen-speeds'] = dict(zip(self.gen_speeds, dummy_var_gen_speed))
        dummy_dict['mapping']['gen-merged'] = {**dummy_dict['mapping']['gen-coords'], **dummy_dict['mapping']['gen-speeds']}
        
        
        simple_dict = dict()
        simple_dict['original gen-coords'] = self.gen_coord
        simple_dict['original gen-speeds'] = self.gen_speeds
        simple_dict['original merged'] = self.gen_coord + self.gen_speeds
        simple_dict['subs gen-coords'] = simple_var_gen_coord
        simple_dict['subs gen-speeds'] = simple_var_gen_speed
        simple_dict['subs merged'] = simple_var_gen_coord + simple_var_gen_speed
        
        
        simple_dict['mapping'] = dict()
        
        simple_dict['mapping']['gen-coords'] = dict(zip(self.gen_coord, simple_var_gen_coord))
        simple_dict['mapping']['gen-speeds'] = dict(zip(self.gen_speeds, simple_var_gen_speed))
        simple_dict['mapping']['gen-merged'] = {**simple_dict['mapping']['gen-coords'], **simple_dict['mapping']['gen-speeds']}
        
        
        self.__map['dummyfied'] = dummy_dict
        self.__map['simplyfied'] = simple_dict
        
        
        
    def getMapping(self):
        """
        Method to be called in order to get the variables mapping dictionary

        Returns
        -------
        dict
            Dictionary of the mapped variables.

        """
        return self.__map
        
    
    def getCurrentState(self):
        tmp = np.array([self.current_gen_coord, self.current_gen_speed]).T
        return tmp
        
        
        
        
        
    def animateMotion(self, t = None, states = None):
        """
        This method is used in order to create a video animation of the siomulation
        conducted on the robot

        Parameters
        ----------
        t : nd-array, optional
            Time array of the simulation. If no argument is passed it will use the internal stored information. The default is None.
        states : nd-array, optional
            nd-array describing the state of the lander during the simulation. If no argument is passed it will use the internal stored information . The default is None.

        Raises
        ------
        RuntimeError
            Raise when no arguments are passed and the class can't find internal informations about the simulation.

        Returns
        -------
        None.
            

        """
        if t is None:
            try:
                t = self.ode_sol['t']
            except:
                    raise RuntimeError('Calling animateMotion without arguments uses class simulation stored data.')
                    
        if states is None:
            try:
                states = self.ode_sol['y']
            except:
                    raise RuntimeError('Calling animateMotion without arguments uses class simulation stored data.')
                    
        
        #Initilize figure objects
        fig = plt.figure()
        ax = Axes3D(fig)
        
        x_max = states[0,:].max() + 0.5
        x_min = states[0,:].min() - 0.5
        
        y_offset = 3*self.config.getfloat('Model Description','link_lenght') + self.config.getfloat('Model Description','wheel_radius') + .1

        y_max = states[1,:].max() + y_offset
        y_min = states[1,:].min() - y_offset
        
        z_max = states[2,:].max() + .2
        z_min = states[2,:].min() - 1
        
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        #ax.set_aspect('equal')
        ax.set_ylim3d(y_min,y_max)
        ax.set_xlim3d(x_min,x_max)
        ax.set_zlim3d(z_min,z_max)
        ax.set_title('Rover motion')
        ax.view_init(0, 180)
        ax.dist = 7
        
# =============================================================================
#         surf = None
#         if self.__ground['coeffs']:
#             X,Y = np.meshgrid([x_min, x_max],[y_min, y_max])
#             Z = self.__ground['coeffs']['a']*X + self.__ground['coeffs']['b']*Y + self.__ground['coeffs']['c']
#             
#             surf = ax.plot_surface(X,Y,Z, alpha=0.5)
# =============================================================================
        surf = None
        if self.__ground_new is not None:
            # x_min = self.__ground_new[:,0].min()
            # x_max = self.__ground_new[:,0].max()
            x_domain = np.linspace(x_min, x_max, 1000)
            
            # y_min = self.__ground_new[:,1].min()
            # y_max = self.__ground_new[:,1].max()
            
            y_domain = np.linspace(y_min, y_max, 1000)
            X,Y = np.meshgrid(x_domain, y_domain)
            
            Z = interpolate.griddata(self.__ground_new[:,:2], self.__ground_new[:,2], (X,Y),method='cubic')

            surf = ax.plot_surface(X,Y,Z, alpha=0.5)
            
            # dt = tri.Triangulation(self.__ground_new[:,0],self.__ground_new[:,1])
            # surf = ax.plot_trisurf(self.__ground_new[:,0], self.__ground_new[:,1], self.__ground_new[:,2], triangles = dt.triangles, alpha=0.1)
        
        
        time_text = ax.text2D(0.04,0.9,'', transform=ax.transAxes)
        
        line1 = ax.plot([],[],[],'-or')[0]
        line2 = ax.plot([],[],[],'-or')[0]
        line3 = ax.plot([],[],[],'-g')[0]
        line4 = ax.plot([],[],[],'-b')[0]
        line5 = ax.plot([],[],[],'-b')[0]
        line6 = ax.plot([],[],[],'-b')[0]
        line7 = ax.plot([],[],[],'-b')[0]
        
        
        
        
        
        
        def init():
            time_text.set_text('')
            line1.set_data(np.array([]), np.array([]))
            line1.set_3d_properties(np.array([]))
            
            line2.set_data(np.array([]), np.array([]))
            line2.set_3d_properties(np.array([]))
            
            line3.set_data(np.array([]), np.array([]))
            line3.set_3d_properties(np.array([]))
            
            line4.set_data(np.array([]), np.array([]))
            line4.set_3d_properties(np.array([]))
            
            line5.set_data(np.array([]), np.array([]))
            line5.set_3d_properties(np.array([]))
            
            line6.set_data(np.array([]), np.array([]))
            line6.set_3d_properties(np.array([]))
            
            line7.set_data(np.array([]), np.array([]))
            line7.set_3d_properties(np.array([]))
            
            
            if surf is not None:
                return time_text, line1, line2, line3, line4, line5, line6, line7, surf
            else:
                return time_text, line1, line2, line3, line4, line5, line6, line7
        
        
        def animate(i):
            time_text.set_text('t = {:2.2f} s'.format(t[i]))
            
            q = states[:len(self.gen_coord),i]
            
            
            leg_right = ['P_{wbr}',
                     'P_{end_leg_br}',
                     'P_{rbh2}',
                     'P_{rbh1}',
                     'P_{SA_r}',
                     'P_{rfh1}',
                     'P_{rfh2}',
                     'P_{end_leg_fr}',
                     'P_{wfr}',]
        
        
            leg_left = ['P_{wbl}',
                         'P_{end_leg_bl}',
                         'P_{lbh_2}',
                         'P_{h_l_b1}',
                         'P_{SA_l}',
                         'P_{h_l_f1}',
                         'P_{h_lf_2}',
                         'P_{end_leg_fl}',
                         'P_{wfl}',]
        
            
            


            leg_right_to_plot = np.empty((0,3))
            leg_left_to_plot = np.empty((0,3))
        
            for pt_name in leg_right:
                for i in range(len(self.points)):
                    if self.points[i].name == pt_name:
                        break
                
                leg_right_to_plot = np.vstack((leg_right_to_plot, self.lambdified_points[i](*q)))
        
        
            for pt_name in leg_left:
                for i in range(len(self.points)):
                    if self.points[i].name == pt_name:
                        break
                
                leg_left_to_plot = np.vstack((leg_left_to_plot, self.lambdified_points[i](*q)))
            
        
        
            rod = np.empty((0,3))
            rod_names = ['P_{SA_l}','P_{SA_r}']
            for pt_name in rod_names:
                for i in range(len(self.points)):
                    if self.points[i].name == pt_name:
                        break
                
                rod = np.vstack((rod, self.lambdified_points[i](*q)))
                
                
            line1.set_data(leg_right_to_plot[:,0],leg_right_to_plot[:,1])
            line1.set_3d_properties(leg_right_to_plot[:,2])
            
            
            line2.set_data(leg_left_to_plot[:,0],leg_left_to_plot[:,1])
            line2.set_3d_properties(leg_left_to_plot[:,2])
            
            line3.set_data(rod[:,0],rod[:,1])
            line3.set_3d_properties(rod[:,2])
            
            
            theta_points = np.linspace(0,2*np.pi)
    
            wheel_to_plot = np.zeros((theta_points.shape[0],3), dtype='float64')
            
            wheel_radius_test = self.config.getfloat('Model Description', 'wheel_radius')
            
            wheel_to_plot[:,0] = 0
            wheel_to_plot[:,1] = wheel_radius_test*np.cos(theta_points)
            wheel_to_plot[:,2] = wheel_radius_test*np.sin(theta_points)
            
            
            frame_names = ['FRS_f','BRS_f', 'FLS_f', 'BLS_f']
            wheel_centre_pt_names = ['P_{wfr}', 'P_{wbr}', 'P_{wfl}', 'P_{wbl}']
            
            
            iterable_var = list(zip(frame_names, wheel_centre_pt_names))
            
            for i in range(len(iterable_var)):
                for n, frame in enumerate(self.frames):
                    if frame.name == iterable_var[i][0]:
                        break
                
                rot_matrix = self.lambdified_frames[n](*q)
                RW_rotated = np.dot(rot_matrix.T, wheel_to_plot.T)
                
                for n, point in enumerate(self.points):
                    if point.name == iterable_var[i][1]:
                        break
                
                translation = self.lambdified_points[n](*q)
                
                RW_to_plot = RW_rotated.T + translation
                #ax.plot(RW_to_plot[:,0],RW_to_plot[:,1],RW_to_plot[:,2],'-b')
                if i==0:
                    line4.set_data(RW_to_plot[:,0], RW_to_plot[:,1])
                    line4.set_3d_properties(RW_to_plot[:,2])
                elif i==1:
                    line5.set_data(RW_to_plot[:,0], RW_to_plot[:,1])
                    line5.set_3d_properties(RW_to_plot[:,2])
                elif i==2:
                    line6.set_data(RW_to_plot[:,0], RW_to_plot[:,1])
                    line6.set_3d_properties(RW_to_plot[:,2])
                elif i==3:
                    line7.set_data(RW_to_plot[:,0], RW_to_plot[:,1])
                    line7.set_3d_properties(RW_to_plot[:,2])
                    
            
                    
                    
            if surf is not None:
                return time_text, line1, line2, line3, line4, line5, line6, line7, surf
            else:
                return time_text, line1, line2, line3, line4, line5, line6, line7
        
        
        anim = animation.FuncAnimation(fig, animate, frames=len(t), init_func=init, blit = True)
        ani_fps = int(self.config.getfloat('Plotting','animation fps'))
        
        if self.config.getboolean('Plotting','live animation'):
            plt.show()
        
        
        
        if self.config.getboolean('Plotting', 'save animation'):
            if self.config.get('Plotting','animation name') == 'auto':
                cwd = os.getcwd()
                
                tmp_str = datetime.datetime.now().strftime("%m-%d-%Y_%H-%M-%S")
                
                if not os.path.isdir(os.path.join(cwd,'data','video')):
                    os.makedirs(os.path.join(cwd, 'data', 'video'))
                    
                video_name = 'video_rover-'+tmp_str+'.mp4'
                anim.save(os.path.join(cwd,'data', 'video', video_name), fps=ani_fps)
            else:
                try:
                    anim.save(self.config.get('Plotting','animation name'), fps=ani_fps)
                except:
                    pass
                
        
        
        
        

    def exportToMatlabFunctions(self):
        print('TO - DO!')
        
        
    def loadModel(self, filename):
        print('TO DO')
        
    
    
    def changeTerrain(self, terrain_type = 'flat_ground'):
        data = self.loadGroundPoints(default_type = terrain_type)
        
        ground_pts = np.array((data['x'],data['y'], data['z'])).T
        self.__ground_new = ground_pts
        
        
        #self.__ground_interpolator = interpolate.SmoothBivariateSpline(ground_pts[0:-1:10,0],
        #                                                               ground_pts[0:-1:10,1],ground_pts[0:-1:10,2])
        
        X = ground_pts[:,0].reshape((data['rows'], data['columns']))
        Y = ground_pts[:,1].reshape((data['rows'], data['columns']))
        Z = ground_pts[:,2].reshape((data['rows'], data['columns']))
        
        
        self.__ground_interpolator = interpolate.RectBivariateSpline(Y[:,0], X[0,:], Z)
        
    def showTerrain(self, samples = 1000):
        x_min = self.__ground_new[:,0].min()
        x_max = self.__ground_new[:,0].max()
        
        y_min = self.__ground_new[:,1].min()
        y_max = self.__ground_new[:,1].max()
        
        x = np.linspace(x_min, x_max, samples)
        y = np.linspace(y_min, y_max, samples)
        
        X,Y = np.meshgrid(x,y)
        
        Z = interpolate.griddata(self.__ground_new[:,:2], self.__ground_new[:,2], (X,Y),method='cubic')
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot_surface(X,Y,Z, alpha = 0.5)
        plt.show()
        
        
    
    def loadGroundPoints(self, filename = None, default_type = 'flat_ground'):
        if filename is not None:
            try:
                with open(filename,'r') as f:
                    data = json.load(f)
            except:
                print('Couldnt load specified data file')
        else:
            cur_file_path = os.path.dirname(__file__)
            filename = os.path.join(cur_file_path,'data','ground',default_type+'.json')
            with open(filename,'r') as f:
                data = json.load(f)
        
        return data
        
        
    def __createWheel(self):
        radius = self.config.getfloat('Model Description','wheel_radius')
        angles = np.arange(-np.pi/2,np.pi/2,np.deg2rad(1))
        
        wheel = generateWheelsPoints(radius, angles)
        wheel[:,-1] = -wheel[:,-1]
        
        self.__wheel_points = wheel
        
    
    def getWheelPoints(self):
        return self.__wheel_points
    
    def getGroundPoints(self):
        return self.__ground_new
    
    
    def setSteerTorques(self, t_fr, t_br, t_fl, t_bl):
        t = [t_fr, t_br, t_fl, t_bl]
        self.__current_driving_torque[:4] = t
        
    def setDrivingTorque(self, t_fr, t_br, t_fl, t_bl):
        t = [t_fr, t_br, t_fl, t_bl]
        self.__current_driving_torque[4:] = t
        
    def getSteerTorques(self):
        return self.__current_driving_torque[:4]
    
    def getDriveTorques(self):
        return self.__current_driving_torque[4:]
    
    
    def loadFrictionModel(self):
        
        self.friction = Friction()
        g = 6*[0]
        if self.config.getboolean('Friction Properties','load_model'):
            g[0] = self.config.getfloat('Friction Properties','g1')
            g[1] = self.config.getfloat('Friction Properties','g2')
            g[2] = self.config.getfloat('Friction Properties','g3')
            g[3] = self.config.getfloat('Friction Properties','g4')
            g[4] = self.config.getfloat('Friction Properties','g5')
            g[5] = self.config.getfloat('Friction Properties','g6')
            
        self.friction.setCoeffiecients(*g)
        
        
    def getFrictionModel(self):
        return self.friction
    
    
    
        
        
        
        
    
    
        
        
        
        
        
    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result
    
    # def __deepcopy__(self, memo):
    #     cls = self.__class__
    #     result = cls.__new__(cls)
    #     memo[id(self)] = result
    #     for k, v in self.__dict__.items():
    #         setattr(result, k, copy.deepcopy(v, memo))
    #     return result
        
        





if __name__ == '__main__':
    sim = RoverSimulator()
    sim.initialize()
    #sim.lambdifyAll('lambdify')
    #print(os.path.dirname(__file__))
    # sim.formEquationsOfMotion('autowrap - f2py')
    # gen_coord = 15*[0]
    # gen_coord[:3] = 3*[0]
    # gen_coord[2] = 1
    # gen_coord[0] = 0
    # sim.createGound(0.0,-0.2,0)
    # sim.setInitialConditions(gen_coord,15*[0])
    # sim.simulate(0.02,3)
    # sim.animateMotion()