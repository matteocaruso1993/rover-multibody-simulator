# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:05:03 2021

@author: Matteo
"""
import numpy as np
import matplotlib.pyplot as plt
from rover_multibody_simulator.utilities import postprocess
import os
from matplotlib import rc
import pandas as pd


rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('font')

if __name__ == '__main__':
    
    sim_current = sim
    sim_folder = '7s_1'
# =============================================================================
#     sim_folder = '1ms_1'
#     sim_folder = '2ms_1'
#    sim_folder = '1ms_new'
#    sim_folder = '2ms_new'
# =============================================================================
    
    
    
    sim_current = sim03_new
# =============================================================================
#    sim_current = sim1ms_v1
#    sim_current = sim2ms_new
# =============================================================================
    
    cm = 1/2.54
    
    show_yaw_and_x_speed = False
    save_for_tikz = True
    single_plot = True
    fig_title = False
    print_errors = True
    tikz_path = r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Paper\tikz\data\Drive\Ramp'
    folder_tikz = 'test'
    
    
    wheel_radius = sim_current.config.getfloat('Model Description','wheel_radius')
    rover_speed = sim_current.q_d_initial_condition[-1]*wheel_radius
    data = np.genfromtxt(os.path.join(r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\Adams Comparison\Drive Test\Ramp\export', sim_folder,'y_pin_and_cm.dat'), skip_header=9)
    
    data1 = data - data[0,:]
    datarpy = np.genfromtxt(os.path.join(r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\Adams Comparison\Drive Test\Ramp\export', sim_folder,'RPY.dat'), skip_header=9)
    dataspeed = np.genfromtxt(os.path.join(r'C:\Users\Matteo\Google Drive\Shared\RoverMatlab\Experiment\Adams Comparison\Drive Test\Ramp\export', sim_folder,'payload_speed.dat'), skip_header=9)
    
    
    fig = plt.figure()
    if not single_plot:
        ax = fig.add_subplot(221)
    else:
        fig.set_size_inches(5*cm,4*cm)
        ax = fig.add_subplot(111)
        
    colors = ['b','r','g']
    ax.plot(data1[:,0],data1[:,1], '-' + colors[0] ,label='Adams-LP')
    ax.plot(data1[:,0],data1[:,2], '-' + colors[1], label='Adams-RP')
    ax.plot(data1[:,0],data1[:,3], '-' + colors[2],label='Adams-CM') #OK!
    
    ax.plot(sim_current.ode_sol['t'],(sim_current.ode_sol['y'][2,:]-sim_current.ode_sol['y'][2,0])*1000,'--' + colors[0],label='Python-CM')
    ax.plot(sim_current.ode_sol['t'],(postprocess.extractPointData(sim_current, 'P_{SA_r}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000,'--' + colors[1],label='Python-RP')
    ax.plot(sim_current.ode_sol['t'],(postprocess.extractPointData(sim_current, 'P_{SA_l}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000,'--' + colors[2],label='Python-LP')
    
    ax.set_xlabel('$t \ [s]$',fontsize=20)
    ax.set_ylabel('$z \ [mm]$',fontsize=20)
    if fig_title:
        ax.set_title('Adams-Python rover points height trend',fontsize=20)
    
    ax.legend()
    
    
    if not single_plot:
        ax1 = fig.add_subplot(222)
    else:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        
        
    if show_yaw_and_x_speed:
        ax1.plot(datarpy[:,0], -np.deg2rad(datarpy[:,2]), label='Adams-Y')
    ax1.plot(datarpy[:,0], -np.deg2rad(datarpy[:,1]), label='Adams-P')
    
    ax1.plot(datarpy[:,0], -np.deg2rad(datarpy[:,3]), label='Adams-R')
    
    if show_yaw_and_x_speed:
        ax1.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'theta', 'displacement'),'--', label ='Python-Y')
    ax1.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'chi', 'displacement'), '--',label ='Python-P')
    ax1.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'psi', 'displacement'), '--',label ='Python-R')
    ax1.set_xlabel('$t \ [s]$',fontsize=20)
    ax1.set_ylabel('$RPY \ [rad]$',fontsize=20)
    ax1.set_title('Adams-Python rover RPY angles trend',fontsize=20)
    ax1.legend()
    
    
    if not single_plot:
        ax2 = fig.add_subplot(223)
    else:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        
        
    ax2.plot(dataspeed[:,0], -dataspeed[:,1]/1000, label='Adams $v_y$')
    ax2.plot(dataspeed[:,0], dataspeed[:,2]/1000, label='Adams $v_z$')
    if show_yaw_and_x_speed:
        ax2.plot(dataspeed[:,0], -dataspeed[:,3]/1000, label='Adams-v_x')
    ax2.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'y', 'speed'),'--', label ='Python $v_y$')
    ax2.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'z', 'speed'), '--',label ='Python $v_z$')
    if show_yaw_and_x_speed:
        ax2.plot(sim_current.ode_sol['t'], postprocess.extractGeneralizedCoord(sim_current, 'x', 'speed'), '--',label ='Python $v_x$')
    
    ax2.set_xlabel('$t \ [s]$',fontsize=20)
    ax2.set_ylabel('$v \ [m/s]$',fontsize=20)
    ax2.set_title('Adams-Python rover centre mass speed trend',fontsize=20)
    
    ax2.legend()
    
    if not single_plot:
        ax3 = fig.add_subplot(224)
    else:
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
    
    names = ['P_{wfr}','P_{wbr}','P_{wfl}','P_{wbl}']
    
    for name in names:
        s,a = postprocess.computeWheelSlipRatio(sim_current, name)
        
        ax3.plot(sim_current.ode_sol['t'], s, label = name.split('P_{')[1].split('}')[0])
        
    ax3.set_xlabel('$t \ [s]$',fontsize=20)
    ax3.set_ylabel('$s$',fontsize=20)
    ax3.legend()
    
    
    
    
    
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(direction="in",right=True, top=True)
    ax.tick_params(which='major',length = 7)
    
    ax1.tick_params(axis='x', labelsize=20)
    ax1.tick_params(axis='y', labelsize=20)
    ax1.tick_params(direction="in",right=True, top=True)
    ax1.tick_params(which='major',length = 7)
    
    ax2.tick_params(axis='x', labelsize=20)
    ax2.tick_params(axis='y', labelsize=20)
    ax2.tick_params(direction="in",right=True, top=True)
    ax2.tick_params(which='major',length = 7)
    
    ax3.tick_params(axis='x', labelsize=20)
    ax3.tick_params(axis='y', labelsize=20)
    ax3.tick_params(direction="in",right=True, top=True)
    ax3.tick_params(which='major',length = 7)
    plt.tight_layout()
    
    
    if print_errors:
        #right_pin
        print('=========================================')
        print('Drive Test. Rover speed:\t%4.4f  [m/s]' % (rover_speed))
        print('=========================================')
        tmp = np.zeros((sim_current.ode_sol['t'].shape[0], 4))
        tmp[:,0] = sim_current.ode_sol['t']
        tmp[:,1] = (sim_current.ode_sol['y'][2,:]-sim_current.ode_sol['y'][2,0])*1000
        tmp[:,2] = (postprocess.extractPointData(sim_current, 'P_{SA_r}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000
        tmp[:,3] = (postprocess.extractPointData(sim_current, 'P_{SA_l}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000
        
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(tmp[:,[0,2]],data1[:,[0,2]], debug=False)
        string_to_print = 'Right Pin:\n---------------\nRMS:\t%4.4f; MSE:\t%4.4f; R2:\t%4.4f\n------------------------------' % (rms,mse,r2)
        print(string_to_print)
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(tmp[:,[0,3]],data1[:,[0,1]], debug=False)
        string_to_print = 'Left Pin:\n---------------\nRMS:\t%4.4f; MSE:\t%4.4f; R2:\t%4.4f\n------------------------------' % (rms,mse,r2)
        print(string_to_print)
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(tmp[:,[0,1]],data1[:,[0,3]], debug=False)
        string_to_print = 'Centre of Mass:\n---------------\nRMS:\t%4.4f; MSE:\t%4.4f; R2:\t%4.4f\n------------------------------' % (rms,mse,r2)
        print(string_to_print)
        
        
        tmp[:,0] = sim_current.ode_sol['t']
        tmp[:,1] = postprocess.extractGeneralizedCoord(sim_current, 'theta', 'displacement')
        tmp[:,2] = postprocess.extractGeneralizedCoord(sim_current, 'chi', 'displacement')
        tmp[:,3] = postprocess.extractGeneralizedCoord(sim_current, 'psi', 'displacement')
        
        tmp1 = np.zeros((datarpy.shape))
        tmp1[:,0] = datarpy[:,0]
        tmp1[:,1] = -np.deg2rad(datarpy[:,1])
        tmp1[:,3] = -np.deg2rad(datarpy[:,3])
        
        
    
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(tmp[:,[0,2]],tmp1[:,[0,1]], debug=False)
        string_to_print = 'Pitch:\n---------------\nRMS:\t%4.4f; MSE:\t%4.4f; R2:\t%4.4f\n------------------------------' % (rms,mse,r2)
        print(string_to_print)
        
        rms, mse, r2 = postprocess.computeErrorBetweenCurve(tmp[:,[0,3]],tmp1[:,[0,3]], debug=False)
        string_to_print = 'Roll:\n---------------\nRMS:\t%4.4f; MSE:\t%4.4f; R2:\t%4.4f\n------------------------------' % (rms,mse,r2)
        print(string_to_print)
        
        
        
    if save_for_tikz:
        if not os.path.exists(os.path.join(tikz_path,'Python',folder_tikz)):
            os.makedirs(os.path.join(tikz_path,'Python',folder_tikz))
            
        if not os.path.exists(os.path.join(tikz_path,'Adams',folder_tikz)):
            os.makedirs(os.path.join(tikz_path,'Adams',folder_tikz))
            
        
        data1_pd = pd.DataFrame(data1, columns=['Time', 'zl','zr','zcm'])
        data1_pd.to_csv(os.path.join(tikz_path,'Adams',folder_tikz,'z_pins.csv'))
        
        tmp = np.zeros((sim_current.ode_sol['t'].shape[0], 4))
        tmp[:,0] = sim_current.ode_sol['t']
        tmp[:,1] = (sim_current.ode_sol['y'][2,:]-sim_current.ode_sol['y'][2,0])*1000
        tmp[:,2] = (postprocess.extractPointData(sim_current, 'P_{SA_r}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000
        tmp[:,3] = (postprocess.extractPointData(sim_current, 'P_{SA_l}', 'displacement')[:,2]-sim_current.ode_sol['y'][2,0])*1000
        
        tmp1 = pd.DataFrame(tmp,columns = ['Time', 'zcm', 'zr','zl'])
        tmp1.to_csv(os.path.join(tikz_path,'Python',folder_tikz,'z_pins.csv'))
        
        tmp = np.vstack((datarpy[:,0], -np.deg2rad(datarpy[:,2]), -np.deg2rad(datarpy[:,1]), -np.deg2rad(datarpy[:,3])))
        #print(tmp.shape)
        tmp_pd = pd.DataFrame(tmp.T, columns=['Time','Yaw','Pitch','Roll'])
        tmp_pd.to_csv(os.path.join(tikz_path,'Adams',folder_tikz,'RPY.csv'))
        
        
        tmp = np.zeros((sim_current.ode_sol['t'].shape[0], 4))
        
        
        tmp[:,0] = sim_current.ode_sol['t']
        tmp[:,1] = postprocess.extractGeneralizedCoord(sim_current, 'theta', 'displacement')
        tmp[:,2] = postprocess.extractGeneralizedCoord(sim_current, 'chi', 'displacement')
        tmp[:,3] = postprocess.extractGeneralizedCoord(sim_current, 'psi', 'displacement')
        
        tmp_pd = pd.DataFrame(tmp, columns = ['Time','Yaw','Pitch','Roll'])
        tmp_pd.to_csv(os.path.join(tikz_path,'Python',folder_tikz,'RPY.csv'))
        
        