# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 16:50:02 2021

@author: Matteo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile, shutil, os
import scipy.interpolate as interpolator
import configparser
from . import functions
from copy import deepcopy
import sklearn.metrics
import pickle
from scipy.io import loadmat
from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('font')


def showCorrelationMatrix(in_dataframe,method='pearson'):
    corr_matrix = in_dataframe.corr(method=method)
    
    sns.heatmap(corr_matrix, annot=True)
    
    
def create_temporary_copy_csv(path, autoremove = True):
    file = os.path.split(path)
    #Create a copy of the file
    shutil.copy2(path,os.path.join(file[0], 'tmp.csv'))
    
    data = loadSummaryData(os.path.join(file[0], 'tmp.csv'))
    if autoremove:
        os.remove(os.path.join(file[0], 'tmp.csv'))
    

    return data


def loadSummaryData(path_to_file):
    try:
        data_frame = pd.read_csv(path_to_file)
    except:
        data_frame = None
        
    return data_frame


def computeErrorBetweenCurve(src, target, metric='least_squares', debug=True, scale=1, mean = True):
    interp = interpolator.interp1d(src[:,0],src[:,1], kind='cubic', bounds_error=False)
    target_x_valid = target[target[:,0] < np.max(src[:,0]),0]
    #error = np.sum((interp(target[:,0]) - target[:,1])**2)
    if mean:
        den = len(target_x_valid)
        #print(den)
    else:
        den = 1
    rms = sklearn.metrics.mean_squared_error(target[:len(target_x_valid),1]*scale, interp(target_x_valid)*scale, squared=False)
    mse = sklearn.metrics.mean_squared_error(target[:len(target_x_valid),1]*scale, interp(target_x_valid)*scale, squared=True)
    r2 = sklearn.metrics.r2_score(target[:len(target_x_valid),1]*scale, interp(target_x_valid)*scale)
    #error = np.sum((interp(target_x_valid)*scale - target[:len(target_x_valid),1]*scale)**2)/den
    if debug:
        plt.plot(target[:,0],target[:,1], label = 'Experimental')
        plt.plot(src[:,0], src[:,1], '--g',label='Numeric')
        plt.plot(target[:,0], interp(target[:,0]),':k', label='Numeric Interpolated')
        plt.legend()
        plt.text(2,.16, 'RMS: ' + '%2.4f'%(mse))
    return rms, mse, r2
    

def preprocessExperimentalData(data, path_to_config = None, save_file=False, save_location = None):
    #data must be a dictionary
    data1 = deepcopy(data)
    #Try to load from default position
    config = configparser.ConfigParser()
    current_folder = os.path.dirname(functions.__file__)
    if path_to_config is None:
        filename = os.path.join(current_folder,'..','data','config','config_experimental.ini')
    else:
        filename = path_to_config
        
    config.read(filename)
    
    for section in config.sections():
        if section == 'General':
            pass
        else:
            pin_name = section[:-2] #key to access data dictionary
            axe = section[-1].lower() #Value for dict
            new_axe = axe + '_proc'
            data_extracted = data1[pin_name][axe].to_numpy() #It's a table
            header = data1[pin_name][axe].columns.to_list()
            data_new = data_extracted[config.getint(section, 'n_initial_data_to_skip'):-config.getint(section, 'n_final_data_to_skip'),:]
            
            data_new[:,0] -= config.getfloat(section, 'time_offset')
            data_new[:,1:] *= config.getfloat('General','scale')/1000
            zero = data_new[0,1]
            data_new[:,1:] += config.getfloat(section,'y_offset') - zero
            
            
            new_database = pd.DataFrame(data_new, columns=header)
            data1[pin_name][new_axe] = new_database
            m = np.where(data_new[:,0] >= 0)[0]
            data_x = data_new[m,:]
            data1[pin_name][new_axe+'_valid'] = pd.DataFrame(data_x, columns=header)
            
            v = list()
            v_smooth = list()
            t = list()
            for i in range(data1[pin_name][new_axe+'_valid'].shape[0]-1):
                t.append(data1[pin_name][new_axe+'_valid'].iloc[i,0])
                v.append((data1[pin_name][new_axe+'_valid'].iloc[i+1,1] - data1[pin_name][new_axe+'_valid'].iloc[i,1])/
                         (data1[pin_name][new_axe+'_valid'].iloc[i+1,0] - data1[pin_name][new_axe+'_valid'].iloc[i,0]))
                v_smooth.append((data1[pin_name][new_axe+'_valid'].iloc[i+1,2] - data1[pin_name][new_axe+'_valid'].iloc[i,2])/
                         (data1[pin_name][new_axe+'_valid'].iloc[i+1,0] - data1[pin_name][new_axe+'_valid'].iloc[i,0]))
            
            vel_data = np.hstack((np.array(t).reshape(-1,1), np.array(v).reshape(-1,1), np.array(v_smooth).reshape(-1,1)))
            data1[pin_name][new_axe+'_speed_valid'] = pd.DataFrame(vel_data, columns=header)
            
            
            if save_file:
                if save_location is None:
                    save_location = os.getcwd()
                
                with open(os.path.join(save_location,'Experimental_data_processed.pickle'),'wb') as f:
                    pickle.dump(data1, f)
                
            
    
    return data1
            
            
            
def loadSavedExperiment(filename):
    with open(filename,'rb') as f:
        data = pickle.load(f)
    
    return data


def saveTable(tab, filename):
    tab.to_csv(filename)
    
                       
            
        
    
    
    
def showDistributionData(in_dataset, var_to_show, n_bins, binwidth):
    if not var_to_show[1] == ' ':
        var_to_show1 = ' ' + var_to_show
    else:
        var_to_show1 = var_to_show
        
    sns.displot(in_dataset,x=var_to_show1, binwidth=binwidth, kde=True)
    
def showScatterData(in_dataset, vars_to_show):
    if not len(vars_to_show) == 2:
        raise ValueError('Input dimensions must be 2')
    vars_to_show1 = list()
    for var in vars_to_show:
        if not var[1] == ' ':
            vars_to_show1.append(' ' + var)
        else:
            vars_to_show1.append(var)
    plt.scatter(in_dataset[vars_to_show1[0]], in_dataset[vars_to_show1[1]])
        
    
    
def computeVerticalOffset(sim_data, exp_data, key, axe):
    y_exp = exp_data[key][axe + '_proc_valid']['Amplitude'].to_numpy()
    y_sim = sim_data.z.to_numpy()
    
    return -y_exp[0] + y_sim[0]
    
    #print('to.do')
    
    
def generateDoeVariationSummary(path_to_folder):
    valid_dir = [directory for directory in os.listdir(path_to_folder) if os.path.isdir(os.path.join(path_to_folder,directory))]
    for n, folder in enumerate(valid_dir):
        summary_to_open = os.path.join(path_to_folder, folder, 'Summary.csv')
        summ = create_temporary_copy_csv(summary_to_open)
        if n == 0:
            #Create header
            sum_keys = list(summ.keys())
            key_to_remove = ['#filename', ' config_name', ' date']
            for k_r in key_to_remove:
                sum_keys.remove(k_r)
            data = pd.DataFrame(columns=sum_keys)
        
        current_row = list()
        current_row.append(folder)
        for key in sum_keys:
            min_sum = np.min(summ[key])
            max_sum = np.max(summ[key])
            if min_sum == max_sum:
                current_row.append(max_sum)
            else:
                current_row.append([min_sum, max_sum])
        
        data = data.append(dict(zip(['DOE ID'] + sum_keys, current_row)), ignore_index = True)
        
    data.to_csv(os.path.join(path_to_folder, 'Variation_Summary.csv'))
    
    

def extractPointData(in_sim, point_name, data_type):
    '''
    

    Parameters
    ----------
    in_sim : Input rover simulator
    
    point_name : name of the point we want to extract data
    
    data_type : either displacement or speed

    Returns
    -------
    None.

    '''
    
    if not data_type == 'displacement' and not data_type == 'speed':
        raise ValueError('Invalid Input argument')
        
        
    
    
    point_names = list()
    for pt in in_sim.points:
        point_names.append(pt.name)
    
    if point_name not in point_names:
        raise ValueError('Input point name is not valid')
        
        
    out_data = np.empty((0,3), dtype=float)
    
    for n,pt in enumerate(in_sim.points):
        if pt.name == point_name:
            break
    
    
    
    for i in range(in_sim.ode_sol['y'].shape[1]):
        q = in_sim.ode_sol['y'][0:len(in_sim.gen_coord),i]
        q_d = in_sim.ode_sol['y'][len(in_sim.gen_coord):,i]
        
        
        
        par = in_sim._map_parameters['vars-to-sub']
        
        if data_type == 'displacement':
            out_data = np.vstack((out_data, in_sim.lambdified_points[n](*np.hstack((q,par)))))
        else:
            out_data = np.vstack((out_data, in_sim.lambdified_points_vel[n](*np.hstack((q,q_d,par)))))
    
    
    return out_data



def extractGeneralizedCoord(in_sim, gen_coord_name, data_type):
    '''
    

    Parameters
    ----------
    in_sim : TYPE
        DESCRIPTION.
    gen_coord_name : TYPE
        DESCRIPTION.
    data_type : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    '''
    gen_coord_names = list()
    for coord in in_sim.gen_coord:
        gen_coord_names.append(coord.name)
        
    if gen_coord_name not in gen_coord_names:
        raise ValueError('Invalid input generalized coordinate name')
        
    data = np.empty((0,),dtype=float)
    if not data_type == 'displacement' and not data_type == 'speed':
        raise ValueError('Invalid input data type argument')
        
    
    if data_type == 'displacement':
        offset = 0
    else:
        offset = len(in_sim.gen_coord)
        
    for n, coord in enumerate(in_sim.gen_coord):
        if coord.name == gen_coord_name:
            break
        
        
        
    for i in range(in_sim.ode_sol['y'].shape[1]):
        data = np.hstack((data, in_sim.ode_sol['y'][n+offset, i]))
        
    return data



def computeWheelSlipRatio(in_sim, wheel_name):
    #Wheel name must be the axle point name
    r = in_sim.config.getfloat('Model Description','wheel_radius')
    header = ['s','alpha']
    
    mapp = {'P_{wfr}':{'drive_variable': 'omega1', 'steer_variable': 'delta1'},
            'P_{wbr}':{'drive_variable': 'omega2', 'steer_variable': 'delta2'},
            'P_{wfl}':{'drive_variable': 'omega3', 'steer_variable': 'delta3'},
            'P_{wbl}':{'drive_variable': 'omega4', 'steer_variable': 'delta4'}}
    
    
    drive = mapp[wheel_name]['drive_variable']
    steer = mapp[wheel_name]['steer_variable']
    
    omega = extractGeneralizedCoord(in_sim,drive,'speed')
    angles = extractGeneralizedCoord(in_sim, steer,'displacement')
    
    axle_speed = extractPointData(in_sim, wheel_name,'speed')
    axle_speed_norm = np.linalg.norm(axle_speed, axis=1)
    
    vel_angle = np.arctan(axle_speed[:,1]/axle_speed[:,0])
    
    slip_L = np.zeros(omega.shape)
    alphas = np.zeros(omega.shape)
    
    for i in range(slip_L.shape[0]):
        if axle_speed_norm[i] > r*abs(omega[i]):
            slip_L[i] = (axle_speed_norm[i] - r*abs(omega[i]))/axle_speed_norm[i]
        else:
            slip_L[i] = -(axle_speed_norm[i] - r*abs(omega[i]))/axle_speed_norm[i]
            
        
        delta = vel_angle[i] - angles[i]
        v = axle_speed_norm[i]*np.sin(delta)
        alphas[i] = np.arctan(v/(r*abs(omega[i])))
        
        
    return slip_L, alphas




def loadSpringMeasures(path, plot = False, measure_name = None, save_location = None):
    names = ['br1','br2','fr1','fr2','bl1','bl2', 'fl1','fl2']
    
    data = dict()
    
    for name_ in names:
        data[name_] = dict()
        tmp = loadmat(os.path.join(path, name_+'.mat'))
        data[name_]['t'] = tmp['data'][:,0]
        data[name_]['y'] = tmp['data'][:,1]
    
    
    if plot and measure_name is not None:
        if measure_name not in names:
            raise ValueError('Measure name must be of the one listed!')
        
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(data[measure_name]['t'],data[measure_name]['y'])
        ax.set_xlabel('$t \ [s]$',fontsize=20)
        ax.set_ylabel('$F [N]$',fontsize=20)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.tick_params(direction="in",right=True, top=True)
        ax.tick_params(which='major',length = 7)
        plt.tight_layout()
        if save_location is not None:
            
            plt.savefig(save_location, dpi=600)
    
    return data
            
            
    
         
    
    
        
        
            
        
        
        
        
        
    
    
    
    
    
        
            
        
        
        
        
    
        
    
    
    
    

if __name__ == '__main__':
    path_to_filename = r'C:\Users\Matteo\Desktop\appuntamento.txt'
    print(path_to_filename)