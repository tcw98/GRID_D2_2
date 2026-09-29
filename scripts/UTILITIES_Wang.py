# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:37:28 2023

@author: Witty
"""
import io
import subprocess
import pandas as pd
import numpy as np
import pyvista as pv
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from tqdm import tqdm
from scipy.spatial import cKDTree
from scipy.stats import norm
from sklearn.preprocessing import QuantileTransformer
import time
from subprocess import Popen
import rasterio
import pickle
from copy import deepcopy
import json
import os
from scipy.optimize import curve_fit
import concurrent.futures
from shapely.geometry import Point,LineString, Polygon
from concurrent.futures import ProcessPoolExecutor
import shapely.affinity
from GSLIB_parfl_011 import (GSLIB_SGSIM, GSLIB_SISIM_LM, GSLIB_GAMV, GSLIB_GAM,
                             GSLIB_VARMAP, GSLIB_SISIM)
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pykrige.ok import OrdinaryKriging



def get_coords_from_tif(tif_file):
    with rasterio.open(tif_file) as src:
        band = src.read(1)  # Read the first band

        x_coords = src.bounds.left + np.arange(src.width) * src.res[0]+0.5*src.res[0]
        x_coords = np.tile(x_coords,src.height)
        y_coords = src.bounds.top - np.arange(src.height) * src.res[1]-0.5*src.res[1]
        y_coords = y_coords.repeat(src.width)
        z_coords = band.flatten()
    return np.c_[x_coords,y_coords,z_coords]

def list_files(directory):
    file_list = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            file_list.append(filename)
    return file_list


def get_dgm1(directory, xmn, xmx, ymn, ymx):
    files = list_files(directory)
    extent = []
    for tif in tqdm(files):
        with rasterio.open(str(directory) + str(tif)) as src:
            extent.append([tif,src.bounds.left,src.bounds.right,src.bounds.top,src.bounds.bottom])
    
    
    ex = np.array(extent)[:,1:].astype(float)
    
    files = np.array(extent)[:,0][(ex[:,0] < xmx) & (ex[:,1] > xmn) & (ex[:,2] > ymn) & (ex[:,3] < ymx)]
    
    coords = []
    for tif in tqdm(files):
        one_tif = get_coords_from_tif(str(directory) + str(tif))
        coords.append(one_tif)
    
    return np.vstack(coords)

def find_closest_data(points, coordinates, data):
    tree = cKDTree(coordinates)
    _, indices = tree.query(points)
    closest_data = data[indices]
    return closest_data

def get_grid_elev(directory,grid_params,origin=0,angle=0):
    
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    x = np.tile(np.arange(xmn,xmn+nx*xsiz,xsiz),ny)
    y = np.repeat(np.arange(ymn,ymn+ny*ysiz,ysiz),nx)
    x,y = rotate_point((x,y), origin, -angle)
    points = np.c_[x,y]
    xmn, xmx, ymn, ymx = min(x), max(x), min(y), max(y)
    dgm1 = get_dgm1(directory,xmn, xmx, ymn, ymx)
    coordinates = dgm1[:,:-1]
    
    data = dgm1[:,-1]
    grid_data = find_closest_data(points, coordinates, data)
    return np.c_[points,grid_data]

    
def len_wei_compositing(SURVEY, FROM_TO, compositing_interval, verticaloption):
    """Function to composite and desurvey drillhole paths

    Parameters
    ----------
    compositing_interval : float
        DESCRIPTION.
    IDS : array
        Array with holeID in survey table.
    x : array
        X coordinates of drilling locations.
    y : array
        Y coordinates of drilling locations.
    z : array
        Z coordinates of drilling locations.
    azm : array
        Azimut of drilling paths in degrees. Set to 0° if vertical.
    dip : array
        dip of drilling path in degrees. Set to 90° if vertical.
    verticaloption : boolean
        if true, all drillings are assumed vertical.
    IDF : TYPE
        Array with holeID in from-to-table df.
    soiltype : array
        Description of the soil.
    lowerdep : array
        lower depth of the layers.

    Returns
    -------
    points:
        composited drillholes
    """

    nbor = len(SURVEY)  # number of Boreholes
    IDS = np.array(SURVEY["ObjektIDObjektID"]).flatten()  # column with holeID in survey df
    x = np.array(SURVEY["xcoord"]).flatten()  # X-location (make sure to use one CRS)
    y = np.array(SURVEY["ycoord"]).flatten()  # Y-location (make sure to use one CRS)
    z = np.array(SURVEY["Ansatzhoeh"]).flatten()  # elevation
    azm = [0] * nbor  # SURVEY.iloc[:,9] # 0° azimuth
    dip = [90] * nbor  # SURVEY.iloc[:,8] # 90° vertical path
    IDF = np.array(FROM_TO["ObjektID"]).flatten()
    lowerdep = np.array(FROM_TO["Untergrenz"].astype(float)).flatten()
    soiltype = np.array(FROM_TO["DIN"]).flatten().astype(str)
    soiltype = np.char.replace(soiltype, '/', '*')
    faulty_H_Id = []
    nbor = len(x)

    for i in tqdm(np.arange(0, nbor, 1)):
        H_Id = IDS[i]  # column with holeID in survey df
        borFT_soiltype = np.array(soiltype[IDF == H_Id])
        borFT_lowerdep = np.array(lowerdep[IDF == H_Id])
        HBA = []
        if borFT_soiltype.size == 0:                                            # profile description missing?
            continue
        silt_option = False
        if ~silt_option:
            for DIN in borFT_soiltype:
                DIN_temp = DIN.split(sep=",")
                HBA_temp = 10  # Default value
                if 'G' in DIN_temp[0]: HBA_temp = 3
                if 'A' in DIN_temp[0]: HBA_temp = 10
                if 'S' in DIN_temp[0]: HBA_temp = 2
                if 'U' in DIN_temp[0]: HBA_temp = 1
                if 'T' in DIN_temp[0]: HBA_temp = 1
                if 'F' in DIN_temp[0]: HBA_temp = 1
                if 'Z' in DIN_temp[0]: HBA_temp = 4
                HBA.append(HBA_temp)
        if silt_option:
            for DIN in borFT_soiltype:
                DIN_temp = DIN.split(sep=",")
                HBA_temp = 10  # Default value
                if 'G' in DIN_temp[0]: HBA_temp = 3
                if 'A' in DIN_temp[0]: HBA_temp = 10
                if 'S' in DIN_temp[0]: 
                    if 't*' in DIN_temp: # more options !!!!!!
                        HBA_temp = 7
                    else:
                        HBA_temp = 2  
                if 'U' in DIN_temp[0]: 
                    if 's*' in DIN_temp: # more options !!!!!!
                        HBA_temp = 7
                    else:
                        HBA_temp = 1 
                if 'T' in DIN_temp[0]: 
                    if 's*' in DIN_temp: # more options !!!!!!
                        HBA_temp = 7
                    else:
                        HBA_temp = 1 
                if 'F' in DIN_temp[0]: HBA_temp = 1
                if 'Z' in DIN_temp[0]: HBA_temp = 4
                HBA.append(HBA_temp)

        DEG2RAD = np.pi / 180.0
        if verticaloption:
            azm = 0
            dip = 90
            razm = azm * DEG2RAD
            rdip = -dip * DEG2RAD
        else:
            razm = azm[i] * DEG2RAD
            rdip = -dip[i] * DEG2RAD
            
        xn = np.sin(razm) * np.cos(rdip)
        yn = np.cos(razm) * np.cos(rdip)
        zn = np.sin(rdip)
        
        xi = x[i]
        yi = y[i]
        zi = z[i]
        
        lenborFT = len(HBA)
        
        xp = xi + xn * borFT_lowerdep
        yp = yi + yn * borFT_lowerdep
        zp = zi + zn * borFT_lowerdep
    
        # classify in interval
        cmin = np.floor(
            zp[lenborFT - 1] / compositing_interval) * compositing_interval     # round interval bottom
        cmax = np.ceil(zi / compositing_interval) * compositing_interval       # round interval top
        zcmps = np.arange(cmin, cmax, compositing_interval)[::-1]               # make z - array
        xcmps = xi + (xn * (zcmps - min(zcmps)) + zi - max(zcmps))
        ycmps = yi + (yn * (zcmps - min(zcmps)) + zi - max(zcmps))
    
        # create composite
        HBA = np.array(HBA)
        borcmps = []
        zcmps0 = -zcmps + zcmps[0] + compositing_interval                       # cmps depth fromp 1. cmps top
        zp0 = -zp + zcmps[0] + compositing_interval                             # interface depth fromp 1. cmps top
        count = 0
        for j in range(len(zcmps0)):
            if zcmps0[j] < zp0[count]:                                          # same interval?
                borcmps.append(HBA[count])
            else:
                weights = []
                HBA_temp = []
                HBA_temp.append(HBA[count])                                     # first interval
                weights.append(zp0[count] - zcmps0[j - 1])                      # interface to top
                count += 1
                for sec in range(
                        sum((zp0 > zcmps0[j - 1]) & (zp0 < zcmps0[j]))):        # loop interfaces in composite
                    if count == len(HBA):                                    # final depth?
                        break
                    if zp0[count] - zcmps0[j] < 0:                              # not last interval?
                        weights.append(
                            zp0[count] - zcmps0[j - 1] - weights[-1])           # interface to interface
                        HBA_temp.append(HBA[count])  
                        count += 1
                    else:                                                       # last interval
                        weights.append(-zp0[count - 1] + zcmps0[j])             # interface to bottom
                        HBA_temp.append(HBA[count])  

                unique_categories, indices = np.unique(                         # aggregate same soil types
                    HBA_temp, return_inverse=True)
                aggregated_weights = np.zeros(unique_categories.shape)
                np.add.at(aggregated_weights, indices, weights)
                max_weight_category = unique_categories[
                    np.argmax(aggregated_weights)]
                borcmps.append(max_weight_category)                             # choose soil type with maximum length
        zcmps += compositing_interval / 2                                       # relocate to center
        if i == 0:
            points = np.c_[xcmps, ycmps, zcmps, borcmps]
        else:
            points1 = np.c_[xcmps, ycmps, zcmps, borcmps]
            points = np.append(points, points1, axis=0)
    return points


def assign_points_2_grid(points_cs, grid_params):
    """
    #######################################################################
    # Assign nearst data (see GSLIB-SISIM: sstrat = 1)
    ######################################################################

    Parameters
    ----------
    points_cs : TYPE
        DESCRIPTION.
    grid_params : TYPE
        DESCRIPTION.

    Returns
    -------
    final_array : TYPE
        DESCRIPTION.

    """
    x,y,z, categories = points_cs[:,0], points_cs[:,1], points_cs[:,2], points_cs[:,3]
    loc, mask = checkgrid(x, y, z, grid_params, XY_option= False)
    print(f"Warning: {len(loc) - sum(mask)} of {len(loc)} points not in grid!")
    valid_locs = loc[mask]
    points_cs= points_cs[mask]
    x,y,z, valid_categories = points_cs[:,0], points_cs[:,1], points_cs[:,2], points_cs[:,3]
    xsiz, ysiz = grid_params["xsiz"], grid_params["ysiz"]
    valid_dist = ((x % xsiz - xsiz/2)**2 + (y % ysiz - ysiz/2)**2)**0.5
    values, indices, counts = np.unique(valid_locs, return_index=True,
                                                    return_counts=True)      # Pick the unique voxel locations, and the first apparent locations in valid_locs
    unique_loc_values = values[counts==1]
    unique_loc_indices = indices[counts==1]
    duplicate_loc_values = values[counts > 1]
    nx,ny,nz = grid_params["nx"], grid_params["ny"], grid_params["nz"]
    final_array = np.full(nx * ny * nz, -99)                                    # Initialize the final array with a default value
    final_array[unique_loc_values] =  valid_categories[unique_loc_indices]      # Assign the unique value to the result array

    # Process each duplicate grid location
    for d_loc_val in tqdm(duplicate_loc_values):
        # Find all indices of points in this duplicate grid location
        indices = np.where(valid_locs == d_loc_val)[0]
        duplicate_categories = valid_categories[indices]
        duplicate_dist = valid_dist[indices]
        # Find the index of the nearest point
        nearest_index = indices[np.argmin(duplicate_dist)]
        # Assign the category value of the nearest point to the corresponding position in the final array
        final_array[d_loc_val] = valid_categories[nearest_index]
           
    return final_array


def fit_2_curve(ax, x_data, y_data, color, label= "None", p0 = [0.05,50,0.01]):
    # Determine if there is a nugget effect
    # if (y_data[1] - y_data[0]) > nugget_threshold:
    #     nugget_effect = y_data[1]
    # else:
    #     nugget_effect = 0
    
    # nugget_effect_adj = nugget_effect
    # # Subtract the nugget effect for fitting
    # y_data_adj = deepcopy(np.array(y_data))
    # y_data_adj[0] = nugget_effect_adj
    
    def exponential_model(x, vsill, vrange, vnugget):
        return vnugget + (vsill - vnugget) * (1 - np.exp(-x / (vrange / 3)))
    
    def sperical_model(x, vsill, vrange, vnugget):
        #vnugget = 0
        x1 = x[x <= vrange]
        y1 = vnugget + (vsill-vnugget) * ((x1/vrange)*(1.5-.5*((x1/vrange)**2))) #
        x2 = x[x > vrange]
        y2 = (x2/x2) * vnugget + (vsill-vnugget)
        y = np.append(y1,y2,axis=0)
        return y
    
    popt_ex, pcov_ex = curve_fit(exponential_model, x_data, y_data)  # Fit the model to the data
    perr_ex = np.sqrt(np.diag(pcov_ex))
    y_pred = exponential_model(x_data, *popt_ex)
    rmse_ex = np.sqrt(np.mean((y_pred - y_data)**2))  # Calculate the RMSE



    #popt_sp, pcov_sp = curve_fit(sperical_model, x_data, y_data, p0)
    popt_sp, pcov_sp = curve_fit(sperical_model, x_data, y_data, p0)# Fit the model to the data
    perr_sp = np.sqrt(np.diag(pcov_sp))
    y_pred = exponential_model(x_data, *popt_sp)
    rmse_sp = np.sqrt(np.mean((y_pred - y_data)**2))  # Calculate the RMSE
    if rmse_ex <= rmse_sp: 
        ax.plot(x_data, y_pred, color, zorder=-1, linestyle="dashed", label=label + "(Exponential)")
        print ("Exponential model is better!")
        return [*popt_ex, *perr_ex, rmse_ex, "Exponential"]
    else:
        ax.plot(x_data, y_pred, color, zorder=-1, linestyle="dotted", label=label+ "(Sperical)")
        print ("Sperical model is better!")
        return [*popt_sp, *perr_sp, rmse_sp, "Sperical"]
        
   


def fit_2_exponential(ax, x_data, y_data, color, label = None):
    def exponential_model(x, vsill, vrange):
        return vsill * (1 - np.exp(-x / (vrange / 3)))
    popt, pcov = curve_fit(exponential_model, x_data, y_data) # Fit the model to the data
    perr = np.sqrt(np.diag(pcov))
    y_pred = exponential_model(x_data, *popt)
    rmse = np.sqrt(np.mean((y_pred - y_data)**2)) # Calculate the RMSE
    # print([*popt, *perr,rmse])
    theovargm = MAKETHEORVARGM(0,*popt,2)#*t_var[i][k])
    ax.plot(theovargm['x'], theovargm['y'], color, zorder=-1,
                linestyle="dotted", label = label)
    return [*popt, *perr,rmse]


def summarize_variograms(variogram_t_north, variogram_t_east, variogram_t_vertical, print_option = True):
    # Determine major and minor directions
    major_direction = variogram_t_north if variogram_t_north[1] >= variogram_t_east[1] else variogram_t_east
    minor_direction = variogram_t_east if variogram_t_north[1] >= variogram_t_east[1] else variogram_t_north
    
    
    # Calculate average sill with emphasis on major direction: 2D 75%, vertical 25%
    average_sill = ((major_direction[0] + minor_direction[0]) / 2) * 0.9 + (variogram_t_vertical[0] * 0.1)
    average_nugget = ((major_direction[2] + minor_direction[2]) / 2) * 0.9 + (variogram_t_vertical[2] * 0.1)
    # Create the summary dictionary
    variogram_t = {
        "type": major_direction[-1],
        "ranges": (major_direction[1], minor_direction[1], variogram_t_vertical[1]),
        "angles": (0, 0, 0),
        "sill": average_sill,
        "nugget": average_nugget
    }
    
    if print_option == True:
        # Print the summary
        print("-------------------Summary--------------------")
        print(f"Type: {major_direction[-1]}")
        print(f"Nugget: {average_nugget}")
        print(f"Major direction range: {major_direction[1]}")
        print(f"Minor direction range: {minor_direction[1]}")
        print(f"Vertical direction range: {variogram_t_vertical[1]}")
        print(f"Sill: {average_sill}")
        print("-----------------------------------------------")
        print("!!!!The type numbers for representing Exponential/Sperical are different in GSLIB and in HPGL!!!!")
        print("-----------------------------------------------")

    
    return variogram_t

def summarize_variograms_2D(variogram_t_north, variogram_t_east, print_option = True):
    # Determine major and minor directions
    major_direction = variogram_t_north if variogram_t_north[1] >= variogram_t_east[1] else variogram_t_east
    minor_direction = variogram_t_east if variogram_t_north[1] >= variogram_t_east[1] else variogram_t_north
    
    
    # Calculate average sill with emphasis on major direction: 2D 75%, vertical 25%
    average_sill = ((major_direction[0] + minor_direction[0]) / 2 )
    average_nugget = ((major_direction[2] + minor_direction[2]) /2 )
    # Create the summary dictionary
    variogram_t = {
        "type": major_direction[-1],
        "ranges": (major_direction[1], minor_direction[1]),
        "angles": (0, 0, 0),
        "sill": average_sill,
        "nugget": average_nugget
    }
    
    if print_option == True:
        # Print the summary
        print("-------------------Summary--------------------")
        print(f"Type: {major_direction[-1]}")
        print(f"Nugget: {average_nugget}")
        print(f"Major direction range: {major_direction[1]}")
        print(f"Minor direction range: {minor_direction[1]}")
        print(f"Sill: {average_sill}")
        print("-----------------------------------------------")
        print("!!!!The type numbers for representing Exponential/Sperical are different in GSLIB and in HPGL!!!!")
        print("-----------------------------------------------")

    
    return variogram_t
    
#     params_array = np.array(list(grid.values()))
#     np.save('D:/Wang/U9/WD/save/grid_params.npy', params_array)
    
    
#     file_path_tohpgl = r'D:\Wang\HPGL\data\tohpgl.json'

#     with open(file_path_tohpgl, 'w') as f:
#         tohpgl = variogram
#         tohpgl['inc_path'] = inc_path
#         tohpgl['nsim'] = nsim
#         tohpgl['single_option'] = single_option
#         json.dump(tohpgl, f) 

        
#     file_path_se = r'D:\Wang\HPGL\data\tohpgl_se.json'
    

#     with open(file_path_se, 'w') as f:
#         json.dump(se, f)

#     # points_path = "D:/Wang/U9/WD/save/tohpgl_points.npy"
#     # np.save('D:/Wang/U9/WD/save/tohpgl_points.npy', points)

#     print(f"Starting script: {pyscript_path}")


#     # Record the start time
#     start_time = time.time()

#     python2_interpreter = r'C:\Python27\python.exe'
#     # Construct the command line command
#     cmd = [python2_interpreter, pyscript_path]
    
#     # Run the command using subprocess
#     process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
#     # Periodically check if the process is still running
#     while process.poll() is None:
#         print("Script is still running...")
#         time.sleep(30)  # Wait for 1 second before checking again
        
#     # Capture the command's output and errors after completion
#     stdout, stderr = process.communicate()

#     # Calculate the total running time
#     total_time = time.time() - start_time
#     print(f"Script finished. Total running time: {total_time:.2f} seconds")

#     # Print the output
#     print("STDOUT:\n", stdout.decode())
    
#     # Print errors, if any
#     if stderr:
#         print("STDERR:\n", stderr.decode())

import threading


def generate_scripts(origin_pyscript, numberofscripts, name):
    # Ensure the input file path is a string
    origin_pyscript = str(origin_pyscript)
    # Read the content of the original script
    with open(origin_pyscript, 'r') as file:
        original_content = file.read()

    for i in range(numberofscripts):  # Start the loop from 0
        # Generate new content for each new script, adding the n_value variable definition
        new_content = f"n_value = {i}\n" + original_content
        # Construct the full path and filename for the new script using the custom name
        new_script_path = f"D:\\Wang\\HPGL\\{name}_{i}.py"
        
        # Write the new content into the new script file
        with open(new_script_path, 'w') as new_file:
            new_file.write(new_content)

        print(f"Generated: {new_script_path}")

def enqueue_output(out, print_lock):
    """
    Continuously read lines from the output stream and print them.
    
    :param out: The output stream of the subprocess.
    :param print_lock: A threading lock to synchronize console output.
    """
    for line in iter(out.readline, ''):
        with print_lock:
            print(line.strip())
    out.close()

def run_hpgl(grid, variogram, se, nsim, pyscript_path, inc_path, single_option=False):
    """
    Runs a specified Python script using a Python 2 interpreter.
    
    :param grid: Grid parameters.
    :param variogram: Variogram parameters.
    :param se: Search ellipse parameters.
    :param nsim: Number of simulations to run.
    :param pyscript_path: The path to the Python 2 script to be run.
    :param inc_path: Path to the include file.
    :param single_option: Flag to force single thread execution.
    """
    
    # Save grid parameters as a numpy array
    params_array = np.array(list(grid.values()))
    np.save('D:/Wang/U9/WD/save/grid_params.npy', params_array)
    
    # Write search ellipse parameters to another JSON file
    file_path_se = r'D:\Wang\HPGL\data\tohpgl_se.json'
    with open(file_path_se, 'w') as f:
        json.dump(se, f)
        
    # Prepare HPGL parameters for JSON file
    tohpgl = variogram
    tohpgl['inc_path'] = inc_path
    tohpgl['nsim'] = nsim
    tohpgl['single_option'] = single_option
    file_path_tohpgl = r'D:\Wang\HPGL\data\tohpgl.json'
    with open(file_path_tohpgl, 'w') as f:
        json.dump(tohpgl, f)
        
    # np.save(f'D:/Wang/U9/WD/save/n_{n}.npy', np.array([n]))
    # Not working....
    
    print("Starting script" + str(pyscript_path))

    python2_interpreter = r'C:\Python27\python.exe'
    cmd = [python2_interpreter, pyscript_path]
    
    # Start the subprocess
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True)

    # Lock for synchronizing console output
    print_lock = threading.Lock()

    # Start threads to read subprocess output and errors
    stdout_thread = threading.Thread(target=enqueue_output, args=(process.stdout, print_lock))
    stderr_thread = threading.Thread(target=enqueue_output, args=(process.stderr, print_lock))
    stdout_thread.start()
    stderr_thread.start()

    # Wait for subprocess to complete and output threads to finish
    process.wait()
    stdout_thread.join()
    stderr_thread.join()

    print("Script finished" + str(pyscript_path))

            
        
def run_hpgl_parallel(grid, variogram, se, nsim, pyscript_path, inc_path, single_option=False):
    """
    Runs the 'run_hpgl' function in parallel for 5 different 'n' values.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for i in range(5):  # Assuming you want to run for n = 0 to 4
            pyscript_path_n = pyscript_path + str(i) + '.py'
            futures.append(executor.submit(run_hpgl, grid, variogram, se, nsim, pyscript_path_n, inc_path, single_option))
        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            future.result()  # This will raise exceptions if any occurred
            
            
def run_hpgl_lm(grid, variogram, se, nsim, pyscript_path, inc_path, single_option=False, use_correlogram = False):

    
    # Save grid parameters as a numpy array
    params_array = np.array(list(grid.values()))
    np.save('D:/Wang/U9/WD/save/grid_params.npy', params_array)
    
    # Write search ellipse parameters to another JSON file
    file_path_se = r'D:\Wang\HPGL\data\tohpgl_se.json'
    with open(file_path_se, 'w') as f:
        json.dump(se, f)
        
    # Prepare HPGL parameters for JSON file
    tohpgl = variogram
    tohpgl['inc_path'] = inc_path
    tohpgl['nsim'] = nsim
    tohpgl['single_option'] = single_option
    tohpgl['use_correlogram'] = use_correlogram
    file_path_tohpgl = r'D:\Wang\HPGL\data\tohpgl.json'
    with open(file_path_tohpgl, 'w') as f:
        json.dump(tohpgl, f)
        
    
    print("Starting script" + str(pyscript_path))

    python2_interpreter = r'C:\Python27\python.exe'
    cmd = [python2_interpreter,
           # f"nvalue = {i}\n", 
           # 'with open("D:/Wang/HPGL/HelloWho.py") as hello: exec(hello.read())',
           pyscript_path]
    
    # Start the subprocess
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True)

    # Lock for synchronizing console output
    print_lock = threading.Lock()

    # Start threads to read subprocess output and errors
    stdout_thread = threading.Thread(target=enqueue_output, args=(process.stdout, print_lock))
    stderr_thread = threading.Thread(target=enqueue_output, args=(process.stderr, print_lock))
    stdout_thread.start()
    stderr_thread.start()

    # Wait for subprocess to complete and output threads to finish
    process.wait()
    stdout_thread.join()
    stderr_thread.join()

    print("Script finished" + str(pyscript_path))

def run_SK(grid, data, pyscript_path):

     
     # Save grid parameters as a numpy array
     params_array = np.array(list(grid.values()))
     np.save('D:/Wang/U9/WD/save/grid_params.npy', params_array)
     
     INC = assign_points_2_grid(data, grid)
     
     file_name = f'D:/Wang/HPGL/data/SK.INC'
     np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')
     
     
     print("Starting script" + str(pyscript_path))
     python2_interpreter = r'C:\Python27\python.exe'
     cmd = [python2_interpreter,
            # f"nvalue = {i}\n", 
            # 'with open("D:/Wang/HPGL/HelloWho.py") as hello: exec(hello.read())',
            pyscript_path]
     
     # Start the subprocess
     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True)

     # Lock for synchronizing console output
     print_lock = threading.Lock()

     # Start threads to read subprocess output and errors
     stdout_thread = threading.Thread(target=enqueue_output, args=(process.stdout, print_lock))
     stderr_thread = threading.Thread(target=enqueue_output, args=(process.stderr, print_lock))
     stdout_thread.start()
     stderr_thread.start()

     # Wait for subprocess to complete and output threads to finish
     process.wait()
     stdout_thread.join()
     stderr_thread.join()

     print("Script finished" + str(pyscript_path))           

        
def run_hpgl_lm_parallel(grid, variogram, se, nsim, pyscript_path, inc_path, single_option=False, use_correlogram = False):
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for i in range(5):  # Assuming you want to run for n = 0 to 4
            pyscript_path_n = pyscript_path + str(i) + '.py'
            futures.append(executor.submit(run_hpgl_lm, grid, variogram, se, nsim, pyscript_path_n, inc_path, single_option, use_correlogram))
        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            future.result()  # This will raise exceptions if any occurred
# ----------------------------------------------------------------------------
def command_RUN_GSLIB(program,parfile):
    return print(str(program) + " " + str(parfile))
# ----------------------------------------------------------------------------
def domain_plus_buffer(grid_params,buffer):
    return [grid_params["xmn"]-buffer,
            grid_params["xmn"]+grid_params["xsiz"]*grid_params["nx"]+buffer,
            grid_params["ymn"]-buffer,
            grid_params["ymn"]+grid_params["ysiz"]*grid_params["ny"]+buffer]
# ----------------------------------------------------------------------------
def RUN_GSLIB(program,parfile):
    
    cmd = str(program) + " " + str(parfile)
    
    # returns output as byte string
    returned_output = subprocess.check_output(cmd)
    
    # using decode() function to convert byte string to string
    print(returned_output.decode("utf-8"))
def run_GSLIB_background(program,parfile):
    cmd = str(program) + " " + str(parfile)
    process = subprocess.Popen(cmd)
    return process
def check_GSLIB_background(process, continuously = False):
    process.poll()
    if continuously:
        while process.poll() is None:
            print("Process is still running...")
            time.sleep(1)  # Optional delay between checks
        
        print("Process has finished.")
#|||Wang:The purpose of these functions is to interact with GSLIB programs in Python, and provide different ways to run them (blocking or background execution), as well as monitoring of the process state while executing in the background.
# ----------------------------------------------------------------------------
def LFU_2_YOU(survey1_path,survey2_path,from_to_path,extent,GK4,rotation=False):
    xmin,xmax,ymin,ymax = extent 
    
    print("-- begin reading")
    drillhole = gpd.read_file(survey1_path)
    print("-- reading finished")
    print("-- transforming coordinates")
    drillhole_t = drillhole
    if GK4 == True: 
        #Bohrungen = Bohrungen.set_crs("epsg: 32632")    # UTM 32
        drillhole_t = drillhole.to_crs("epsg: 31468")     # Gauss Krueger 4
    if rotation != False:
        angle, origin = rotation
        def rotate_points_gdf(points, angle, origin):
            return Point(shapely.affinity.rotate(points, angle, origin=origin))
        drillhole['geometry'] = drillhole['geometry'].apply(lambda point: rotate_points_gdf(point, angle, origin))

    drillhole_t["xcoord"] = drillhole_t.geometry.x
    drillhole_t["ycoord"] = drillhole_t.geometry.y
    survey = drillhole_t[["ObjektID","xcoord","ycoord"]]
    print("-- join tables")
    SURVEY_1 = survey[(survey.xcoord > xmin) & (survey.xcoord < xmax) & (survey.ycoord > ymin) & (survey.ycoord < ymax)] 
    SURVEY_2 = pd.read_csv(survey2_path, sep=";", decimal=",")
    FROM_TO = pd.read_csv(from_to_path, sep=";", decimal=",").sort_values(by=["Untergrenz"]) 
    SURVEY = SURVEY_1.join(SURVEY_2, lsuffix = "ObjektID", rsuffix = "2", how = "inner")
    # Ansatzhöhe wird entfernt, da wo Fehler/nicht passt
    SURVEY = SURVEY.dropna(subset=["Ansatzhoeh"])
    SURVEY = SURVEY[SURVEY["Ansatzhoeh"]!=0]
    print("-- job finished")
    
    
    return SURVEY, FROM_TO
# ----------------------------------------------------------------------------

def qnd_compositing(SURVEY,FROM_TO,compositing_interval,vtkoption,
                    verticaloption,qbasisoption,qbasistol):
    """Function to composite and desurvey drillhole paths

    Parameters
    ----------
    compositing_interval : float
        DESCRIPTION.
    IDS : array
        Array with holeID in survey table.
    x : array
        X coordinates of drilling locations.
    y : array
        Y coordinates of drilling locations.
    z : array
        Z coordinates of drilling locations.
    azm : array
        Azimut of drilling paths in degrees. Set to 0° if vertical.
    dip : array
        dip of drilling path in degrees. Set to 90° if vertical.
    verticaloption : boolean
        if true, all drillings are assumed vertical.
    IDF : TYPE
        Array with holeID in from-to-table df.
    soiltype : array
        Description of the soil.
    lowerdep : array
        lower depth of the layers.
    vtkoption : boolean
        if True, a pyvista vtk object will be generated.
    qbasisoption : boolean
        get lowest gravel depth.

    Returns
    -------
    points:
        composited drillholes
    point_cloud:
        pyvista object with composited points
    line:
        pyvista object with drillhole paths

    """
    nbor = len(SURVEY) # number of Boreholes
    IDS = np.array(SURVEY["ObjektIDObjektID"]).flatten() # column with holeID in survey df
    #IDS = np.array(IDS.loc[:,~IDS.columns.duplicated()].copy()).flatten()
    x = np.array(SURVEY["xcoord"]).flatten()  # X-location (make sure to use one CRS)
    #x = x + 32000000
    y = np.array(SURVEY["ycoord"]).flatten()  # Y-location (make sure to use one CRS)
    z = np.array(SURVEY["Ansatzhoeh"]).flatten()  # elevation
    azm= [0]*nbor #SURVEY.iloc[:,9] # 0° azimuth
    dip= [90]*nbor #SURVEY.iloc[:,8] # 90° vertical path 
    IDF = np.array(FROM_TO["ObjektID"]).flatten()
    lowerdep = np.array(FROM_TO["Untergrenz"].astype(float)).flatten()
    soiltype = np.array(FROM_TO["DIN"]).flatten().astype(str)
    soiltype = np.char.replace(soiltype, '/', '*')
    
    faulty_H_Id = []
    nbor = len(x)
    if vtkoption:
        import pyvista as pv
        #check= 0
        HB_array = []
        r = 0
        p = np.array([[0,0,0]])
        cells = np.array([[0,0,0]])  
    ###########################################################################
    # if qbasisoption:
    #     qbasis = []
    # else:
    #     qbasis = False
    ###########################################################################
    qbasis = []

    for i in tqdm(np.arange(0,nbor,1)):
        #print(f"\r{round((i+1)/nbor*100,1)} %, borehole {i}", flush=True, end="")
        H_Id = IDS[i]                     # column with holeID in survey df
        #borFT = np.array(soiltype[IDF==H_Id])    # column with holeID in FROM_TO df
        borFT_soiltype = np.array(soiltype[IDF==H_Id])
        borFT_lowerdep = np.array(lowerdep[IDF==H_Id])
        
        HBA = []
        for DIN in borFT_soiltype:
            #print(ba)
            DIN_temp = DIN.split(sep=",")
            HBA_temp = 10
            if 'G' in DIN_temp[0]: HBA_temp = 3
            if 'A' in DIN_temp[0]: HBA_temp = 10
            if 'S' in DIN_temp[0]: HBA_temp = 2
            if 'U' in DIN_temp[0]: HBA_temp = 1
            if 'T' in DIN_temp[0]: HBA_temp = 1
            if 'F' in DIN_temp[0]: HBA_temp = 1
            if 'Z' in DIN_temp[0]: HBA_temp = 4
            HBA.append(HBA_temp)
 
        
        
        DEG2RAD = np.pi/180.0
        
        if verticaloption:
            azm = 0
            dip = 90
            razm = azm * DEG2RAD
            rdip = -dip * DEG2RAD 
        else:
            razm = azm[i] * DEG2RAD
            rdip = -dip[i] * DEG2RAD 
        # do the conversion    
        xn = np.sin(razm) * np.cos(rdip)
        yn = np.cos(razm) * np.cos(rdip)
        zn = np.sin(rdip)
        xi = x[i] 
        yi = y[i]
        zi = z[i]

        lenborFT = len(HBA)

        xp = xi + xn * borFT_lowerdep 
        yp = yi + yn * borFT_lowerdep 
        zp = zi + zn * borFT_lowerdep
        
        if qbasisoption:
            temp_list = []
            temp_list2 = []
            for DIN in borFT_soiltype:
                DIN_temp = DIN.split(sep=",")
                temp_list.append('G' in DIN_temp[0] or 'A' in DIN_temp[0])
                if len(DIN_temp) > 1:
                    temp_list2.append('g*' in DIN_temp[1])
                if len(DIN_temp) == 1:
                    temp_list2.append(False)
            count = 9999
            count2 = 7777
            if True in temp_list:
                
                for var in range(len(temp_list)):
                    if temp_list[var] == True:
                        count = var
                        count3 = 0
                    if temp_list2[var] == True:
                        count2 = var
                    if temp_list[var] == False:
                        count3 = count + 1
                        if count3 == qbasistol:
                            break
            
            if count2 == count + 1:
                count = count2

            if len(temp_list) > count+1:
                xii = x[i]
                yii = y[i]
                zii = z[i] - borFT_lowerdep[count]
                if np.isfinite(zi):
                    qbasis.append([H_Id,xii,yii,zii])
                else:
                    qbasis.append([H_Id,xii,yii,np.NaN])
            else:
                qbasis.append([H_Id,xii,yii,np.NaN])

        try:    
    # classify in interval
            cmin = np.ceil(zp[lenborFT-1] / compositing_interval) * compositing_interval    # round interval min
            cmax = np.ceil(zi / compositing_interval) * compositing_interval                # round interval min
            zcmps = np.arange(cmin,cmax,compositing_interval)[::-1]                         # make z - array 
            xcmps = xi + (xn * (zcmps - min(zcmps)) + zi - max(zcmps))
            ycmps = yi + (yn * (zcmps - min(zcmps)) + zi - max(zcmps)) 
    # create composite 

            cmps = np.digitize(zcmps,zp, right=False)
            borcmps = np.array(HBA)[cmps]

            if i == 0:
                points = np.c_[xcmps,ycmps,zcmps,borcmps]
            else:
                points1 = np.c_[xcmps,ycmps,zcmps,borcmps]
                points = np.append(points,points1,axis=0)
        # make vtk
            if vtkoption:
                    p0 = np.array([xi,yi,zi])
                    p1 = np.append([p0],np.array([xp,yp,zp]).transpose((1,0)), axis=0)         
                    p = np.append(p, p1, axis=0)
                    lenborFT = len(borFT_soiltype)
                        
                    a = np.full((1, lenborFT), 2)[0]
                    b = np.arange(r, r + lenborFT)
                    c = np.arange(r + 1, r + lenborFT+1)
                    c1 = np.array([a,b,c]).transpose((1,0))
                    cells = np.append(cells, c1, axis=0)
                    r = len(p)
                    if i == 0:
                        p = p[1:,:]
                        cells = cells[1:,:]
                        r = r-1
                       

                    HB_array = np.append(HB_array,HBA,axis=0)
        except:          
            faulty_H_Id.append(H_Id)

    if vtkoption:
        poly = pv.PolyData()
        poly.points = p  
        poly.lines = cells
        poly.lines = cells
        line = poly
        line
        line["scalars"] = np.arange(len(cells))
    
        line['litho'] = HB_array #1 + FROM_TO['Tertiär'] + FROM_TO['Bindig']
        #tube = line.tube(radius=1)
        #tube.plot(scalars='litho',smooth_shading=True)
        point_cloud = pv.PolyData(np.array(points[:,0:3],dtype=float))
        point_cloud["soiltype"] = np.array(points[:,-1],dtype=float)
        print("\n VTK - Error in:\n ")
        print(faulty_H_Id)
    else:
        line = False
        point_cloud = False
    return points,point_cloud,line,qbasis


# ----------------------------------------------------------------------------
# def moving_average_3d_kdtree(points,values,grid_params,radius_xy,radius_z):
    
#     nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()

#     scaling_factor_z = radius_xy / radius_z
#     scale_points = np.c_[points[:,0],points[:,1],points[:,2]*scaling_factor_z]
#     tree = cKDTree(scale_points)

#     x_coordinates = [xmn + x * xsiz for x in range(nx)]
#     y_coordinates = [ymn + y * ysiz for y in range(ny)]
#     z_coordinates = [zmn*scaling_factor_z + z * zsiz*scaling_factor_z for z in range(nz)]

#     coordinates_grid = [(x, y, z) for z in z_coordinates for y in y_coordinates for x in x_coordinates  ]

#     grid_n_neighbours = []
#     grid_data = []

#     for k in tqdm(range(len(coordinates_grid))):

#         x,y,z_scaled = coordinates_grid[k]
#         distances, indices = tree.query([(x, y, z_scaled)], k=1000, distance_upper_bound=radius_xy)
#         mask = indices != tree.n
#         distances = distances[mask]
#         indices = indices[mask]
#         #weights = np.exp(-distances ** 2 / radius_xy ** 2)
#         # 1. gaussian
#         # weights = np.exp((-(distances/(radius_xy*1.5/3))**2))
#         # linear
#         weights = 1.001 - distances / radius_xy
#         threshold = 10
#         factor = 1
#         grid_n_neighbours.append(len(distances))
#         # if sum(weights) < threshold:
#         #     factor = sum(weights) / threshold
            
#         if len(weights) > 0: 
#             grid_data.append(
#                 factor * np.average(values[indices], weights=weights
#                                     ) + (1-factor) * 0.5)            
#         else: grid_data.append(0.5)
#     return grid_data, grid_n_neighbours


def moving_average_3d_kdtree(points,values,grid_params,radius_xy,radius_z):
    
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()

    scaling_factor_z = radius_xy / radius_z
    scale_points = np.c_[points[:,0],points[:,1],points[:,2]*scaling_factor_z]
    tree = cKDTree(scale_points)

    x_coordinates = [xmn + x * xsiz for x in range(nx)]
    y_coordinates = [ymn + y * ysiz for y in range(ny)]
    z_coordinates = [zmn*scaling_factor_z + z * zsiz*scaling_factor_z for z in range(nz)]

    coordinates_grid = [(x, y, z) for z in z_coordinates for y in y_coordinates for x in x_coordinates  ]

    grid_n_neighbours = []
    grid_data = []

    for k in tqdm(range(len(coordinates_grid))):

        x,y,z_scaled = coordinates_grid[k]
        distances, indices = tree.query([(x, y, z_scaled)], k=1000, distance_upper_bound=radius_xy)
        mask = indices != tree.n
        distances = distances[mask]
        indices = indices[mask]
        #weights = np.exp(-distances ** 2 / radius_xy ** 2)
        # # 1. gaussian
        # weights = np.exp((-(distances/(radius_xy*1.5/3))**2))
        # 2. IDW
        power = 1.8
        weights = 1 / (distances ** power)

        grid_n_neighbours.append(len(indices))

        if len(weights) > 0:
            grid_data.append(np.average(values[indices], weights=weights))

        else:
            grid_data.append(np.nan)

    return grid_data, grid_n_neighbours
        
from concurrent.futures import ThreadPoolExecutor, as_completed 
from multiprocessing import Pool  
           
def moving_average_3d_kdtree_parallel(points, values, grid_params, radius_xy, radius_z):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn, zmn = grid_params.values()

    scaling_factor_z = radius_xy / radius_z
    scale_points = np.c_[points[:, 0], points[:, 1], points[:, 2] * scaling_factor_z]
    tree = cKDTree(scale_points)

    x_coordinates = [xmn + x * xsiz for x in range(nx)]
    y_coordinates = [ymn + y * ysiz for y in range(ny)]
    z_coordinates = [zmn * scaling_factor_z + z * zsiz * scaling_factor_z for z in range(nz)]

    coordinates_grid = [(x, y, z) for z in z_coordinates for y in y_coordinates for x in x_coordinates]

    def process_point(coord):
        x, y, z_scaled = coord
        distances, indices = tree.query([(x, y, z_scaled)], k=1000, distance_upper_bound=radius_xy)
        mask = indices != tree.n
        distances = distances[mask]
        indices = indices[mask]
        weights = 1.001 - distances / radius_xy
        threshold = 10
        factor = 1
        if sum(weights) < threshold:
            factor = sum(weights) / threshold
        if len(weights) > 0:
            return factor * np.average(values[indices], weights=weights) + (1 - factor) * 0.5
        else:
            return 0.5

    grid_data = []
    
    with Pool() as pool:
        grid_data = list(tqdm(pool.imap(process_point, coordinates_grid), total=len(coordinates_grid)))

    return grid_data
        # # 2. IDW
        # power = 1.8
        # weights = 1 / (distances ** power)

        # grid_n_neighbours.append(len(indices))

        # if len(weights) > 20:
        #     if len(weights) < 100:
        #         factor =  (len(weights)-20)/80
        #         grid_data.append(
        #             factor * np.average(values[indices], weights=weights
        #                                 ) + (1-factor) * 0.5)
        #     else:
        #         grid_data.append(np.average(values[indices], weights=weights))

        # else:
        #     grid_data.append(0.5)

    return grid_data, grid_n_neighbours
# ----------------------------------------------------------------------------
def moving_average_2d_kdtree(qbasis,grid_params,radius_x,radius_y,tolerance = 10):
    data = np.array(qbasis)[:,1:].astype(float)
    data = data[~np.isnan(data[:,2])]
    data = data[data[:,2] > np.mean(data[:,2])-tolerance]
    points,values = data[:,:2],data[:,2]
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()

    scaling_factor_y = radius_x / radius_y
    scale_points = np.c_[points[:,0],points[:,1]*scaling_factor_y]
    tree = cKDTree(scale_points)

    x_coordinates = [xmn + x * xsiz for x in range(nx)]
    y_coordinates = [ymn*scaling_factor_y + y * ysiz*scaling_factor_y for y in range(ny)]

    coordinates_grid = [(x, y) for y in y_coordinates for x in x_coordinates]

    grid_n_neighbours = []
    grid_data = []

    for k in tqdm(range(len(coordinates_grid))):
        x,y_scaled = coordinates_grid[k]
        distances, indices = tree.query([(x, y_scaled)], k=100, distance_upper_bound=radius_x)
        mask = indices != tree.n
        distances = distances[mask]
        indices = indices[mask]
        weights = np.exp(-distances / (radius_x/3) )
        grid_n_neighbours.append(len(indices))
        if len(weights) > 0:
            grid_data.append(np.average(values[indices], weights=weights))
        else:
            grid_data.append(np.nan)
    return np.array(grid_data), np.array(grid_n_neighbours), points, values




    
def checkgrid(xlist,ylist,zlist,grid_params,XY_option= True):
    # check if out of bounds
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()

    #  CELL SEARCH
    xsc = np.arange(0, nx+1, 1)*xsiz + xmn - 0.5 * xsiz
    ysc = np.arange(0, ny+1, 1)*ysiz + ymn - 0.5 * ysiz
    zsc = np.arange(0, nz+1, 1)*zsiz + zmn - 0.5 * zsiz
    ixlist = np.digitize(xlist,xsc,right=True)
    iylist = np.digitize(ylist,ysc,right=True)
    izlist = np.digitize(zlist,zsc,right=True)
    #np.searchsorted(np.sort(ysc), vym, side='left')
    # PARAMETER SEARCH
    # search for values in grid
    # search grid cell (index location)
    loc = (izlist-1)*nx*ny + (iylist-1)*nx + ixlist - 1                         # Python counts from 0 !!!!
    if XY_option:
        loc = (iylist-1)*nx + ixlist - 1                                        # Python counts from 0 !!!!
    
    check = (ixlist>0)&(ixlist<=nx)&(iylist>0)&(iylist<=ny)&(izlist>0)&(izlist<=nz)
    return(loc,check)
# =====================================================
# Statistics
# =====================================================


@np.errstate(invalid='ignore')
def Soil_to_Y(Soil, can_be_negative):
    """
    Calculation of Y from the soil parameters;
    Application of ln() if soil parameter cannot be negative
    """

    Y = np.empty(Soil.shape)  # Presetting Y as an array
    for i in range(can_be_negative.size):
        if can_be_negative[i]:
            # If soil parameter can be negative, it is simply taken over as such
            Y[:, i] = Soil[:, i]
        else:  # If soil parameter cannot be negative, calculate ln() from this
            Y[:, i] = np.log(Soil[:, i])
    return Y


@np.errstate(invalid='ignore')
def calc_dist_params(Y):
    """
    Calculation of the Johnson distribution parameters ax, bx, ay, by, bx_star 
    and the Johnson distribution type dist_type (SB, SU or SL) accordingly: 
    'Ching & Phoon (2014) - Correlations among some clay parameters — the 
    multivariate distribution'
    """

    z = 0.7
    quantiles = norm.cdf((-3 * z, -z, z, 3 * z))  # Calculation of the percentiles
    # Calculation of the corresponding y-values
    ya, yb, yc, yd = np.nanquantile(Y, quantiles, axis=0)
    m = yd - yc
    n = yb - ya
    p = yc - yb
    D = m * n / p ** 2

    # Set the distribution types according to the value of D
    dist_type = np.select(
        [D > 1, D < 1, D == 1],
        ['SU', 'SB', 'SL'])

    # Calculation of the ax parameter according to the distribution type
    ax = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [np.maximum(0, 2 * z / np.arccosh(0.5 * (m / p + n / p))),
         np.maximum(0, z / np.arccosh(0.5 * ((1 + p / m) * (1 + p / n)) ** 0.5)),
         2 * z / np.log(m / p)])

    # Calculation of the bx parameter according to the distribution type
    # (Not present for SL)
    bx = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [ax * np.arcsinh((n / p - m / p) / (2 * (D - 1) ** 0.5)),
         ax * np.arcsinh(
             (p / n - p / m) * ((1 + p / m) * (1 + p / n) - 4) ** 0.5 / (
                         2 * (1 / D - 1))),
         np.nan])

    # Calculation of the bx* parameter according to the distribution type
    # (Not present for SU and SB)
    bx_star = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [np.nan,
         np.nan,
         ax * np.log((m / p - 1) / (p * (m / p) ** 0.5))])

    # Calculation of the ay parameter according to the distribution type
    # (Not present for SL)
    ay = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [np.maximum(0, 2 * p * (D - 1) ** 0.5 / (
                    (m / p + n / p - 2) * (m / p + n / p + 2) ** 0.5)),
         np.maximum(0, p * (((1 + p / m) * (1 + p / n) - 2) ** 2 - 4) ** 0.5 / (
                     1 / D - 1)),
         np.nan])

    # Calculation of the by parameter according to the distribution type
    by = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [(yc + yb) / 2 + p * (n / p - m / p) / (2 * (m / p + n / p - 2)),
         (yc + yb) / 2 - ay / 2 + p * (p / n - p / m) / (2 * (1 / D - 1)),
         (yc + yb) / 2 - 0.5 * p * (m / p + 1) / (m / p - 1)])

    return dist_type, ax, bx, ay, by, bx_star


@np.errstate(invalid='ignore')
def Y_to_X(Y, dist_type, ax, bx, ay, by, bx_star):
    """
    Calculation of X from Y according to the respective Johnson distribution type
    (SU, SB or SL)
    """

    # Select conversion according to distribution type (SU, SB, SL)
    X = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [np.arcsinh((Y - by) / ay) * ax + bx,
         np.log((Y - by) / (ay + by - Y)) * ax + bx,
         np.log(Y - by) * ax + bx_star])
    return X


@np.errstate(over='ignore', invalid='ignore')
def X_to_Y(X, dist_type, ax, bx, ay, by, bx_star):
    """
    Calculation of Y from X according to the respective Johnson distribution type
    (SU, SB or SL)
    """

    # Select conversion according to distribution type (SU, SB, SL)
    Y = np.select(
        [dist_type == 'SU', dist_type == 'SB', dist_type == 'SL'],
        [np.sinh((X - bx) / ax) * ay + by,
         (np.exp((X - bx) / ax) * (ay + by) + by) / (1 + np.exp((X - bx) / ax)),
         np.exp((X - bx_star) / ax) + by])
    return Y


@np.errstate(over='ignore')
def Y_to_Soil(Y, can_be_negative):
    """
    Calculation of soil parameters from Y;
    Application of exp() if soil parameter cannot be negative
    """

    Soil = np.empty(Y.shape)
    for i in range(can_be_negative.size):
        if can_be_negative[i]:
            # If soil parameter can be negative, Y is simply taken over as such
            Soil[:, i] = Y[:, i]
        else:  # If soil parameter cannot be negative, calculate exp() of Y
            Soil[:, i] = np.exp(Y[:, i])
    return Soil


def cumdist(array):
    x = np.sort(array)
    y = (np.arange(len(array))+1)/len(array)
    return x,y


def probdist(array):
    x = np.sort(array)
    mu = np.mean(array)
    sigma = np.std(array)
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *np.exp(-0.5 * (1 / sigma * (x - mu))**2))
    return x,y

def plane_from_points(x, y, z):
    from numpy.linalg import inv
    # Create this matrix correctly without transposing it later?
    points = np.array([x, y, z])
    A = np.array([
        points[0,:],
        points[1,:],
        np.ones(points.shape[1])
    ]).T
    b = np.array([points[2, :]]).T

    # fit = (A.T * A).I * A.T * b
    fit = np.dot(np.dot(inv(np.dot(A.T, A)), A.T), b)
    errors = b - np.dot(A, fit)
    residual = np.linalg.norm(errors)
    return fit, errors, residual

def remove_trend(x,y,z):
    fit_params, errors, residual = plane_from_points(x, y, z)
    np.mean(abs(errors))
    # detrending
    zd = z - (fit_params[0] * x + fit_params[1] * y + fit_params[2])
    return(zd,fit_params)

def add_trend_single(zd,x,y,fit_params,nsim):
    z = zd + (fit_params[0] * x + fit_params[1] * y + fit_params[2])
    return(z)

def add_trend_grid(zd,fit_params,grid_params,nsim):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    x,y = mgrid_2D(grid_params)
    z = zd + (fit_params[0] * np.tile(x, nsim) + fit_params[1] * np.tile(y, nsim) + fit_params[2])
    return(z)

def transform_to_gaussian(zd, trans_method):
    if trans_method == "Johnson":
        dist_type, ax, bx, ay, by, bx_star = calc_dist_params(zd)
        trans_params = dist_type, ax, bx, ay, by, bx_star
        zdn = Y_to_X(zd, dist_type, ax, bx, ay, by, bx_star)
        return(zdn,trans_method,trans_params)
    
def backtransform_from_gaussian(df, trans_method, trans_params):
    if trans_method == "Johnson":
        dist_type, ax, bx, ay, by, bx_star = trans_params
        Y = X_to_Y(df, dist_type, ax, bx, ay, by, bx_star)
        return(Y)
    
def mgrid_2D(grid_params):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    nz = 1
    grid = np.mgrid[0:nz,0:ny,0:nx]
    x = grid[2].flatten() * xsiz + xmn
    y = grid[1].flatten() * ysiz + ymn
    #z = grid[0].flatten() * zsiz + zmn
    return x,y

def buffer_grid_params(grid_params,buffer):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    buffered_grid_params = {"nx" : int(nx + 2*buffer/xsiz),
                            "ny" : int(ny + 2*buffer/ysiz),
                            "nz" : int(nz + 2*buffer/zsiz),
                            "xsiz" : xsiz,"ysiz" : ysiz,"zsiz" : zsiz,
                            "xmn" : xmn-buffer,
                            "ymn" : ymn-buffer,
                            "zmn" : zmn-buffer}
    return buffered_grid_params


def unit_vector(v):
    # Calculate the magnitude
    magnitude = np.linalg.norm(v,axis=2)
    # Calculate the unit vector
    unit_vector = np.divide(v , np.expand_dims(magnitude, axis=2))
    return unit_vector


def slope(grid, raster_size=[25,25]):
    """Calculates the slope and azimuth of a grid.
      
    Args:
      grid: A NumPy array of elevation values.
      raster_size: The size of each raster cell in meters.
      
    Returns:
      A tuple of two NumPy arrays, one with the slopes in degrees and one with
      the azimuths in degrees from 0 to 360, with 0 degrees being north.
    """
      
    # Calculate the gradients in the x and y directions.
    dx = np.gradient(grid, axis=1)
    dy = np.gradient(grid, axis=0)
    
    ones = np.ones(np.shape(grid))
    
    vx = unit_vector(np.dstack([raster_size[0]*ones,ones-1,-dx]))
    vy = unit_vector(np.dstack([ones-1,raster_size[1]*ones,-dy]))
    
    vxy = np.cross(vy,vx)

    vxy[vxy[:,:,2] > 0] = vxy[vxy[:,:,2] > 0] * -1 # always dip
    azm = np.arctan2(vxy[:,:,1],vxy[:,:,0])*180/np.pi

    azm = (90-azm)%360
    hor = (vxy[:,:,1]**2+vxy[:,:,0]**2)**0.5
    dip = np.arctan2(vxy[:,:,2],hor)*180/np.pi + 90
    
    return azm, dip


def detrend_2D(qbasis,grid_params,moving_window=False,
                           radius_xy=[1000,1000], tolerance = 10,
                           outlier_tolerance = 5):
    
    buffer=max(radius_xy)
    buffered_grid_params = buffer_grid_params(grid_params,buffer)
    data = np.array(qbasis)[:,1:].astype(float)
    data = data[~np.isnan(data[:,2])]
    data = data[data[:,2] > np.mean(data[:,2])-tolerance]
    x,y,z = data[:,0],data[:,1],data[:,2]
    xd,yd = x,y
    if moving_window:
        radius_x,radius_y = radius_xy
        result = moving_average_2d_kdtree(np.c_[np.ones(len(data)),data],buffered_grid_params,
                                          radius_x,radius_y,tolerance=tolerance)
        grid_data_buf, grid_n_neighbours, points, values = result
        loc, check = checkgrid(x, y, z, buffered_grid_params)
        zd = z[check] - grid_data_buf[loc[check]]
        azm, dip = slope(grid_data_buf.reshape((buffered_grid_params["ny"],
                                                buffered_grid_params["nx"])),           
                         raster_size=[buffered_grid_params["xsiz"], 
                                      buffered_grid_params["ysiz"]])
        z_dip = dip.flatten()[loc[check]]
        data2 = data[check]
        data2[(z_dip > 0.6) & (abs(zd) > outlier_tolerance)] = np.nan
        data2[(z_dip <= 0.6) & (abs(zd) > outlier_tolerance)] = np.nan
        data3 = data2[~np.isnan(data2[:,-1])]
        x,y = mgrid_2D(grid_params)
        z1 = np.ones(len(x))
        loc1, _ = checkgrid(x, y, z1, buffered_grid_params)
        grid_data = grid_data_buf[loc1]
        coords_d = np.c_[data3[:,0],data3[:,1],zd[~np.isnan(data2[:,-1])]]
    else:
        nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = buffered_grid_params.values()
        zd, fit_params = remove_trend(x,y,z)
        x,y =mgrid_2D(buffered_grid_params)
        grid_data_buf = fit_params[0] * x + fit_params[1] * y + fit_params[2]
        x,y =mgrid_2D(grid_params)
        grid_data = fit_params[0] * x + fit_params[1] * y + fit_params[2]
        loc, check = checkgrid(xd, yd, z, buffered_grid_params)
        zd = z[check] - grid_data_buf[loc[check]]
        coords_d = np.c_[xd[check],yd[check],zd]
        data3 = data
    return data3, coords_d, grid_data
#Wang：  Data Output: The function returns the detrended data (data3), the coordinates of detrended data (coords_d), and the grid data after detrending (grid_data)


from scipy.stats import rankdata, norm

def transform(coords, trans_method="Quantile"):
    zd = coords[:,2]
    if trans_method == "Johnson":
    # transform to gaussian
        zt, trans_method, trans_params = transform_to_gaussian(zd, trans_method)  
    if trans_method == "Quantile":
        rng = np.random.RandomState(304)
        trans_params = QuantileTransformer(n_quantiles=100, output_distribution="normal",
                         random_state=rng)
        zt = trans_params.fit_transform(zd.reshape(-1, 1))
    if trans_method == "Normal_Score":
        zt = norm.ppf(rankdata(zd)/(len(zd) + 1))
        trans_params = zd
    return np.c_[coords[:,0],coords[:,1],zt], trans_method, trans_params
#wang:The purpose of this code is to convert z-values in geologic data to a specified distribution for further analysis.

def backtransform(zdn,trans_method, trans_params):
    if trans_method == "Johnson":
        zd = backtransform_from_gaussian(zdn, trans_method, trans_params)  
    if trans_method == "Quantile":
        zd = trans_params.inverse_transform(zdn.reshape(-1, 1)).flatten()
    if trans_method == "Normal_Score":
        sorted_indices = np.argsort(trans_params)
        rank = np.digitize(np.argsort(zdn) * len(zdn) / len(trans_params),
                           np.arange(len(trans_params)))
        rank[rank == 0] = 1
        zd = trans_params[rank-1]
    return zd





def post_processing_2D(outfl,grid_data,GSLIB_SGSIM_params,
                    trans_method, trans_params):
    nsim = GSLIB_SGSIM_params["nsim"]
    sgsim = readsgsim(outfl)
    sgsim_dt = backtransform(sgsim, trans_method, trans_params)
    sgsim_dtbt = sgsim_dt + np.tile(grid_data,nsim)
    return sgsim, sgsim_dt, sgsim_dtbt
    
def plot_grid_2D(grid_data,grid_params,points=[[]], title = "unknown"):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    x = xmn + xsiz * (np.arange(nx))
    y = ymn + ysiz * (np.arange(ny))

    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(8,12),dpi=600)

    im = ax.pcolormesh(X, Y, grid_data.reshape(ny,nx,order="C"), cmap='viridis', shading='auto')
    cbar = plt.colorbar(im, ax=ax)
    plt.scatter(points[:,0],points[:,1],c="Black",marker="x")
    plt.xlim(xmn-0.5*xsiz, xmn-0.5*xsiz + xsiz * nx)
    plt.ylim(ymn-0.5*ysiz, ymn-0.5*ysiz + ysiz * ny)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(title, pad=20)
    plt.show()
    
def plot_contours_2D(grid_data,grid_params,points):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    x = xmn + xsiz * (np.arange(nx))
    y = ymn + ysiz * (np.arange(ny))

    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(8,12),dpi=600)
    im = ax.contourf(X,Y,grid_data.reshape(ny, nx,order="C"))

    #im = ax.pcolormesh(X, Y, grid_data.reshape(ny,nx,order="C"), cmap='viridis', shading='auto')
    cbar = plt.colorbar(im, ax=ax)
    plt.scatter(points[:,0],points[:,1],c="Black",marker="x")
    plt.xlim(xmn-0.5*xsiz, xmn-0.5*xsiz + xsiz * nx)
    plt.ylim(ymn-0.5*ysiz, ymn-0.5*ysiz + ysiz * ny)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D Grid')
    plt.show()


def plot_result_2D(sgsim, sgsim_dt, sgsim_dtbt,
                   GSLIB_SGSIM_params,coords):
    nsim = GSLIB_SGSIM_params["nsim"]
    grid_params_2D = GSLIB_SGSIM_params["grid_params"]
    plt.figure(figsize=(8, 6), dpi=350) 
    plt.hist(sgsim,bins=np.arange(-6,6,0.3),density=True)
    plt.title('Histogram of SGSIM Results')
    plt.hist(coords[:,2],bins=np.arange(-6,6,0.3),density=True)
    

    sgsim_mean = np.mean(np.array_split(sgsim_dtbt, nsim), axis=0)
    plot_grid_2D(sgsim_mean, grid_params_2D, coords, "Mean of SGSIM Results")
    
    sgsim_min = np.min(np.array_split(sgsim_dtbt, nsim), axis=0)
    plot_grid_2D(sgsim_min, grid_params_2D, coords, "Minimum of SGSIM Results")
    
    sgsim_max = np.max(np.array_split(sgsim_dtbt, nsim), axis=0)
    plot_grid_2D(sgsim_max, grid_params_2D, coords, "Maximum of SGSIM Results")
    
    plot_grid_2D(sgsim_max - sgsim_min, grid_params_2D, coords, "Range of SGSIM Results")
    
    sgsim_std = np.std(np.array_split(sgsim_dt, nsim), axis=0)
    plot_grid_2D(sgsim_std, grid_params_2D, coords, "Standard Deviation of SGSIM Results")


def PANDASDF2GSLIBGeoEAS(df,gslib_output):
    # Creating the output file containing the composited points in compliant with GSLIB format
    f = open(gslib_output,'w',encoding='utf8') #Wang: create a Object name f and use open function to open it with mode write and with encoding.... and assign it to f.

    # Writing the title or name of the output file at the first row
    f.write(gslib_output)
    f.write('\n')
    # Extracting 'n' number of properties to be modeled (number of columns) at the second row
    #used?# property_number =                    # Defining quantity of properties: number of columns - HoleID, from, to, and Geology
    f.write(str(len(df.columns)))                         # Writing number of variables: coordinates X,Y,Z plus quantity of properties
    f.write('\n')

    # Writing 'n' lines with the names of the properties at each row
    for n in range(len(df.columns)):
        f.write(df.columns[n])
        f.write('\n')
    # Writing the dataframe properties
    f.write(df.to_string(header=False, index=False, index_names=False))
    f.close()
    
def readsgsim(data_file):
    with open(data_file) as file:
        data = file.read()
    test =data[data.find("value\n")+6:]   
    array = np.array(test.split(), dtype=float)
    return array


def split_array(k, coords, coords_t=None):
    train = []
    train_t = []
    test = []
    test_t = []
    start_idx = 0
    
    if coords_t != None:
        parallel = np.c_[coords,coords_t]
        group_size = len(coords_t) // k
        np.random.shuffle(parallel)
        coords,coords_t = parallel[:,:3],parallel[:,3:]

        for i in range(k):
            if i == 0:
                group_t = coords_t[ : start_idx + group_size]
                remaining_t = coords_t[start_idx + group_size:] 
            if i == k-1:
                group = coords_t[start_idx : ]
                remaining_t = coords_t[:start_idx]
                remaining = coords[:start_idx]
            else:
                group_t = coords_t[start_idx : start_idx + group_size]
                remaining_t = np.concatenate((coords_t[:start_idx] , coords_t[start_idx + group_size:]), axis=0)

            start_idx += group_size
            train.append(remaining_t)
            test_t.append(group_t)
    else:
        group_size = len(coords) // k
        np.random.shuffle(coords)
    for i in range(k):
        if i == 0:
            group = coords[ : start_idx + group_size]
            remaining = coords[start_idx + group_size:] 
        if i == k-1:
            group = coords[start_idx : ]
            remaining = coords[:start_idx]
        else:
            group = coords[start_idx : start_idx + group_size]
            remaining = np.concatenate((coords[:start_idx] , coords[start_idx + group_size:]), axis=0)
        start_idx += group_size
        train.append(remaining)
        test.append(group)
    return train, train_t, test, test_t
def split_array_3D(k, Cl_Sa_points):
    coords = Cl_Sa_points
    train = []

    test = []

    start_idx = 0
    
    _, indices = np.unique(coords[:,:2], axis=0, return_inverse=True)

    group_size = (max(indices)+1) // k
    index = np.arange(0,max(indices)+1,1)
    np.random.shuffle(index)
    for i in range(k):
        print(i)
        if i == 0:
            states = index[ : group_size]
            remaining = coords[~np.in1d(indices, states)]
            group = coords[np.in1d(indices, states)]
            
 
        if i == k-1:
            states = index[start_idx : ]
            remaining = coords[~np.in1d(indices, states)]
            group = coords[np.in1d(indices, states)]
            
        else:
            states = index[start_idx : start_idx + group_size]
            group = coords[np.in1d(indices, states)]
            remaining = coords[~np.in1d(indices, states)]
        start_idx += group_size
        train.append(remaining)
        test.append(group)
    return train, test




def run_parallel(WD_exe,WD_par,k):
    # run file
    commands = []
    for i in range(k):
        parfl = str(WD_par) + str(i) + ".par"  
        cmd = str(WD_exe) + " " + str(parfl)
        commands.append(cmd)
        
    from subprocess import Popen
    procs = [ Popen(i) for i in commands ]
    for p in procs:
        p.wait()
        
def readsisim(data_file, lines2ignore = 3):
    with open(data_file) as file:
        count = 0
        lines = []
        for line in file:
            count += 1
            if count <= lines2ignore:
                continue  # Skip the line
            else:
                lines.append(int(line[2]))
    return lines

def run_sisim_lm_parallel(sisim_lmexe,WD_par,WD_out,WD_dbg,datafl,priormfl,
                          GSLIB_SISIM_LM_params):
    print("-- generate parfiles")
    GSLIB_SISIM_LM_params["nsim"] = 5 
    commands = []   
    for part in tqdm(range(20)):
        parfl = str(WD_par) + str(part) + ".par"
        output = str(WD_out)+str(part)+".out"
        dbgfl = str(WD_dbg) + str(part)+".dbg"
        GSLIB_SISIM_LM_params["seed"] = 69069 + 11 * part
        
        GSLIB_SISIM_LM(sisim_lmexe,parfl,datafl,priormfl,
                       dbgfl,output,GSLIB_SISIM_LM_params)   
        cmd = str(sisim_lmexe) + " " + str(parfl)
        commands.append(cmd)
        
    print("-- running simulations")        
    
    for batch in tqdm(range(4)):
        b_start = batch * 5
        b_end = 5 + batch * 5
        procs = [ Popen(i) for i in commands[b_start:b_end] ]
        for p in procs:
            p.wait()
            
def run_sisim_parallel(sisimexe,WD_par,WD_out,WD_dbg,datafl,
                          GSLIB_SISIM_params):
    print("-- generate parfiles")
    GSLIB_SISIM_params["nsim"] = 5 
    commands = []   
    for part in tqdm(range(20)):
        parfl = str(WD_par) + str(part) + ".par"
        output = str(WD_out)+str(part)+".out"
        dbgfl = str(WD_dbg) + str(part)+".dbg"
        GSLIB_SISIM_params["seed"] = 69069 + 11 * part
        
        GSLIB_SISIM(sisimexe,parfl,datafl,
                       dbgfl,output,GSLIB_SISIM_params)   
        cmd = str(sisimexe) + " " + str(parfl)
        commands.append(cmd)
        
    print("-- running simulations")        
    
    for batch in tqdm(range(4)):
        b_start = batch * 5
        b_end = 5 + batch * 5
        procs = [ Popen(i) for i in commands[b_start:b_end] ]
        for p in procs:
            p.wait()           
            

def read_sisim_lm(WD_out):
    print("-- reading simulations")
    all_sim = []
    for part in tqdm(range(20)):
        output = str(WD_out)+str(part)+".out"
        int_list = readsisim(output)
        chunk_size = len(int_list) / 5
        start = 0
        for i in range(5):
            end = start + int(chunk_size)
            all_sim.append(int_list[start:end])
            start = end
    return all_sim
    
def read_sisim(WD_out):
    print("-- reading simulations")
    all_sim = []
    for part in tqdm(range(20)):
        output = str(WD_out)+str(part)+".out"
        int_list = readsisim(output)
        chunk_size = len(int_list) / 5
        start = 0
        for i in range(5):
            end = start + int(chunk_size)
            all_sim.append(int_list[start:end])
            start = end
    return all_sim
    
        
def save_train(train,WD_data,k,
                       grid_params, moving_window=False,
                       radius_xy=[600,600], tolerance = 10,
                       trans_method="Quantile"):
    detrend_transform_params = []
    for i in range(k):
        data = np.c_[np.ones(len(train[i])),train[i]]
        result = detrend_2D(data,grid_params, moving_window,
                                   radius_xy, tolerance)
        coords, coords_d, grid_data = result
        coords_d = coords_d[(coords_d[:,2]>-10)&(coords_d[:,2]<10)]
        coords_t,trans_method, trans_params = transform(coords_d,"Quantile")

        file_path = str(WD_data) + str(i) + ".out"
        PANDASDF2GSLIBGeoEAS(pd.DataFrame(coords_t,columns=["X","Y","Z"]),file_path)
        detrend_transform_params.append([trans_method, trans_params,grid_data])
    return detrend_transform_params

def save_train_no_transofrm(train,WD_data,k,
                       grid_params, moving_window=False,
                       radius_xy=[600,600], tolerance = 10):
    detrend_transform_params = []
    for i in range(k):
        data = np.c_[np.ones(len(train[i])),train[i]]
        result = detrend_2D(data,grid_params, moving_window,
                                   radius_xy, tolerance)
        coords, coords_d, grid_data = result
        coords_d = coords_d[(coords_d[:,2]>-10)&(coords_d[:,2]<10)]

        file_path = str(WD_data) + str(i) + ".out"
        PANDASDF2GSLIBGeoEAS(pd.DataFrame(coords_d,columns=["X","Y","Z"]),file_path)
        detrend_transform_params.append([grid_data])
    return detrend_transform_params

def k_fold_cross_validation(k,test,WD_outfl,GSLIB_SGSIM_params,
                            detrend_transform_params):
    
    
    grid_params = GSLIB_SGSIM_params["grid_params"]
    nsim = GSLIB_SGSIM_params["nsim"]
    zmn = grid_params["zmn"]    
    error_list = []
    quantile_list = []
    width_list = []
    for i in range(k):
        trans_method, trans_params,grid_data = detrend_transform_params[i]
        outfl = str(WD_outfl) + str(i) + ".out" 
        sgsim = readsgsim(outfl)
        sgsim_dt = backtransform(sgsim, trans_method, trans_params)
        sgsim_bt = sgsim_dt + np.tile(grid_data, nsim)
        #sgsim_mean = np.mean(np.array_split(sgsim_bt,nsim),axis=0)

        xd,yd = test[i][:,0],test[i][:,1]
        loc, check = checkgrid(xd, yd, zmn, grid_params)
        truevals = test[i][:,2][check]


        simar = np.array(np.array_split(sgsim_bt,nsim))
        simvals = simar[:,loc[check]]
        validation = np.mean(abs(truevals[:,None] - simvals.T),axis=1)
        
        error_list.append(validation)
        print("at k = " + str(i) + ", mean absolute error: " + str(np.mean(validation)))
        
        sortedvals = np.sort(simvals,axis=0)
        width = sortedvals[1:,:]-sortedvals[:-1,:]
        width_list.append(np.mean(width,axis=1))
        for i in range(0,len(simvals[0])):
            quantile_list.append(np.digitize(truevals[i],sortedvals[:,i]))
    width_list = np.mean(np.array(width_list),axis=0)
    return error_list, quantile_list, width_list

def plot_KCV(width_list,quantile_list,nsim):
    fig, ax = plt.subplots(figsize=(5,3),dpi=600)
    x,y = np.arange(nsim-1)+1,np.cumsum(width_list)
    y -= y[len(y)//2]
    plt.plot(x/100,y,label="SGSIM")
    plt.legend()
    ax.set_xlim((0,1))
    ylim = np.ceil(max(abs(y)))
    ax.set_ylim((-ylim,ylim))
    ax.set_xlabel('Probability Interval')  
    ax.set_ylabel('Mean Width of\nProbability Interval') 
    plt.show()
    
    fig, ax = plt.subplots(figsize=(6,3),dpi=600)
    x,y = cumdist(quantile_list)
    plt.plot(x/100,y,label="SGSIM")
    plt.plot([0,1],[0,1])
    plt.legend()
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    ax.set_xlabel('Probability Interval')  
    ax.set_ylabel('\nProportion in this Interval') 
    plt.savefig('result SGSIM.svg', format='svg')
    plt.show()

def save_Cl_Sa_points(points,WD_data = None):
    cp = points[:,[0,1,2,-1]].astype(float)
    cp = cp[cp[:,3] < 3]
    cp = cp[cp[:,3] > 0]
    cp[:,-1] =  cp[:,-1] - 1 # san is 1, clay is 0
    if WD_data != None:
        cpdf = pd.DataFrame(cp,columns=["X","Y","Z","ST"])
        PANDASDF2GSLIBGeoEAS(cpdf,WD_data)
    print("Sa : 1, Cl : 0")
    return cp

    
def plot_2_variograms_q(variograms):
    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(8.8*cm,5*cm), dpi=600)
    sill = 1
    ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval']/sill, 'o',
            c='Orange',zorder=-1,label="North", alpha=1)
    ax.plot(variograms[1]['avsepdist'],variograms[1]['svargval']/sill, 'o',
            c='Blue',zorder=-1,label="East", alpha=1)
    theovargm = MAKETHEORVARGM(0/sill,1.2/sill,50,2)
    ax.plot(theovargm['x'],theovargm['y'],c='Black',zorder=-1,label="var", alpha=1) 
    ax.legend()
    ax.set_xlim(0,200)
    # ax.set_ylim(0,2)
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Semivariance [-]')
    plt.show()
    
def plot_2_variograms_t(variograms,t_var = [[0,1,50,2]]):
    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(8.8*cm,5*cm), dpi=600)
    ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval'], 'o',
            c='Orange',zorder=-1,label="North", alpha=1)
    ax.plot(variograms[1]['avsepdist'],variograms[1]['svargval'], 'o',
            c='Blue',zorder=-1,label="East", alpha=1)
    for i in range(len(t_var)):
        vnugget,vsill,vrange,vtype = t_var[i]
        theovargm = MAKETHEORVARGM(vnugget,vsill,vrange,vtype)
        ax.plot(theovargm['x'],theovargm['y'],c='Black',zorder=-1,label="var", alpha=1) 
    ax.legend()
    ax.set_xlim(0,100)
    # ax.set_ylim(0,2)
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Semivariance [-]')
    plt.show()   
def plot_1_variogram_t(variograms,t_var = [[0,1,50,2]]):
    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(8.8*cm,5*cm), dpi=600)
    ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval'], 'o',
            c='Orange',zorder=-1,label="North", alpha=1)
    for i in range(len(t_var)):
        vnugget,vsill,vrange,vtype = t_var[i]
        theovargm = MAKETHEORVARGM(vnugget,vsill,vrange,vtype)
        ax.plot(theovargm['x'],theovargm['y'],c='Black',zorder=-1,label="var", alpha=1) 
    ax.legend()
    ax.set_xlim(0,100)
    # ax.set_ylim(0,2)
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Semivariance [-]')
    plt.show()   


def visualize_grid(grid_data, points, values, alldata_grid_params):

    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = alldata_grid_params.values()
    
    # Create the structured grid
    point_cloud = pv.PolyData(points)

    # Add the scalars array as point data
    point_cloud['scalars'] = values

    x = np.arange(xmn, xmn + (nx + 1) * xsiz, xsiz)
    y = np.arange(ymn, ymn + (ny + 1) * ysiz, ysiz)
    z = np.arange(zmn, zmn + (nz + 1) * zsiz, zsiz)
    
    x, y, z = np.meshgrid(x, y, z)
    structured_grid = pv.StructuredGrid(x, y, z)
 
    # Add the scalar data to the structured grid
    structured_grid.cell_data['scalars'] = grid_data
    
    slices = structured_grid.slice_orthogonal()

    # Create a plotter and add the structured grid
    plotter = pv.Plotter()
    plotter.add_mesh(slices)
    plotter.add_mesh(point_cloud, render_points_as_spheres=True,
                     point_size=10, scalars='scalars', cmap='viridis')

    # Show the plot
    plotter.show()

# def save_grid(grid_data_list, name_list, alldata_grid_params, file_name):

#     nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = alldata_grid_params.values()
    
#     # Create the structured grid
#     point_cloud = pv.PolyData(points)

#     # Add the scalars array as point data
#     point_cloud['scalars'] = values

#     x = np.arange(xmn, xmn + (nx + 1) * xsiz, xsiz)
#     y = np.arange(ymn, ymn + (ny + 1) * ysiz, ysiz)
#     z = np.arange(zmn, zmn + (nz + 1) * zsiz, zsiz)
    
#     x, y, z = np.meshgrid(x, y, z)
#     structured_grid = pv.StructuredGrid(x, y, z)
 
#     # Add the scalar data to the structured grid
#     for i in range(len(grid_data_list)):
#         structured_grid.cell_data[name_list[i]] = grid_data_list[i]
    
#     structured_grid.save(file_name)
    
def moving_window_depth_average_fusion(grid_data, grid_n_neighbours,depth_average,vrange=200):
    weights = 1-np.exp((-(np.array(grid_n_neighbours)/(vrange/3))**2))
    # weights[grid_n_neighbours==0] = 0
    grid_data_temp = np.array(grid_data)
    grid_data_temp[np.isnan(grid_data_temp)] = 0
    fusion_grid = grid_data_temp * weights + depth_average * (1-weights)
    return fusion_grid

def save_prior(fusion_grid,priorfl):
    prior = np.c_[np.ones(len(fusion_grid))-fusion_grid,fusion_grid]
    prior = pd.DataFrame(prior,columns=["p0","p1"])
    #PANDASDF2GSLIBGeoEAS(pd.DataFrame(prior,columns=["p0","p1"]),priorfl)
    # Creating the output file containing the composited points in compliant with GSLIB format
    f = open(priorfl,'w',encoding='utf8')

    # Writing the title or name of the output file at the first row
    f.write(priorfl)
    f.write('\n')
    # Extracting 'n' number of properties to be modeled (number of columns) at the second row
    #used?# property_number =                    # Defining quantity of properties: number of columns - HoleID, from, to, and Geology
    f.write(str(len(prior.columns)))                         # Writing number of variables: coordinates X,Y,Z plus quantity of properties
    f.write('\n')

    # Writing 'n' lines with the names of the properties at each row
    for n in range(len(prior.columns)):
        f.write(prior.columns[n])
        f.write('\n')
    # Writing the dataframe properties
    f.write(prior.to_string(header=False, index=False, index_names=False, float_format="{:.2f}".format))
    f.close()
    

def make_prior(Cl_Sa_points,grid_params,no_data_value=0.5,correctiom = 0.25, nbh_thres = 20):
    cpdf = pd.DataFrame(Cl_Sa_points,columns=["X","Y","Z","ST"])
    meanst = cpdf.astype(float).groupby(by="Z").mean() 
    lenst = cpdf.astype(float).groupby(by="Z").size() 
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    # meanst.plot(y="priorp")
    # gererate array over depth
    ppnew = np.arange(zmn,zmn+nz*zsiz,zsiz)
    prior_mean = np.array(meanst["ST"])
    prior_index = np.array(meanst.index)
    prior = []
    for i in ppnew:
        mean_value = prior_mean[prior_index == i+correctiom]
        if mean_value.size == 0:
            mean_value = no_data_value
            print("Warning, check z interval and grid definition zmn!")
        prior.append(mean_value)
    depth_average = np.ones(nx*ny*nz)*np.repeat(np.array(prior),nx*ny)
    return depth_average


def save_train_prior_3D(train,k,WD_prior,WD_data,
                        grid_params,radius_xy,radius_z):
    for i in tqdm(range(k)):
        print("-- calculate depth dependent average")
        depth_average = make_prior(train[i],grid_params)
        print("-- calculate moving window")
        ma_temp = moving_average_3d_kdtree(train[i][:,:3],train[i][:,3],grid_params,radius_xy,radius_z)
        grid_data, grid_n_neighbours = ma_temp
        print("-- join depth dependent average and moving window")
        fusion_grid = moving_window_depth_average_fusion(grid_data, grid_n_neighbours,
                                                        depth_average)
        print("-- save train file")
        priormfl = str(WD_prior)+str(i) + ".gslib"
        save_prior(fusion_grid,priormfl)
        print("-- save prior file")#???
        file_path = str(WD_data) +str(i)+".gslib"
        PANDASDF2GSLIBGeoEAS(pd.DataFrame(train[i],columns=["X","Y","Z","VR"]),file_path)
        
def adjust_grid_values(grid_data, grid_n_neighbours, threshold=200):
    """
    Adjusts values in grid_data based on the number of neighbours in grid_n_neighbours.

    Args:
    grid_data (list): List of averaged values from moving_average_3d_kdtree.
    grid_n_neighbours (list): List of neighbour counts from moving_average_3d_kdtree.
    threshold (int, optional): Minimum number of neighbours required to keep the original value in grid_data.

    Returns:
    list: Modified grid_data where entries with neighbours below the threshold are set to 0.5.
    """
    adjusted_data = grid_data.copy()  # Make a copy to avoid modifying the original data
    for index, value in enumerate(grid_n_neighbours):
        if value < threshold:
            adjusted_data[index] = 0.5
    return adjusted_data
        
        
def save_train_prior_3D_HPGL(train,
                        grid_params,radius_xy,radius_z):
    
        print("-- calculate depth dependent average")
        # depth_average = np.ones(123*692*210) * 0.5
        depth_average = make_prior(train,grid_params)
        print("-- calculate moving window")
        
        prior, prior_neighbours = moving_average_3d_kdtree(train[:,:3],train[:,3],grid_params,radius_xy,radius_z)
        
        WD = "D:/Wang/U9/WD/"
        
        priormfl = str(WD) + "save/prior_k.gslib"
        
        depth_average = np.ones(123*692*210) * 0.5
        
        fusion_grid = moving_window_depth_average_fusion(prior, prior_neighbours,
                                                          depth_average)
        
        save_prior(fusion_grid, priormfl)

        prior = GSLIBGeoEAS2PANDASDF(priormfl)
        
        prior.p0[prior.p0 < 0.1] = 0.1
        prior.p1[prior.p1 < 0.1] = 0.1
        prior.p0[prior.p0 > 0.9] = 0.9
        prior.p1[prior.p1 > 0.9] = 0.9
        
        
        
        prior_tohpgl_p0 = prior.p0.values
        
        final_p0 = adjust_grid_values(prior_tohpgl_p0, prior_neighbours, threshold=200)
        
        file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p0.npy'
        # file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p0_500.npy'
        np.save(file_path, final_p0)
        prior_tohpgl_p1 = prior.p1.values
        
        final_p1 = adjust_grid_values(prior_tohpgl_p1, prior_neighbours, threshold=200)
        
        file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p1.npy'
        # file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p1_500.npy'
        np.save(file_path, final_p1)

    
        
def save_train_3D(train,k,WD_data):
    for i in tqdm(range(k)):
        # print("-- calculate moving window")
        # ma_temp = moving_average_3d_kdtree(train[i][:,:3],train[i][:,3],grid_params,radius_xy,radius_z)
        # grid_data,_ = ma_temp
        print("-- save train file")
        file_path = str(WD_data) +str(i)+".gslib"
        PANDASDF2GSLIBGeoEAS(pd.DataFrame(train[i],columns=["X","Y","Z","VR"]),file_path)



def readsisim(data_file):
    with open(data_file) as file:
        data = file.read()
    test =data[data.find("Simulated Value\n")+17:]   
    array = np.array(test.split(), dtype=float)
    return array

def save_npy(WD_outfl,k):
    for i in range(k):
        output = str(WD_outfl) + str(i) + ".out"
        arr = readsisim(output)
        file = str(WD_outfl) + str(i) + ".npy"
        np.save(file,arr)
        print(i)
        
def GSLIBGeoEAS2PANDASDF(data):

    """Convert GSLIB Geo-EAS file to a 1D or 2D numpy ndarray for use with
    Python methods
    :param data_file: file name
    :param kcol: TODO
    :param nx: shape along x dimension
    :param ny: shape along y dimension
    :return: ndarray, column name
    """

    #data_file = open(data).read()

    with open(data) as f:
        head = [next(f) for _ in range(2)]                  # read first two lines
        line2 = head[1].split()                             # get number of columns
        ncol = int(line2[0])                                # get the number of columns
        col_name = [next(f).strip() for _ in range(ncol)]   # read over the column names
        data = np.loadtxt(f, skiprows=0)                    # read table
        df = pd.DataFrame(data)                             # create data frame
        df.columns = col_name                               # assign column names

    return(df)
        
# def save_train_3D(train,k,nsim, WD_data):
#         print("-- save train file")
#         for i in range(k):
#             print(str(i+1)+"/"+str(k))
#             for j in tqdm(range(nsim)):
#                 file_path = str(WD_data) + str(i) + "_" + str(j) + ".gslib"
#                 # Get the total number of rows in the array
#                 num_rows = train[i].shape[0]
#                 # Calculate the number of rows to save (80% of the total rows)
#                 num_rows_to_save = int(0.7 * num_rows)
#                 # Generate random indices to select rows
#                 random_indices = np.random.choice(num_rows, size=num_rows_to_save, replace=False)
#                 # Select the random rows from the array
#                 random_coords_t = train[i][random_indices]
#                 PANDASDF2GSLIBGeoEAS(pd.DataFrame(random_coords_t,columns=["X","Y","Z","VR"]),file_path)
def prior_2_gslib(df,gslib_output):
    # Creating the output file containing the composited points in compliant with GSLIB format
    f = open(gslib_output,'w',encoding='utf8')

    # Writing the title or name of the output file at the first row
    f.write(gslib_output)
    f.write('\n')
    # Extracting 'n' number of properties to be modeled (number of columns) at the second row
    #used?# property_number =                    # Defining quantity of properties: number of columns - HoleID, from, to, and Geology
    f.write(str(len(df.columns)))                         # Writing number of variables: coordinates X,Y,Z plus quantity of properties
    f.write('\n')

    # Writing 'n' lines with the names of the properties at each row
    for n in range(len(df.columns)):
        f.write(df.columns[n])
        f.write('\n')
    # Writing the dataframe properties
    for i in tqdm(range(len(df))):
        f.write(' '.join([str(df.iloc[i,0]),str(df.iloc[i,1]), '\n']))
    f.close()
    
def rotate_point(point, origin, angle_degrees):
    """Rotate a point around another point by a given angle."""
    angle_radians = np.radians(angle_degrees)
    x, y = point
    ox, oy = origin
    qx = ox + np.cos(angle_radians) * (x - ox) - np.sin(angle_radians) * (y - oy)
    qy = oy + np.sin(angle_radians) * (x - ox) + np.cos(angle_radians) * (y - oy)
    return [qx, qy]


# def write_boolean_arrays_to_binary(arrays, file_path):
#     with open(file_path, 'wb') as f:
#         for array in arrays:
#             # Convert the boolean array to a bytearray
#             byte_data = bytearray(len(array) // 8 + 1)
#             for i, value in enumerate(array):
#                 if value:
#                     byte_data[i // 8] |= 1 << (7 - i % 8)

#             # Write the bytearray to the binary file
#             f.write(byte_data)
        
def write_boolean_array_to_binary(data, file_path):
    # Convert the 2D list to a NumPy array
    numpy_array = np.array(data)

    # Convert the NumPy array to a boolean array
    boolean_array = (numpy_array != 0)

    # Pack the boolean array into a uint8 array
    byte_data = np.packbits(boolean_array, axis=-1)

    # Write the uint8 array to the binary file
    with open(file_path, 'wb') as f:
        f.write(byte_data.tobytes())

# def read_boolean_array_from_binary(file_path):
#     # Read the binary data from the file
#     with open(file_path, 'rb') as f:
#         byte_data = f.read()

#     # Convert the bytearray back to a boolean array
#     data = []
#     for byte in byte_data:
#         for i in range(8):
#             data.append((byte >> (7 - i)) & 1)

#     return data
#Unse unpactbits is more effecient 
def makecolorvariogram(nst, df):
    for i in range(nst+1):
        if i == 0:  plt.fill_between(df.x, 0, df["c"+str(i)], alpha=0.8)
        else:      plt.fill_between(df.x, df["c"+str(i-1)], df["c"+str(i)], alpha=0.5)
    plt.ylim(0,1.2)


def read_boolean_array_from_binary(file_path):
    # Read the binary file into a uint8 array
    with open(file_path, 'rb') as f:
        byte_data = np.frombuffer(f.read(), dtype=np.uint8)

    # Unpack the uint8 array into a boolean array
    boolean_array = np.unpackbits(byte_data)

    return boolean_array

def k_fold_cross_validation_3D(k,test,WD_outfl,GSLIB_SISIM_params):

    grid_params = GSLIB_SISIM_params["grid_params"]
    # nsim1 = GSLIB_SISIM_params["nsim"]
    nsim1 = 100
    correctlist = []
    wronglist = []
    completelist =[]
    for i in range(k):
        # arr_sis = np.load(str(WD_outfl) + str(i) + ".npy")
        arr_sis = read_boolean_array_from_binary(str(WD_outfl) + str(i) + "binary.binary")
        arr_sis = np.array_split(arr_sis,nsim1)
        sis_mean = np.mean(arr_sis,axis=0)
        sis_etype = np.digitize(sis_mean,(0.5,1),right= True)
        loc1, check = checkgrid(test[i][:,0], test[i][:,1], test[i][:,2], 
                                grid_params,XY_option=False)
        stlist = test[i][:,3][check].astype(int)
        loc = loc1[check]
        for j in range(len(loc)):
            completelist.append(float(sis_mean[loc[j]]))
            if sis_etype[loc[j]] == stlist[j]:
                correctlist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
            if sis_etype[loc[j]] != stlist[j]:
                wronglist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
        
    return correctlist, wronglist, completelist

def calculate_classification_metrics_k_fold(k, test_hpgl, hpgl_validation_file_path, grid_params):
    nsim1 = 100
    Cor = 0  # Correct predictions
    Wrong = 0  # Incorrect predictions
    TP = 0  # True Positives
    TN = 0  # True Negatives
    FP = 0  # False Positives
    FN = 0  # False Negatives

    for i in range(k):
        arr_sis = read_boolean_array_from_binary(str(hpgl_validation_file_path) + str(i) + "binary.binary")
        arr_sis = np.array_split(arr_sis, nsim1)
        sis_mean = np.mean(arr_sis, axis=0)
        sis_etype = np.digitize(sis_mean, (0.5, 1), right=True)
        
        loc1, check = checkgrid(test_hpgl[i][:,0], test_hpgl[i][:,1], test_hpgl[i][:,2], 
                                grid_params, XY_option=False)
        stlist = test_hpgl[i][:,3][check].astype(int)
        loc = loc1[check]

        for j in range(len(loc)):
            if stlist[j] == sis_etype[loc[j]]:
                Cor += 1
                if stlist[j] == 1:
                    TP += 1
                elif stlist[j] == 0:
                    TN += 1
            if sis_etype[loc[j]] != stlist[j]:
                Wrong += 1
                if sis_etype[loc[j]] == 1:
                    FP += 1
                elif sis_etype[loc[j]] == 0:
                    FN += 1

    precision_sand = TP / (TP + FP) if (TP + FP) > 0 else 0
    precision_clay = TN / (TN + FN) if (TN + FN) > 0 else 0
    accuracy = (TP + TN) / (TP + TN + FP + FN) if (TP + TN + FP + FN) > 0 else 0

    return {
        'TP': TP,
        'FP': FP,
        'TN': TN,
        'FN': FN,      
        'precision_clay': precision_clay,
        'precision_sand': precision_sand,
        'precision_clay': precision_clay,
        'accuracy': accuracy

    }

def k_fold_cross_validation_3D_hpgl(k,test_hpgl,hpgl_validation_file_path,grid_params):

    nsim1 = 100
    correctlist = []
    wronglist = []
    completelist =[]
    for i in range(k):
        # arr_sis = np.load(str(WD_outfl) + str(i) + ".npy")
        arr_sis = read_boolean_array_from_binary(str(hpgl_validation_file_path) + str(i) + "binary.binary")
        arr_sis = np.array_split(arr_sis, nsim1)[:-1]
        sis_mean = np.mean(arr_sis,axis=0)
        sis_etype = np.digitize(sis_mean,(0.5,1),right= True)
        loc1, check = checkgrid(test_hpgl[i][:,0], test_hpgl[i][:,1], test_hpgl[i][:,2], 
                                grid_params,XY_option=False)
        stlist = test_hpgl[i][:,3][check].astype(int)
        loc = loc1[check]
        for j in range(len(loc)):
            completelist.append(float(sis_mean[loc[j]]))
            if sis_etype[loc[j]] == stlist[j]:
                correctlist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
            if sis_etype[loc[j]] != stlist[j]:
                wronglist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
        
    return correctlist, wronglist, completelist

def k_fold_cross_validation_3D_HPGL(k, test_hpgl, WD_outfl, grid_params):
    nsim1 = 100  # 
    correctlist_sand = []
    wronglist_sand = []
    correctlist_ton = []
    wronglist_ton = []
    completelist = []

    for i in range(k):
        #
        arr_sis = read_boolean_array_from_binary(str(WD_outfl) + str(i) + "binary.binary")
        arr_sis = np.array_split(arr_sis, nsim1)[:-1]
        sis_mean = np.mean(arr_sis, axis=0)
        sis_etype = np.digitize(sis_mean, (0.5, 1), right=True)

        loc1, check = checkgrid(test_hpgl[i][:,0], test_hpgl[i][:,1], test_hpgl[i][:,2], 
                                grid_params, XY_option=False)
        stlist = test_hpgl[i][:,3][check].astype(int)
        loc = loc1[check]

        for j in range(len(loc)):
            completelist.append(float(sis_mean[loc[j]]))
            # 
            if stlist[j] == 1:  # 
                if sis_etype[loc[j]] == stlist[j]:
                    correctlist_sand.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
                else:
                    wronglist_sand.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
            elif stlist[j] == 0:  # 
                if sis_etype[loc[j]] == stlist[j]:
                    correctlist_ton.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
                else:
                    wronglist_ton.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
        
    return (correctlist_sand, wronglist_sand, correctlist_ton, wronglist_ton, completelist)


def k_fold_cross_validation_3D_moving_average(k,test,train_prior ,GSLIB_SISIM_params):

    grid_params = GSLIB_SISIM_params["grid_params"]
    # nsim1 = GSLIB_SISIM_params["nsim"]
    correctlist = []
    wronglist = []
    completelist =[]
    for i in range(k):
        # read_gs GSLIB2
        priormfl = str(train_prior)+str(i) + ".gslib"
        df = GSLIBGeoEAS2PANDASDF(priormfl)
        sis_mean = np.array(df.p1)
        sis_etype = np.digitize(sis_mean,(0.5,1),right= True)
        loc1, check = checkgrid(test[i][:,0], test[i][:,1], test[i][:,2], 
                                grid_params,XY_option=False)
        stlist = test[i][:,3][check].astype(int)
        loc = loc1[check]
        for j in range(len(loc)):
            completelist.append(float(sis_mean[loc[j]]))
            if sis_etype[loc[j]] == stlist[j]:
                correctlist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
            if sis_etype[loc[j]] != stlist[j]:
                wronglist.append(float(abs(sis_mean[loc[j]]-0.5)+0.5))
        
    return correctlist, wronglist, completelist


def plot_combined_accuracy(correctlist_sand, wronglist_sand, correctlist_ton, wronglist_ton, correctlist, wronglist, combined=True):
    bins = np.linspace(0.5, 1.0, 26)
    cm = 1/2.54  # centimeters to inches conversion
    dpi = 300  # dots per inch for image resolution

    if combined:
        # Only plot combined calibration curves
        fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)
        datasets = [
            (correctlist_sand, wronglist_sand, "Sand", 'green', 'dotted'),
            (correctlist_ton, wronglist_ton, "Clay", 'red','dotted'),
            (correctlist, wronglist, "Combined", 'blue', '-')
        ]
        for correct, wrong, label, color, linestyle in datasets:
            correct_hist, _ = np.histogram(correct, bins)
            wrong_hist, _ = np.histogram(wrong, bins)
            probability_forecast = correct_hist / (correct_hist + wrong_hist)
            ax.step(bins[:-1]+0.02, probability_forecast, where='mid', label=f"{label} Forecast", color=color, linestyle = linestyle)
        
        ax.plot([0.5, 1], [0.5, 1], label="Reference Line", color='black', linestyle='--')
        ax.set_xlabel('Predicted Probability')
        ax.set_ylabel('Accuracy')
        ax.set_ylim(0.5, 1)
        ax.set_xlim(0.5, 1)
        ax.legend(loc='upper left')
        ax.set_title("Calibration Curves Comparison")
        plt.show()

    else:
        # Plot separate histograms and calibration curves for each dataset
        fig, axs = plt.subplots(3, 2, figsize=(25*cm, 28*cm), dpi=dpi)
        datasets = [
            (correctlist_sand, wronglist_sand, "Sand Dataset"),
            (correctlist_ton, wronglist_ton, "Clay Dataset"),
            (correctlist, wronglist, "Combined Dataset")
        ]
        for i, (correct, wrong, label) in enumerate(datasets):
            # Histogram plot
            axs[i, 0].hist([correct, wrong], bins, stacked=True, label=['Correct Predictions', 'Incorrect Predictions'], color=['green', 'red'])
            axs[i, 0].legend(loc='upper right')
            axs[i, 0].set_title(f'{label} - Predictions Histogram')
            axs[i, 0].set_xlabel('Predicted Probability')
            axs[i, 0].set_ylabel('Frequency')

            # Calibration curve
            correct_hist, _ = np.histogram(correct, bins)
            wrong_hist, _ = np.histogram(wrong, bins)
            probability_forecast = correct_hist / (correct_hist + wrong_hist)
            axs[i, 1].step(bins[:-1]+0.02, probability_forecast, where='mid', label="Probability Forecast", color='blue')
            axs[i, 1].plot([0.5, 1], [0.5, 1], label="Reference Line", color='black', linestyle='--')
            axs[i, 1].set_xlabel('Predicted Probability')
            axs[i, 1].set_ylabel('Accuracy')
            axs[i, 1].set_ylim(0.5, 1)
            axs[i, 1].set_xlim(0.5, 1)
            axs[i, 1].legend(loc='upper left')
            axs[i, 1].set_title(f'{label} - Accuracy')

        # Adjust layout to prevent overlap
        plt.tight_layout()
        plt.savefig('calibration_curves.svg', format='svg')
        
        plt.show()


def plot_accuracy_3D(correctlist,wronglist):
    bins = np.linspace(0.5, 1.0, 26)
    # correctlist_length = len(correctlist)
    # wronglist_length = len(wronglist)
    # min(correctlist_length, wronglist_length)
    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(15*cm, 8*cm), dpi=600)
    # Stacked

    plt.hist([correctlist, wronglist], bins, stacked=True, label=['Number of Correct Predictions', 'Number of Incorrect Predictions'])
    plt.legend()
    # plt.hist(correctlist,bins)
    # plt.hist(wronglist,bins)
    plt.show()
    
    c = np.histogram(correctlist,bins)[0]
    w = np.histogram(wronglist,bins)[0]

    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(8.8*cm,8.8*cm), dpi=600)
    # plt.scatter(bins[:-1]+0.025,c/(w+c), 
    #             label="Wahrscheinlichkeit Prognose", c="Black",marker="x")
    plt.step(bins, np.append(c/(w+c), c[-1]/(w[-1]+c[-1])), where='post', label="Probability forecast", color="Blue")
    
    plt.plot([0.5,1],[0.5,1], label="p(Refrence line)", c="Black")
    plt.xlabel('Predicted probability')
    plt.ylabel('Accuracy')
    plt.ylim(0.5,1)
    plt.xlim(0.5,1)
    plt.legend()
    plt.show()
    #plt.savefig(str(WD) + "save/Fig_7.svg")



def plot_accuracy_3D_sperate(correctlist, wronglist, title):
    bins = np.linspace(0.5, 1.0, 26)
    cm = 1/2.54  # centimeters to inches conversion
    dpi = 300  # dots per inch for image resolution

    # Plotting the histogram of correct and wrong predictions
    fig, ax = plt.subplots(figsize=(20*cm, 12*cm), dpi=dpi)
    ax.hist([correctlist, wronglist], bins, stacked=True, label=['Correct Predictions', 'Incorrect Predictions'], color=['green', 'red'])
    ax.legend(loc='upper right')
    ax.set_title(title)
    ax.set_xlabel('Predicted Probability')
    ax.set_ylabel('Frequency')
    plt.show()

    # Calculating probabilities for step plot
    correct_hist, _ = np.histogram(correctlist, bins)
    wrong_hist, _ = np.histogram(wronglist, bins)
    probability_forecast = correct_hist / (correct_hist + wrong_hist)

    # Plotting the step plot for probability forecast
    fig, ax = plt.subplots(figsize=(12*cm, 12*cm), dpi=dpi)
    ax.step(bins[:-1]+0.02, probability_forecast, where='mid', label="Probability Forecast", color='blue')
    ax.plot([0.5, 1], [0.5, 1], label="Reference Line", color='black', linestyle='--')
    ax.set_xlabel('Predicted Probability')
    ax.set_ylabel('Accuracy')
    ax.set_ylim(0.5, 1)
    ax.set_xlim(0.5, 1)
    ax.legend(loc='upper left')
    ax.set_title(title)
    plt.show()



def make_uniform_grid(grid_params):
    # create vtk
    grid = pv.UniformGrid()
    # Set the grid dimensions: shape because we want to inject our values on the
    #   POINT data
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    grid.dimensions = np.array((nx,ny,nz)) + 1 #values.shape
    # Edit the spatial reference
    grid.origin = (xmn-xsiz/2,ymn-ysiz/2,zmn-zsiz/2)  # The bottom left corner of the data set
    grid.spacing = (xsiz, ysiz, zsiz)  # These are the cell sizes along each axis
    return grid

# def to_3D(data_2D, grid_params ,nsim):
#     nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
#     loc = np.arange(0,nx*ny*nz,1)+1
#     iz = 1 + np.trunc((loc-1)/(nx*ny))
#     # iy = 1 + np.trunc( (loc-1-(iz-1)*nx*ny)/nx) 
#     # ix = loc - (iz-1)*nx*ny - (iy-1)*nx 
    
#     # ax = (ix-1)*xsiz+xmn-xsiz/2
#     # ay = (iy-1)*ysiz+ymn+ysiz/2
#     az = (iz-1)*zsiz+zmn+zsiz/2
    
#     # Model = np.c_[loc,ax,ay,az]
#     # Defining variables required for this calculation as arrays
#     ncell = nx*ny*nz
#     result = np.zeros((ncell, nsim))                                            # Variable to store binary coding
#     rs = data_2D.reshape(nsim,int(len(data_2D)/nsim))
#     # Defining binary coding for each model: 0 if elevation of block is higher than simulated elevation and 1 if lower
#     for iloc in loc-1:
#         result[iloc] = rs[:,iloc % (nx*ny)] > az[iloc]
#     return result
# ----------------------------------------------------------------------------
def to_3D(data_2D, grid_params ,nsim):
    nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
    az = np.arange(zmn,zmn+nz*zsiz,zsiz)
    rs = np.split(data_2D,nsim)
    result = []
    for i in tqdm(range(nsim)):
        result_temp = [] 
        for iloc in range(nx*ny):
            result_temp.append(rs[i][iloc] > az)
        arr = np.array(result_temp).flatten(order="F")
        result.append(arr)
    return result

def combine_sim(sisim_3D, sgsim_3D,nsim):
    nmodel = []
    for sim in tqdm(range(nsim)):        # assign sand and clay to tertiary
        sisim_3D
        qt = sgsim_3D[sim]*1
        one_sim = (qt - 1) * (-3)  # assign Quaternary (erosion layer)
        #print(sim)    
        one_sim[one_sim == 0] = sisim_3D[sim][one_sim == 0]+1
        nmodel.append(one_sim)
    return nmodel

def combine_sim_ohnelm(sisim_3D_ohnelm, sgsim_3D,nsim):
    nmodel = []
    for sim in tqdm(range(nsim)):        # assign sand and clay to tertiary
        sisim_3D_ohnelm
        qt = sgsim_3D[sim]*1
        one_sim = (qt - 1) * (-3)  # assign Quaternary (erosion layer)
        #print(sim)    
        one_sim[one_sim == 0] = sisim_3D_ohnelm[sim][one_sim == 0]+1
        nmodel.append(one_sim)
    return nmodel


def save_data(file, data):
    with open(file, 'wb') as file:
        pickle.dump(data, file)
def load_data(file_path):
    with open(file_path, 'rb') as file:
        loaded_data = pickle.load(file)
    return loaded_data


def evaluate_nmodel(nmodel, nsim):
    # probabilities
    print("-- calculate probabilities")
    nmodel = np.array(nmodel)
    p_clay = np.count_nonzero((nmodel == 1),axis=0) / nsim 
    p_sand = np.count_nonzero((nmodel == 2),axis=0) / nsim 
    p_gravel = np.count_nonzero((nmodel == 3),axis=0) / nsim
    print("-- calculate etype")
    c = ((p_clay >= p_sand) & (p_clay > p_gravel)) * 1
    s = ((p_sand > p_clay) & (p_sand >= p_gravel)) * 2
    g = ((s == 0) & (c == 0)) * 3
    s_type = np.vstack([c,s,g]).max(axis=0)
    fs_type = np.vstack([1 + p_sand,g]).max(axis=0)
    print("-- calculate variability and entropy")
    p_min = 0.0001
    p_c = p_clay.copy()
    p_c[p_c == 0] = p_min
    p_s = p_sand.copy()
    p_s[p_s == 0] = p_min
    p_g = p_gravel.copy()
    p_g[p_g == 0] = p_min

    var = p_clay.copy()
    var[s_type == 2] = p_sand[s_type == 2] 
    var[s_type == 3] = p_gravel[s_type == 3] 
    var = 1 - var
    
    entropy = -((p_c*np.log(p_c))+(p_s*np.log(p_s))+(p_g*np.log(p_g)))
    
    print("total_ent = " + str(np.mean(entropy)))
    return s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var

def make_model(grid_params,s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var, elev_3D):
    grid = make_uniform_grid(grid_params)
    grid.cell_data["s_type"] = s_type
    grid.cell_data["fs_type"] = fs_type
    grid.cell_data["p_clay"] = p_clay
    grid.cell_data["p_sand"] = p_sand
    grid.cell_data["p_gravel"] = p_gravel
    grid.cell_data["entropy"] = entropy
    grid.cell_data["var"] = var
    grid.cell_data["elev_3D"] = elev_3D
    return grid

def plot_VARMAP(varmap, VARMAP_params, plane="XY", range_xyz=[50, 50, 4], levels=None):
    (dxlag, dylag, dzlag, nxlag, nylag, nzlag) = (VARMAP_params["dxlag"], VARMAP_params["dylag"], VARMAP_params["dzlag"],
                                                  VARMAP_params["nxlag"], VARMAP_params["nylag"], VARMAP_params["nzlag"])
    
    x_axis_label = "X-Distance [m]"
    y_axis_label = "Y-Distance [m]"
    
    nx = nxlag*2 + 1
    ny = nylag*2 + 1
    nz = nzlag*2 + 1
    
    # Reshape varmap based on selected plane
    if plane == "XY":
        vslice = varmap[nzlag * ny * nx:(nzlag + 1) * ny * nx].reshape(ny, nx)
    elif plane == "XZ":
        vslice = varmap.reshape(nz, ny, nx)[nylag, :, :].reshape(nx, nz).T
        x_axis_label = "X-Distance [m]"
        y_axis_label = "Z-Distance [m]"
    elif plane == "YZ":
        vslice = varmap.reshape(nz, ny, nx)[:, nxlag, :].reshape(ny, nz).T
        x_axis_label = "Y-Distance [m]"
        y_axis_label = "Z-Distance [m]"

    fig, ax = plt.subplots(dpi=600, figsize=(6, 5))
    X = np.arange(-(nxlag) * dxlag, (nxlag + 1) * dxlag, dxlag)
    Y = np.arange(-(nylag) * dylag, (nylag + 1) * dylag, dylag)
    XX, YY = np.meshgrid(X, Y)
    
    cmap = 'plasma'  # Setting the colormap to 'viridis'
    if levels is not None:
        cs = ax.contourf(XX, YY, vslice, levels=levels, cmap=cmap)
    else:
        cs = ax.contourf(XX, YY, vslice, cmap=cmap)
    plt.colorbar(cs, ax=ax, label='Variogram value')

    ax.set_xlabel(x_axis_label, size=14)
    ax.set_ylabel(y_axis_label, size=14)
    plt.title(f'Variogram Map ({plane}-plane)')
    plt.show()


def plot_VARMAP_3d(varmap, VARMAP_params, plane="XY", range_xyz=[50,50,4], levels=False, view_option=False, view_elev=90, view_azim=-90):
    (dxlag, dylag, dzlag, nxlag, nylag, nzlag) = (VARMAP_params["dxlag"], VARMAP_params["dylag"], VARMAP_params["dzlag"],
                                                  VARMAP_params["nxlag"], VARMAP_params["nylag"], VARMAP_params["nzlag"])
    
    nx = nxlag*2+1
    ny = nylag*2+1

    # Assuming varmap is a 1D array, reshape it to extract the specific slice
    varmap_reshaped = varmap.reshape((nzlag*2+1, ny, nx))
    vslice = varmap_reshaped[nzlag, :, :]  # Middle slice in Z direction
    
    X = np.arange(-(nxlag) * dxlag, (nxlag+1) * dxlag, dxlag)
    Y = np.arange(-(nylag) * dylag, (nylag+1) * dylag, dylag)
    XX, YY = np.meshgrid(X, Y)
    
    window_size = 3
    smooth_vslice = np.convolve(vslice.flatten(), np.ones(window_size)/window_size, mode='same').reshape(vslice.shape)
    
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    surface = ax.plot_surface(XX, YY, smooth_vslice, cmap='plasma')
    
    if view_option:
        ax.view_init(elev=view_elev, azim=view_azim)

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    
    plt.title('3D Variomap')

    
    # Adding a color bar
    cbar = fig.colorbar(surface, ax=ax, pad=0.1)
    cbar.set_label('Value')

    plt.savefig('varmap_3d_plot.png', dpi=300)
    plt.show()

            

# # plots all 12 variograms separately but can not plot the all at once !!!V 
# def plot_vertical_variogram_variants(list_of_vertical_variograms):
    
#     (figsize,nxlag,nylag,nx,ny) = vario_params.values()
       
#     for i in range(len(list_of_vertical_variograms)):
        
#         varmap_vario = list_of_vertical_variograms[i]      
        
#         slice_X = varmap_vario[np.repeat(np.arange(nx*nxlag,nx*nx*nx,nx*nx),nx)+np.tile(np.arange(0,nx,1),nx)]         
#         slice_Y = varmap_vario[np.repeat(np.arange(nxlag,nx*nx*nx,nx*nx),nx)+np.tile(np.arange(0,nx,1)*nx,nx)]
#         slice_Z = varmap_vario[(nx*nx*nxlag):(nx*nx*(nxlag+1))]
        
#         for j, (slice_, slice_name) in enumerate(zip([slice_X, slice_Y, slice_Z], ["Slice X", "Slice Y", "Slice Z"])):
            
#             #plt.sca(axes[i, j])
#             plot_varmap_2D(vario_params, slice_, f"Level {i+1} Variogram 2D: {slice_name}")
    
#%%
# ----------------------------------------------------------------------------
# MAKETHEORVARGM
# ----------------------------------------------------------------------------
#%%   
def MAKETHEORVARGM(vnugget,vsill,vrange,vtype):         # vtype 1: exponential, 2: spherical, 3: gaussian
    """
    Function to make the variogram

    Parameters
    ----------
    vnugget : float
        discontinuity in the variogram at distaances less than the minimum data 
        spacing.
    vsill : float
        the upper limit or bound of the variogram where there is no more linear
        correlation.
    vrange : float
        a quantitative measure of the geological zone of influence of the sample.
    vtype : int
        type of variogram:
            1 - spherical   
            2 - exponential
            3 - gaussian

    Returns
    -------
    None.

    """
    vstep = vrange/100
    vmax = 3*vrange
    x = np.arange(0,vmax,vstep)
    if vtype == 1:                                  # spherical variogram
        x1 = x[x <= vrange]
        y1 = vnugget + (vsill-vnugget) * ((x1/vrange)*(1.5-.5*((x1/vrange)**2))) #
        x2 = x[x > vrange]
        y2 = (x2/x2) * vnugget + (vsill-vnugget)
        y = np.append(y1,y2,axis=0)
    if vtype == 2:                                  # exponential variogram
        y = vnugget + (vsill-vnugget) * (1-np.exp(-x/(vrange/3)))
    if vtype == 3:                                  # gaussian variogram
        y = vnugget + (vsill-vnugget) * (1-np.exp((-(x/(vrange/3))**2)))
    theovargm = pd.DataFrame(data = {'x' : list(x), 'y' : list(y)})
    return(theovargm)
#%%
# ----------------------------------------------------------------------------
# MAKEDHVARGM
# ----------------------------------------------------------------------------
#%%   
def MAKEDHVARGM(nugget, v1sill, v1range, v1type, v2sill, v2cyl, d):
    """
    Function to make the variogram

    Parameters
    ----------
    nugget : float
        discontinuity in the variogram at distaances less than the minimum data 
        spacing.
    v1sill : float
        the upper limit or bound of the variogram where there is no more linear
        correlation.
    v1range : float
        a quantitative measure of the geological zone of influence of the sample.
    v1type : int
        type of variogram:
            1 - spherical   
            2 - exponential
            3 - gaussian
 
    v2sill : float
    v2cyl : float
       model dened by a length a to the first peak (size of the underlying cyclic features)
     h : float
       d is the distance at which 95% of the hole effect is dampened out (the variance magnitude of the periodic component is then 5% of c)
         
      

    Returns
    -------
    pd.DataFrame
        Theoretical variogram.

    """
    
    vstep = v1range / 100
    vmax = 3 * v1range
    x = np.arange(0, vmax, vstep)

    # Base variogram
    if v1type == 1:  # spherical variogram
        x1 = x[x <= v1range]
        y1 = nugget + (v1sill-nugget) * ((x1/v1range)*(1.5-.5*((x1/v1range)**2))) #
        x2 = x[x > v1range]
        y2 = (x2/x2) * nugget + (v1sill-nugget)
        y_base = np.append(y1,y2,axis=0)
    elif v1type == 2:  # exponential variogram
        y_base = nugget + (v1sill - nugget) * (1 - np.exp(-x / (v1range / 3)))
    elif v1type == 3:  # gaussian variogram
        y_base = nugget + (v1sill - nugget) * (1 - np.exp(-(x / (v1range / 3)) ** 2))

    # Additional variogram for hole effect

    y_hole = v2sill * (1.0 - np.exp((-3 * x)/d) * np.cos( x * np.pi/ v2cyl))
    y_base += y_hole  # Combine base and hole effect

    theovargm = pd.DataFrame(data={'x': list(x), 'y': list(y_base)})

    return (theovargm)

# ----------------------------------------------------------------------------
#%%
# ----------------------------------------------------------------------------
# MAKEHEVARGM
# ----------------------------------------------------------------------------
#%%   
def MAKEHEVARGM(nugget, v1sill, v1range, v1type, v2sill, v2cyl):
    """
    Function to make the variogram

    Parameters
    ----------
    nugget : float
        discontinuity in the variogram at distaances less than the minimum data 
        spacing.
    v1sill : float
        the upper limit or bound of the variogram where there is no more linear
        correlation.
    v1range : float
        a quantitative measure of the geological zone of influence of the sample.
    v1type : int
        type of variogram:
            1 - spherical   
            2 - exponential
            3 - gaussian
 
    v2sill : float
    v2cyl : float
       model dened by a length a to the first peak (size of the underlying cyclic features)
      
    Returns
    -------
    pd.DataFrame
        Theoretical variogram.

    """
    
    vstep = v1range / 100
    vmax = 3 * v1range
    x = np.arange(0, vmax, vstep)

    # Base variogram
    if v1type == 1:  # spherical variogram
        x1 = x[x <= v1range]
        y1 = nugget + (v1sill-nugget) * ((x1/v1range)*(1.5-.5*((x1/v1range)**2))) #
        x2 = x[x > v1range]
        y2 = (x2/x2) * nugget + (v1sill-nugget)
        y_base = np.append(y1,y2,axis=0)
    elif v1type == 2:  # exponential variogram
        y_base = nugget + (v1sill - nugget) * (1 - np.exp(-x / (v1range / 3)))
    elif v1type == 3:  # gaussian variogram
        y_base = nugget + (v1sill - nugget) * (1 - np.exp(-(x / (v1range / 3)) ** 2))

    # Additional variogram for hole effect

    y_hole = v2sill * (1.0 - np.cos( x * np.pi/ v2cyl))
    y_base += y_hole  # Combine base and hole effect

    theovargm = pd.DataFrame(data={'x': list(x), 'y': list(y_base)})

    return (theovargm)

# ----------------------------------------------------------------------------
#%%
# ----------------------------------------------------------------------------
# MAKENESTEDVARGM_Major
# ----------------------------------------------------------------------------
#%%   

    
def MAKE_THEOR_VARGM(vnugget,vsill,vrange,vtype,vmax,vstep):         # vtype 1: exponential, 2: spherical, 3: gaussian
    """
    Function to make the variogram

    Parameters
    ----------
    vnugget : float
        discontinuity in the variogram at distaances less than the minimum data 
        spacing.
    vsill : float
        the upper limit or bound of the variogram where there is no more linear
        correlation.
    vrange : float
        a quantitative measure of the geological zone of influence of the sample.
    vtype : int
        type of variogram:
            1 - spherical   
            2 - exponential
            3 - gaussian

    Returns
    -------
    None.

    """
    # vstep = vrange/100
    # vmax = 3*vrange
    x = np.arange(0,vmax,vstep)
    if vtype == 1:                                  # spherical variogram
        x1 = x[x <= vrange]
        y1 = vnugget + (vsill-vnugget) * ((x1/vrange)*(1.5-.5*((x1/vrange)**2))) #
        x2 = x[x > vrange]
        y2 = (x2/x2) * vnugget + (vsill-vnugget)
        y = np.append(y1,y2,axis=0)
    if vtype == 2:                                  # exponential variogram
        y = vnugget + (vsill-vnugget) * (1-np.exp(-x/(vrange/3)))
    if vtype == 3:                                  # gaussian variogram
        y = vnugget + (vsill-vnugget) * (1-np.exp((-(x/(vrange/3))**2)))
    theovargm = pd.DataFrame(data = {'x' : list(x), 'y' : list(y)})
    return(theovargm) 

def MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep):
    
    vargm_list = []
    # Search hole effect, is there a hole effect?
    has_hole_effect = any(var[0] == 5 for var in nest)
    
    if not has_hole_effect:
        for i in range (nst):
            vargm_list.append(MAKE_THEOR_VARGM(c0, nest[i][1], nest[i][5], nest[i][0], vmax, vstep))
    else:
        for i in range (nst - 1):
            vargm_list.append(MAKE_THEOR_VARGM(c0, nest[i][1], nest[i][5], nest[i][0], vmax, vstep))
            
        x = np.arange(0, vmax, vstep)
        y_hole = nest[nst - 1][1] * (1.0 - np.cos( x * np.pi/ nest[nst - 1][5]))
        hole = pd.DataFrame(data={'x': list(x), 'y': list(y_hole)})
        vargm_list.append(hole)
            
    df = vargm_list[0]
    # if nst > 1:
    #     for i in range (nst-1):
    #         df.y = df.y + vargm_list[i+1].y
    # df.y = df.y + c0
    ###################################################
    df["c0"] = df.y * 0 + c0
    if nst > 1:
        for i in range (nst):
            df["c"+str(i+1)] = df["c"+str(i)] + vargm_list[i].y  
            df.y = df["c"+str(i+1)] 

    return df  
   
    
# # Example:
# nst = 3
# c0 = 0.15
# vmax = 500
# vstep = 1
# nest     = [[1, 0.15, 0, 0, 0, 30, 10, 6],
#             [2, 0.358, 0, 0, 0, 450, 150, 6], 
#             [5, 0.34, 0, 0, 0, 100, 160, 6]]

# # Plot with matplotlib using separate colors
# plt.plot(df['x'], df['y'], label='DataFrame 1', color='blue')


# # Add labels and a title
# plt.xlabel('X-axis Label')
# plt.ylabel('Y-axis Label')
# plt.title('Two DataFrames Plotted Together')

# # Show the legend
# plt.legend()

# # Display the plot
# plt.show()

# ----------------------------------------------------------------------------
#%%
# ----------------------------------------------------------------------------
# MAKENESTEDVARGM_Minor
# ----------------------------------------------------------------------------
#%%   

    
def MAKE_THEOR_VARGM(vnugget,vsill,vrange,vtype,vmax,vstep):         # vtype 1: exponential, 2: spherical, 3: gaussian
    """
    Function to make the variogram

    Parameters
    ----------
    vnugget : float
        discontinuity in the variogram at distaances less than the minimum data 
        spacing.
    vsill : float
        the upper limit or bound of the variogram where there is no more linear
        correlation.
    vrange : float
        a quantitative measure of the geological zone of influence of the sample.
    vtype : int
        type of variogram:
            1 - spherical   
            2 - exponential
            3 - gaussian

    Returns
    -------
    None.

    """
    # vstep = vrange/100
    # vmax = 3*vrange
    x = np.arange(0,vmax,vstep)
    if vtype == 1:                                  # spherical variogram
        x1 = x[x <= vrange]
        y1 = vnugget + (vsill-vnugget) * ((x1/vrange)*(1.5-.5*((x1/vrange)**2))) #
        x2 = x[x > vrange]
        y2 = (x2/x2) * vnugget + (vsill-vnugget)
        y = np.append(y1,y2,axis=0)
    if vtype == 2:                                  # exponential variogram
        y = vnugget + (vsill-vnugget) * (1-np.exp(-x/(vrange/3)))
    if vtype == 3:                                  # gaussian variogram
        y = vnugget + (vsill - vnugget) * (1 - np.exp(-(3 * x)**2 / vrange**2))
    theovargm = pd.DataFrame(data = {'x' : list(x), 'y' : list(y)})
    return(theovargm) 

def MAKENESTEDVARGM_Minor(nst, c0, nest, vmax,vstep):
    
    vargm_list = []
    # Search hole effect, is there a hole effect?
    has_hole_effect = any(var[0] == 5 for var in nest)
    
    if not has_hole_effect:
        for i in range (nst):
            vargm_list.append(MAKE_THEOR_VARGM(c0, nest[i][1], nest[i][6], nest[i][0], vmax, vstep))
    else:
        for i in range (nst - 1):
            vargm_list.append(MAKE_THEOR_VARGM(c0, nest[i][1], nest[i][6], nest[i][0], vmax, vstep))
            
        x = np.arange(0, vmax, vstep)
        y_hole = nest[nst - 1][1] * (1.0 - np.cos( x * np.pi/ nest[nst - 1][6]))
        hole = pd.DataFrame(data={'x': list(x), 'y': list(y_hole)})
        vargm_list.append(hole)
            
    df = vargm_list[0]
    # if nst > 1:
    #     for i in range (nst-1):
    #         df.y = df.y + vargm_list[i+1].y
    # df.y = df.y + c0
    ###################################################
    df["c0"] = df.y * 0 + c0
    if nst > 1:
        for i in range (nst):
            df["c"+str(i+1)] = df["c"+str(i)] + vargm_list[i].y  
            df.y = df["c"+str(i+1)] 

    return df  
   
   
    
# Example:
# nst = 3
# c0 = 0.15
# nest     = [[1, 0.15, 0, 0, 0, 30, 10, 6],
#             [2, 0.358, 0, 0, 0, 450, 150, 6], 
#             [2, 0.34, 0, 0, 0, 8000, 160, 6]]

# # Plot with matplotlib using separate colors
# plt.plot(df1['x'], df1['y'], label='DataFrame 1', color='blue')
# plt.plot(df2['x'], df2['y'], label='DataFrame 2', color='orange')

# # Add labels and a title
# plt.xlabel('X-axis Label')
# plt.ylabel('Y-axis Label')
# plt.title('Two DataFrames Plotted Together')

# # Show the legend
# plt.legend()

# # Display the plot
# plt.show()

# ----------------------------------------------------------------------------
#%%
# ----------------------------------------------------------------------------
# read_GAMV_GAM_VARMAP
# ----------------------------------------------------------------------------
#%%
def read_GAMV_GAM(filename,params):

    col_names = ["lag","avsepdist","svargval","npairs","meantail","meanhead","tailheadvar"]
    num_tables = params["ndir"]


    with open(filename) as f:
        content = f.readlines()
    n_lines = len(content)
    rows_per_table = (n_lines // num_tables) - 1
    
    tables = []
    for i in range(num_tables):
        start_line = i * rows_per_table + i 
        table = pd.read_csv(io.StringIO(''.join(content[start_line + 1:start_line + rows_per_table + 1])),
                            delim_whitespace=True, names=col_names)
        tables.append(table)
        
    return tables

    f.close()
    
def GSLIBGeoEAS2PANDASDF(data):

    """Convert GSLIB Geo-EAS file to a 1D or 2D numpy ndarray for use with
    Python methods
    :param data_file: file name
    :param kcol: TODO
    :param nx: shape along x dimension
    :param ny: shape along y dimension
    :return: ndarray, column name
    """

    #data_file = open(data).read()

    with open(data) as f:
        head = [next(f) for _ in range(2)]                  # read first two lines
        line2 = head[1].split()                             # get number of columns
        ncol = int(line2[0])                                # get the number of columns
        col_name = [next(f).strip() for _ in range(ncol)]   # read over the column names
        data = np.loadtxt(f, skiprows=0)                    # read table
        df = pd.DataFrame(data)                             # create data frame
        df.columns = col_name                               # assign column names

    return(df)

def read_VARMAP(outfl):
    
    # open varmap.out file and replace ********** with 0.
    with open(outfl, 'r') as f:
        lines = f.readlines()
    with open(outfl, 'w') as f:
        for line in lines:
            if '**********' in line:
                line = line.replace('**********', '        0.')
            f.write(line)
    
    varmap = GSLIBGeoEAS2PANDASDF(outfl)
    varmap_vario = np.array(varmap.variogram.replace(-999,np.nan))
    
    return varmap_vario
