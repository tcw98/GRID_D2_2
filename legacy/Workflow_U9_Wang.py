# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:17:32 2023

@author: Witty
"""
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
WD = "YOUR FOLDER/WD/"
HOME = YOUR FOLDER/WD/PYTHON/"
os.chdir(str(HOME))
sys.path.append(str(HOME))

from GSLIB_parfl_011 import (GSLIB_SGSIM, GSLIB_SISIM_LM, GSLIB_GAMV, GSLIB_GAM,
                             GSLIB_VARMAP, GSLIB_SISIM)
from UTILITIES_013_Wa import (MAKETHEORVARGM, MAKEDHVARGM, MAKEHEVARGM, read_GAMV_GAM, RUN_GSLIB, LFU_2_YOU, qnd_compositing,
                           domain_plus_buffer, detrend_2D, PANDASDF2GSLIBGeoEAS,
                           transform, readsisim, read_VARMAP, plot_VARMAP,
                           plot_2_variograms_q, post_processing_2D,
                           fit_2_exponential,summarize_variograms, fit_2_curve,
                           plot_result_2D, save_Cl_Sa_points, save_train,
                           plot_2_variograms_t, plot_1_variogram_t, plot_VARMAP_3d, 
                           split_array_3D, save_train_prior_3D, save_npy,
                           k_fold_cross_validation_3D, plot_accuracy_3D,
                           to_3D, combine_sim, evaluate_nmodel, make_model,
                           moving_average_3d_kdtree, make_prior, run_parallel,
                           moving_window_depth_average_fusion, save_prior,
                           split_array, k_fold_cross_validation, plot_KCV,
                           read_sisim_lm,rotate_point,write_boolean_array_to_binary, run_sisim_parallel, read_sisim, read_boolean_array_from_binary, MAKENESTEDVARGM_Major, moving_average_2d_kdtree
                           , MAKENESTEDVARGM_Minor, makecolorvariogram, GSLIBGeoEAS2PANDASDF, prior_2_gslib,
                           k_fold_cross_validation_3D_moving_average, run_hpgl, assign_points_2_grid, run_hpgl_parallel, k_fold_cross_validation_3D_HPGL, generate_scripts, run_hpgl_lm_parallel, len_wei_compositing) 


grid_params = {  'nx': 123, #+20
                 'ny': 692, #+20
                 'nz': 210,
                 'xsiz': 10,
                 'ysiz': 10,
                 'zsiz': 0.5,
                 'xmn': 688440.0, #-100
                 'ymn': 5332620.0, #-100
                 'zmn': 475.25}                                                 # Compositing interval/2


origin = [688000,5332000]
#point = [692657.3,5338082.9]
angle = 25.67

# ############################################################################
# 0 SETUP
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


survey1_path = str(WD) + "data/bohrungen_epsg25832_shp/bohrungen.shp"
survey2_path = str(WD) + "data/bohrungen_epsg25832_shp/bohrungen_stammdaten.csv"
from_to_path = str(WD) + "data/bohrungen_epsg25832_shp/bohrungen_schichten.csv"
# grid_params = {"nx": 225, "ny": 324, "nz": 210,
#                "xsiz": 25, "ysiz": 25, "zsiz": 0.5,
#                "xmn": 4464876.708008, "ymn": 5329868.822266, "zmn": 429.75}
# U9 UTM 32


# grid_params = {"nx": 50,"ny": 40,"nz": 90,
#                 "xsiz": 5,"ysiz": 5,"zsiz": 0.5,
#                 "xmn": 691500,"ymn": 5334700, "zmn": 475}

# ############################################################################
# 1 COMPOSITING
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

buffer = 500  # m

extent = domain_plus_buffer(grid_params, buffer)
rotation = angle, origin
SURVEY, FROM_TO = LFU_2_YOU(survey1_path, survey2_path,
                            from_to_path, extent, False, rotation)

comp_temp = qnd_compositing(SURVEY, FROM_TO, compositing_interval=0.5,
                            vtkoption=True, verticaloption=True,
                            qbasisoption=True, qbasistol=3)

tertiary_points = len_wei_compositing (SURVEY, FROM_TO, compositing_interval=0.5, verticaloption = True )

#qbasis from function
points, point_cloud, line, qbasis = comp_temp 

np.save(str(WD)+"save/points_rotated.npy", points)
np.save(str(WD)+"save/qbasis.npy", qbasis)

#--------------------- No manually check--------------
qbasis_no_check = np.load(str(WD)+"save/no_check_qbasis.npy",allow_pickle=True)

#################################start form here!!!############################(No compositing)

#result of compositing (version 08.12.2023)
qbasis = np.load(str(WD)+"save/qbasis.npy",allow_pickle=True)
points = np.load(str(WD)+"save/points_rotated.npy")


# save and load point_cloud
line.save(str(WD)+"line.npy")
point_cloud_array = np.array(point_cloud)
np.save(str(WD) + "point_cloud.npy", point_cloud_array)

loaded_point_cloud_array = np.load(str(WD) + "point_cloud.npy")
point_cloud = loaded_point_cloud_array.astype(bool)

# ############################################################################
# 2 PLEISTOCENE
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# ############################################################################
# 2.1 Detrending and Transformation
# ############################################################################

result = detrend_2D(qbasis_no_check, grid_params, moving_window=True,
                    radius_xy=[300, 300], tolerance=20,
                    outlier_tolerance = 5 )


coords_no_check, coords_d_no_check, grid_data_no_check = result

np.save(str(WD) + "save/detrend_2D.npy", result)


mask = ~np.isnan(coords_d[:,2].astype(float))
coords, coords_d = coords[mask], coords_d[mask]

mask = ~np.isnan(coords_d_no_check[:,2].astype(float))
coords_no_check, coords_d_no_check = coords_no_check[mask], coords_d_no_check[mask]

#--------------------start form here!!!---------------------(No compositing/Detrending)

#result of detrending/compositing (version 08.12.2023) outlierr 3
#result of detrending/compositing (version 18.12.2023) outlier 99

# coords, coords_d, grid_data = np.load(str(WD)+"save/detrend_2D.npy",allow_pickle=True)
# coords, coords_d, grid_data = np.load(str(WD)+"save/detrend_2D_no_check.npy",allow_pickle=True)
# np.save(str(WD)+"save/no_check_qbasis_out99.npy", result)
# np.save(str(WD)+"save/no_check_qbasis_out5.npy", result)
# coords_no_check, coords_d_no_check, grid_data_no_check = np.load(str(WD)+"save/no_check_qbasis_out99.npy",allow_pickle=True)
# coords_no_check, coords_d_no_check, grid_data_no_check = np.load(str(WD)+"save/no_check_qbasis_out5.npy",allow_pickle=True)

coords_t, trans_method, trans_params = transform(coords_d, "Johnson")
coords_t_no_check, trans_method, trans_params = transform(coords_d_no_check, "Johnson")
coords_t_no_check_tol, trans_method, trans_params = transform(coords_d_no_check_tol, "Johnson")
np.save(str(WD) + "save/trans_params.npy",
        trans_params, allow_pickle=True)
# np.save(str(WD) + "save/no_check_trans_params.npy",        trans_params, allow_pickle=True)

file_path = str(WD) + "save/qbasis.gslib"
# file_path = str(WD) + "save/qbasis_no_check.gslib"


PANDASDF2GSLIBGeoEAS(pd.DataFrame(coords_t_no_check, columns=["X", "Y", "Z"]), str(WD) + "save/qbasis_no_check_300.gslib")


# ############################################################################
# 2.2 Variography
# ############################################################################
# ----------------------------------------------------------------------------
# 2.2.1 VARMAP
# ----------------------------------------------------------------------------

varmapexe = str(WD) + "gslib90/varmap.exe"
varmappar = str(WD) + "parfiles/varmap.par"
data = str(WD) + "save/qbasis_no_check_300.gslib"  # file with data
outfl = str(WD) + "save/varmap.out"  # file for variogram output

# Input parameters for running the GSLIB_VARMAP function for Tertiary variograms
VARMAP_params = {
    "nx": 20, "ny": 20, "nz": 1,
    "xsiz": 5.0, "ysiz": 5.0, "zsiz": 1,
    "icolx": 1, "icoly": 2, "icolz": 0,
    # if =0: columns for x,y, z coordinates
    "nvar": 1, "ivar1": 3,
    "tmin": -1.0e21, "tmax": 1.0e21, "igrid": 0,
    "nxlag": 50, "nylag": 50, "nzlag": 0,
    "dxlag": 5, "dylag": 5, "dzlag": 1,
    "minpairs": 2, "standardize": 1, "nvarg": 1,
    "ivtail": 1, "ivhead": 1, "ivtype": 1, "cut": 0.5, }

GSLIB_VARMAP(varmapexe, varmappar, data, outfl, VARMAP_params)
RUN_GSLIB(varmapexe, varmappar)
# 3D darstellen
varmap = read_VARMAP(outfl)
levels = np.arange(0, 1.5, 0.1)
levels = 3-np.log(np.arange(20, 1, -1))
plot_VARMAP(varmap, VARMAP_params, range_xyz=[102, 54, 4], levels=levels)



plot_VARMAP_3d(varmap, VARMAP_params, range_xyz=[102, 54, 4], levels=levels, view_option = True, view_elev=90, view_azim=-90)


# ----------------------------------------------------------------------------
# 2.2.2 GAMV
# ----------------------------------------------------------------------------



gamvexe = str(WD) + "gslib90/gamv.exe"  # GAMV executable
gamvpar = str(WD) + "parfiles/gamv.par"  # file for GAMV input
data = str(WD) + "save/qbasis_no_check_300.gslib"  # file with data

outfl = str(WD) + "save/gamv_Q.out"  # file for variogram output

GSLIB_GAMV_params = {"icolx": 1, "icoly": 2, "icolz": 0,
                     "nvar": 1, "ivar1": 3,
                     "tmin": -1.0e21, "tmax": 1.0e21,
                     "nlag": 100, "xlag": 5, "lagtol": 2.5,
                     "ndir": 2,  # azm,atol,bandwh,dip,dtol,bandwd
                     "xdir": [[90, 45, 100, 0, 0, 2], [0, 45, 100, 0, 0, 2]],
                     "standardize": 0,
                     "nvarg": 1, "ivtail": 1, "ivhead": 1, "ivtype": 1}

GSLIB_GAMV(gamvexe, gamvpar, data, outfl, GSLIB_GAMV_params)
# command_RUN_GSLIB(gamvexe,gamvpar)
RUN_GSLIB(gamvexe, gamvpar)
variograms = read_GAMV_GAM(outfl, GSLIB_GAMV_params)


var_list = []

for xlag in tqdm(range(5, 30, 5)):
    for lagtol in range(3, 18, 3):
        GSLIB_GAMV_params["xlag"] = xlag
        GSLIB_GAMV_params["lagtol"] = lagtol
        GSLIB_GAMV_params["xdir"][0][2] = xlag *3
        GSLIB_GAMV_params["xdir"][1][2] = xlag *3
        GSLIB_GAMV(gamvexe, gamvpar, data, outfl, GSLIB_GAMV_params)
        RUN_GSLIB(gamvexe, gamvpar)
        variograms = read_GAMV_GAM(outfl, GSLIB_GAMV_params)
        var_list.append(variograms)

# ---------------------------------
# ---------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
for i in range(len(var_list)):
    x = var_list[i][0]['avsepdist']
    y = var_list[i][0]['svargval']
    plt.scatter(x, y, c="Orange")
    x = var_list[i][1]['avsepdist']
    y = var_list[i][1]['svargval']
    plt.scatter(x, y, c="Blue")
    
plt.xlim(0, 500)
plt.ylim(0, 1.2)
plt.show()

# theovargm = MAKEDHVARGM(0.25, 0.65, 110, 2, 0.25, 110, 200)

#############################################################
# Plot vaiogram for q with auto fit
#############################################################
cm = 1/2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval'], 'o',
                    c='Orange',zorder=-1,label="E-W Experimental Variogram", alpha=1)
ax.plot(variograms[1]['avsepdist'],variograms[1]['svargval'], 'o',
                    c='Blue',zorder=-1,label="N-S Experimental Variogram", alpha=1)

color_north = 'Blue'
color_east = 'Orange'

#Auto fit
variogram_q_east = fit_2_curve(ax, variograms[0]['avsepdist'], variograms[0]['svargval'], color_east, "E-W Auto Fitted Curve")
variogram_q_north = fit_2_curve(ax, variograms[1]['avsepdist'], variograms[1]['svargval'], color_north, "N-S Auto Fitted Curve")

ax.set_xlabel('Distance [m]')
ax.set_ylabel('Horizontal Semivariance [-]')
# Move the legend to the right bottom corner
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
plt.savefig('variogram_q_auto_1000.png', dpi=300)
ax.set_yticks(np.arange(0, 1.3, 0.1))
plt.xlim(0, 500)
plt.ylim(0, 1.2)
plt.savefig('exp auto Q.svg', format='svg')
plt.show()

# variogram_q = summarize_variograms_2D(variogram_q_north, variogram_q_east)

# nst = 1
# c0 = 0.3
# nest = [[1, 0.977, 0, 0, 0, 186.786, 65.181, 6]]



#############################################################
# Plot vaiogram for q with manual fit
#############################################################
cm = 1/2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval'], 'o',
                    c='Orange',zorder=-1,label="E-W Experimental Variogram", alpha=1)
ax.plot(variograms[1]['avsepdist'],variograms[1]['svargval'], 'o',
                    c='Blue',zorder=-1,label="N-S Experimental Variogram", alpha=1)


vmax = 500
vstep = 1
nst = 4
c0 = 0
# ---------------------------------
# Test Variogram here (Without Hole Effect)
# --------------------------------

nest = [[1, 0.35, 0, 0, 0, 10, 10, 6], 
        [2, 0.3, 0, 0, 0, 40, 30, 6], 
        [2, 0.28, 0, 0, 0, 300, 110, 6],
        [5, 0.05, 0, 0, 0, 1.0e21, 80, 6]]
# # Outlier 3
# nest = [[1, 0.3, 0, 0, 0, 10, 10, 6], 
#         [2, 0.2, 0, 0, 0, 40, 45, 6], 
#         [2, 0.38, 0, 0, 0, 250, 140, 6],
#         [5, 0.07, 0, 0, 0, 1.0e21, 95, 6]]


theovargm1 = MAKENESTEDVARGM_Minor(nst, c0, nest, vmax, vstep)
plt.plot(theovargm1['x'], theovargm1['y'], c='red', zorder=-1, alpha=0.5, linewidth=2)
theovargm2 = MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep)
plt.plot(theovargm2['x'], theovargm2['y'], c='green', zorder=-1, alpha=0.5, linewidth=2)

legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10),
                  plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10),
                  plt.Line2D([0], [0], color='Red'),
                  plt.Line2D([0], [0], color='Green')]

plt.legend(title="Legend", handles=legend_handles, labels=["N-S", "E-W", "Fitted East", "Fitted Nord"])

ax.set_xlim(0, 500)
ax.set_ylim(0, 1.2)
plt.xticks(np.arange(0, 501, 50)) 
# ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
#ax.set_ylim(0, 4) #for raw data...
plt.savefig('exp Q.svg', format='svg')
plt.show()


# ############################################################################
# 2.3 Simulation
# ############################################################################
data = str(WD) + "save/qbasis_no_check_300.gslib"  # file with data
sgsimexe = str(WD) + "gslib90/sgsim.exe"  # GAMV executable
parfl = str(WD) + "parfiles/sgsim.par"  # file for variogram output
dbgfl = str(WD) + "save/sgsim.dbg"  # file for debugging output
outfl = str(WD) + "save/sgsim.out"  # the output grid is written to this file.
# data = str(WD) + "save/qbasis.gslib"  # file with data
# data = str(WD) + "save/qbasis_no_check.gslib"  # file with data
grid_params_2D = grid_params.copy()
grid_data = grid_data_no_check.copy()
grid_params_2D["nz"] = 1
GSLIB_SGSIM_params = {"icolx": 1, "icoly": 2, "icolz": 0,
                      "icolvr": 3, "icolwt": 0, "icolsec": 0,
                      "tmin": -3, "tmax": 3,
                      "itrans": 0, "transfl": "none.trn",
                      "ismooth": 0, "smthfl": "none.dat",
                      "icolvrsmthfl": 3, "icolwtsmthfl": 0,
                      "zmin": 0, "zmax": 40,
                      "ltail": 1, "ltpar": 0, "utail": 1, "utpar": 1,
                      "idbg": 1, "dbgfl": dbgfl,
                      "nsim": 100, "grid_params": grid_params_2D,
                      "seed": 12345, "ndmin": 1, "ndmax": 12, "ncnode": 12,
                      "sstrat": 1, "multgrid": 1, "nmult": 3, "noct": 0,
                      "radiushmax": 300, "radiushmin": 200, "radiusvert": 10,
                      "sang1": 0, "sang2": 0, "sang3": 0,
                      "covtab1": 60, "covtab2": 60, "covtab3": 60,
                      "ktype": 0, "rho": 0.6, "varred": 1,
                      "secfl": "none.dat", "icolsecfl": 0,
                      # it,cc,ang1,ang2,ang3, a_hmax, a_hmin, a_vert
                      "nst": nst, "c0": c0, "variogram": nest}


GSLIB_SGSIM(sgsimexe, data, parfl, outfl, GSLIB_SGSIM_params)
RUN_GSLIB(sgsimexe, parfl)


pp_temp = post_processing_2D(outfl, grid_data, GSLIB_SGSIM_params,
                             trans_method, trans_params)
sgsim, sgsim_dt, sgsim_dtbt = pp_temp


plot_result_2D(sgsim, sgsim_dt, sgsim_dtbt,
               GSLIB_SGSIM_params, coords_no_check)


# ############################################################################
# 2.4 Validation
# ############################################################################

# ----------------------------------------------------------------------------
# 2.4.1 CHECK GRID VARIOGRAMS
# ----------------------------------------------------------------------------
gamexe = str(WD) + "/gslib90/gam.exe"  # GAMV executable
gampar = str(WD) + "/parfiles/gam.par"  # file for GAM input
data = str(WD) + "save/sgsim.out"  # file with data
outfl = str(WD) + "save/gam.out"  # file with data

GSLIB_GAM_params = {"nvar": 1, "ivar1": 1,
                    "tmin": -1.0e21, "tmax": 1.0e21, "nsim": 0,
                    "grid_params": grid_params_2D,
                    "ndir": 2, "nlag": 50,
                    "xdir": [[0, 1, 0], [1, 0, 0]],
                    "standardize": 0, "nvarg": 1,
                    "ivtail": 1, "ivhead": 1, "ivtype": 1}
GSLIB_GAM(gamexe, gampar, data, outfl, GSLIB_GAM_params)
RUN_GSLIB(gamexe, gampar)

# ----------------------------------------------------------------------------
# To check the simulation result with experimental result
# ----------------------------------------------------------------------------
##################################################
# THIS PART MUST DO IN THE PART 2.2.2 GAMV!!!!!
##################################################
# var_list = []

# for xlag in tqdm(range(10, 50, 5)):
#     for lagtol in range(10, 20, 5):
#         GSLIB_GAMV_params["xlag"] = xlag
#         GSLIB_GAMV_params["lagtol"] = lagtol
#         GSLIB_GAMV(gamvexe, gamvpar, data, outfl, GSLIB_GAMV_params)
#         RUN_GSLIB(gamvexe, gamvpar)
#         variograms = read_GAMV_GAM(outfl, GSLIB_GAMV_params)
#         var_list.append(variograms)


# fig, ax = plt.subplots(figsize=(8, 5), dpi=330)

# # Result of GAMV (The real result)
# for i in range(len(var_list)):
#     x = var_list[i][0]['avsepdist']
#     y = var_list[i][0]['svargval']
#     ax.scatter(x, y, c="Blue", alpha=0.5, label="N-S Experimental Variogram" if i == 0 else "")

#     x = var_list[i][1]['avsepdist']
#     y = var_list[i][1]['svargval']
#     ax.scatter(x, y, c="Orange", alpha=0.5, label="E-W Experimental Variogram" if i == 0 else "")
    


fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
ax.plot(variograms[0]['avsepdist'],variograms[0]['svargval'], 'o',
                    c='Orange',zorder=-1,label="E-W Experimental Variogram", alpha=1)
ax.plot(variograms[1]['avsepdist'],variograms[1]['svargval'], 'o',
                    c='Blue',zorder=-1,label="N-S Experimental Variogram", alpha=1)

# Result of GAM (Reproduced Variogram)
for nsim in tqdm(range(100)):
    GSLIB_GAM_params["nsim"] = nsim + 1
    GSLIB_GAM(gamexe, gampar, data, outfl, GSLIB_GAM_params)
    RUN_GSLIB(gamexe, gampar)
    
    tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
    ax.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='blue', zorder=-1, alpha=0.5, linewidth=1, label="Reproduced N-S Variogram" if nsim == 0 else "")
    ax.plot(tables[1]['avsepdist'], tables[1]['svargval'], c='orange', zorder=-1, alpha=0.5, linewidth=1, label="Reproduced E-W Variogram" if nsim == 0 else "")

# Theoretical Variograms
theovargm1 = MAKENESTEDVARGM_Minor(nst, c0, nest, vmax, vstep)
plt.plot(theovargm1['x'], theovargm1['y'], c='red', zorder=-1, alpha=1, linewidth=3, label="Theoretical E-W")

theovargm2 = MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep)
plt.plot(theovargm2['x'], theovargm2['y'], c='green', zorder=-1, alpha=1, linewidth=3, label="Theoretical N-S")

plt.legend(title="Legend")

ax.set_ylim(0, 1.2)
plt.xticks(np.arange(0, 501, 50))
ax.set_xlabel('Distance')
ax.set_ylabel('Semivariance [-]')

plt.title('Variogram Analysis')
plt.savefig('Q exp mit re.svg', format='svg')
plt.show()

# ----------------------------------------------------------------------------
# 2.4.2 CROSS VALIDATION
# ----------------------------------------------------------------------------
import copy
# ---------copy form no check----
coords = copy.deepcopy(coords_no_check)
# ------------------------------
k = 5

train, train_t, test, test_t = split_array(k, coords,
                                           coords_t=None)

WD_data = str(WD) + "save/qbasis_"
dt_params = save_train(train, WD_data, k,
                       grid_params, moving_window=True,
                       radius_xy=[300, 300], tolerance=20,
                       trans_method="Johnson")


WD_outfl = str(WD) + "save/sgsim_"
WD_parfl = str(WD) + "parfiles/sgsim_"
for i in range(k):
    data = str(WD_data) + str(i) + ".gslib"
    outfl = str(WD_outfl) + str(i) + ".out"
    parfl = str(WD_parfl) + str(i) + ".par"
    GSLIB_SGSIM(sgsimexe, data, parfl, outfl, GSLIB_SGSIM_params)

run_parallel(sgsimexe, WD_parfl, k)

result = k_fold_cross_validation(k, test, WD_outfl,
                                 GSLIB_SGSIM_params, dt_params)
error_list, quantile_list, width_list = result

plot_KCV(width_list, quantile_list, GSLIB_SGSIM_params["nsim"])

np.save(str(WD) + "save/sg_error.npy", error_list)
np.save(str(WD) + "save/sg_quantile.npy", quantile_list)
np.save(str(WD) + "save/sg_width.npy", width_list)

#---------------------------------------------------
def calculate_overall_MAE(error_list):
    all_errors = np.concatenate(error_list)
    overall_MAE = np.mean(all_errors)
    return overall_MAE

print(calculate_overall_MAE(error_list))



# ############################################################################
# 3 MIOCENE
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ############################################################################
# 3.1 Variography
# ############################################################################



gamvexe = str(WD) + "gslib90/gamv.exe"  # GAMV executable
gamvpar = str(WD) + "parfiles/gamv.par"  # file for GAMV input
data = str(WD) + "save/tertiary.gslib"  # file with data
#data = str(WD) + "save/tertiarytest.gslib"  # file with data
outfl = str(WD) + "save/gamv.out"  # file for variogram output

# Cl_Sa_points = save_Cl_Sa_points(points,data)
Cl_Sa_points = save_Cl_Sa_points(tertiary_points,data)

# ############################################################################
# 3.1.2 GAMV
# ############################################################################

GSLIB_GAMV_params = {"icolx": 1, "icoly": 2, "icolz": 3,
                     "nvar": 1, "ivar1": 4,
                     "tmin": -1.0e21, "tmax": 1.0e21,
                     "nlag": 30, "xlag": 5, "lagtol": 2.5,
                     "ndir": 2,  # azm,atol,bandwh,dip,dtol,bandwd
                     "xdir": [[90, 30, 30, 0, 2, 1], [0, 30, 30, 0, 2, 1]],
                     "standardize": 0,
                     "nvarg": 1, "ivtail": 1, "ivhead": 1, "ivtype": 1}


GSLIB_GAMV(gamvexe, gamvpar, data, outfl, GSLIB_GAMV_params)
# command_RUN_GSLIB(gamvexe,gamvpar)

RUN_GSLIB(gamvexe, gamvpar)
variograms_t = read_GAMV_GAM(outfl, GSLIB_GAMV_params)

#--------------------------------------------------------------
# horizontal Auto fitted
#--------------------------------------------------------------



# cm = 1/2.54
# fig, ax = plt.subplots(figsize=(8, 5))
# ax.plot(variograms_t[0]['avsepdist'],variograms_t[0]['svargval'], 'o',
#                     c='Orange',zorder=-1,label="N-S Experimental Variogram", alpha=1)
# ax.plot(variograms_t[1]['avsepdist'],variograms_t[1]['svargval'], 'o',
#                     c='Blue',zorder=-1,label="E-W Experimental Variogram", alpha=1)
# # Define color
# color_north = 'Orange'
# color_east = 'Blue'
# variogram_t_north = fit_2_curve(ax, variograms_t[0]['avsepdist'], variograms_t[0]['svargval'], color_north, "N-S Auto Fitted Curve")
# variogram_t_east = fit_2_curve(ax, variograms_t[1]['avsepdist'], variograms_t[1]['svargval'], color_east, "E-W Auto Fitted Curve")
# ax.set_xlabel('Distance [m]')
# ax.set_ylabel('Horizontal Semivariance [-]')
# # Move the legend to the right bottom corner
# ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
# plt.show()

#--------------------------------------------------------------
# Manual fitting
#--------------------------------------------------------------

cm = 1/2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
ax.plot(variograms_t[0]['avsepdist'],variograms_t[0]['svargval'], 'o',
                    c='Orange',zorder=-1,label="E-W Experimental Variogram", alpha=1)
ax.plot(variograms_t[1]['avsepdist'],variograms_t[1]['svargval'], 'o',
                    c='Blue',zorder=-1,label="N-S Experimental Variogram", alpha=1)
# Define color
color_north = 'Orange'
color_east = 'Blue'

ax.set_xlabel('Distance [m]')
ax.set_ylabel('Horizontal Semivariance [-]')
# Move the legend to the right bottom corner
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
# plt.show()

#--------------------------------------------------------------
# MANUALLY Test here!!!
#--------------------------------------------------------------
vmax = 500
vstep = 1
nst = 1
c0 = 0.068


# --------------------------------
# nest = [[3, 0.14, 0, 0, 0, 100, 120, 6]]
# nest = [[2, 0.135, 0, 0, 0, 80, 90, 6]]

nest1 = [[2, 0.147, 0, 0, 0, 125, 150, 1.55]]
vmax = 500
vstep = 1
nst = 1
c01 = 0

# nest2 = [[3, 0.14, 0, 0, 0, 100, 120, 6]]
# vmax = 500
# vstep = 1
# nst = 1
# c02 = 0.07


#-----------------Gaussian with Exp.-------------------------------
theovargm_minor1 = MAKENESTEDVARGM_Minor(nst, c01, nest1, vmax, vstep)
theovargm_major1 = MAKENESTEDVARGM_Major(nst, c01, nest1, vmax, vstep)

ax.plot(theovargm_minor1['x'], theovargm_minor1['y'], c='green', zorder=-1, alpha=0.5, linewidth=2, label="Exp. E-W Fitted Variogram")
ax.plot(theovargm_major1['x'], theovargm_major1['y'], c='red', zorder=-1, alpha=0.5, linewidth=2, label="Exp. N-S Fitted Variogram")

# theovargm_minor2 = MAKENESTEDVARGM_Minor(nst, c02, nest2, vmax, vstep)
# theovargm_major2 = MAKENESTEDVARGM_Major(nst, c02, nest2, vmax, vstep)

# ax.plot(theovargm_minor2['x'], theovargm_minor2['y'], c='red', zorder=-1, alpha=0.5, linewidth=1, label="Gaussian N-S Fitted Variogram")
# ax.plot(theovargm_major2['x'], theovargm_major2['y'], c='green', zorder=-1, alpha=0.5, linewidth=1, label="Gaussian E-W Fitted Variogram")
# ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
# ax.set_xlim(0, 150)
# plt.xticks(np.arange(0, 300, 25)) 
# #ax.set_ylim(0, 4) #for raw data...
# plt.show()


# theovargm_minor = MAKENESTEDVARGM_Minor(nst, c0, nest, vmax, vstep)
# ax.plot(theovargm_minor['x'], theovargm_minor['y'], c='red', zorder=-1, alpha=0.5, linewidth=2, label="N-S Fitted Variogram")
# theovargm_major = MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep)
# ax.plot(theovargm_major['x'], theovargm_major['y'], c='green', zorder=-1, alpha=0.5, linewidth=2, label="E-W Fitted Variogram")

ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
ax.set_xlim(0, 150)
ax.set_ylim(0, 0.15)
# plt.xticks(np.arange(0, 300, 25)) 
# #ax.set_ylim(0, 4) #for raw data...
plt.savefig('exp Ter.svg', format='svg')
plt.show()

#--------------------------------------------------------------
# MANUALLY Test here!!!
#--------------------------------------------------------------

nest = [[2, 0.125, 0, 0, 0, 50, 50, 6]]
vmax = 500
vstep = 1
nst = 1


#-----------------Gaussian compare-------------------------------
# ax.plot(theovargm_minor2['x'], theovargm_minor2['y'], c='red', zorder=-1, alpha=0.5,linestyle= 'dashed', linewidth=1, label="Gaussian N-S Fitted Variogram")
# ax.plot(theovargm_major2['x'], theovargm_major2['y'], c='green', zorder=-1, alpha=0.5, linewidth=1,linestyle= 'dashed', label="Gaussian E-W Fitted Variogram")

#-----------------pure exp.-------------------------------

theovargm_minor = MAKENESTEDVARGM_Minor(nst, c0, nest, vmax, vstep)
ax.plot(theovargm_minor['x'], theovargm_minor['y'], c='red', zorder=-1, alpha=0.5, linewidth=2, label="N-S Fitted Variogram")
theovargm_major = MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep)
ax.plot(theovargm_major['x'], theovargm_major['y'], c='green', zorder=-1, alpha=0.5, linewidth=2, label="E-W Fitted Variogram")

ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
ax.set_xlim(0, 60)
# plt.xticks(np.arange(0, 300, 25)) 
# #ax.set_ylim(0, 4) #for raw data...
plt.show()

#--------------------------------------------------------------
# Vertical
#--------------------------------------------------------------

outfl_v = str(WD) + "save/gamv_v.out"
  # file for variogram output
GSLIB_GAMV_params["xlag"] = 0.5
GSLIB_GAMV_params["lagtol"] = 0.25
GSLIB_GAMV_params["ndir"] = 1
GSLIB_GAMV_params["xdir"] = [[0, 30, 30, 90, 2, 5]]
GSLIB_GAMV_params["standardize"] = 1
GSLIB_GAMV(gamvexe, gamvpar, data, outfl_v, GSLIB_GAMV_params)
RUN_GSLIB(gamvexe, gamvpar)
#direction is vertical
# run_GSLIB_background(gamvexe, gamvpar)
variogram_t_v = read_GAMV_GAM(outfl_v, GSLIB_GAMV_params)

#--------------------------------------------------------------
# Autofitted vertical....
#--------------------------------------------------------------
cm = 1/2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
color_vertical = 'green'
ax.plot(variogram_t_v[0]['avsepdist'],variogram_t_v[0]['svargval'], 'o', c='Green',zorder=-1,label="V Experimental Variogram", alpha=1)
variogram_t_vertical = fit_2_curve(ax, variogram_t_v[0]['avsepdist'], variogram_t_v[0]['svargval'], color_vertical, "V Auto Fitted Curve")
ax.set_xlabel('Distance [m]')
ax.set_ylabel('Vertical Semivariance [-]')
# Move the legend to the right bottom corner
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
plt.show()

#--------------------------------------------------------------
# Summarize the variogram and out put
#--------------------------------------------------------------

# variogram_t = summarize_variograms(variogram_t_north, variogram_t_east, variogram_t_vertical)

# ---------------------------------
#   MANUALLY Test Variogram here (Without Hole Effect)
# --------------------------------
cm = 1/2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi = 300)
color_vertical = 'green'
ax.plot(variogram_t_v[0]['avsepdist'],variogram_t_v[0]['svargval'], 'o', c='Green',zorder=-1,label="V Experimental Variogram", alpha=1)

vmax = 500
vstep = 0.1
nst = 1

c01v = 0

# ---------------------------------
# Test Variogram here (Without Hole Effect)
# --------------------------------
nest_v_E = [[2, 1, 0, 0, 0, 7, 100, 6]]

theovargm_v_E = MAKENESTEDVARGM_Major(nst, c01v, nest_v_E, vmax, vstep)
plt.plot(theovargm_v_E['x'], theovargm_v_E['y'], c='green', zorder=-1, alpha=0.5, linewidth=2, label="V. Fitted Variogram")


ax.set_xlabel('Distance [m]')
ax.set_ylabel('Vertical Semivariance [-]')
# Move the legend to the right bottom corner
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
ax.set_xlim(0, 15)
plt.savefig('exp Ter.svg', format='svg')
plt.show()

# ----------------------------------------------------------------------------
varmapexe = str(WD) + "gslib90/varmap.exe"
parfl = str(WD) + "parfiles/varmap_tertiary.par"
outfl = str(WD) + "save/varmap_tertiary.out"  # file for variogram output
datfl = data

GSLIB_VARMAP_params = {
    "nx": 50, "ny": 50, "nz": 30,
    "xsiz": 10, "ysiz": 10, "zsiz": 0.5,
    "icolx": 1, "icoly": 2, "icolz": 3,
    "nvar": 1, "ivar1": 4,
    "tmin": -1.0e21, "tmax": 1.0e21,
    "igrid": 0,
    "nxlag": 20, "nylag": 20, "nzlag": 20,
    "dxlag": 5.0, "dylag": 5.0, "dzlag": 0.5,
    "minpairs": 2,
    "standardize": 0,
    "nvarg": 1,
    "ivtail": 1, "ivhead": 1, "ivtype": 1,
    "cut": 0.5, }
GSLIB_VARMAP(varmapexe, parfl, data, outfl, GSLIB_VARMAP_params)
RUN_GSLIB(varmapexe, parfl)
# varmap = read_VARMAP(outfl, varmapexe, parfl)
# function to execute and read varmap 
varmap = read_VARMAP(outfl)                                                                                                                                                              
range_xyz = [30, 40, 5]
plot_VARMAP(varmap, GSLIB_VARMAP_params, plane="XY", range_xyz=range_xyz, levels=np.arange(0, 0.25, 0.01))
plot_VARMAP(varmap, GSLIB_VARMAP_params, plane="XZ", range_xyz=range_xyz, levels=np.arange(0, 0.25, 0.01))
plot_VARMAP(varmap, GSLIB_VARMAP_params, plane="YZ", range_xyz=range_xyz, levels=np.arange(0, 0.25, 0.01))

levels = np.arange(0, 1.5, 0.1)
levels = 3-np.log(np.arange(20, 1, -1))
plot_VARMAP_3d(varmap, GSLIB_VARMAP_params, range_xyz=[102, 54, 4], levels=levels, view_option = True, view_elev=21, view_azim=-70)

# ############################################################################
# 3.2 Generate Prior Information
# ############################################################################

from scipy.interpolate import RBFInterpolator

scaling_factor_z = 100
y,d = np.c_[Cl_Sa_points[:,0],Cl_Sa_points[:,1],Cl_Sa_points[:,2]*scaling_factor_z],Cl_Sa_points[:,3]
nx, ny, nz, xsiz, ysiz, zsiz, xmn, ymn,zmn = grid_params.values()
x_coordinates = [xmn + x * xsiz for x in range(nx)]
y_coordinates = [ymn + y * ysiz for y in range(ny)]
z_coordinates = [zmn*scaling_factor_z + z * zsiz*scaling_factor_z for z in range(nz)]

coordinates_grid = [(x, y, z) for z in z_coordinates for y in y_coordinates for x in x_coordinates  ]

rbfi = RBFInterpolator(y,d,neighbors = 100, kernel = "thin_plate_spline", smoothing = 10 )(coordinates_grid)

##############################################################################


prior = moving_average_3d_kdtree(Cl_Sa_points[:, :3], Cl_Sa_points[:, 3],
                                    grid_params, 500, 10)


# ma_temp = moving_average_3d_kdtree_parallel(Cl_Sa_points[:, :3], Cl_Sa_points[:, 3],
#                                     grid_params, 500, 10)

prior = np.array(prior)

# Set values less than 0.1 to 0.1
prior[prior < 0.1] = 0.1

# Set values greater than 0.9 to 0.9
prior[prior > 0.9] = 0.9


grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['SIM1'] = prior
grid_proior.plot()
grid_proior.save(str(WD) + 'save/0423.vtk')


prior_sand = prior.copy()
prior_clay = 1 - prior_sand


# Define the file paths
file_path_sand = r'D:\Wang\HPGL\data\prior_sand.npy'
file_path_clay = r'D:\Wang\HPGL\data\prior_clay.npy'

# Save prior_sand
np.save(file_path_sand, prior_sand)

# Save prior_clay
np.save(file_path_clay, prior_clay)





# depth_average = np.ones(103*672*210) * 0.5 # make_prior(Cl_Sa_points, grid_params)
depth_average = np.ones(123*692*210) * 0.5

fusion_grid = moving_window_depth_average_fusion(grid_data, grid_n_neighbours,
                                                  depth_average)
priormfl = str(WD) + "save/prior.gslib"
save_prior(fusion_grid, priormfl)

prior = GSLIBGeoEAS2PANDASDF(priormfl)

prior.p0[prior.p0 < 0.1] = 0.1
prior.p1[prior.p1 < 0.1] = 0.1
prior.p0[prior.p0 > 0.9] = 0.9
prior.p1[prior.p1 > 0.9] = 0.9



# ############################################################################
# TO HPGL (already done!!!!!!!)
# ############################################################################
prior_tohpgl_p0 = prior.p0.values
file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p0.npy'
# file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p0_500.npy'
np.save(file_path, prior_tohpgl_p0)
prior_tohpgl_p1 = prior.p1.values
file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p1.npy'
# file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p1_500.npy'
np.save(file_path, prior_tohpgl_p1)
# ############################################################################

priormfl = str(WD) + "save/prior12.gslib"
priormfl = str(WD) + "save/prior_tohpgl.gslib"
# Prior10 21.11.2023
# Prior11 14.12.2023
PANDASDF2GSLIBGeoEAS(prior, priormfl)

prior_2_gslib(prior, priormfl)


# ############################################################################
# 3.3.1 Simulation SISIM_LM
# ############################################################################

from UTILITIES_013_Wa import (make_uniform_grid,
                            run_sisim_lm_parallel)

sisim_lmexe = str(WD) + "gslib90/sisim_lm.exe"  # sisim executable
# priormfl = str(WD) + "save/prior11.gslib" # !!!! prior10 is from Nov. prior 11 form 15.12.2023
priormfl = str(WD) + "save/prior20.gslib" # !!!! IDW 1.8
parfl = str(WD) + "parfiles/sisim_lm.par"  # parameter file
datafl = str(WD) + "save/tertiary.gslib"  # file with data
dbgfl = str(WD) + "save/sisim_lm.dbg"  # file for debugging output
output = str(WD) + "save/sisim_lm.out"  # file for kriging output
GSLIB_SISIM_LM_params = {
    "vartype": 0, "ncat": 2, "catx": [0, 1, "", "", ""], "gpdfx": [0.5, 0.5, "", "", ""],
    "icolx": 1, "icoly": 2, "icolz": 3, "icolvr": 4,
    "tmin": 0, "tmax": 1.0e21, "zmin": 0, "zmax": 1,
    "ltail": 1, "ltpar": 1, "middle": 3, "midpar": 1, "utail": 1, "utpar": 1,
    "tabfl": datafl, "icolvrt": 0, "icolwtt": 0,
    "idbg": 1, "nsim": 1, "grid_params": grid_params, "seed": 12345,
    "ndmax":  12, "ncnode": 12,
    "sstrat": 1, "multgrid": 1, "nmult": 2, "noct": 0,
    "radiushmax": 100, "radiushmin": 100, "radiusvert": 10,
    "sang1": 90, "sang2": 0, "sang3": 0,
    "covtab1": 50, "covtab2": 50, "covtab3": 50, "mik": 0, "mikcat": 1,
    # nst,c0,it,cc,ang1,ang2,ang3,a_hmax,a_hmin,a_vert
    "variogram": [[nst, c0, nest],[nst, c0, nest]]
}
#In the source code, ndmax is a single variable int, not an array.

# GSLIB_SISIM_LM(sisim_lmexe,parfl,datafl,priormfl,
#                dbgfl,output,GSLIB_SISIM_LM_params)
# run_GSLIB_background(sisim_lmexe,parfl)
# RUN_GSLIB(sisim_lmexe,parfl)
# sisim_lm = readsisim(output)
# grid = make_uniform_grid(grid_params)
# grid.cell_data["sim"] = sisim_lm
# grid.plot()

WD_par = str(WD) + "parfiles/sisim_lm_"
WD_out = str(WD) + "save/sisim_lm_"
WD_dbg = str(WD) + "save/sisim_lm_"

# run_sisim_lm_parallel(sisim_lmexe, WD_par, WD_out, WD_dbg, datafl, priormfl,
#                       GSLIB_SISIM_LM_params,nsim_single= 5, n_parallel = 10,
#                       n_loops= 2)
run_sisim_lm_parallel(sisim_lmexe, WD_par, WD_out, WD_dbg, datafl, priormfl,
                      GSLIB_SISIM_LM_params
                      )
# -- generate parfiles
# 100%|██████████| 20/20 [00:00<00:00, 1250.58it/s]
# -- running simulations
# 100%|██████████| 4/4 [17:11:29<00:00, 15472.49s/it]  


sisim_3D = read_sisim_lm(WD_out)
write_boolean_array_to_binary(sisim_3D, str(WD_outfl) + 'binary.binary')
# np.save(str(WD) + 'save/sisim_3D.npy', sisim_3D)
# sisim_3D = np.load(str(WD)+'save/sisim_3D.npy')


# grid = make_uniform_grid(grid_params)
# grid.cell_data["sim"] = np.mean(sisim_3D, axis=0)
# grid.cell_data["fusion_grid"] = fusion_grid
# grid.cell_data["sim0"] = sisim_3D[0]
# grid.cell_data["grid_data"] = grid_data
# grid.cell_data["grid_n_neighbours"] = grid_n_neighbours
# grid.plot()
# grid.save(str(WD) + "save/sisim_lm.vtk")

# ############################################################################
# 3.3.2.1 HPGL Simulation SISIM_LM
# ############################################################################


# # ----------------- already done ----------------------
# # Prior information
# prior_tohpgl_p0 = prior.p0.values
# file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p0.npy'
# np.save(file_path, prior_tohpgl_p0)
# prior_tohpgl_p1 = prior.p1.values
# file_path = r'D:\Wang\HPGL\data\prior_tohpgl_p1.npy'
# np.save(file_path, prior_tohpgl_p1)
# # In hpgl_lm scripts it will load those file_paths as prior information...
# ----------------- make INK file -----------------------
INC, loc = assign_points_2_grid(Cl_Sa_points, grid_params)
# unique_elements, counts = np.unique(loc, return_counts=True)
# Filter to include only elements that appear more than once
# duplicates = unique_elements[counts > 1]

# Define the file name
file_name = 'D:/Wang/HPGL/data/SISIMINC.INC'
# Save the data to a text file with a header
np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')

# ----------------- define se and variogram -----------------------
# variogram = {
#     "type": 1, # 1 for exponentail
#     "ranges": variogram_t["ranges"],
#     "angles": variogram_t["angles"],
#     "sill": variogram_t["sill"],
#     "nugget": variogram_t["nugget"]
# }

# se = {
#     "radiuses": (250, 200, 15),
#     "max_neighbours": 12 # ncnode
# }
  

# ----------------- make py scripts -----------------------
origin_pyscript = "D:/Wang/HPGL/SISIM_LM.py"
numberofscripts = 5 # number of prozeses, numberofscripts x nsim_hpgl should equal to 100
name = "SISIM_LM"

pyscript_path = generate_scripts(origin_pyscript, numberofscripts, name)
# It will change the n_value, so will create result_from_hpgl_{n_value}_20.bin
nsim_hpgl = 20
pyscript_path= "D:/Wang/HPGL/SISIM_LM_"
# ----------------- run simulations -----------------------
run_hpgl_lm_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path = pyscript_path,   
                  inc_path = file_name, single_option=False, use_correlogram = True)
# run_hpgl_lm_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path = pyscript_path,   
#                   inc_path = file_name, single_option=False, use_correlogram = False)
# ----------------- read simulations -----------------------
file_path_tohpgl_result = r'D:\Wang\HPGL\data\result_from_hpgl.bin'
file_path_tohpgl_sim = "D:/Wang/HPGL/data/sim/result_sim"
sisim_3D_ohnelm=[]
for j in range(5):
    for i in range(nsim_hpgl):
         file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{j}_{i}.bin'
         sisim_3D_ohnelm.append(read_boolean_array_from_binary(file_path))


        
write_boolean_array_to_binary(sisim_3D_ohnelm, file_path_tohpgl_sim + 'binary.binary')
     


# ############################################################################
# 3.3.2.1 GSLIB Simulation SISIM
# ############################################################################

sisimexe = str(WD) + "gslib90/sisim.exe"  # sisim executable
parfl = str(WD) + "parfiles/sisim.par"  # parameter file
datafl = str(WD) + "save/tertiary.gslib"  # file with data
dbgfl = str(WD) + "save/sisim.dbg"  # file for debugging output
output = str(WD) + "save/sisim.out"  # file for kriging output
# GSLIB_SISIM_params = {
#     "vartype": 0, "ncat": 2, "catx": [0, 1, "", "", ""], "gpdfx": [0.5, 0.5, "", "", ""],
#     "icolx": 1, "icoly": 2, "icolz": 3, "icolvr": 4, "directik": "direct.ik", "icolsx": 1,
#     "icolsy": 2, "icolsz": 3, "icol": [4,5,6,7], "imbsim": 5, "bz": [1,2,3,4,5],
#     "tmin": 0, "tmax": 1.0e21, "zmin": 0, "zmax": 30, "ltail": 1, "ltbar": 0,
#     "middle": 1, "midpar": 1, "utail": 1, "utpar": 1, "tabfl": "cluster.dat", "icolvrt": 3,
#     "icolwtt": 0, "idbg": 1, "nsim": 1, "grid_params": grid_params, "seed": 12345,
#     "ndmax": 12, "ncnode": 12, "maxsec": 0, "sstrat": 1, "multgrid": 0, "nmult": 0,
#     "noct": 0, "radiushmax": 100, "radiushmin": 100, "radiusvert": 10, "sang1": 90,
#     "sang2": 0, "sang3": 0, "covtab1": 50, "covtab2": 50, "covtab3": 50, "mik": 0,
#     "mikcat": 1, "ktype": 1, "variogram": [[2, 0.05,[ 2, 0.1, 90, 0.0, 0.0, 70, 70, 7],
#                                            [2, 0.85, 90, 0.0, 0.0, 999, 999, 7],
#                                            [1, 0.3, 2, 1, 90, 0.0, 0.0, 70, 70, 7]]
# }
GSLIB_SISIM_params = {
    "vartype": 0, "ncat": 2, "catx": [0, 1, "", "", ""], "gpdfx": [0.5, 0.5, "", "", ""],
    "icolx": 1, "icoly": 2, "icolz": 3, "icolvr": 4, "directik": "direct.ik", "icolsx": 1,
    "icolsy": 2, "icolsz": 3, "icol": [4,5,6,7], "imbsim": 0, "bz": [1,2,3,4,5],
    "tmin": 0, "tmax": 1.0e21, "zmin": 0, "zmax": 1, "ltail": 1, "ltbar": 0,
    "middle": 1, "midpar": 1, "utail": 1, "utpar": 1, "tabfl": "cluster.dat", "icolvrt": 3,
    "icolwtt": 0, "idbg": 1, "nsim": 1, "grid_params": grid_params, "seed": 12345,
    "ndmax": 12, "ncnode": 12, "maxsec": 0, "sstrat": 1, "multgrid": 0, "nmult": 0,
    "noct": 5, "radiushmax": 100, "radiushmin": 100, "radiusvert": 10, "sang1": 0,
    "sang2": 0, "sang3": 0, "covtab1": 50, "covtab2": 50, "covtab3": 50, "mik": 0,
    "mikcat": 1, "ktype": 1, "variogram": [[nst, c0, nest],[nst, c0, nest]]
}


nst = 1
c0 = 0.5
nest = [[3, 1, 0, 0.0, 0.0, 100, 120, 7]]


############################ MIK ################################
GSLIB_SISIM_params = {
    "vartype": 0, "ncat": 2, "catx": [0, 1, "", "", ""], "gpdfx": [0.5, 0.5, "", "", ""],
    "icolx": 1, "icoly": 2, "icolz": 3, "icolvr": 4, "directik": "direct.ik", "icolsx": 1,
    "icolsy": 2, "icolsz": 3, "icol": [4,5,6,7], "imbsim": 0, "bz": [1,2,3,4,5],
    "tmin": 0, "tmax": 1.0e21, "zmin": 0, "zmax": 1, "ltail": 1, "ltbar": 0,
    "middle": 1, "midpar": 1, "utail": 1, "utpar": 1, "tabfl": "cluster.dat", "icolvrt": 3,
    "icolwtt": 0, "idbg": 3, "nsim": 1, "grid_params": grid_params, "seed": 12345,
    "ndmax": 12, "ncnode": 12, "maxsec": 0, "sstrat": 1, "multgrid": 0, "nmult": 0,
    "noct": 5, "radiushmax": 100, "radiushmin": 100, "radiusvert": 10, "sang1": 0,
    "sang2": 0, "sang3": 0, "covtab1": 50, "covtab2": 50, "covtab3": 50, "mik":1,
    "mikcat": 0, "ktype": 1, "variogram": [[nst, c0, nest],[nst, c0, nest]]
}

nst = 1
c0 = 0.5
nest = [[3, 1, 0, 0.0, 0.0, 100, 120, 7]] 

    # "icolwtt": 0, "idbg": 1, "nsim": 1, "grid_params": grid_params, "seed": 12345,
    # "ndmax": 40, "ncnode": 40, "maxsec": 0, "sstrat": 1, "multgrid": 0, "nmult": 0,
    # "noct": 0, "radiushmax": 100, "radiushmin": 100, "radiusvert": 10, "sang1": 90,
# GSLIB_SISIM(sisimexe,parfl,datafl,dbgfl,output,GSLIB_SISIM_params) 

WD_outfl = str(WD) + "save/sisim_test_"
WD_parfl = str(WD) + "parfiles/sisim_test_"
WD_dbg = str(WD) + "save/sisim_test_"

def run_sisim_parallel(sisimexe,WD_par,WD_out,WD_dbg,datafl,
                          GSLIB_SISIM_params):
    print("-- generate parfiles")
    GSLIB_SISIM_params["nsim"] = 1
    commands = []   
    for part in tqdm(range(5)):
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
            
nst = 1
c0 = 0.3
nest = [[3, 1, 0, 0.0, 0.0, 100, 120, 7]] 

WD_parfl = str(WD) + "parfiles/sisim_test_1.par"
RUN_GSLIB(sisimexe,WD_parfl)

run_sisim_parallel(sisimexe,WD_parfl,WD_outfl,WD_dbg,datafl,
                          GSLIB_SISIM_params)

# run_sisim_parallel_10(sisimexe,WD_parfl,WD_outfl,WD_dbg,datafl,
#                           GSLIB_SISIM_params)

#---------------Runtime SISIM----------------
# sstrat= 1, 10:54 start
#---------------Runtime SISIM----------------

# from datetime import datetime
# now = datetime.now()

# current_time = now.strftime("%H:%M:%S")

# print("Current Time =", current_time)
# RUN_GSLIB(sisimexe,parfl)

# now_after = datetime.now()

# current_time = now.strftime("%H:%M:%S")
# print("Current Time =", current_time)
#sstrat = 0: For 1 sim: 16:22:34 --- 18:34:13, about 135mins.
#sstrat = 1: For 1 sim: 12:43:56 --- 13:18:56, about 35mins.

sisim_3D_ohnelm = read_sisim(WD_outfl)
    # np.save(str(WD) + 'save/sisim_3D_ohnelm.npy', sisim_3D_ohnelm)
write_boolean_array_to_binary(sisim_3D_ohnelm, str(WD_outfl) + 'binary.binary')
# np.save(str(WD) + 'save/sisim_3D_ohnelm.npy', np.array(sisim_3D_ohnelm, dtype=object)) 


# restored_array = read_boolean_array_from_binary('D:/Wang/Result/SISIM_sstrat1_ndmax40_ncnode40/simulation/sisim_binary.binary')
# sisim_3D_ohnelm_restored = np.array([np.array(item) for item in np.split(restored_array, 100)])
# sisim_3D_ohnelm = sisim_3D_ohnelm_restored.astype(object)

# ############################################################################
# 3.3.2.1 HPGL Simulation SISIM
# ############################################################################
# for HPGL only single structure variogram
# type ranges(10,10,10) angles(0,0,0) sill nugget
# import json

# params_array = np.array(list(grid_params.values()))
# np.save('D:/Wang/U9/WD/save/grid_params.npy', params_array)

# tohpgl = {
#     "type": 0,
#     "ranges": (100, 100, 6),
#     "angles": (0, 0, 0),
#     "sill": 0.35,
#     "nugget": 0.0
# }


# ## (for nested struccture,but not possible....
# # nst_c0_array = np.array([nst, c0])
# # np.save('nst_c0.npy', nst_c0_array) 
# # nest_array = np.array(nest)
# # np.save('nest.npy', nest_array)

# # define search ellipsoid
# tohpgl_se = {
#     "radiuses": (100, 100, 20),
#     "max_neighbours": 12 # ncnode
# }

#

variogram = {
    "type": 1, # 1 for exponentail
    "ranges": [50, 50, 30],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0
}

se = {
    "radiuses": (100, 100, 7),
    "max_neighbours": 12 # ncnode
}
  

se = {
    "radiuses": (20, 20, 20),
    "max_neighbours": 12 # ncnode
}

INC, loc = assign_points_2_grid(Cl_Sa_points, grid_params)
# unique_elements, counts = np.unique(loc, return_counts=True)
# Filter to include only elements that appear more than once
# duplicates = unique_elements[counts > 1]

# Define the file name
file_name = 'D:/Wang/HPGL/data/SISIMINC.INC'
#file_name = 'D:/Wang/HPGL/data/SISIMINCtest.INC'
# Save the data to a text file with a header
np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')


#run_hpgl(grid_params, variogram, se, nsim, Cl_Sa_points)
nsim_hpgl = 1

pyscript_path= "D:/Wang/HPGL/SISIM_"


run_hpgl_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path,   # 5 py scripts "D:/Wang/HPGL/SISIM_i.py"
                  inc_path = file_name, single_option=False)


file_path_tohpgl_result = r'D:\Wang\HPGL\data\result_from_hpgl.bin'
file_path_tohpgl_sim = "D:/Wang/HPGL/data/result_sim"


sisim_3D_ohnelm=[]
for j in range(5):
    for i in range(nsim_hpgl):
         file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{j}_{i}.bin'
         sisim_3D_ohnelm.append(read_boolean_array_from_binary(file_path))
         
write_boolean_array_to_binary(sisim_3D_ohnelm, file_path_tohpgl_sim + 'binary.binary')


# ############################################################################
# 3.4 Validation
# ############################################################################

# ----------------------------------------------------------------------------
# 3.4.1 CHECK GRID VARIOGRAMS（SISIM_LM）
# ----------------------------------------------------------------------------
gamexe = str(WD) + "/gslib90/gam.exe"  # GAMV executable
gampar = str(WD) + "/parfiles/gam.par"  # file for GAM input
data = str(WD) + "save/sisim.out"  # file with data
outfl = str(WD) + "save/gam.out"  # file for variogram output

GSLIB_GAM_params = {"nvar": 1, "ivar1": 1,
                    "tmin": -1.0e21, "tmax": 1.0e21, "nsim": 1,
                    "grid_params": grid_params,
                    "ndir": 2, "nlag": 50,
                    "xdir": [[0, 1, 0], [1, 0, 0]],
                    "standardize": 0, "nvarg": 1,
                    "ivtail": 1, "ivhead": 1, "ivtype": 1}
tables = GSLIB_GAM(gamexe, gampar, data, outfl, GSLIB_GAM_params)

cm = 1 / 2.54
fig, ax = plt.subplots(figsize=(8.8 * cm, 5 * cm), dpi=600)


for i in range (20):
    data_temp = str(WD) + "save/sisim_lm_" + str(i) + ".out"
    for nsim in tqdm(range(5)):
        GSLIB_GAM_params["nsim"]= nsim+1
        GSLIB_GAM(gamexe, gampar, data_temp, outfl, GSLIB_GAM_params)
        RUN_GSLIB(gamexe, gampar)
        
        tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
        plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='Orange', zorder=-1, alpha=0.5)
        plt.plot(tables[1]['avsepdist'], tables[1]['svargval'], c='blue', zorder=-1, alpha=0.5)
    # tables[1].plot(kind='line',x='avsepdist',y='svargval',ax=ax,c='Blue',zorder=-1, alpha=1)
plt.show()

# ----------------------------------------------------------------------------
# 3.4.1 CHECK GRID VARIOGRAMS（SISIM）
# ----------------------------------------------------------------------------
gamexe = str(WD) + "/gslib90/gam.exe"  # GAMV executable
gampar = str(WD) + "/parfiles/gam.par"  # file for GAM input
data = str(WD) + "save/sisim_.out"  # file with data





GSLIB_GAM_params = {"nvar": 1, "ivar1": 1,
                    "tmin": -1.0e21, "tmax": 1.0e21, "nsim": 1,
                    "grid_params": grid_params,
                    "ndir": 1, "nlag": 50,
                    "xdir": [[0, 0, 1]],
                    "standardize": 0, "nvarg": 1,
                    "ivtail": 1, "ivhead": 1, "ivtype": 1}
tables = GSLIB_GAM(gamexe, gampar, data, outfl, GSLIB_GAM_params)



cm = 1 / 2.54
fig, ax = plt.subplots(figsize=(8.8 * cm, 5 * cm), dpi=600)

for i in range (1):
    # data_temp = str(WD) + "save/sisim_" + str(i) + ".out"
    data_temp = str(WD) + "save/sisim_test_0.out"
    for nsim in tqdm(range(1)):
        GSLIB_GAM(gamexe, gampar, data_temp, outfl, GSLIB_GAM_params)
        RUN_GSLIB(gamexe, gampar)
        
        tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
        # plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='Orange', zorder=-1, alpha=0.5)
        # plt.plot(tables[1]['avsepdist'], tables[1]['svargval'], c='blue', zorder=-1, alpha=0.5)
        plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='blue', zorder=-1, alpha=0.5)
    # tables[1].plot(kind='line',x='avsepdist',y='svargval',ax=ax,c='Blue',zorder=-1, alpha=1)
ax.set_xlim(0, 30)
ax.set_ylim(-0, 0.3)
plt.show()



data_temp = str(WD) + "save/sisim_0_11.out"      
GSLIB_GAM(gamexe, gampar, data_temp, outfl, GSLIB_GAM_params)
RUN_GSLIB(gamexe, gampar)
tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='blue', zorder=-1, alpha=0.5)


cm = 1 / 2.54
fig, ax = plt.subplots(figsize=(8.8 * cm, 5 * cm), dpi=600)

for i in range (1):
    # data_temp = str(WD) + "save/sisim_" + str(i) + ".out"
    data_temp = str(WD) + "save/sisim_test_5.out"
    for nsim in tqdm(range(5)):
        GSLIB_GAM_params["nsim"]= nsim+1
        GSLIB_GAM(gamexe, gampar, data_temp, outfl, GSLIB_GAM_params)
        RUN_GSLIB(gamexe, gampar)
        
        tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
        # plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='Orange', zorder=-1, alpha=0.5)
        # plt.plot(tables[1]['avsepdist'], tables[1]['svargval'], c='blue', zorder=-1, alpha=0.5)
        plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='blue', zorder=-1, alpha=0.5)
    # tables[1].plot(kind='line',x='avsepdist',y='svargval',ax=ax,c='Blue',zorder=-1, alpha=1)
ax.set_xlim(0, 30)
ax.set_ylim(-0, 0.3)
plt.show()


# ----------------------------------------------------------------------------
# 3.4.1 CHECK GRID VARIOGRAMS（HPGL）
# ----------------------------------------------------------------------------


gamexe = str(WD) + "/gslib90/gam.exe"  # GAMV executable
gampar = str(WD) + "/parfiles/gam.par"  # file for GAM input



GSLIB_GAM_params = {"nvar": 1, "ivar1": 1,
                    "tmin": -1.0e21, "tmax": 1.0e21, "nsim": 1,
                    "grid_params": grid_params,
                    "ndir": 3, "nlag": 20,
                    "xdir": [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
                    "standardize": 0, "nvarg": 1,
                    "ivtail": 1, "ivhead": 1, "ivtype": 1}

GSLIB_GAM_params = {"nvar": 1, "ivar1": 1,
                    "tmin": -1.0e21, "tmax": 1.0e21, "nsim": 1,
                    "grid_params": grid_params,
                    "ndir": 1, "nlag": 30,
                    "xdir": [[0, 0, 1]],
                    "standardize": 0, "nvarg": 1,
                    "ivtail": 1, "ivhead": 1, "ivtype": 1}

vargm_val_list = []  
variogram_t_lists = []


cm = 1 / 2.54
fig, ax = plt.subplots(figsize=(8, 5), dpi=600)


# ----------------------------------------------------------------------------
# Run and read
# ----------------------------------------------------------------------------
for i in range(4):
    for j in range(0, 20):
        # file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{i}_{j}.bin'
        # df = pd.DataFrame(read_boolean_array_from_binary(file_path), columns = ["simulations"])
        # PANDASDF2GSLIBGeoEAS(df, f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{i}_{j}.out')
        # data_temp = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{i}_{j}.out'
        outfl = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_vargm_{i}_{j}.out'
        # GSLIB_GAM(gamexe, gampar, data_temp, outfl, GSLIB_GAM_params)
        # RUN_GSLIB(gamexe, gampar)
        tables = read_GAMV_GAM(outfl, GSLIB_GAM_params)
        # ------------------------ Horizontal ----------------------------------------
        # plt.plot(tables[0]['avsepdist'], tables[0]['svargval'], c='blue', zorder=-1, alpha=0.5) #0 1 0 means y axe, which means the N-S direction
        # # # # ns_temp = fit_2_curve(ax, tables[0]['avsepdist'], tables[0]['svargval'], color = 'Orange', label = "N-S Auto Fitted Curve")
        # plt.plot(tables[1]['avsepdist'], tables[1]['svargval'], c='orange', zorder=-1, alpha=0.5)
        # # # ew_temp = fit_2_curve(ax, tables[1]['avsepdist'], tables[1]['svargval'], color = 'Blue', label = "E-W Auto Fitted Curve")
        # ------------------------ Vertical ----------------------------------------

        plt.plot(tables[2]['avsepdist'], tables[2]['svargval'], c='Green', zorder=-1, alpha=0.5)
        # vertical_temp = fit_2_curve(ax, tables[2]['avsepdist'], tables[2]['svargval'], color = 'Green', label = "Vertical Auto Fitted Curve")
        # ter_temp = [ns_temp, ew_temp]
        # vargm_val_list.append(ter_temp)
        # variogram_t_temp = summarize_variograms(vargm_val_list[i][0], vargm_val_list[i][1], variogram_t_vertical, print_option = False)
        # variogram_t_lists.append(variogram_t_temp)

# ------------------------ Plot experimental Variogram ----------------------------
# scale_factor = 0.25/0.14

# ax.plot(variograms_t[0]['avsepdist'],variograms_t[0]['svargval'] * scale_factor, 'o',
#                     c='orange',zorder=-1,label="E-W Experimental Variogram", alpha=1)
# ax.plot(variograms_t[1]['avsepdist'],variograms_t[1]['svargval'] * scale_factor, 'o',
#                     c='blue',zorder=-1,label="N-S Experimental Variogram", alpha=1)
ax.plot(variogram_t_v[0]['avsepdist'],variogram_t_v[0]['svargval'] * 0.25, 'o', c='Green',zorder=-1,label="V Experimental Variogram", alpha=1)




ax.set_xlabel('Distance [m]')
ax.set_ylabel('Semivariance [-]')

# nest = [[2, 0.25, 0, 0, 0, 100, 125, 6]]
# c0 = 0
# nst = 1
# vmax = 300
# vstep = 1
# # # ------------------------ Plot fitted variogram horizontal ----------------------------

# theovargm_major = MAKENESTEDVARGM_Major(nst, c0, nest, vmax, vstep)
# ax.plot(theovargm_major['x'], theovargm_major['y'], c='green', zorder=-1, alpha=0.6, linewidth=2,linestyle = '--', label="N-S Theoretical Variogram")
# theovargm_minor = MAKENESTEDVARGM_Minor(nst, c0, nest, vmax, vstep)
# ax.plot(theovargm_minor['x'], theovargm_minor['y'], c='red', zorder=-1, alpha=0.6, linewidth=2, linestyle = '--',label="E-W Theoretical Variogram")

# ------------------------ Plot fitted variogram vertical ----------------------------
c0 = 0
nest_v = [[2, 0.25, 0, 0, 0, 7, 120, 6]]
vstep = 0.1
nst = 1


theovargm_v = MAKENESTEDVARGM_Major(nst, c0, nest_v, vmax, vstep)
plt.plot(theovargm_v['x'], theovargm_v['y'], c='red', zorder=-1, alpha=0.6, linewidth=2, linestyle = '--', label="V Theoretical Variogram")

# ax.plot([], [], c='orange', label="E-W Reproduced Variogram")
# ax.plot([], [], c='blue', label="N-S Reproducted Variogram")
ax.plot([], [], c='Green', label="V Reproduced Variogram")

ax.set_xlim(0, 10)
ax.set_ylim(-0, 0.3)
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), shadow=True, ncol=1)
plt.savefig('reproduction.svg', format='svg')
plt.show()



df = pd.DataFrame(read_boolean_array_from_binary(data_temp), columns = ["simulations"])

grid= make_uniform_grid(grid_params)
grid.cell_data['SIM1'] = df[:17874360]
grid.plot()

grid_1= make_uniform_grid(grid_params)
grid_1.cell_data['SIM1'] = df[:17874360]
grid_1.plot()
# ----------------------------------------------------------------------------
# 3.4.2 K FOLD CROSS VALIDATION (for SISIM_LM)
# ----------------------------------------------------------------------------
k = 3

train, test = split_array_3D(k, Cl_Sa_points)
# save_data(str(WD)+"save/train_test_3D.pkl", [train,test])
# train, test = load_data(str(WD)+"save/train_test_3D.pkl")
WD_prior = str(WD) + "save/prior_"
WD_data = str(WD) + "data/tertiary_"

save_train_prior_3D(train, k, WD_prior, WD_data,
                    grid_params, radius_xy=500, radius_z=10)


WD_outfl = str(WD) + "save/sisim_lm_"
WD_parfl = str(WD) + "parfiles/sisim_lm_"
for i in range(k):
    datafl = str(WD_data) + str(i) + ".gslib"

    priormfl = str(WD_prior) + str(i) + ".gslib"

    WD_par = str(WD_parfl) + str(i) + "_"
    WD_out = str(WD_outfl) + str(i) + "_"
    WD_dbg = str(WD_outfl) + str(i) + "_"
    run_sisim_lm_parallel(sisim_lmexe,WD_par,WD_out,WD_dbg,datafl,priormfl,
                              GSLIB_SISIM_LM_params)

for i in range(k):
    WD_out = str(WD_outfl) + str(i) + "_"
    sisim_3D = read_sisim_lm(WD_out)

    write_boolean_array_to_binary(sisim_3D, str(WD_outfl) + str(i) + 'binary.binary') 


WD_outfl = str(WD) + 'save/sisim_lm_'
kcv_temp = k_fold_cross_validation_3D(k, test, WD_outfl, GSLIB_SISIM_LM_params)
correctlist, wronglist, completelist = kcv_temp
plot_accuracy_3D(correctlist, wronglist)

correctlist_save_path = str(WD) + "save/correctlist.npy"
wronglist_save_path = str(WD) + "save/wronglist.npy"
completelist_save_path = str(WD) + "save/completelist.npy"
np.save(correctlist_save_path, correctlist)
np.save(wronglist_save_path, wronglist)
np.save(completelist_save_path, completelist)

correctlist = np.load(str(WD) + "save/correctlist.npy" ,allow_pickle=True)
wronglist = np.load(str(WD) + "save/wronglist.npy" ,allow_pickle=True)

# ############################################################################
# 3.3.2.1 HPGL validation SISIM_LM
# ############################################################################

k = 3
train_hpgl, test_hpgl = split_array_3D(k, Cl_Sa_points)

variogram = {
    "type": 2, # 1 for exponentail
    "ranges": [100, 120, 7],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0.5
}


variogram = {
    "type": 1, # 1 for exponentail
    "ranges": [10, 12.5, 3],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0
}


variogram = {
    "type": 1, # 1 for exponentail
    "ranges": [10, 12.5, 14],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0
}


variogram = { # E-w  N-S V
    "type": 1, # 1 for exponentail
    "ranges": [14, 12.5, 14],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0
    }


variogram = { # E-w  N-S V
    "type": 1, # 1 for exponentail
    "ranges": [14, 12.5, 14],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0.5
    }
 

# ############################################################################
# SK prior
# ############################################################################


grid_data, grid_n_neighbours = moving_average_3d_kdtree(Cl_Sa_points[:, :3], Cl_Sa_points[:, 3],
                                    grid_params, 500, 10)


file_path = 'D:/Wang/Result/FInal SISIM test/Final model/SIS 50/SIS/prior_tohpgl_p1.npy'

grid_data = np.load(file_path)

grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['SIM'] = grid_data

grid_proior.save(str(WD) + 'save/IDW_new.vtk')



file_path = str(WD) + 'save/IDW neighbours.vtk' 

import pyvista as pv
grid = pv.read(file_path)

grid_n_neighbours = grid.cell_data['neighbours']

grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['neighbours'] = grid_n_neighbours
grid_proior.cell_data['infos']  = grid_data


grid_proior.save(str(WD) + 'save/IDW neighbours.vtk')



threshold = 200
# Option 1: "Hard Cut"
grid_data_opt1 = grid_data.copy()

for index, value in enumerate(grid_n_neighbours):

    if value < threshold:

        grid_data_opt1[index] = 0.5

grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['SIM'] = grid_data_opt1
grid_proior.save(str(WD) + 'save/IDW 100threshold.vtk')



# save "Hard Cut" SIM
file_path_p0 = 'D:/Wang/Result/FInal SISIM test/Final model/SIS 50/SIS/prior_tohpgl_p0.npy'
file_path_p1 = 'D:/Wang/Result/FInal SISIM test/Final model/SIS 50/SIS/prior_tohpgl_p1.npy'

IDW_p0 = np.load(file_path_p0)
IDW_p1 = np.load(file_path_p1)

p0 = IDW_p0.copy()
p1 = IDW_p1.copy()

for index, value in enumerate(grid_n_neighbours):

    if value < threshold:

        p0[index] = 0.5
        p1[index] = 0.5

np.save('D:/Wang/Result/FINAL FINAL/threeshold 200 IDW 1.8/prior_p0.npy', p0)
np.save('D:/Wang/Result/FINAL FINAL/threeshold 200 IDW 1.8/prior_p1.npy', p1)

# save "Hard Cut" VAL
file_path_p0 = ['D:/Wang/HPGL/Prior Val/prior_tohpgl_p0_a.npy', 'D:/Wang/HPGL/Prior Val/prior_tohpgl_p0_b.npy', 'D:/Wang/HPGL/Prior Val/prior_tohpgl_p0_c.npy']
file_path_p1 = ['D:/Wang/HPGL/Prior Val/prior_tohpgl_p1_a.npy', 'D:/Wang/HPGL/Prior Val/prior_tohpgl_p1_b.npy', 'D:/Wang/HPGL/Prior Val/prior_tohpgl_p1_c.npy']


for i in range (3):
    IDW_p0 = np.load(file_path_p0[i])
    IDW_p1 = np.load(file_path_p1[i])
    
    p0 = IDW_p0.copy()
    p1 = IDW_p1.copy()
    
    for index, value in enumerate(grid_n_neighbours):
    
        if value < threshold:
    
            p0[index] = 0.5
            p1[index] = 0.5

np.save('D:/Wang/Result/FINAL FINAL/threeshold 200 IDW 1.8/prior_p0.npy', p0)
np.save('D:/Wang/Result/FINAL FINAL/threeshold 200 IDW 1.8/prior_p1.npy', p1)


# Option 2: "Smooth Transition"
grid_data_opt2 = grid_data.copy()
factor = np.ones(np.shape(grid_n_neighbours))

# factor = np. exp((factor[grid_n_neighbours < threshold] / threshold)**2)
# grid_data_opt2 = factor * grid_data + (1-factor) * 0.5

for index, value in enumerate(grid_n_neighbours):

    if value < threshold:
        factor[index] = np.exp((value / threshold)**2)
        grid_data_opt2[index] = factor[index] * grid_data[index] + (1-factor[index]) * 0.5
        
grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['SIM'] = grid_data_opt2
grid_proior.save(str(WD) + 'save/IDWsmooth.vtk')




for i in range (3):
    run_SK(grid_params, train[i], "D:/Wang/HPGL/SK.py")




run_SK(grid_params, train_hpgl[0], "D:/Wang/HPGL/SK.py")

SK_pro = np.load("D:/Wang/HPGL/data/SK_RESULT.npy")

# Set values less than 0.1 to 0.1
SK_pro[SK_pro < 0.1] = 0.1

# Set values greater than 0.9 to 0.9
SK_pro[SK_pro > 0.9] = 0.9


grid_proior= make_uniform_grid(grid_params)
grid_proior.cell_data['SIM1'] = SK_pro
grid_proior.plot()

grid_proior.save(str(WD) + 'save/SK100.vtk')



Li_0_a = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_a.npy')
Li_0_b = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_b.npy')
Li_0_c = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_c.npy')
Li_1_a = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_a.npy')
Li_1_b = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_b.npy')
Li_1_c = np.load(r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_c.npy')




file_path_0_a = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_a.npy'
file_path_0_b = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_b.npy'
file_path_0_c = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_c.npy'
file_path_1_a = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_a.npy'
file_path_1_b = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_b.npy'
file_path_1_c = r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_c.npy'


np.save(file_path_0_a, Li_0_a)
np.save(file_path_0_b, Li_0_b)
np.save(file_path_0_c, Li_0_c)
np.save(file_path_1_a, Li_1_a)
np.save(file_path_1_b, Li_1_b)
np.save(file_path_1_c, Li_1_c)



file_paths = [
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_a.npy',
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_b.npy',
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p0_c.npy',
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_a.npy',
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_b.npy',
    r'D:\Wang\HPGL\Prior Val\Lineal\prior_tohpgl_p1_c.npy'
]


file_paths = [
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p0_a.npy',
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p0_b.npy',
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p0_c.npy',
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p1_a.npy',
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p1_b.npy',
    r'D:\Wang\HPGL\Prior Val\prior_tohpgl_p1_c.npy'
]


for file_path in file_paths:
    data = np.load(file_path)
    data = (data - 0.5) * 0.8 + 0.5
    new_file_path = file_path.replace("prior_tohpgl", "adjusted_prior_tohpgl")
    np.save(new_file_path, data)
# ----------------------------------------------------------------------------
# generate py scripts
# ----------------------------------------------------------------------------
origin_pyscript = 'D:/Wang/HPGL/Val/LM_test_aj_a.py'
generate_scripts(origin_pyscript, 5, 'LM_test_aj_a')
origin_pyscript = 'D:/Wang/HPGL/Val/LM_test_aj_b.py'
generate_scripts(origin_pyscript, 5, 'LM_test_aj_b')
origin_pyscript = 'D:/Wang/HPGL/Val/LM_test_aj_c.py'
generate_scripts(origin_pyscript, 5, 'LM_test_aj_c')

# ----------------------------------------------------------------------------
# Run HPGL with alrady exsits Prior a/b/c
# ----------------------------------------------------------------------------


se = {
    "radiuses": (20, 20, 20),
    "max_neighbours": 12 # ncnode
}



hpgl_validation_file_path = "D:/Wang/HPGL/data/validation/result_validation_"

# ------------------IDW 1.8 500-10 ------------------
pyscript_path_val= ["D:/Wang/HPGL/SISIM_LM_a_", "D:/Wang/HPGL/SISIM_LM_b_", "D:/Wang/HPGL/SISIM_LM_c_"]

# ------------------Gaussian 500-10 ------------------
pyscript_path_val = ["D:/Wang/HPGL/SISIM_LM_G_a_", "D:/Wang/HPGL/SISIM_LM_G_b_", "D:/Wang/HPGL/SISIM_LM_G_c_"] 

# ------------------Lineal advanced logic ------------------
pyscript_path_val = ["D:/Wang/HPGL/LM_test_a_", "D:/Wang/HPGL/LM_test_b_", "D:/Wang/HPGL/LM_test_c_"] 

# ------------------Lineal advanced logic with adjusted prior ------------------
pyscript_path_val = ["D:/Wang/HPGL/LM_test_aj_a_", "D:/Wang/HPGL/LM_test_aj_b_", "D:/Wang/HPGL/LM_test_aj_c_"] 


nsim_hpgl = 20

for i in range(k):
    INC, loc = assign_points_2_grid(train_hpgl[i], grid_params)
    file_name = f'D:/Wang/HPGL/data/SISIMINC_{i}.INC'
    # Save the data to a text file with a header
    np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')
    
    run_hpgl_lm_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path = pyscript_path_val[i],   
                      inc_path = file_name, single_option=False, use_correlogram = True)
    
    sisim_3D_lm_hpgl=[]
    for j in range(5):
        for z in range(nsim_hpgl):
             file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{j}_{z}.bin'
             sisim_3D_lm_hpgl.append(read_boolean_array_from_binary(file_path))
    write_boolean_array_to_binary(sisim_3D_lm_hpgl, str(hpgl_validation_file_path) + str(i) + 'binary.binary')
  
    
kcv_temp = k_fold_cross_validation_3D_hpgl(k, test_hpgl, hpgl_validation_file_path, grid_params)
correctlist, wronglist, completelist = kcv_temp

kcv_temp = k_fold_cross_validation_3D_HPGL(k, test_hpgl, hpgl_validation_file_path, grid_params)
correctlist_sand, wronglist_sand, correctlist_ton, wronglist_ton, completelist = kcv_temp



correctlist_save_path = str(WD) + "save/correctlist.npy"
wronglist_save_path = str(WD) + "save/wronglist.npy"
completelist_save_path = str(WD) + "save/completelist.npy"
np.save(correctlist_save_path, correctlist)
np.save(wronglist_save_path, wronglist)
np.save(completelist_save_path, completelist)

correctlist = np.load(correctlist_save_path ,allow_pickle=True)
wronglist = np.load(wronglist_save_path,allow_pickle=True)

plot_accuracy_3D(correctlist, wronglist)


############################# TEsting 13.05.2024 #############################
# from 100 simulations delete the last one and perform the calculation with 99.....
def k_fold_cross_validation_3D_hpgl(k, test_hpgl, train_hpgl, hpgl_validation_file_path, grid_params):
    nsim1 = 100
    correctlist = []
    wronglist = []
    completelist = []
    for i in range(k):
        arr_sis = read_boolean_array_from_binary(str(hpgl_validation_file_path) + str(i) + "binary.binary")
        arr_sis = np.array_split(arr_sis, nsim1)[:-1]
        sis_mean = np.mean(np.stack(arr_sis), axis=0)
        sis_etype = np.digitize(sis_mean, (0.5, 1), right=True)

        # Check grid for both test and train datasets
        loc1, check1 = checkgrid(test_hpgl[i][:,0], test_hpgl[i][:,1], test_hpgl[i][:,2], grid_params, XY_option=False)
        loc2, check2 = checkgrid(train_hpgl[i][:,0], train_hpgl[i][:,1], train_hpgl[i][:,2], grid_params, XY_option=False)

        test_locations = set(loc1[check1])  # Only locations within the grid
        train_locations = set(loc2[check2])

        stlist = test_hpgl[i][:,3][check1].astype(int)
        loc = loc1[check1]

        for idx in range(len(loc)):
            if loc[idx] not in train_locations:  # Check if the location is not in the train locations
                completelist.append(float(sis_mean[loc[idx]]))
                if sis_etype[loc[idx]] == stlist[idx]:
                    correctlist.append(float(abs(sis_mean[loc[idx]] - 0.5) + 0.5))
                else:
                    wronglist.append(float(abs(sis_mean[loc[idx]] - 0.5) + 0.5))

    return correctlist, wronglist, completelist

kcv_temp = k_fold_cross_validation_3D_hpgl(k, test_hpgl, train_hpgl, hpgl_validation_file_path, grid_params)
correctlist, wronglist, completelist = kcv_temp

correctlist_tem_save_path = str(WD) + "save/99correctlist.npy"
wronglist_tem__save_path = str(WD) + "save/99wronglist.npy"
completelist_save_path = str(WD) + "save/99completelist.npy"
np.save(correctlist_tem_save_path, correctlist)
np.save(wronglist_tem__save_path, wronglist)



correctlist = np.load(correctlist_save_path ,allow_pickle=True)
wronglist = np.load(wronglist_save_path,allow_pickle=True)

plot_accuracy_3D(correctlist, wronglist)



# ----------------------------------------------------------------------------
# Run HPGL with no Prior
# ----------------------------------------------------------------------------
origin_pyscript = 'D:/Wang/HPGL/LM_test.py'
generate_scripts(origin_pyscript, 5, 'LM_test_')



import time

hpgl_validation_file_path = "D:/Wang/HPGL/data/validation/result_validation_"

pyscript_path="D:/Wang/HPGL/SISIM_LM_"

# pyscript_path = 'D:/Wang/HPGL/LM_test__'

nsim_hpgl = 20
k = 3
train_hpgl, test_hpgl = split_array_3D(k, Cl_Sa_points)

start_time = time.time()

for i in range(k):
    INC = assign_points_2_grid(train_hpgl[i], grid_params)
    file_name = f'D:/Wang/HPGL/data/SISIMINC_{i}.INC'
    # Save the data to a text file with a header
    np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')
    save_train_prior_3D_HPGL(train_hpgl[i],
                        grid_params, radius_xy= 500, radius_z=10)
    # pyscript_path = r'D:\Wang\HPGL\SISIM_' + str(i) + '.py'
    # run_hpgl(grid_params, variogram, se, train[i], pyscript_path)
    run_hpgl_lm_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path = pyscript_path,   
                      inc_path = file_name, single_option=False, use_correlogram = True)
    
    sisim_3D_lm_hpgl=[]
    for j in range(5):
        for z in range(nsim_hpgl):
             file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{j}_{z}.bin'
             sisim_3D_lm_hpgl.append(read_boolean_array_from_binary(file_path))
    write_boolean_array_to_binary(sisim_3D_lm_hpgl, str(hpgl_validation_file_path) + str(i) + 'binary.binary')
         
end_time = time.time()


# ----------------------------------------------------------------------------
# correctlist/wronglist for sand and clay
# ----------------------------------------------------------------------------
readed_filepath = f'D:\\Wang\\Result\\FINAL FINAL\\12 12.5 14 0 Nugget\\Val\\result_validation_'
readed_filepath = f'D:\\Wang\\Result\\FINAL FINAL\\ncnode 24\\Val\\result_validation_'



kcv_temp = k_fold_cross_validation_3D_HPGL(k, test_hpgl, readed_filepath, grid_params)
correctlist_sand, wronglist_sand, correctlist_ton, wronglist_ton, completelist = kcv_temp


correctlist_save_path = str(WD) + "save/correctlist_sand.npy"
wronglist_save_path = str(WD) + "save/wronglist_sand.npy"
completelist_save_path = str(WD) + "save/completelist.npy"
np.save(correctlist_save_path, correctlist_sand)
np.save(wronglist_save_path, wronglist_sand)
np.save(completelist_save_path, completelist)

correctlist_sand = np.load(correctlist_save_path ,allow_pickle=True)
wronglist_sand = np.load(wronglist_save_path,allow_pickle=True)



plot_accuracy_3D_sperate(correctlist_sand, wronglist_sand, "Sand Calibration Curve")

correctlist_save_path = str(WD) + "save/correctlist_clay.npy"
wronglist_save_path = str(WD) + "save/wronglist_clay.npy"
np.save(correctlist_save_path, correctlist_ton)
np.save(wronglist_save_path, wronglist_ton)

# sperate probability
plot_accuracy_3D_sperate(correctlist_ton, wronglist_ton, "Clay Calibration Curve")

# together
plot_accuracy_3D(correctlist, wronglist)



plot_combined_accuracy(correctlist_sand, wronglist_sand, correctlist_ton, wronglist_ton, correctlist, wronglist, combined=False)



bins = np.linspace(0.5, 1.0, 26)
c = np.histogram(correctlist,bins)[0]
w = np.histogram(wronglist,bins)[0]
# Mean Squared Error, MSE
# print ("Mean Squared Error:", np.mean(((bins[1:] - 0.01) - (c/(w+c))).reshape((-1,1))**2))
mse = np.mean((bins[1:] - (c / (w + c)))**2)
print("Mean Squared Error:", mse)

# Root Mean Squared Error, RMSE
rmse = np.sqrt(mse)
print("Root Mean Squared Error:", rmse)
# print ("Root Mean Squared Error:", np.mean(((bins[1:] - 0.02) - (c/(w+c))).reshape((-1,1))**2)**0.5)

acc = len(correctlist)/ (len(correctlist) + len(wronglist))
print("acc:", acc)


print(f"RUN_GSLIB Runtime：{start_time - end_time} s")

# ----------------------------------------------------------------------------
# prescision
# ----------------------------------------------------------------------------

result_prescision = calculate_classification_metrics_k_fold(k, test_hpgl, hpgl_validation_file_path ,grid_params)

k = 3

precision_sand = result_prescision['precision_sand']
precision_clay = result_prescision['precision_clay']
accuracy = result_prescision['accuracy']
# ----------------------------------------------------------------------------
# 3.4.2 K FOLD CROSS VALIDATION (for SISIM)
# ----------------------------------------------------------------------------
from UTILITIES_013_Wa import save_train_3D
k = 3

train, test = split_array_3D(k, Cl_Sa_points)
WD_data = str(WD) + "data/tertiary_"
save_train_3D(train, k, WD_data)


#temp_file_location  = "D:/Wang/Result/SISIM_sstrat1_ndmax12_ncnode12_newnest/correctlist_wronglist/sisim_"

# Use train[] directly


WD_outfl = str(WD) + "save/sisim_"
WD_parfl = str(WD) + "parfiles/sisim_"
for i in range(k):
    datafl = str(WD_data) + str(i) + ".gslib"

    WD_par = str(WD_parfl) + str(i) + "_"
    WD_out = str(WD_outfl) + str(i) + "_"
    WD_dbg = str(WD_outfl) + str(i) + "_"

    run_sisim_parallel(sisimexe,WD_par,WD_out,WD_dbg,datafl,
                              GSLIB_SISIM_params)

for i in range(k):
    WD_out = str(WD_outfl) + str(i) + "_"
    sisim_3D_ohnelm = read_sisim_lm(WD_out)
    # np.save(str(WD) + 'save/sisim_3D_ohnelm' + str(i) + '.npy', sisim_3D_ohnelm)   
    write_boolean_array_to_binary(sisim_3D_ohnelm, str(WD_outfl) + str(i) + 'binary.binary') ################################################### PASTE HERE !!!
    
  
# k=2
WD_outfl = str(WD) + 'save/sisim_'
WD_outfl = "D:/Wang/Result/SISIM_sstrat1_ndmax12_ncnode12_newnest/correctlist_wronglist/sisim_"

kcv_temp = k_fold_cross_validation_3D(k, test, WD_outfl, GSLIB_SISIM_params)
correctlist, wronglist, completelist = kcv_temp

correctlist_save_path = str(WD) + "save/correctlist.npy"
wronglist_save_path = str(WD) + "save/wronglist.npy"
completelist_save_path = str(WD) + "save/completelist.npy"
np.save(correctlist_save_path, correctlist)
np.save(wronglist_save_path, wronglist)
np.save(completelist_save_path, completelist)

correctlist = np.load(str(WD) + "save/correctlist.npy" ,allow_pickle=True)
wronglist = np.load(str(WD) + "save/wronglist.npy" ,allow_pickle=True)

plot_accuracy_3D(correctlist, wronglist)




# ############################################################################
# 3.3.2.1 HPGL validation SISIM
# ############################################################################

variogram = { # E-w  N-S V
    "type": 1, # 1 for exponentail
    "ranges": [14, 12.5, 14],
    "angles": [0, 0, 0],
    "sill": 1,
    "nugget": 0
    }


se = {
    "radiuses": (20, 20, 20),
    "max_neighbours": 12 # ncnode
}



k = 3
from UTILITIES_013_Wa import save_train_3D

train_hpgl, test_hpgl = split_array_3D(k, Cl_Sa_points)


hpgl_validation_file_path = "D:/Wang/HPGL/data/validation/result_validation_"
pyscript_path="D:/Wang/HPGL/SISIM_"

nsim_hpgl = 20

for i in range(k):
    INC = assign_points_2_grid(train_hpgl[i], grid_params)
    file_name = f'D:/Wang/HPGL/data/SISIMINC_{i}.INC'
    # Save the data to a text file with a header
    np.savetxt(file_name, INC, fmt='%d', header='ind_data', comments='')
    # pyscript_path = r'D:\Wang\HPGL\SISIM_' + str(i) + '.py'
    # run_hpgl(grid_params, variogram, se, train[i], pyscript_path)
    run_hpgl_parallel(grid_params, variogram, se, nsim_hpgl, pyscript_path,
                      inc_path = file_name, single_option=False)
    
    sisim_3D_ohnelm_hpgl=[]
    for j in range(5):
        for z in range(nsim_hpgl):
             file_path = f'D:\\Wang\\HPGL\\data\\result_from_hpgl_{j}_{z}.bin'
             sisim_3D_ohnelm_hpgl.append(read_boolean_array_from_binary(file_path))
    write_boolean_array_to_binary(sisim_3D_ohnelm_hpgl, str(hpgl_validation_file_path) + str(i) + 'binary.binary')
         

kcv_temp = k_fold_cross_validation_3D_hpgl(k, test_hpgl, hpgl_validation_file_path, grid_params)
correctlist, wronglist, completelist = kcv_temp


correctlist_save_path = str(WD) + "save/correctlist_.npy"
wronglist_save_path = str(WD) + "save/wronglist_.npy"
completelist_save_path = str(WD) + "save/completelist_.npy"
np.save(correctlist_save_path, correctlist)
np.save(wronglist_save_path, wronglist)
np.save(completelist_save_path, completelist)

correctlist = np.load(str(WD) + "save/correctlist_.npy" ,allow_pickle=True)
wronglist = np.load(str(WD) + "save/wronglist_.npy" ,allow_pickle=True)

plot_accuracy_3D(correctlist, wronglist)



bins = np.linspace(0.5, 1.0, 26)
c = np.histogram(correctlist,bins)[0]
w = np.histogram(wronglist,bins)[0]

mse = np.mean((bins[1:] - (c / (w + c)))**2)
print("Mean Squared Error:", mse)

# Root Mean Squared Error, RMSE
rmse = np.sqrt(mse)
print("Root Mean Squared Error:", rmse)

acc = len(correctlist)/ (len(correctlist) + len(wronglist))
print("acc:", acc)

rmse12 = rmse

rmse40 = rmse

print(rmse12, rmse40)





# ############################################################################
# 3.4.3 Cross validation for prior information
# ############################################################################

train = WD_prior

kcv_temp_prior = k_fold_cross_validation_3D_moving_average(k, test, train, GSLIB_SISIM_LM_params)

correctlist_p, wronglist_p, completelist_p = kcv_temp

plot_accuracy_3D(correctlist, wronglist)


prior_plot = np.load("D:/Wang/HPGL/Prior Val/prior_tohpgl_p0_a.npy")
grid_2= make_uniform_grid(grid_params)
grid_2.cell_data['SIM1'] = prior_plot
grid_2.plot()
grid_2.save(str(WD) + 'save/prior_plot.vtk')
# ----------------------------------------------------------------------------
# ############################################################################
# 4.1 Combine and calculate Probability SISIM_LM
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
from UTILITIES_013_Wa import (get_grid_elev)
 
directory = 'D:/Wang/U9/WD/data/geo_tif_muc/'
elev_grid = get_grid_elev(directory, grid_params, origin, angle)
np.save(str(WD) + "save/elev_grid.npy", elev_grid)
data_2D = sgsim_dtbt
 
np.save(str(WD) + "save/data_2D.npy", data_2D)
elev_grid = np.load(str(WD) + "save/elev_grid.npy")
elev_3D = to_3D(elev_grid[:, -1].flatten(), grid_params, 1)[0] * 1
mean_qb = to_3D(grid_data, grid_params, 1)[0] * 1
data_2D = np.load(str(WD) + "save/data_2D.npy")

# sisim_3D = np.load(str(WD) + 'save/sisim_3D.npy')
# sisim_3D_temp = read_boolean_array_from_binary(str(WD) + "save/sisim_lm_binary.binary")
# sisim_3D = np.array([np.array(item) for item in np.split(sisim_3D_temp, 100)])
# sisim_3D = sisim_3D.astype(object)
# del sisim_3D_temp

nsim = GSLIB_SGSIM_params["nsim"]
grid_params = GSLIB_SISIM_LM_params["grid_params"]
sgsim_3D = to_3D(data_2D, grid_params, nsim)
# nmodel = combine_sim(sisim_3D_lm_hpgl, sgsim_3D, nsim)
nmodel = combine_sim(sisim_3D, sgsim_3D, nsim)
stat_temp = evaluate_nmodel(nmodel, nsim)
np.save(str(WD) + "save/stat_temp.npy", stat_temp)
s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var = stat_temp
grid = make_model(grid_params, s_type, fs_type, p_clay,
                  p_sand, p_gravel, entropy, var, elev_3D)
grid.set_active_scalars("s_type")
grid = grid.rotate_z(-angle, point = origin+[0], inplace =False)
grid.save(str(WD) + 'save/model_rotated.vtk')
# grid.save(str(WD) + 'save/model.vtk')

# grid = make_uniform_grid(grid_params)
# grid.cell_data["s_type"] = mean_qb
# grid.cell_data["s_type"] = np.mean(np.array(sgsim_3D), axis=0)
# grid.cell_data["s_type"] = nmodel[0]
# grid.set_active_scalars("s_type")
# grid.plot()



# ############################################################################
# 4.2 Combine and calculate Probability SISIM
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
from UTILITIES_013_Wa import (get_grid_elev)
from UTILITIES_013_Wa import (combine_sim_ohnelm)

directory = 'D:/Wang/U9/WD/data/geo_tif_muc/'
elev_grid = get_grid_elev(directory, grid_params, origin, angle)
np.save(str(WD) + "save/elev_grid.npy", elev_grid)
data_2D = sgsim_dtbt

# load result of simulation in SGSIM
# data_2D = np.load(str(WD) + "save/data_2D.npy", allow_pickle=True)
# data_2D = np.load(str(WD) + "save/data_2D.npy")

np.save(str(WD) + "save/data_2D.npy", data_2D)
elev_grid = np.load(str(WD) + "save/elev_grid.npy")
elev_3D = to_3D(elev_grid[:, -1].flatten(), grid_params, 1)[0] * 1
mean_qb = to_3D(grid_data, grid_params, 1)[0] * 1


# sisim_3D_ohnelm = np.load(str(WD) + 'save/sisim_3D_ohnelm.npy', allow_pickle=True)
# sisim_3D_ohnelm_temp = read_boolean_array_from_binary(str(WD) + "save/sisim_binary.binary")
# sisim_3D = np.array([np.array(item) for item in np.split(sisim_3D_ohnelm_temp, 100)])
# sisim_3D = sisim_3D.astype(object)
# del sisim_3D_temp

nsim = 100
sgsim_3D = to_3D(data_2D, grid_params, nsim)
nmodel = combine_sim_ohnelm(sisim_3D_ohnelm, sgsim_3D, nsim)
stat_temp = evaluate_nmodel(nmodel, nsim)
np.save(str(WD) + "save/stat_temp.npy", stat_temp)
s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var = stat_temp
grid = make_model(grid_params, s_type, fs_type, p_clay,
                  p_sand, p_gravel, entropy, var, elev_3D)
grid.set_active_scalars("s_type")
# grid = pv.read(str(WD) + 'save/model.vtk')

grid = grid.rotate_z(-angle, point = origin+[0], inplace =False)
grid.save(str(WD) + 'save/model_rotated.vtk')
# grid.plot()

grid = make_uniform_grid(grid_params)
grid.cell_data["s_type"] = mean_qb
grid.cell_data["s_type"] = np.mean(np.array(sgsim_3D), axis=0)
grid.cell_data["s_type"] = nmodel[0]
grid.set_active_scalars("var")
grid.plot()

# ############################################################################
# 4.2 Combine and calculate Probability HPGL
# ############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
from UTILITIES_013_Wa import (get_grid_elev)
from UTILITIES_013_Wa import (combine_sim_ohnelm)

directory = 'D:/Wang/U9/WD/data/geo_tif_muc/'
elev_grid = get_grid_elev(directory, grid_params, origin, angle)
np.save(str(WD) + "save/elev_grid.npy", elev_grid)
data_2D = sgsim_dtbt

# load result of simulation in SGSIM
data_2D = np.load(str(WD) + "save/data_2D.npy", allow_pickle=True)
data_2D = np.load(str(WD) + "save/data_2D.npy")

np.save(str(WD) + "save/data_2D.npy", data_2D)
elev_grid = np.load(str(WD) + "save/elev_grid.npy")
elev_3D = to_3D(elev_grid[:, -1].flatten(), grid_params, 1)[0] * 1
mean_qb = to_3D(grid_data, grid_params, 1)[0] * 1

# sisim_3D_ohnelm = np.load(str(WD) + 'save/sisim_3D_ohnelm.npy', allow_pickle=True)
# sisim_3D_ohnelm_temp = read_boolean_array_from_binary(str(WD) + "save/sisim_binary.binary")
# sisim_3D = np.array([np.array(item) for item in np.split(sisim_3D_ohnelm_temp, 100)])
# sisim_3D = sisim_3D.astype(object)
# del sisim_3D_temp

nsim = 100
sgsim_3D = to_3D(data_2D, grid_params, nsim)
nmodel = combine_sim_ohnelm(sisim_3D_ohnelm, sgsim_3D, nsim)
stat_temp = evaluate_nmodel(nmodel, nsim)
np.save(str(WD) + "save/stat_temp.npy", stat_temp)
s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var = stat_temp
grid = make_model(grid_params, s_type, fs_type, p_clay,
                  p_sand, p_gravel, entropy, var, elev_3D)
grid.set_active_scalars("s_type")
# grid = pv.read(str(WD) + 'save/model.vtk')

grid = grid.rotate_z(-angle, point = origin+[0], inplace =False)
grid.save(str(WD) + 'save/model_rotated.vtk')
# grid.plot()

grid = make_uniform_grid(grid_params)
grid.cell_data["s_type"] = mean_qb
grid.cell_data["s_type"] = np.mean(np.array(sgsim_3D), axis=0)
grid.cell_data["s_type"] = nmodel[0]
grid.set_active_scalars("var")
grid.plot()

# -------------------------------- Read binary --------------------------------
# sisim_3D_ohnelm = np.load(str(WD) + 'save/sisim_3D_ohnelm.npy', allow_pickle=True)
# sisim_3D_ohnelm_temp = read_boolean_array_from_binary(str(WD) + "save/sisim_binary.binary")
# sisim_3D = np.array([np.array(item) for item in np.split(sisim_3D_ohnelm_temp, 100)])
# sisim_3D = sisim_3D.astype(object)
# del sisim_3D_temp
# ---------------------------------------------------------------------------
nsim = GSLIB_SGSIM_params["nsim"]
sgsim_3D = to_3D(data_2D, grid_params, nsim)
nmodel = combine_sim_ohnelm(sisim_3D_ohnelm, sgsim_3D, nsim)
stat_temp = evaluate_nmodel(nmodel, nsim)
np.save(str(WD) + "save/stat_temp.npy", stat_temp)
s_type, fs_type, p_clay, p_sand, p_gravel, entropy, var = stat_temp
grid = make_model(grid_params, s_type, fs_type, p_clay,
                  p_sand, p_gravel, entropy, var, elev_3D)
grid.set_active_scalars("s_type")


grid = grid.rotate_z(-angle, point = origin+[0], inplace =False)
grid.save(str(WD) + 'save/model_rotated.vtk')


grid = make_uniform_grid(grid_params)
grid.cell_data["s_type"] = mean_qb
grid.cell_data["s_type"] = np.mean(np.array(sgsim_3D), axis=0)
grid.cell_data["s_type"] = nmodel[0]
grid.set_active_scalars("var")
grid.plot()



# SEE OTHER PY-DATA!!!!
# # ############################################################################
# # 4.3.1 Make 3D Tunnel
# # ############################################################################
# import pyvista as pv
# import geopandas as gpd
# from shapely.geometry import Point
# from make_tunnel_001 import (lines_from_points)

# tunnel_geometry = pd.read_excel("D:/Wang/U9/WD/data/tunnel_geometry.xlsx")

# tunnel_geometry.columns
# tunnel_geometry["MA_Tunnel"] = tunnel_geometry.OK_Tunnel - ((tunnel_geometry.OK_Tunnel - tunnel_geometry.UK_Tunnel)/2)
# tunnel_geometry_sel = tunnel_geometry[tunnel_geometry["Gleis"] == 'U9 Gleis 2']
# points = np.array(tunnel_geometry_sel[['Rechtswert', 'Hochwert',"MA_Tunnel"]])

# gdf_gk4 = gpd.GeoDataFrame(geometry=[Point(xyz) for xyz in points], crs="31468")
# gdf_utm32 = gdf_gk4.to_crs("25832")

# points = np.c_[np.array(gdf_utm32.geometry.apply(lambda geom: (geom.x, geom.y)).to_list()),points[:,2]]

# spline = pv.Spline(points, 1000)
# tube = spline.tube(radius=5)

# line = lines_from_points(points)
# line.save("D:/Wang/U9/WD/data/tunnel_geometry.vtk")
# tube.save("D:/Wang/U9/WD/data/tunnel_geometry.vtk")
# # ############################################################################
# # 4.3.2 3Dmap
# # ############################################################################
# import pyvista as pv
# import geopandas as gpd
# from shapely.geometry import Point
# import rasterio
# from matplotlib.colors import Normalize
# from matplotlib.cm import ScalarMappable
# from scipy.ndimage import map_coordinates



# # Replace with the paths to your elevation and color TIFF files
# elevation_tif_path = 'D:/Wang/U9/WD/data/geo_tif_muc/639_5358.tif'
# color_tif_path = 'D:/Wang/U9/WD/data/geo_tif_muc/32639_5358.tif'

# # Open the elevation TIFF file using rasterio
# with rasterio.open(elevation_tif_path) as elevation_src:
#     # Read the elevation raster data
#     elevation_data = elevation_src.read(1) # Assuming it's a single band
#     elevation_data = elevation_data[::-1,:]
#     origin = (elevation_src.bounds.left, elevation_src.bounds.bottom, 0)
#     spacing = elevation_src.res + (0,)
#     dimension = (elevation_src.width, elevation_src.height, 1)
#     # Get the spatial transform information
#     transform = elevation_src.transform

# # Create a PyVista mesh from the elevation data
# mesh = pv.UniformGrid()
# mesh.dimensions = dimension
# mesh.origin = origin
# mesh.spacing = spacing
# mesh['Elevation'] = elevation_data.flatten()
# mesh = mesh.warp_by_scalar('Elevation')


# # Open the color TIFF file using rasterio
# with rasterio.open(color_tif_path) as color_src:
#     # Read the color raster data
#     #print(color_src.read().shape)
#     color_data = color_src.read().transpose([1, 2, 0])  # Transpose to (height, width, bands)
#     color_data = color_data[::-1,:,:]

# reduced = np.ceil(np.arange(0,1000,1)*2.5).astype(int)
# interpolated_colors =color_data[np.repeat(reduced,1000),np.tile(reduced,1000),:]

# # Convert RGB to grayscale
# grayscale_values = 0.299 * interpolated_colors[:, 0] + 0.587 * interpolated_colors[:, 1] + 0.114 * interpolated_colors[:, 2]

# # Set grayscale values as scalar to the mesh
# mesh['Colors'] = grayscale_values

# mesh.save('D:/Wy/U9/WD/data/639_5358.vtk')

