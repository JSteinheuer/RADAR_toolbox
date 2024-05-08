#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_PPI_with_ML.py                                                         #
#                                                                             #
# plot PPI scans of moments with ERA5 temperatur rings of bottom and top temp #
# per range bin.                                                              #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import xarray as xr
import datatree as dttree
import numpy as np
import matplotlib.pyplot as plt
import wradlib as wrl
import glob
import pandas as pd
import os
import matplotlib as mpl
from PLOT_PPI import plot_PPI, plot_PPI_temp_ring

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20221222"
location = 'oft'
time_i = 19
pdf_or_png = 'png'

include_sweep = np.array([
    True,
    True, False, True, False, True, True,
    False, True, False, True
])
# include_sweep = np.array([
#     True,
#     True, True, True, True, True, True,
#     True, True, True, True
# ])
elevation_degs = np.array([
    5.5,
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8., 12., 17., 25.
])
range_maxes = np.array([
    None,
    25, 30, 40, 55, 150, None,
    20, 10, 10, 10
])
# range_maxes = np.array([
#     None,
#     None, None, None, None, None, None,
#     None, None, None, None
# ])
temp_thicknesses = np.array([
    .2,
    .2, .2, .2, .2, .2, .2,
    .3, .4, .4, .5
])
modes = np.array([
    'pcp',
    'vol', 'vol', 'vol', 'vol', 'vol', 'vol',
    'vol', 'vol', 'vol', 'vol'
])
elevation_degs = elevation_degs[include_sweep]
range_maxes = range_maxes[include_sweep]
temp_thicknesses = temp_thicknesses[include_sweep]
modes = modes[include_sweep]
# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# plot parameters
n_rows = 4
n_cols = elevation_degs.size
n_i_zh = 0
n_i_zdr = n_cols*1
n_i_rho = n_cols*2
n_i_kdp = n_cols*3
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
# --------------------------------------------------------------------------- #
# loop over all elevations:
# --------------------------------------------------------------------------- #
# plot parameters
elevation_degs_2_sort = elevation_degs.copy()
elevation_degs_2_sort[modes == 'pcp'] = 0
sort_i = elevation_degs_2_sort.argsort()
range_maxes = np.array(range_maxes)[sort_i]
temp_thicknesses = np.array(temp_thicknesses)[sort_i]
elevation_degs = elevation_degs[sort_i]
modes = modes[sort_i]
for elevation_deg, range_max, temp_thickness, mode in \
        zip(elevation_degs, range_maxes, temp_thicknesses, modes):
    print(elevation_deg)
    # ------------------------------------------------------------------------#
    # folder and file search
    folder_plot = header.folder_ppi_plot
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    folder_in = "/".join([header.dir_data_obs + '*',
                          year, year + '-' + mon,
                          year + '-' + mon + '-' + day,
                          location, mode + '*', sweep])
    nc_file_mom = glob.glob(folder_in + '/*allmoms*')
    if len(nc_file_mom) > 1:
        print('mom: too many files')
    elif len(nc_file_mom) == 0:
        print('mom: no files')
    else:
        nc_file_mom = nc_file_mom[0]
    nc_file_temp = glob.glob(folder_in + '/*temp*')
    if len(nc_file_temp) > 1:
        print('temp: too many files')
    elif len(nc_file_temp) == 0:
        print('temp: no files')
    else:
        nc_file_temp = nc_file_temp[0]
    nc_file_rho = glob.glob(folder_in + '/*rhohv_nc*')
    if len(nc_file_rho) > 1:
        print('rho: too many files')
    elif len(nc_file_rho) == 0:
        print('rho: no files')
        nc_file_rho = nc_file_mom
    else:
        nc_file_rho = nc_file_rho[0]
    nc_file_kdp = glob.glob(folder_in + '/*kdp_nc*')
    if len(nc_file_kdp) > 1:
        print('kdp: too many files')
        nc_file_kdp = nc_file_mom
    elif len(nc_file_kdp) == 0:
        print('kdp: no files')
        nc_file_kdp = nc_file_mom
    else:
        nc_file_kdp = nc_file_kdp[0]
    # ----------------------------------------------------------------------- #
    # time
    time_start = date + nc_file_mom.split('-')[-4][-4:]
    time_end = date + nc_file_mom.split('-')[-3][-4:]
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
    time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
    time_i_era = int(round(time_i * 24 / 288, 0))
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # plot moments
    moment = 'DBZH'
    n_i_zh = n_i_zh + 1
    ax = plt.subplot(n_rows, n_cols, n_i_zh)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_mom, ax, time_i, moment, title=title, range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=277.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beamtop', range_max=range_max, title=title)
    # ----------------------------------------------------------------------- #
    moment = 'ZDR'
    n_i_zdr = n_i_zdr + 1
    ax = plt.subplot(n_rows, n_cols, n_i_zdr)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_mom, ax, time_i, moment, title=title, range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=277.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beamtop', range_max=range_max, title=title)
    # ----------------------------------------------------------------------- #
    moment = ['RHOHV_NC2P', 'RHOHV']
    n_i_rho = n_i_rho + 1
    ax = plt.subplot(n_rows, n_cols, n_i_rho)
    if mode == 'vol':
        title = 'RHO at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'RHO at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_rho, ax, time_i, moment, title=title, range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=277.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beamtop', range_max=range_max, title=title)
    # ------------------------------------------------------------------------#
    moment = ['PHI_C', 'PHI_NC', 'UPHIDP', ]
    n_i_kdp = n_i_kdp + 1
    ax = plt.subplot(n_rows, n_cols, n_i_kdp)
    if mode == 'vol':
        title = 'PHIDP at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'PHIDP at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_kdp, ax, time_i, moment, title=title, range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)
    plot_PPI_temp_ring(nc_file_temp, ax, time_i_era, temp=277.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beamtop', range_max=range_max, title=title)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_in = nc_file_mom.split('/')[-1]
file_out = file_in.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
                           '_' + str(n_rows) + 'x' + str(n_cols) +
                           '.' + pdf_or_png
                           ).replace('_' + sweep, '_all')
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'VOL_' + str(n_rows) + 'x' + str(n_cols) + '_' +
            file_out, format=pdf_or_png, transparent=True)
plt.close()
