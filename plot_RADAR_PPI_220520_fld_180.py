#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_RADAR_PPI.py                                                           #
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
from PLOT_RADAR import plot_PPI, plot_PPI_temp_ring

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20220520"
# location = 'hnr'
# location = 'ess'
location = 'fld'
time_i = 180
# time_i = 179
pdf_or_png = 'png'

include_sweep = np.array([
    False,
    False, False, False, True, False, False,
    False, False, False, False
])
# include_sweep = np.array([
#     False,
#     False, True, False, False, False, False,
#     False, False, False, False
# ])
elevation_degs = np.array([
    5.5,
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8., 12., 17., 25.
])
range_maxes = np.array([
    None,
    60, 60, 70, 100, 150, None,
    60, 60, 60, 60
])
# range_maxes = np.array([
#     None,
#     None, None, None, None, None, None,
#     None, None, None, None
# ])
temp_thicknesses = np.array([
    .3,
    .3, .3, .3, .3, .3, .7,
    .3, .5, .4, .4
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
n_rows = elevation_degs.size
n_cols = 5
n_cols = 4
n_i_zh = 0
n_i_zdr = 1*1
n_i_rho = 1*2
n_i_kdp = 1*3
# n_i_vrad = 1*4
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
    nc_file_comb = glob.glob(folder_in + '/*polmom*')
    if len(nc_file_comb) > 1:
        print('mom: too many files')
    elif len(nc_file_comb) == 0:
        print('mom: no files')
    else:
        nc_file_comb = nc_file_comb[0]

    # ----------------------------------------------------------------------- #
    # time
    time_start = date + nc_file_comb.split('-')[-4][-4:]
    time_end = date + nc_file_comb.split('-')[-3][-4:]
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
    time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
    # time_i_era = int(round(time_i * 24 / 288, 0))
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # plot moments
    moment = 'ZH_AC'
    n_i_zh = n_i_zh + 1
    ax = plt.subplot(n_rows, n_cols, n_i_zh)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)

    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beamtop', range_max=range_max, title=title)
    # ----------------------------------------------------------------------- #
    moment = 'ZDR_AC_OC'
    n_i_zdr = n_i_zdr + 1
    ax = plt.subplot(n_rows, n_cols, n_i_zdr)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)

    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
                       temp_thickness=temp_thickness, moment='temp_beamtop',
                       range_max=range_max, title=title)
    # ----------------------------------------------------------------------- #
    moment = ['RHOHV_NC2P']
    n_i_rho = n_i_rho + 1
    ax = plt.subplot(n_rows, n_cols, n_i_rho)
    if mode == 'vol':
        title = 'RHO at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'RHO at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)
    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
                       temp_thickness=temp_thickness, moment='temp_beamtop',
                       range_max=range_max, title=title)
    # ------------------------------------------------------------------------#
    moment = ['KDP_NC']
    n_i_kdp = n_i_kdp + 1
    ax = plt.subplot(n_rows, n_cols, n_i_kdp)
    if mode == 'vol':
        title = 'KDP at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'KDP at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
                       temp_thickness=temp_thickness,
                       moment='temp_beambottom', range_max=range_max)

    plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
                       temp_thickness=temp_thickness, moment='temp_beamtop',
                       range_max=range_max, title=title)
    # ------------------------------------------------------------------------#
    # n_i_vrad = n_i_vrad + 1
    # moment = ['VRADH']
    # ax = plt.subplot(n_rows, n_cols, n_i_vrad)
    # if mode == 'vol':
    #     title = 'VRADH at ' + str(elevation_deg) + '° ' + \
    #             location.upper() + ' ' + time_UTC
    # else:
    #     title = 'VRADH at pcp ' + \
    #             location.upper() + ' ' + time_UTC
    #
    # plot_PPI(nc_file_comb, ax, time_i, moment, title=title,
    #          range_max=range_max)
    # plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
    #                    temp_thickness=temp_thickness,
    #                    moment='temp_beambottom', range_max=range_max)
    # plot_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
    #                    temp_thickness=temp_thickness, moment='temp_beamtop',
    #                    range_max=range_max, title=title)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_in = nc_file_comb.split('/')[-1]
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

if sum(include_sweep) == 1:
    file_out = file_in.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
                               '_' + str(n_rows) + 'x' + str(n_cols) +
                               '.' + pdf_or_png
                               )
    plt.savefig(folder_plot + 'PPI_' + str(elevation_deg) + '_' +
                file_out, format=pdf_or_png, transparent=True)
else:
    file_out = file_in.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
                               '_' + str(n_rows) + 'x' + str(n_cols) +
                               '.' + pdf_or_png
                               ).replace('_' + sweep, '_all')
    plt.savefig(folder_plot + 'VOL_' + str(n_rows) + 'x' + str(n_cols) + '_' +
                file_out, format=pdf_or_png, transparent=True)
plt.close()
