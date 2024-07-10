#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_syn_RADAR_PPI.py                                                       #
#                                                                             #
# plot PPI scans of moments                                                   #
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
from PLOT_SYN_RADAR import plot_syn_PPI, plot_syn_PPI_temp_ring

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20220520"
location = 'fld'
# time_i = 180
time_utc = 1455
time_utc = 1500
# time_utc = 1505
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
da_run = 'ASS_2405'
icon_emvorado_run = 'MAIN_2405.1/EMVO_00400000.2'
spin_up_mm = '60'
# --------------------------------------------------------------------------- #
# folder and file search
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_utc - time_utc % 600).zfill(4)
hhmm_end = str(time_utc - time_utc % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
# --------------------------------------------------------------------------- #

# include_sweep = np.array([
#     False,
#     False, False, False, True, False, False,
#     False, False, False, False
# ])
include_sweep = np.array([
    False,
    False, True, False, False, False, False,
    False, False, False, False
])

# include_sweep = np.array([
#     False,
#     False, False, False, False, False, True,
#     False, False, False, False
# ])
#
# include_sweep = np.array([
#     False,
#     False, False, False, False, True, False,
#     False, False, False, False
# ])
elevation_degs = np.array([
    5.5,
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8., 12., 17., 25.
])
range_maxes = np.array([
    None,
    80, 80, 80, 100, 150, None,
    60, 60, 60, 60
])
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
n_cols = 4
n_i_zh = 0
n_i_zdr = 1 * 1
n_i_rho = 1 * 2
n_i_kdp = 1 * 3
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
# ------------------------------------------------------------------------#

# for elevation_deg, range_max, temp_thickness, mode in \
#         zip(elevation_degs, range_maxes, temp_thicknesses, modes):
elevation_deg = elevation_degs[0]
range_max = range_maxes[0]
temp_thickness = temp_thicknesses[0]
mode = modes[0]

# ------------------------------------------------------------------------#
print(elevation_deg)
# ------------------------------------------------------------------------#
# sweep
nc_file_comb = syn_nc.sel(elevation=elevation_deg)
# time
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i = np.where(dti == pd.to_datetime(date + str(time_utc).zfill(4)))[0][0]
time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# plot moments
moment = 'zrsim'
n_i_zh = n_i_zh + 1
ax = plt.subplot(n_rows, n_cols, n_i_zh)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC

plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
#                    temp_thickness=temp_thickness,
#                    moment='temp_beambottom', range_max=range_max)
#
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
#                    temp_thickness=temp_thickness,
#                    moment='temp_beamtop', range_max=range_max, title=title)
# ----------------------------------------------------------------------- #
moment = 'zdrsim'
n_i_zdr = n_i_zdr + 1
ax = plt.subplot(n_rows, n_cols, n_i_zdr)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC

plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
#                    temp_thickness=temp_thickness,
#                    moment='temp_beambottom', range_max=range_max)
#
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
#                    temp_thickness=temp_thickness, moment='temp_beamtop',
#                    range_max=range_max, title=title)
# ----------------------------------------------------------------------- #
moment = ['rhvsim']
n_i_rho = n_i_rho + 1
ax = plt.subplot(n_rows, n_cols, n_i_rho)
if mode == 'vol':
    title = 'RHO at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = 'RHO at pcp ' + \
            location.upper() + ' ' + time_UTC

plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
#                    temp_thickness=temp_thickness,
#                    moment='temp_beambottom', range_max=range_max)
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
#                    temp_thickness=temp_thickness, moment='temp_beamtop',
#                    range_max=range_max, title=title)
# ------------------------------------------------------------------------#
moment = ['kdpsim']
n_i_kdp = n_i_kdp + 1
ax = plt.subplot(n_rows, n_cols, n_i_kdp)
if mode == 'vol':
    title = 'KDP at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = 'KDP at pcp ' + \
            location.upper() + ' ' + time_UTC

plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
             range_max=range_max)
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=273.15,
#                    temp_thickness=temp_thickness,
#                    moment='temp_beambottom', range_max=range_max)
#
# plot_syn_PPI_temp_ring(nc_file_comb, ax, time_i, temp=277.15,
#                    temp_thickness=temp_thickness, moment='temp_beamtop',
#                    range_max=range_max, title=title)

###

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm) + 'min.'
str_location = '_'.join([location.upper(), date + str(time_utc).zfill(4)])
if sum(include_sweep) == 1:
    file_out = folder_plot + 'SYN_PPI_' + str(elevation_deg) + \
               '°_' + str_location + '_' + str_mod + pdf_or_png
else:
    file_out = folder_plot + 'SYN_VOL_' + str(n_rows) + 'x' + \
               str(n_cols) + '_' + str_location + '_' + str_mod + pdf_or_png

plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(file_out, format=pdf_or_png, transparent=True)
plt.close()
