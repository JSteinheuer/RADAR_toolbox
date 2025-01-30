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
# Colors NINJO DWD JM                                                         #
# --------------------------------------------------------------------------- #

header.colors_radar = np.array(
    [[1.00, 1.00, 1.00, 1],
     [0.6, 1.0, 1.0, 1.],
     [0.2, 1.0, 1.0, 1.],
     [0.0, 0.7921569, 0.7921569, 1.],
     [0.0, 0.6, 0.20392157, 1.],
     [0.3019608, 0.7490196, 0.101960786, 1.],
     [0.6, 0.8, 0.0, 1.],
     [0.8, 0.9019608, 0.0, 1.],
     [1.0, 1.0, 0.0, 1.],
     [1.0, 0.76862746, 0.0, 1.],
     [1.0, 0.5372549, 0.0, 1.],
     [1.0, 0.0, 0.0, 1.],
     [0.7058824, 0.0, 0.0, 1.],
     [0.28235295, 0.28235295, 1.0, 1.],
     [0.0, 0.0, 0.7921569, 1.],
     [0.6, 0.0, 0.6, 1.],
     [1.0, 0.2, 1.0, 1.],
     [1.0, 0.8, 1.0, 1.],
     ])
header.cmap_radar = mpl.colors.ListedColormap(header.colors_radar)

# Zh
header.levels_zh = np.append(np.arange(1., 56., 4.5),
                             np.array([60., 65., 75., 85.]))
header.norm_zh = mpl.colors.BoundaryNorm(header.levels_zh,
                                         len(header.levels_zh) - 1)

# ZDR
header.levels_zdr = np.append(np.append([-10, -1], np.arange(-0.1, 0.5, 0.1)),
                              [0.6, 0.8, 1., 1.5, 2., 3., 4.5, 6, 8])
header.norm_zdr = mpl.colors.BoundaryNorm(header.levels_zdr,
                                          len(header.levels_zdr) - 1)

# KDP
header.levels_kdp = [-.2, -.1, 0,
                     .05, .1, .15,
                     .2, .3, .4,
                     .6, .8, 1,
                     2, 3, 4, 5, 6]
header.norm_kdp = mpl.colors.BoundaryNorm(header.levels_kdp,
                                          len(header.levels_kdp) - 1)

# RHOHV
header.levels_rhohv = [.7, .8, .9,
                       .92, .94,
                       .95, .96, .97, .98,
                       .985, .99,
                       .9925, .995,
                       .997, .999,
                       .9995, 1]
header.norm_rhohv = mpl.colors.BoundaryNorm(header.levels_rhohv,
                                            len(header.levels_rhohv) - 1)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20210714"
location = 'ess'
# time_i = 180
time_utc = 1700
# time_utc = 1500
# time_utc = 1505
pdf_or_png = 'png'
# pdf_or_png = 'pdf'
folder_plot = header.folder_ppi_plot
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.1'
icon_emvorado_run = 'MAIN_2411.1/EMVO_00510000.2'
spin_up_mm = '120'
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
path2 = '/'.join([header.dir_data_vol[:-1], date, da_run,
                  icon_run, '/ICONdata', str(spin_up_mm) +
                  'min_spinup/ICON_Vol_' + location.upper() + '_' +
                  date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
syn_nc2 = xr.open_dataset(path2)
# --------------------------------------------------------------------------- #

elevation_deg = 2.5
range_max = 100
temp_thickness = .3
mode = 'vol'

elevation_deg = 1.5
range_max = 150
temp_thickness = .3
mode = 'vol'

elevation_deg = 4.5
range_max = 75
temp_thickness = .3
mode = 'vol'

elevation_deg = 0.5
range_max = 180
temp_thickness = .3
mode = 'vol'

elevation_deg = 12
range_max = 50
temp_thickness = .3
mode = 'vol'

temp_thickness = .3
mode = 'vol'
for elevation_deg, range_max in zip([0.5, 1.5, 2.5, 3.5, 4.5, 12, ],
                                    [180, 150, 150, 100, 75, 50, ]):

    # case (adjust!):
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    # plot parameters
    n_rows = 4
    n_cols = 6
    fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))

    # ------------------------------------------------------------------------#
    print(elevation_deg)
    # ------------------------------------------------------------------------#
    # sweep
    nc_file_comb = syn_nc.sel(elevation=elevation_deg)
    nc_file_comb2 = syn_nc2.sel(elevation=elevation_deg)
    # time
    dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                        freq="5min", inclusive='both')
    time_i = np.where(dti == pd.to_datetime(date +
                                            str(time_utc).zfill(4)))[0][0]
    time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # plot moments
    moment = 'zrsim'
    n_i = 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
                 range_max=range_max)
    # ----------------------------------------------------------------------- #
    moment = 'zdrsim'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    if mode == 'vol':
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = moment + ' at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
                 range_max=range_max)
    # ----------------------------------------------------------------------- #
    moment = ['rhvsim']
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    if mode == 'vol':
        title = 'RHO at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'RHO at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
                 range_max=range_max)
    # ------------------------------------------------------------------------#
    moment = ['kdpsim']
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    if mode == 'vol':
        title = 'KDP at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
    else:
        title = 'KDP at pcp ' + \
                location.upper() + ' ' + time_UTC

    plot_syn_PPI(nc_file_comb, ax, time_i, moment, title=title,
                 range_max=range_max)

    # ------------------------------------------------------------------------#
    moment = 'w'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = 'w at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_PPI(nc_file_comb2, ax, time_i, moment, title=title,
                 range_max=range_max)
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key, i in zip(hm, range(len(hm))):
        moment = 'q' + key
        n_i = n_cols * 2 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC

        plot_syn_PPI(nc_file_comb2, ax, time_i, moment, title=title,
                     range_max=range_max)

    # ----------------------------------------------------------------------- #
    from PROCESS_SYN_RADAR import adjust_icon_fields, \
        mgdparams, calc_moments, calc_multimoments, calc_Dmean

    q_dens, qn_dens = adjust_icon_fields(nc_file_comb2)
    multi_params = mgdparams(q_dens, qn_dens)
    moments = calc_moments(mgd=multi_params)
    multimoments = calc_multimoments(moments)
    mean_volume_diameter = calc_Dmean(multimoments)
    for key in mean_volume_diameter.keys():
        print(key)
        mean_volume_diameter[key][mean_volume_diameter[key] <= 0] = np.nan
        nc_file_comb2 = nc_file_comb2.assign({'Dm_' + key: (
            ('time', 'range', 'azimuth'), mean_volume_diameter[key] * 1000,)})
        nc_file_comb2['Dm_' + key].attrs = {'units': 'mm'}

    hm = ['rain', 'cloud', 'graupel', 'hail', 'ice', 'snow']
    for key, i in zip(hm, range(len(hm))):
        moment = 'Dm_' + key
        n_i = n_cols * 1 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC

        plot_syn_PPI(nc_file_comb2, ax, time_i, moment, title=title,
                     range_max=range_max)

    # ------------------------------------------------------------------------#
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key, i in zip(hm, range(len(hm))):
        moment = 'qn' + key
        n_i = n_cols * 3 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC

        plot_syn_PPI(nc_file_comb2, ax, time_i, moment, title=title,
                     range_max=range_max)

    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # SAVE                                                                    #
    # ----------------------------------------------------------------------- #
    str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + \
              str(spin_up_mm) + 'min_q_N.'
    str_location = '_'.join([location.upper(), date + str(time_utc).zfill(4)])
    file_out = folder_plot + 'SYN_PPI_' + str(elevation_deg) + \
               '°_' + str_location + '_' + str_mod + pdf_or_png
    plt.tight_layout()
    if not os.path.exists(folder_plot):
        os.makedirs(folder_plot)

    plt.savefig(file_out, format=pdf_or_png, transparent=True)
    plt.close()
