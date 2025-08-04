#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_syn_RADAR_pseudoRHI.py                                                 #
#                                                                             #
# plot pseudoRHI scans of moments                                             #
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
from PLOT_SYN_RADAR import plot_syn_pseudoRHI

# --------------------------------------------------------------------------- #
# Colors DWD JM                                                               #
# --------------------------------------------------------------------------- #

# header.colors_radar = np.array(
#     [[1.00, 1.00, 1.00], [0.70, 1.00, 1.00],  # white, light cyan,
#      [0.00, 1.00, 1.00],  # cyan
#      [0.50, 1.00, 0.00], [0.40, 0.80, 0.00], [0.27, 0.55, 0.00],  # greens
#      [1.00, 1.00, 0.00], [0.80, 0.80, 0.00], [1.00, 0.65, 0.00],  # yellows
#      [1.00, 0.27, 0.00], [0.80, 0.22, 0.00], [0.55, 0.15, 0.00],  # reds
#      [0.00, 0.70, 0.93], [0.00, 0.00, 1.00],  # blues
#      [1.00, 0.00, 1.00], [0.58, 0.44, 0.86]])  # pinks
# header.cmap_radar = mpl.colors.ListedColormap(header.colors_radar)
#
# # Zh
# header.levels_zh = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80]
# header.norm_zh = mpl.colors.BoundaryNorm(header.levels_zh,
#                                          len(header.levels_zh) - 1)
#
# # ZDR
# header.levels_zdr = [-.2, -.1, 0, .1, .2, .3, .4,
#                      .6, .8, 1.2, 1.6, 2.3, 3, 4.5, 6]
# header.norm_zdr = mpl.colors.BoundaryNorm(header.levels_zdr,
#                                           len(header.levels_zdr) - 1)
#
# # KDP
# header.levels_kdp = [-.4, -.2, 0, .05, .1, .2, .3, .45, .6, .8, 1, 2, 3, 4, 7]
# header.norm_kdp = mpl.colors.BoundaryNorm(header.levels_kdp,
#                                           len(header.levels_kdp) - 1)
#
# # RHOHV
# header.levels_rhohv = [.7, .8, .9, .92, .94, .95,
#                 .96, .97, .98, .985, .99, .9925, .995, .9975, 1]
# header.norm_rhohv = mpl.colors.BoundaryNorm(header.levels_rhohv,
#                                             len(header.levels_rhohv) - 1)

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
date = "20220520"
location = 'fld'
# time_i = 180
# time_utc = 1455
time_utc = 1500
azimuth_deg = 229.5  # c2
azimuth_deg = 230.5  # c2
azimuth_deg = 231.5  # c2
azimuth_deg = 301.5  # c1
azimuth_deg = 302.5  # c1
azimuth_deg = 303.5  # c1
azimuth_deg = 304.5  # c1
for azimuth_deg in [229.5, 230.5, 231.5, 301.5, 302.5, 303.5, 304.5]:
    pdf_or_png = 'png'
    folder_plot = header.folder_rhi_plot
    da_run = 'ASS_2405'
    icon_run = 'MAIN_2405.1'
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
    path2 = '/'.join([header.dir_data_vol[:-1], date, da_run,
                      icon_run, '/ICONdata', str(spin_up_mm) +
                      'min_spinup/ICON_Vol_' + location.upper() + '_' +
                      date + hhmm_start + '_' + date + hhmm_end + '.nc'])
    syn_nc = xr.open_dataset(path)
    syn_nc2 = xr.open_dataset(path2)
    # ----------------------------------------------------------------------- #

    elevation_degs = np.array([
        0.5, 1.5, 2.5, 3.5, 4.5, 5.5,
        8., 12., 17., 25.
    ])
    range_max = 200
    # case (adjust!):
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # plot parameters
    n_rows = 4  # elevation_degs.size
    n_cols = 6
    n_i = 0
    fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
    # ----------------------------------------------------------------------- #
    # loop over all elevations:
    # ----------------------------------------------------------------------- #

    # ------------------------------------------------------------------------#
    print(azimuth_deg)
    # ------------------------------------------------------------------------#
    # sweep
    # nc_file_comb = syn_nc.sel(azimuth=azimuth_deg, method='nearest')
    # time
    dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                        freq="5min", inclusive='both')
    time_i = np.where(dti == pd.to_datetime(date + str(time_utc).zfill(4)))[0][
        0]
    time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
    nc_file_comb = syn_nc.isel(time=time_i)
    nc_file_comb2 = syn_nc2.isel(time=time_i)
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # plot moments
    moment = 'zrsim'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = moment + ' at ' + str(azimuth_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
                       range_max=range_max)
    # ----------------------------------------------------------------------- #
    moment = 'zdrsim'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = moment + ' at ' + str(azimuth_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
                       range_max=range_max)
    # ----------------------------------------------------------------------- #
    moment = 'rhvsim'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = moment + ' at ' + str(azimuth_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
                       range_max=range_max)
    # ------------------------------------------------------------------------#
    moment = 'kdpsim'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = moment + ' at ' + str(azimuth_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
                       range_max=range_max)
    # ------------------------------------------------------------------------#
    moment = 'w'
    n_i = n_i + 1
    ax = plt.subplot(n_rows, n_cols, n_i)
    title = 'w at ' + str(azimuth_deg) + '° ' + \
            location.upper() + ' ' + time_UTC

    plot_syn_pseudoRHI(nc_file_comb2, ax, azimuth_deg, moment,
                       title=title, range_max=range_max)
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    from PROCESS_SYN_RADAR import adjust_icon_fields, \
        mgdparams, calc_moments, calc_multimoments, calc_Dmean

    # for key in ['qr', 'qc', 'qg', 'qh', 'qi', 'qs']:
    #     nc_file_comb2[key] = xr.where(nc_file_comb2[key] < 1e-40,
    #                                   10, nc_file_comb2[key])

    q_dens, qn_dens = adjust_icon_fields(nc_file_comb2)
    multi_params = mgdparams(q_dens, qn_dens)
    moments = calc_moments(mgd=multi_params)
    multimoments = calc_multimoments(moments)
    mean_volume_diameter = calc_Dmean(multimoments)
    for key in mean_volume_diameter.keys():
        print(key)
        mean_volume_diameter[key][mean_volume_diameter[key] <= 0] = np.nan
        nc_file_comb2 = nc_file_comb2.assign({'Dm_' + key: (
            ('range', 'azimuth', 'elevation'),
            mean_volume_diameter[key] * 1000,)})
        nc_file_comb2['Dm_' + key].attrs = {'units': 'mm'}

    hm = ['rain', 'cloud', 'graupel', 'hail', 'ice', 'snow']
    for key, i in zip(hm, range(len(hm))):
        moment = 'Dm_' + key
        n_i = n_cols * 1 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(azimuth_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
        plot_syn_pseudoRHI(nc_file_comb2, ax, azimuth_deg, moment,
                           title=title, range_max=range_max)

    # ------------------------------------------------------------------------#
    # ----------------------------------------------------------------------- #
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key, i in zip(hm, range(len(hm))):
        moment = 'q' + key
        n_i = n_cols * 2 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(azimuth_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
        plot_syn_pseudoRHI(nc_file_comb2, ax, azimuth_deg, moment,
                           title=title, range_max=range_max)

    # ----------------------------------------------------------------------- #
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key, i in zip(hm, range(len(hm))):
        moment = 'qn' + key
        n_i = n_cols * 3 + 1 + i
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(azimuth_deg) + '° ' + \
                location.upper() + ' ' + time_UTC
        plot_syn_pseudoRHI(nc_file_comb2, ax, azimuth_deg, moment,
                           title=title, range_max=range_max)

    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # SAVE                                                                        #
    # --------------------------------------------------------------------------- #
    str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(
        spin_up_mm) + 'min_JM.'
    str_location = '_'.join([location.upper(), date + str(time_utc).zfill(4)])
    # if sum(include_sweep) == 1:
    file_out = folder_plot + 'SYN_pseudoRHI_' + str(azimuth_deg) + \
               '°_' + str_location + '_Dm_' + str_mod + pdf_or_png
    # else:
    #     file_out = folder_plot + 'SYN_VOL_' + str(n_rows) + 'x' + \
    #                str(n_cols) + '_' + str_location + '_' + str_mod + pdf_or_png

    plt.tight_layout()
    if not os.path.exists(folder_plot):
        os.makedirs(folder_plot)

    plt.savefig(file_out, format=pdf_or_png, transparent=True)
    plt.close()
