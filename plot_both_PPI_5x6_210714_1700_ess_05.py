#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 05.02.25                                                 #
# plot_both_PPI_5x6.py                                                        #
#                                                                             #
# [...] Description here [...]                                                #
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
from PLOT_RADAR import plot_PPI
from PLOT_SYN_RADAR import plot_syn_PPI, plot_syn_PPI_temp_ring, calc_vradh
from PROCESS_SYN_RADAR import adjust_icon_fields, \
    mgdparams, calc_moments, calc_multimoments, calc_Dmean

# --------------------------------------------------------------------------- #
# initialize case                                                             #
# --------------------------------------------------------------------------- #
date = "20210714"
time_hhmm = 1700
location = 'ess'  # for obs
elevation_degs = [0.5, ]
range_maxes = [None, ]
elevation_degs = [0.5, 1.5, 2.5, 3.5, 4.5, 12, ]
range_maxes = [180, 150, 150, 100, 75, 50, ]
mode = 'vol'
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.1'
icon_emvorado_run = 'MAIN_2411.1/EMVO_00510000.2'
# icon_emvorado_run = 'MAIN_2411.1/EMVO_00410000.2'
spin_up_mm = '120'
# --------------------------------------------------------------------------- #
date = "20220520"
time_hhmm = 1500
location = 'fld'  # for obs
elevation_degs = [0.5, ]
range_maxes = [None, ]
elevation_degs = [0.5, 1.5, 2.5, 3.5, 4.5, 12, ]
range_maxes = [180, 150, 150, 100, 75, 50, ]
mode = 'vol'
da_run = 'ASS_2405'
icon_run = 'MAIN_2405.1'
icon_emvorado_run = 'MAIN_2405.1/EMVO_00400000.2'
spin_up_mm = '120'
spin_up_mm = '60'
# --------------------------------------------------------------------------- #
# initialize mod/obs                                                          #
# --------------------------------------------------------------------------- #
year = date[0:4]
mon = date[4:6]
day = date[6:8]
# --------------------------------------------------------------------------- #
# initialize mod                                                              #
# --------------------------------------------------------------------------- #
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
                      icon_run, '/ICONdata', str(spin_up_mm) +
                      'min_spinup/ICON_Vol_' + location.upper() + '_' +
                      date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
syn_nc2 = xr.open_dataset(path_icon)
syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')

# --------------------------------------------------------------------------- #
# loop over elevations                                                        #
# --------------------------------------------------------------------------- #
for elevation_deg, range_max in zip(elevation_degs, range_maxes):
    # ----------------------------------------------------------------------- #
    # initialize plot                                                         #
    # ----------------------------------------------------------------------- #
    n_rows = 5
    n_cols = 6
    fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
    pdf_or_png = 'png'
    folder_plot = header.folder_ppi_plot
    n_i = 0
    # ----------------------------------------------------------------------- #
    # initialize mod                                                          #
    # ----------------------------------------------------------------------- #
    nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
                             syn_nc2.sel(elevation=elevation_deg)])
    # ----------------------------------------------------------------------- #
    # initialize obs                                                          #
    # ----------------------------------------------------------------------- #
    sweep = '0' + str(np.where(
        header.ELEVATIONS_ALL == float(elevation_deg))[0][0])
    time_i_obs = int(time_hhmm * 0.12)
    folder_in = "/".join([header.dir_data_obs + '*', year,
                          year + '-' + mon,
                          year + '-' + mon + '-' + day,
                          location, mode + '*', sweep])
    nc_file = glob.glob(folder_in + '/*polmoms*')
    if len(nc_file) > 1:
        print('polmom: too many files')
    elif len(nc_file) == 0:
        print('polmom: no files')
    else:
        nc_file = nc_file[0]

    time_start = date + nc_file.split('-')[-4][-4:]
    time_end = date + nc_file.split('-')[-3][-4:]
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
    time_UTC_obs = dti[time_i_obs].strftime('%Y-%m-%d %H:%M')
    # ----------------------------------------------------------------------- #
    # row 1: obs                                                              #
    # ----------------------------------------------------------------------- #
    moments = ['ZH_AC', 'ZDR_AC_OC', 'RHOHV_NC2P', 'PHI_NC', 'KDP_NC', 'VRADH']
    for moment in moments:
        n_i = n_i + 1
        ax = plt.subplot(n_rows, n_cols, n_i)
        if mode == 'vol':
            title = moment + ' at ' + str(elevation_deg) + '° ' + \
                    location.upper() + ' ' + time_UTC_obs
        else:
            title = moment + ' at pcp ' + \
                    location.upper() + ' ' + time_UTC_obs
        plot_PPI(nc_file, ax, time_i_obs, moment, title=title,
                     range_max=range_max)
    # ----------------------------------------------------------------------- #
    # row 2: mod pol moms                                                     #
    # ----------------------------------------------------------------------- #
    moments = ['zrsim', 'zdrsim', 'rhvsim', 'temp', 'kdpsim', 'vradh']
    # moments = ['zrsim', 'u', 'v', 'w', 'temp', 'vradh']
    for moment in moments:
        n_i = n_i + 1
        ax = plt.subplot(n_rows, n_cols, n_i)
        if mode == 'vol':
            title = moment + ' at ' + str(elevation_deg) + '° ' + \
                    location.upper() + ' ' + time_UTC_mod
        else:
            title = moment + ' at pcp ' + \
                    location.upper() + ' ' + time_UTC_mod
        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                     range_max=range_max)
    # ----------------------------------------------------------------------- #
    # row 3: hm Dm                                                            #
    # ----------------------------------------------------------------------- #
    hm = ['rain', 'cloud', 'graupel', 'hail', 'ice', 'snow']
    q_dens, qn_dens = adjust_icon_fields(nc_file_comb)
    multi_params = mgdparams(q_dens, qn_dens)
    moments = calc_moments(mgd=multi_params)
    multimoments = calc_multimoments(moments)
    mean_volume_diameter = calc_Dmean(multimoments)
    for key in mean_volume_diameter.keys():
        print(key)
        mean_volume_diameter[key][mean_volume_diameter[key] <= 0] = np.nan
        nc_file_comb = nc_file_comb.assign({'Dm_' + key: (
            ('time', 'range', 'azimuth'), mean_volume_diameter[key] * 1000,)})
        nc_file_comb['Dm_' + key].attrs = {'units': 'mm'}

    for key in hm:
        moment = 'Dm_' + key
        n_i = n_i + 1
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC_mod

        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                     range_max=range_max)
    # ----------------------------------------------------------------------- #
    # row 4: hm q                                                             #
    # ----------------------------------------------------------------------- #
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key in hm:
        moment = 'q' + key
        n_i = n_i + 1
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC_mod
        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                     range_max=range_max)
    # ----------------------------------------------------------------------- #
    # row 5: hm qn                                                            #
    # ----------------------------------------------------------------------- #
    hm = ['r', 'c', 'g', 'h', 'i', 's']
    for key in hm:
        moment = 'qn' + key
        n_i = n_i + 1
        ax = plt.subplot(n_rows, n_cols, n_i)
        title = moment + ' at ' + str(elevation_deg) + '° ' + \
                location.upper() + ' ' + time_UTC_mod
        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                     range_max=range_max)
    # ----------------------------------------------------------------------- #
    # SAVE                                                                    #
    # ----------------------------------------------------------------------- #
    str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm)
    str_location = '_'.join([location.upper(), date, str(time_hhmm)])
    file_out = 'BOTH_PPI_' + str(elevation_deg) + \
               '°_' + str_location + '_' + str_mod + '.' + pdf_or_png
    plt.tight_layout()
    if not os.path.exists(folder_plot):
        os.makedirs(folder_plot)

    plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
                transparent=True)
    plt.close()
    # ----------------------------------------------------------------------- #
