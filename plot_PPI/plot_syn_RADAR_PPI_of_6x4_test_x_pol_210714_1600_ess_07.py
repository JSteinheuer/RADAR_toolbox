#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 30.10.25                                                 #
# plot_syn_RADAR_PPI_6x5.py                                                   #
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
time_hhmm = 1600
location = 'ess'  # for obs
elevation_deg = 12
range_max = 25
elevation_deg = 0.5
range_max = 150
elevation_deg = 5.5
range_max = 50
elevation_deg = 1.5
range_max = 100
elevation_deg = 3.5
range_max = 100
mode = 'vol'
# ----------------------------------------------------------------------- #
# initialize syn                                                          #
# ----------------------------------------------------------------------- #
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.3'
icon_emvorado_runs = [
    'MAIN_2411.3/EMVO_00510000.2',
    'MAIN_2411.3/EMVO_005100NN.2',
    'MAIN_2411.3/EMVO_00510010.2',
    'MAIN_2411.3/EMVO_00510020.2',
    'MAIN_2411.3/EMVO_00510030.2',
    'MAIN_2411.3/EMVO_00510040.2',
]
spin_up_mm = '120'
short_names = [ i[-8:] for i in icon_emvorado_runs]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
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
# initialize plot                                                         #
# ----------------------------------------------------------------------- #
n_rows = icon_emvorado_runs.__len__()+1
n_cols = 4
factor=1
factor=1.1
fig = plt.figure(figsize=(factor*5 * n_cols, factor*4 * n_rows))
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
letters=('abcdefghijklmnopqrstuvwxyz'
         '\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03B9'
         '\u03BA\u03BB\u03BC\u03BD\u03BE\u03BF'
         '\u03C1\u03C2\u03C3\u03C4\u03C5\u03C6\u03C7\u03C8\u03C9')
col=0
row=0


# --------------------------------------------------------------------------- #
# rows mod                                                                    #
# --------------------------------------------------------------------------- #
for icon_emvorado_run in icon_emvorado_runs:
    path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                     icon_emvorado_run, str(spin_up_mm) +
                     'min_spinup/EMV_Vol_' + location.upper() + '_' +
                     date + hhmm_start + '_' + date + hhmm_end + '.nc'])
    # path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
    #                       icon_run, '/ICONdata', str(spin_up_mm) +
    #                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
    #                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
    syn_nc = xr.open_dataset(path)
    # syn_nc2 = xr.open_dataset(path_icon)
    # syn_nc2 = calc_vradh(syn_nc2)
    dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                        freq="5min", inclusive='both')
    time_i_mod = np.where(dti == pd.to_datetime(
        date + str(time_hhmm).zfill(4)))[0][0]
    time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
    nc_file_comb=syn_nc
    # nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
    #                          syn_nc2.sel(elevation=elevation_deg)])
    nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
    nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
    time_i_mod=0

    col=0
    row=row+1
    moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
    cmaps = [header.cmap_radar, header.cmap_radar,
             header.cmap_radar, header.cmap_radar, ]
    norms = [header.norm_zh, header.norm_zdr,
             header.norm_kdp, header.norm_rhohv,]
    levelss = [header.levels_zh, header.levels_zdr,
               header.levels_kdp, header.levels_rhohv,]
    for moment, cmap, norm, levels in zip(
            moments, cmaps, norms, levelss):
        col=col+1
        n_i=(row-1)*n_cols+col
        ax = plt.subplot(n_rows, n_cols, n_i)
        letters_i = col + (row-1)*(n_cols)
        if mode == 'vol':
            title = (letters[letters_i] +') ' + icon_emvorado_run[-4:] +
                     ': ' + moment[:-3] + ' ' +
                     str(elevation_deg) + '° ' +
                     location.upper() + ' ' + time_UTC_mod)
        else:
            title = (letters[letters_i] +') ' + icon_emvorado_run[-4:] +
                     ': ' + moment[:-3] + ' pcp ' +
                     location.upper() + ' ' + time_UTC_mod)
        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment,
                     cmap=cmap,
                     norm=norm,
                     levels=levels,
                     title=title,range_max=range_max)

# ----------------------------------------------------------------------- #
# row 1: obs                                                              #
# ----------------------------------------------------------------------- #

row = row + 1
col=0
moments = ['ZH_AC', 'ZDR_AC_OC', 'KDP_NC', 'RHOHV_NC2P', ]
for moment in moments:
    col=col+1
    n_i=(row-1)*n_cols+col
    ax = plt.subplot(n_rows, n_cols, n_i)
    letters_i = col + (row-1)*(n_cols)
    if mode == 'vol':
        title = (letters[letters_i] + ') obs: ' + moment[:-3] + ' ' +
                 str(elevation_deg) + '° ' +
                 location.upper() + ' ' + time_UTC_obs)
    else:
        title = (letters[letters_i] + ') obs: ' + moment[:-3] + ' pcp ' +
                     location.upper() + ' ' + time_UTC_obs)
    plot_PPI(nc_file, ax, time_i_obs, moment, title=title,
             range_max=range_max)

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm)
str_location = '_'.join([location.upper(), date, str(time_hhmm)])
file_out = 'threshold_tests_PPI_' + str(elevation_deg) + \
           '°_' + str_location + '_' + str_mod + '.' + pdf_or_png
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=False, dpi=300, bbox_inches='tight')
plt.close()
# ----------------------------------------------------------------------- #

