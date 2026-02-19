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
# Initialize plot                                                             #
# --------------------------------------------------------------------------- #
n_rows = 7
n_cols = 3
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
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=1
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510000.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=1
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510000.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=1
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] +") " +
    "default (RUC 2.0 E3):\n" +
    "q_x > 1e-7\n" +
    "qn_x > 1e-7\n",
    ha="center",
    va="center",
    fontsize=24
)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=2
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510050.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=2
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510050.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=2
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] +") " +
    "december try 50:\n" +
    "q_x > 1e-7\n"+
    "qn_x > 1e-0\n",
    ha="center",
    va="center",
    fontsize=24
)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=3
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510060.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=3
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510060.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=3
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] +") " +
    "december try 60:\n" +
    "q_x > 1e-7\n"+
    "qn_x > 1e-1\n",
    ha="center",
    va="center",
    fontsize=24
)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=4
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510810.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=4
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510810.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=4
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] + ") " +
    "new 81:\n" +
    "q_x > 1e-7\n" +
    "qn_x > 1e-0 and " +
    "qn_h > 1e-3\n",
    ha="center",
    va="center",
    fontsize=24
)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=5
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510820.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=5
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510820.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=5
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] + ") " +
    "new 82:\n" +
    "q_x > 1e-5\n" +
    "qn_x > 1e-0 and " +
    "qn_h > 1e-3\n",
    ha="center",
    va="center",
    fontsize=24
)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=6
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510830.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=6
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510830.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=6
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] + ") " +
    "new 83:\n" +
    "q_x > 1e-6\n" +
    "qn_x > 1e-0 and " +
    "qn_h > 1e-3\n",
    ha="center",
    va="center",
    fontsize=24
)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=1
row=7
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510840.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=0
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
col=2
row=7
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 12
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_20510840.2'
spin_up_mm = '120'
# ------------------------------------ #
short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# ------------------------------------- #
# ------------------------------------- #
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
nc_file_comb=syn_nc
# path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                       icon_run, '/ICONdata', str(spin_up_mm) +
#                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
#                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
# nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
#                          syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0
# ------------------------------------- #
moment_i=1
moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
cmaps = [header.cmap_radar, header.cmap_radar,
         header.cmap_radar, header.cmap_radar, ]
norms = [header.norm_zh, header.norm_zdr,
         header.norm_kdp, header.norm_rhohv,]
levelss = [header.levels_zh, header.levels_zdr,
           header.levels_kdp, header.levels_rhohv,]
moment = moments[moment_i]
cmap = cmaps[moment_i]
norm = norms[moment_i]
levels =levelss[moment_i]
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
letters_i = - 1 + col + (row-1)*(n_cols)
if mode == 'vol':
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' ' +
             str(elevation_deg) + '° ' +
             location.upper() + ' ' + time_UTC_mod)
else:
    title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
             icon_emvorado_run[-10:-3] +
             ': ' + moment[:-3] + ' pcp ' +
             location.upper() + ' ' + time_UTC_mod)
plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
             cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# ------------------------------------- #
# text add right
# ------------------------------------- #
col=3
row=7
# ------------------------------------- #
n_i=(row-1)*n_cols+col
ax = plt.subplot(n_rows, n_cols, n_i)
ax.axis('off')
letters_i = - 1 + col + (row-1)*(n_cols)
ax.text(
    0.5, 0.5,
    letters[letters_i] + ") " +
    "new 84:\n" +
    "q_x > 1e-7 and " +
    "q_h > 1e-5\n" +
    "qn_x > 1e-0 and " +
    "qn_h > 1e-3\n",
    ha="center",
    va="center",
    fontsize=24
)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one obs PPI                                     #
# --------------------------------------------------------------------------- #
# row = 1
# col=1
# date = "20210714"
# time_hhmm = 1600
# location = 'ess'
# elevation_deg = 12
# range_max = 50
# mode = 'vol'
# # ------------------------------------- #
# short_name = 'obs'
# year = date[0:4]
# mon = date[4:6]
# day = date[6:8]
# hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
# hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
# date_start = '-'.join([year, mon, day, hhmm_start])
# date_end = '-'.join([year, mon, day, hhmm_end])
# sweep = '0' + str(np.where(
#     header.ELEVATIONS_ALL == float(elevation_deg))[0][0])
# time_i_obs = int(time_hhmm * 0.12)
# folder_in = "/".join([header.dir_data_obs + '*', year,
#                       year + '-' + mon,
#                       year + '-' + mon + '-' + day,
#                       location, mode + '*', sweep])
# nc_file = glob.glob(folder_in + '/*polmoms*')
# if len(nc_file) > 1:
#     print('polmom: too many files')
# elif len(nc_file) == 0:
#     print('polmom: no files')
# else:
#     nc_file = nc_file[0]
#
# time_start = date + nc_file.split('-')[-4][-4:]
# time_end = date + nc_file.split('-')[-3][-4:]
# dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
# dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
# dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
# time_UTC_obs = dti[time_i_obs].strftime('%Y-%m-%d %H:%M')
# # ------------------------------------- #
# # ------------------------------------- #
# moment_i=0
# moments = ['ZH_AC', 'ZDR_AC_OC', 'KDP_NC', 'RHOHV_NC2P', ]
# cmaps = [header.cmap_radar, header.cmap_radar,
#          header.cmap_radar, header.cmap_radar, ]
# norms = [header.norm_zh, header.norm_zdr,
#          header.norm_kdp, header.norm_rhohv,]
# levelss = [header.levels_zh, header.levels_zdr,
#            header.levels_kdp, header.levels_rhohv,]
# moment = moments[moment_i]
# cmap = cmaps[moment_i]
# norm = norms[moment_i]
# levels =levelss[moment_i]
# n_i=(row-1)*n_cols+col
# ax = plt.subplot(n_rows, n_cols, n_i)
# letters_i = - 1 + col + (row-1)*(n_cols)
# if mode == 'vol':
#     title = (letters[letters_i] + ') obs: ' + moment[:-3] + ' ' +
#              str(elevation_deg) + '° ' +
#              location.upper() + ' ' + time_UTC_obs)
# else:
#     title = (letters[letters_i] + ') obs: ' + moment[:-3] + ' pcp ' +
#                  location.upper() + ' ' + time_UTC_obs)
# plot_PPI(nc_file, ax, time_i_obs, moment, title=title,
#          cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
# col=2
# row=2
# # ------------------------------------- #
# date = "20210714"
# time_hhmm = 1600
# location = 'ess'
# elevation_deg = 12
# range_max = 50
# mode = 'vol'
# # ------------------------------------- #
# da_run = 'ASS_2411'
# icon_emvorado_run = 'MAIN_2411.3/EMVO_20510000.2'
# spin_up_mm = '120'
# # ------------------------------------ #
# short_name =  icon_emvorado_run[-17:-15]+ icon_emvorado_run[-10:]
# year = date[0:4]
# mon = date[4:6]
# day = date[6:8]
# hhmm_start = str(time_hhmm - time_hhmm % 600).zfill(4)
# hhmm_end = str(time_hhmm - time_hhmm % 600 + 555).zfill(4)
# date_start = '-'.join([year, mon, day, hhmm_start])
# date_end = '-'.join([year, mon, day, hhmm_end])
# dti = pd.date_range(date + hhmm_start, date + hhmm_end,
#                     freq="5min", inclusive='both')
# time_i_mod = np.where(dti == pd.to_datetime(
#     date + str(time_hhmm).zfill(4)))[0][0]
# time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
# # ------------------------------------- #
# # ------------------------------------- #
# path = '/'.join([header.dir_data_vol[:-1], date, da_run,
#                  icon_emvorado_run, str(spin_up_mm) +
#                  'min_spinup/EMV_Vol_' + location.upper() + '_' +
#                  date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# syn_nc = xr.open_dataset(path)
# nc_file_comb=syn_nc
# # path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
# #                       icon_run, '/ICONdata', str(spin_up_mm) +
# #                       'min_spinup/ICON_Vol_' + location.upper() + '_' +
# #                       date + hhmm_start + '_' + date + hhmm_end + '.nc'])
# # syn_nc2 = xr.open_dataset(path_icon)
# # syn_nc2 = calc_vradh(syn_nc2)
# # nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
# #                          syn_nc2.sel(elevation=elevation_deg)])
# nc_file_comb=nc_file_comb.sel(elevation=elevation_deg)
# nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
# time_i_mod=0
# # ------------------------------------- #
# moment_i=0
# moments =['zrsim','zdrsim','kdpsim', 'rhvsim',]
# cmaps = [header.cmap_radar, header.cmap_radar,
#          header.cmap_radar, header.cmap_radar, ]
# norms = [header.norm_zh, header.norm_zdr,
#          header.norm_kdp, header.norm_rhohv,]
# levelss = [header.levels_zh, header.levels_zdr,
#            header.levels_kdp, header.levels_rhohv,]
# moment = moments[moment_i]
# cmap = cmaps[moment_i]
# norm = norms[moment_i]
# levels =levelss[moment_i]
# n_i=(row-1)*n_cols+col
# ax = plt.subplot(n_rows, n_cols, n_i)
# letters_i = - 1 + col + (row-1)*(n_cols)
# if mode == 'vol':
#     title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
#              icon_emvorado_run[-10:-3] +
#              ': ' + moment[:-3] + ' ' +
#              str(elevation_deg) + '° ' +
#              location.upper() + ' ' + time_UTC_mod)
# else:
#     title = (letters[letters_i] +') ' + icon_emvorado_run[10:12]+
#              icon_emvorado_run[-10:-3] +
#              ': ' + moment[:-3] + ' pcp ' +
#              location.upper() + ' ' + time_UTC_mod)
# plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
#              cmap=cmap, norm=norm, levels=levels, range_max=range_max)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm)
str_location = '_'.join([location.upper(), date, str(time_hhmm)])
# file_out = ('7x2_PPI_zh_05_2last_try_2runs_' + str_location + '_' + str_mod +
#             '.' + pdf_or_png)
file_out = ('7x2_PPI_zh_12_last_try_2runs_' + str_location + '_' + str_mod +
            '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=False, dpi=300, bbox_inches='tight')
plt.close()
# --------------------------------------------------------------------------- #
