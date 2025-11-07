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
# range_max = 55
range_max = 25
mode = 'vol'
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.3'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510000.2'
spin_up_mm = '120'
short_names = 'R2E3'

# --------------------------------------------------------------------------- #
# initialize mod                                                              #
# --------------------------------------------------------------------------- #
year = date[0:4]
mon = date[4:6]
day = date[6:8]
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
dti = pd.date_range(date + hhmm_start, date + hhmm_end,
                    freq="5min", inclusive='both')
time_i_mod = np.where(dti == pd.to_datetime(
    date + str(time_hhmm).zfill(4)))[0][0]
time_UTC_mod = dti[time_i_mod].strftime('%Y-%m-%d %H:%M')
nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
                         syn_nc2.sel(elevation=elevation_deg)])
nc_file_comb=nc_file_comb.isel(time=[time_i_mod,time_i_mod+1])
time_i_mod=0

# --------------------------------------------------------------------------- #
# calculate vol specific q qn                                                 #
# --------------------------------------------------------------------------- #
from PROCESS_SYN_RADAR import adjust_icon_fields, \
    mgdparams, calc_moments, calc_multimoments, calc_Dmean

q_dens, qn_dens = adjust_icon_fields(nc_file_comb)
multi_params = mgdparams(q_dens, qn_dens)
moments = calc_moments(mgd=multi_params)
multimoments = calc_multimoments(moments)
mean_volume_diameter = calc_Dmean(multimoments)
for hm in ['graupel', 'ice', 'rain', 'hail', 'cloud', 'snow']:
    nc_file_comb['vol_q' + hm[0]] = (
        ['time', 'range', 'azimuth', ], q_dens[hm]*1000, dict(
            standard_name='volume ' + nc_file_comb[
                'q' + hm[0]].standard_name,
            units='g m-3'))
    nc_file_comb['vol_q' + hm[0]] = nc_file_comb['vol_q' + hm[0]].where(
        nc_file_comb['vol_q' + hm[0]]>0)  # necessary as we want =0 as white!

    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]]*1000
    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['q' + hm[0]]>0)  # necessary as we want =0 as white!
    nc_file_comb['q' + hm[0]].attrs['units']='g kg-1'

    nc_file_comb['vol_qn' + hm[0]] = (
        ['time', 'range', 'azimuth', ], np.log10((qn_dens[hm]-1e-9)/1000.),
        # necessary as we want =-3 as bright yellow!
        dict(standard_name=nc_file_comb['qn' + hm[0]].standard_name +
                           ' per volume', units='log_10(L-1)'))

    nc_file_comb['qn' + hm[0]] = np.log10(nc_file_comb['qn' + hm[0]]/1000.)
    nc_file_comb['qn' + hm[0]].attrs['units']='log_10(g-1)'

    nc_file_comb['D0_' + hm[0]] = (
        ['time', 'range', 'azimuth', ],
        mean_volume_diameter[hm] * 1000,
        dict(standard_name='mean volume diameter of ' +
                           nc_file_comb['qn' + hm[0]].standard_name[21:],
             units='mm'))

# ----------------------------------------------------------------------- #
# initialize plot                                                         #
# ----------------------------------------------------------------------- #
n_rows = 7
n_cols = 5
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
letters=('abcdefghijklmnopqrstuvwxyz'
         '\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03B9'
         '\u03BA\u03BB\u03BC\u03BD\u03BE\u03BF'
         '\u03C1\u03C2\u03C3\u03C4\u03C5\u03C6\u03C7\u03C8\u03C9')
row=0
for hm in ['i','s','h', 'g', 'c','r']:
    col=0
    row=row+1
    moments = ['q', 'vol_q', 'qn', 'vol_qn', 'D0_', ]
    cmaps = [header.cmap_radar_white, header.cmap_radar_white,
             header.cmap_radar, header.cmap_radar,
             header.cmap_radar_white, ]
    norms = [header.norm_iwc, header.norm_iwc,
             header.norm_nt_iwc, header.norm_nt_iwc,
             header.norm_d0_ice,]
    levelss = [header.levels_iwc, header.levels_iwc,
               header.levels_nt_iwc, header.levels_nt_iwc,
               header.levels_d0_ice,]
    for moment, cmap, norm, levels in zip(
            moments, cmaps, norms, levelss):
        moment=moment+hm
        col=col+1
        n_i=(row-1)*n_cols+col
        ax = plt.subplot(n_rows, n_cols, n_i)
        letters_i = col + (row-1)*(n_cols)
        if mode == 'vol':
            title = (letters[letters_i] +') ' + moment + ' at ' +
                     str(elevation_deg) + '° ' +
                     location.upper() + ' ' + time_UTC_mod)
        else:
            title = (letters[letters_i] +') ' + moment + ' at pcp ' +
                     location.upper() + ' ' + time_UTC_mod)
        plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment,
                     cmap=cmap,
                     norm=norm,
                     levels=levels,
                     title=title,range_max=range_max)

# ----------------------------------------------------------------------- #
# row -1: mod pol moms                                                    #
# ----------------------------------------------------------------------- #
row=n_rows
col=0
moments = ['zrsim', 'zdrsim', 'kdpsim', 'rhvsim', 'temp']
cmaps = [header.cmap_radar,
         header.cmap_radar,
         header.cmap_radar,
         header.cmap_radar,
         None, ]
norms = [header.norm_zh,
         header.norm_zdr,
         header.norm_kdp,
         header.norm_rhohv,
         None,]
levelss = [header.levels_zh,
           header.levels_zdr,
           header.levels_kdp,
           header.levels_rhohv,
           None,]
for moment, cmap, norm, levels in zip(
        moments, cmaps, norms, levelss):
    col=col+1
    n_i=(row-1)*n_cols+col
    ax = plt.subplot(n_rows, n_cols, n_i)
    letters_i = col + (row-1)*(n_cols)
    if mode == 'vol':
        title = (letters[letters_i] +') ' + moment + ' at ' +
                 str(elevation_deg) + '° ' +
                 location.upper() + ' ' + time_UTC_mod)
    else:
        title = (letters[letters_i] +') ' + moment + ' at pcp ' +
                 location.upper() + ' ' + time_UTC_mod)
    plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                 range_max=range_max)

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm)
str_location = '_'.join([location.upper(), date, str(time_hhmm)])
file_out = 'HMs_PPI_' + str(elevation_deg) + \
           '°_' + str_location + '_' + str_mod + '.' + pdf_or_png
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=False, dpi=300, bbox_inches='tight')
plt.close()
# ----------------------------------------------------------------------- #

