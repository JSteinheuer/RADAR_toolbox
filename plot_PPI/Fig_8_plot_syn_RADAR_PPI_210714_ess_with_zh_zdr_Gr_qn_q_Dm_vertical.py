#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 02.03.26                                                 #
# plot_syn_RADAR_PPI_6x5.py                                                   #
#                                                                             #
# [...] Description here [...]                                                #
# --------------------------------------------------------------------------- #

import sys
for entry in sys.path.copy():
    if '/RADAR_toolbox/' in entry:
        entry_folders = entry.split('/')
        index_mother = entry_folders.index('RADAR_toolbox') + 1
        sys.path.extend(['/'.join(entry_folders[:index_mother])])

import HEADER_RADAR_toolbox as header
import ColorBlindFriendlyRadarColorMaps as radar_colors
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
n_rows = 5
n_cols = 2
# factor = .7
# fig = plt.figure(figsize=(factor *n_cols * 3.5, factor *n_rows * 2.65), layout='constrained')
factor = .75
fig = plt.figure(figsize=(factor *n_cols * 3.46, factor *n_rows * 2.65), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols)
axs = gs.subplots(sharex=True, sharey=True)
pdf_or_png = 'png'
# pdf_or_png = 'pdf'
folder_plot = header.folder_ppi_plot
letters = ('abcdefghijklmnopqrstuvwxyz'
           '\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03B9'
           '\u03BA\u03BB\u03BC\u03BD\u03BE\u03BF'
           '\u03C1\u03C2\u03C3\u03C4\u03C5\u03C6\u03C7\u03C8\u03C9')
col = 0
row = 0

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
# col = 0
# row = row + 1
row = 0
col = col + 1
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 0.5
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.3'
emvorado_run = 'EMVO_20510000.2'
icon_emvorado_run = icon_run + '/' + emvorado_run
short_name = 'I2E3'
spin_up_mm = '120'
qn_i = 1e-7
q_i = 1e-7
qn_c = 1e-7
q_c = 1e-7
qn_r = 1e-7
q_r = 1e-7
qn_s = 1e-7
q_s = 1e-7
qn_g = 1e-7
q_g = 1e-7
qn_h = 1e-7
q_h = 1e-7
# ------------------------------------ #
short_name = 'I2E3'
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
# nc_file_comb = syn_nc
path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
                      icon_run, '/ICONdata', str(spin_up_mm) +
                      'min_spinup/ICON_Vol_' + location.upper() + '_' +
                      date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
                         syn_nc2.sel(elevation=elevation_deg)], compat='override')
# nc_file_comb = nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb = nc_file_comb.isel(time=[time_i_mod, time_i_mod + 1])
# ----------------------------------------------------------------------- #
# row 3: hm Dm                                                            #
# ----------------------------------------------------------------------- #
for hm in ['graupel', 'ice', 'rain', 'hail', 'cloud', 'snow']:
    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], 0)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], 0)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], 0)
    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], 0)

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
    nc_file_comb['vol_qn' + hm[0]] = (
        ['time', 'range', 'azimuth', ], np.log10(qn_dens[hm]/1000), dict(
            standard_name=nc_file_comb['qn' + hm[0]].standard_name +
                          ' per volume', units='log_10(L-1)'))

    nc_file_comb['D0_' + hm[0]] = (
        ['time', 'range', 'azimuth', ],
        mean_volume_diameter[hm] * 1000,
        dict(standard_name='mean volume diameter of ' +
                           nc_file_comb['qn' + hm[0]].standard_name[21:],
             units='mm'))
    nc_file_comb['vol_qn' + hm[0]]=nc_file_comb['vol_qn' + hm[0]].where(
        nc_file_comb['D0_' + hm[0]]>0)
    nc_file_comb['vol_q' + hm[0]]=nc_file_comb['vol_q' + hm[0]].where(
        nc_file_comb['D0_' + hm[0]]>0)

    nc_file_comb['D0_' + hm[0]]=nc_file_comb['D0_' + hm[0]].where(
        nc_file_comb['vol_q' + hm[0]]>0)
    nc_file_comb['vol_qn' + hm[0]]=nc_file_comb['vol_qn' + hm[0]].where(
        nc_file_comb['vol_q' + hm[0]]>0)

    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], np.nan)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], np.nan)

    nc_file_comb['q' + hm[0]]=np.log10(nc_file_comb['q' + hm[0]])
    nc_file_comb['q' + hm[0]]=nc_file_comb['q' + hm[0]].assign_attrs(dict(units='log_10(kg kg-1)'))

    nc_file_comb['qn' + hm[0]]=np.log10(nc_file_comb['qn' + hm[0]])
    nc_file_comb['qn' + hm[0]]=nc_file_comb['qn' + hm[0]].assign_attrs(dict(units='log_10(kg-1)'))

time_i_mod = 0
# ------------------------------------- #
# moments = ['zrsim', 'zdrsim', 'D0_g', 'vol_qg', 'vol_qng']
moments = ['zrsim', 'zdrsim', 'D0_g', 'qg', 'qng']
# cmaps = [radar_colors.cmap_radar, radar_colors.cmap_radar, radar_colors.cmap_radar_dm,
#          radar_colors.cmap_radar_cont, radar_colors.cmap_radar_nt2, ]
# norms = [radar_colors.norm_zh, radar_colors.norm_zdr, radar_colors.norm_dm,
#          radar_colors.norm_cont, radar_colors.norm_nt2, ]
# levelss = [radar_colors.levels_zh, radar_colors.levels_zdr,radar_colors.levels_dm,
#            radar_colors.levels_cont, radar_colors.levels_nt2, ]
cmaps = [radar_colors.cmap_radar, radar_colors.cmap_radar, radar_colors.cmap_radar_dm2,
         radar_colors.cmap_radar_nt2, radar_colors.cmap_radar_nt2, ]
norms = [radar_colors.norm_zh, radar_colors.norm_zdr, radar_colors.norm_dm2,
         radar_colors.norm_nt2, radar_colors.norm_qnt2, ]
levelss = [radar_colors.levels_zh, radar_colors.levels_zdr,radar_colors.levels_dm2,
           radar_colors.levels_nt2, radar_colors.levels_qnt2, ]
labels = ['$Z_{H}$ [dBZ]', '$Z_{DR}$ [dB]', '$D_{m,\,g}$ [mm]',
          '$q_g$ [$log_{10}(kg\,\,kg^{-1})$]', '$qn_{g}$ [$log_{10}(kg^{-1})$]', ]
for moment, cmap, norm, levels, label in zip(
            moments, cmaps, norms, levelss, labels):
    # col = col + 1
    row = row + 1
    n_i = (row - 1) * n_cols + col
    # n_i = (col - 1) * n_cols + row
    ax = plt.subplot(n_rows, n_cols, n_i)
    letters_i = - 1 + col + (row - 1) * (n_cols)
    if mode == 'vol':
        title = (letters[letters_i] + ') ' + icon_emvorado_run[10:12] +
                 icon_emvorado_run[-10:-3] +
                 ': ' + moment + ' ' +
                 str(elevation_deg) + '° ' +
                 location.upper() + ' ' + time_UTC_mod)
    else:
        title = (letters[letters_i] + ') ' + icon_emvorado_run[10:12] +
                 icon_emvorado_run[-10:-3] +
                 ': ' + moment + ' pcp ' +
                 location.upper() + ' ' + time_UTC_mod)
    title=''
    panel=letters[letters_i] + ') ' + short_name
    ylabel = False if col > 1 else True
    xlabel=False if row < n_rows else True
    pcm=plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                 cmap=cmap, norm=norm, levels=levels, range_max=range_max,
                     add_colorbar=False,xlabel=xlabel, ylabel=ylabel, panel= panel)



# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Initialize case and plot of one syn PPI                                     #
# --------------------------------------------------------------------------- #
# col = 0
# row = row + 1
row = 0
col = col + 1
# ------------------------------------- #
date = "20210714"
time_hhmm = 1600
location = 'ess'
elevation_deg = 0.5
range_max = 180
mode = 'vol'
# ------------------------------------- #
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.3'
emvorado_run = 'EMVO_20510840.2'
icon_emvorado_run = icon_run + '/' + emvorado_run
short_name = 'I2E4'
spin_up_mm = '120'
qn_i = 1e-0
q_i = 1e-7
qn_c = 1e-0
q_c = 1e-7
qn_r = 1e-0
q_r = 1e-7
qn_s = 1e-0
q_s = 1e-7
qn_g = 1e-0
q_g = 1e-7
qn_h = 1e-3
q_h = 1e-5
# ------------------------------------ #
short_name = 'I2E4'
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
# nc_file_comb = syn_nc
path_icon = '/'.join([header.dir_data_vol[:-1], date, da_run,
                      icon_run, '/ICONdata', str(spin_up_mm) +
                      'min_spinup/ICON_Vol_' + location.upper() + '_' +
                      date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc2 = xr.open_dataset(path_icon)
# syn_nc2 = calc_vradh(syn_nc2)
nc_file_comb = xr.merge([syn_nc.sel(elevation=elevation_deg),
                         syn_nc2.sel(elevation=elevation_deg)], compat='override')
# nc_file_comb = nc_file_comb.sel(elevation=elevation_deg)
nc_file_comb = nc_file_comb.isel(time=[time_i_mod, time_i_mod + 1])
# ----------------------------------------------------------------------- #
# row 3: hm Dm                                                            #
# ----------------------------------------------------------------------- #
for hm in ['graupel', 'ice', 'rain', 'hail', 'cloud', 'snow']:
    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], 0)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], 0)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], 0)
    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], 0)

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
    nc_file_comb['vol_qn' + hm[0]] = (
        ['time', 'range', 'azimuth', ], np.log10(qn_dens[hm]/1000), dict(
            standard_name=nc_file_comb['qn' + hm[0]].standard_name +
                          ' per volume', units='log_10(L-1)'))

    nc_file_comb['D0_' + hm[0]] = (
        ['time', 'range', 'azimuth', ],
        mean_volume_diameter[hm] * 1000,
        dict(standard_name='mean volume diameter of ' +
                           nc_file_comb['qn' + hm[0]].standard_name[21:],
             units='mm'))
    nc_file_comb['vol_qn' + hm[0]]=nc_file_comb['vol_qn' + hm[0]].where(
        nc_file_comb['D0_' + hm[0]]>0)
    nc_file_comb['vol_q' + hm[0]]=nc_file_comb['vol_q' + hm[0]].where(
        nc_file_comb['D0_' + hm[0]]>0)

    nc_file_comb['D0_' + hm[0]]=nc_file_comb['D0_' + hm[0]].where(
        nc_file_comb['vol_q' + hm[0]]>0)
    nc_file_comb['vol_qn' + hm[0]]=nc_file_comb['vol_qn' + hm[0]].where(
        nc_file_comb['vol_q' + hm[0]]>0)

    nc_file_comb['q' + hm[0]] = nc_file_comb['q' + hm[0]].where(
        nc_file_comb['q' + hm[0]] >= locals()['q_' + hm[0]], np.nan)
    nc_file_comb['qn' + hm[0]] = nc_file_comb['qn' + hm[0]].where(
        nc_file_comb['qn' + hm[0]] >= locals()['qn_' + hm[0]], np.nan)

    nc_file_comb['q' + hm[0]]=np.log10(nc_file_comb['q' + hm[0]])
    nc_file_comb['q' + hm[0]]=nc_file_comb['q' + hm[0]].assign_attrs(dict(units='log_10(kg kg-1)'))

    nc_file_comb['qn' + hm[0]]=np.log10(nc_file_comb['qn' + hm[0]])
    nc_file_comb['qn' + hm[0]]=nc_file_comb['qn' + hm[0]].assign_attrs(dict(units='log_10(kg-1)'))

time_i_mod = 0
# ------------------------------------- #
# moments = ['zrsim', 'zdrsim', 'D0_g', 'vol_qg', 'vol_qng']
moments = ['zrsim', 'zdrsim', 'D0_g', 'qg', 'qng']
# cmaps = [radar_colors.cmap_radar, radar_colors.cmap_radar, radar_colors.cmap_radar_dm,
#          radar_colors.cmap_radar_cont, radar_colors.cmap_radar_nt2, ]
# norms = [radar_colors.norm_zh, radar_colors.norm_zdr, radar_colors.norm_dm,
#          radar_colors.norm_cont, radar_colors.norm_nt2, ]
# levelss = [radar_colors.levels_zh, radar_colors.levels_zdr,radar_colors.levels_dm,
#            radar_colors.levels_cont, radar_colors.levels_nt2, ]
cmaps = [radar_colors.cmap_radar, radar_colors.cmap_radar, radar_colors.cmap_radar_dm2,
         radar_colors.cmap_radar_nt2, radar_colors.cmap_radar_nt2, ]
norms = [radar_colors.norm_zh, radar_colors.norm_zdr, radar_colors.norm_dm2,
         radar_colors.norm_nt2, radar_colors.norm_qnt2, ]
levelss = [radar_colors.levels_zh, radar_colors.levels_zdr,radar_colors.levels_dm2,
           radar_colors.levels_nt2, radar_colors.levels_qnt2, ]
labels = ['$Z_{H}$ [dBZ]', '$Z_{DR}$ [dB]', '$D_{m,\,g}$ [mm]',
          '$q_g$ [$log_{10}(kg\,\,kg^{-1})$]', '$qn_{g}$ [$log_{10}(kg^{-1})$]', ]
for moment, cmap, norm, levels, label in zip(
            moments, cmaps, norms, levelss, labels):
    # col = col + 1
    row = row + 1
    n_i = (row - 1) * n_cols + col
    # n_i = (col - 1) * n_cols + row
    ax = plt.subplot(n_rows, n_cols, n_i)
    letters_i = - 1 + col + (row - 1) * (n_cols)
    if mode == 'vol':
        title = (letters[letters_i] + ') ' + icon_emvorado_run[10:12] +
                 icon_emvorado_run[-10:-3] +
                 ': ' + moment + ' ' +
                 str(elevation_deg) + '° ' +
                 location.upper() + ' ' + time_UTC_mod)
    else:
        title = (letters[letters_i] + ') ' + icon_emvorado_run[10:12] +
                 icon_emvorado_run[-10:-3] +
                 ': ' + moment + ' pcp ' +
                 location.upper() + ' ' + time_UTC_mod)

    title=''
    panel = letters[letters_i] + ') ' + short_name
    ylabel = False if col > 1 else True
    xlabel = False if row < n_rows else True
    pcm = plot_syn_PPI(nc_file_comb, ax, time_i_mod, moment, title=title,
                       cmap=cmap, norm=norm, levels=levels,
                       range_max=range_max,
                       add_colorbar=False, xlabel=xlabel, ylabel=ylabel,
                       panel=panel)
    cbar=fig.colorbar(pcm,  location='right', label=label, fraction=.048,
                      ticks=levels)
    # cbar.ax.tick_params(rotation=90)

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm)
str_location = '_'.join([location.upper(), date, str(time_hhmm)])
file_out = 'Fig_8_v_' + str(elevation_deg) + \
           '°_' + str_location + '_' + str_mod + '.' + pdf_or_png

# plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=False, dpi=400, bbox_inches='tight')
plt.close()
# --------------------------------------------------------------------------- #
