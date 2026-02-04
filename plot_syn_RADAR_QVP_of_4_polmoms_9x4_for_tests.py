#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 13.02.25                                                 #
# plot_syn_RADAR_QVP_of_4_polmoms.py                                          #
#                                                                             #
# Run the QVP functions in PLOT_SYN_RADAR.py for generating specific QVP plot.#
# --------------------------------------------------------------------------- #

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list, mom_plot_dict
from SET_SYN_RADAR import rad_dict
from statsmodels.stats.weightstats import DescrStatsW

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #

# ------------------------------------ #
# QVPs                                 #
# ------------------------------------ #
location = 'ESS'
date = '20210714'
hhmm_start_qvp = '16:00'
hhmm_end_qvp = '17:00'
year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start_qvp])
date_end = '-'.join([year, mon, day, hhmm_end_qvp])
top_height = 12
mode = 'vol'
elevation_deg = 12
# elevation_deg = 17
# elevation_deg = 5.5
# elevation_deg = 0.5
# for elevation_deg in [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]:
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
letters=('abcdefghijklmnopqrstuvwxyz'
         '\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03B9'
         '\u03BA\u03BB\u03BC\u03BD\u03BE\u03BF'
         '\u03C1\u03C2\u03C3\u03C4\u03C5\u03C6\u03C7\u03C8\u03C9')
filter_moms = False

# ------------------------------------ #
# MODELS                               #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
# ------------------------------------ #
# SYN data rows 2-8                    #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('00.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510010.2')
spin_up_mms.append('120')
short_names.append('10.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510020.2')
spin_up_mms.append('120')
short_names.append('20.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510030.2')
spin_up_mms.append('120')
short_names.append('30.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510040.2')
spin_up_mms.append('120')
short_names.append('40.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510050.2')
spin_up_mms.append('120')
short_names.append('50.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510060.2')
spin_up_mms.append('120')
short_names.append('60.2')
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510070.2')
spin_up_mms.append('120')
short_names.append('70.2')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# QVPs A: no filter                                                           #
# --------------------------------------------------------------------------- #
filter_entr = False
filter_entr_at = 0

# ------------------------------------ #
# QVPs plot parameters                 #
# ------------------------------------ #
folder_plot = header.folder_qvp_plot
mod_names = ''
letters_i = 0
n_rows = len(da_runs)+1
n_cols = 4
current_row = -1
current_col = -1
scale_font = 1.
scale_numbers = 1.
fig = plt.figure(figsize=(n_cols * 2.7, n_rows * 2.7), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.002,wspace=0.002)
axs = gs.subplots(sharex=True, sharey=True)
# ------------------------------------ #

# --------------------------------------------------------------------------- #
# QVPs OBS row 1                                                              #
# --------------------------------------------------------------------------- #
path_obs = "/".join([header.dir_obs_qvp + '*',
                     year, year + '-' + mon,
                     year + '-' + mon + '-' + day,
                     location.lower(), mode + '*', sweep, ])
path_obs = sorted(glob.glob(path_obs))
if len(path_obs) == 1:
    path_obs = path_obs[0]
    path_obs = glob.glob(path_obs + '_'.join(
        ['/ras*qvp*', '*_polmoms_nc_*' + date + '*' + location.lower() + '*']))
    if len(path_obs) == 1:
        path_obs = path_obs[0]
        print(path_obs)
    else:
        print('no file or too many -> return')
else:
    print('no folder or too many-> return')

obs_nc = xr.open_dataset(path_obs)
if filter_entr:
    obs_nc=obs_nc.where(obs_nc['min_entropy']>filter_entr_at)

# ------------------------------------ #
qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zh_obs = obs_nc['ZH_AC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_temp_obs = obs_nc['temp'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time') - 273.15

# ------------------------------------ #
current_row = 0
# ------------------------------------ #
current_col = 0
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    # mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$Z_{H}$ [dBZ]',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
current_col = current_col + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    # mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$Z_{DR}$ [dB]',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    ylabel='',
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
current_col = current_col + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_obs,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    # mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$K_{DP}$ [°/km]',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    ylabel='',
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
current_col = current_col + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_obs,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    # mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$\u03C1_{HV}$ [1]',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    ylabel='',
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
obs_nc.close()

# --------------------------------------------------------------------------- #
# QVPs SYN row i                                                              #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm, short_name in zip(
        da_runs, icon_emvorado_runs, spin_up_mms, short_names):
    date_start = '-'.join([year, mon, day, hhmm_start_qvp])
    date_end = '-'.join([year, mon, day, hhmm_end_qvp])
    path_mod = '/'.join([header.dir_data_qvp + date, da_run, icon_emvorado_run,
                         str(spin_up_mm) + 'min_spinup', 'QVP_' +
                         str(elevation_deg) + '_Syn_' + location + '_' +
                         date + '0000_' + date + '2355.nc'])
    path_mod = sorted(glob.glob(path_mod))
    if len(path_mod) == 1:
        path_mod = path_mod[0]
        syn_nc = xr.open_dataset(path_mod)
        print(path_mod)
    else:
        path_mod = '/'.join([header.dir_data_qvp + date, da_run,
                             icon_emvorado_run,
                             str(spin_up_mm) + 'min_spinup', 'QVP_' +
                             str(elevation_deg) + '_Syn_' + location + '_' +
                             date + '*_' + date + '*.nc'])
        path_mod = sorted(glob.glob(path_mod))
        if (len(path_mod) < 4) and (len(path_mod) > 0):
            syn_nc = xr.merge([xr.open_dataset(syn_nc_i) for
                               syn_nc_i in path_mod])
            print(path_mod)
        else:
            print('no file or too many -> return')
            continue

    if filter_entr:
        syn_nc = syn_nc.where(syn_nc['min_entropy'] > filter_entr_at)

    model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
                                icon_emvorado_run.split('/')[1][5:]])
    mod_names = '_'.join([mod_names, model_name_file])
    model_name = '-'.join([da_run[4:],
                           icon_emvorado_run.split('/')[0][5:],
                           icon_emvorado_run.split('/')[1][5:],
                           spin_up_mm + 'min'])
    # ------------------------------------ #
    h_syn = syn_nc['height']
    qvp_zh_syn = syn_nc['zrsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_zdr_syn = syn_nc['zdrsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_rho_syn = syn_nc['rhvsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_kdp_syn = syn_nc['kdpsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_d0q_syn = syn_nc['D0_r'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')  # * 1000
    qvp_temp_syn = syn_nc['temp'].sel(
        time=slice(date_start, date_end)).transpose(..., 'time') - 273.15
    # ------------------------------------ #
    current_row = current_row +1
    add_colorbar = False
    xlabel = None
    if current_row == n_rows - 1:
        add_colorbar = True
        xlabel='UTC [mm-dd hh]'

    # ------------------------------------ #
    current_col = 0
    plot_qvp_of_polarimetric_variable(
        mom=qvp_zh_syn,
        cmap=header.cmap_radar,
        norm=header.norm_zh,
        levels=header.levels_zh,
        # mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 5),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$Z_{H}$ [dBZ]',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_zdr_syn,
        cmap=header.cmap_radar,
        norm=header.norm_zdr,
        levels=header.levels_zdr,
        # mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 5),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$Z_{DR}$ [dB]',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_kdp_syn,
        cmap=header.cmap_radar,
        norm=header.norm_kdp,
        levels=header.levels_kdp,
        # mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 5),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$K_{DP}$ [°/km]',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_rho_syn,
        cmap=header.cmap_radar,
        norm=header.norm_rhohv,
        levels=header.levels_rhohv,
        # mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 5),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$\u03C1_{HV}$ [1]',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    syn_nc.close()

# --------------------------------------------------------------------------- #
# QVPs SAVE                                                                   #
# --------------------------------------------------------------------------- #
# hh_at=[2,4,6,8,10,12,14,16,18]
#hh_at=[15,16,17,18]
hh_25=np.linspace(axs[-1,-1].get_xticks()[0],
                   axs[-1,-1].get_xticks()[4],5,endpoint=True)
str_hh_at=['16:10', '16:20', '16:30', '16:40', '16:50']
axs[-1,-1].set_xticks(hh_25,str_hh_at)
axs[-1,-1].set_xlabel('UTC [hh:mm]', fontsize=12 * scale_font)
axs[-1,-2].set_xticks(hh_25,str_hh_at)
axs[-1,-2].set_xlabel('UTC [hh:mm]', fontsize=12 * scale_font)
axs[-1,-3].set_xticks(hh_25,str_hh_at)
axs[-1,-3].set_xlabel('UTC [hh:mm]', fontsize=12 * scale_font)
axs[-1,-4].set_xticks(hh_25,str_hh_at)
axs[-1,-4].set_xlabel('UTC [hh:mm]', fontsize=12 * scale_font)
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot +
    '/QVPs_' + str(n_rows) + 'x4polmoms_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location + '_' +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x4polmoms_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.pdf', format='pdf', transparent=False)
# plt.close()

# # --------------------------------------------------------------------------- #
# # --------------------------------------------------------------------------- #
# # QVPs B: filter at 0.8                                                       #
# # --------------------------------------------------------------------------- #
# filter_entr = True
# filter_entr_at = 0.8
#
# # ------------------------------------ #
# # QVPs plot parameters                 #
# # ------------------------------------ #
# folder_plot = header.folder_qvp_plot
# mod_names = ''
# letters_i = 0
# n_rows = len(da_runs)+1
# n_cols = 4
# current_row = -1
# current_col = -1
# scale_font = 1.
# scale_numbers = 1.
# fig = plt.figure(figsize=(n_cols * 2.7, n_rows * 2.7), layout='constrained')
# gs = fig.add_gridspec(n_rows, n_cols, hspace=0.002,wspace=0.002)
# axs = gs.subplots(sharex=True, sharey=True)
# # ------------------------------------ #
#
# # --------------------------------------------------------------------------- #
# # QVPs OBS row 1                                                              #
# # --------------------------------------------------------------------------- #
# path_obs = "/".join([header.dir_obs_qvp + '*',
#                      year, year + '-' + mon,
#                      year + '-' + mon + '-' + day,
#                      location.lower(), mode + '*', sweep, ])
# path_obs = sorted(glob.glob(path_obs))
# if len(path_obs) == 1:
#     path_obs = path_obs[0]
#     path_obs = glob.glob(path_obs + '_'.join(
#         ['/ras*qvp*', '*_polmoms_nc_*' + date + '*' + location.lower() + '*']))
#     if len(path_obs) == 1:
#         path_obs = path_obs[0]
#         print(path_obs)
#     else:
#         print('no file or too many -> return')
# else:
#     print('no folder or too many-> return')
#
# obs_nc = xr.open_dataset(path_obs)
# if filter_entr:
#     obs_nc=obs_nc.where(obs_nc['min_entropy']>filter_entr_at)
#
# # ------------------------------------ #
# qvp_kdp_obs = obs_nc['KDP_NC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_zh_obs = obs_nc['ZH_AC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_temp_obs = obs_nc['temp'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time') - 273.15
#
# # ------------------------------------ #
# current_row = 0
# # ------------------------------------ #
# current_col = 0
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_zh_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_zh,
#     levels=header.levels_zh,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$Z_{H}$ [dBZ]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_zdr_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_zdr,
#     levels=header.levels_zdr,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$Z_{DR}$ [dB]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_kdp_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_kdp,
#     levels=header.levels_kdp,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$K_{DP}$ [°/km]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_rho_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_rhohv,
#     levels=header.levels_rhohv,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$\u03C1_{HV}$ [1]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# obs_nc.close()
#
# # --------------------------------------------------------------------------- #
# # QVPs SYN row i                                                              #
# # --------------------------------------------------------------------------- #
# for da_run, icon_emvorado_run, spin_up_mm, short_name in zip(
#         da_runs, icon_emvorado_runs, spin_up_mms, short_names):
#     date_start = '-'.join([year, mon, day, hhmm_start_qvp])
#     date_end = '-'.join([year, mon, day, hhmm_end_qvp])
#     path_mod = '/'.join([header.dir_data_qvp + date, da_run, icon_emvorado_run,
#                          str(spin_up_mm) + 'min_spinup', 'QVP_' +
#                          str(elevation_deg) + '_Syn_' + location + '_' +
#                          date + '0000_' + date + '2355.nc'])
#     path_mod = sorted(glob.glob(path_mod))
#     if len(path_mod) == 1:
#         path_mod = path_mod[0]
#         syn_nc = xr.open_dataset(path_mod)
#         print(path_mod)
#     else:
#         path_mod = '/'.join([header.dir_data_qvp + date, da_run,
#                              icon_emvorado_run,
#                              str(spin_up_mm) + 'min_spinup', 'QVP_' +
#                              str(elevation_deg) + '_Syn_' + location + '_' +
#                              date + '*_' + date + '*.nc'])
#         path_mod = sorted(glob.glob(path_mod))
#         if (len(path_mod) < 4) and (len(path_mod) > 0):
#             syn_nc = xr.merge([xr.open_dataset(syn_nc_i) for
#                                syn_nc_i in path_mod])
#             print(path_mod)
#         else:
#             print('no file or too many -> return')
#             continue
#
#     if filter_entr:
#         syn_nc = syn_nc.where(syn_nc['min_entropy'] > filter_entr_at)
#
#     model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
#                                 icon_emvorado_run.split('/')[1][5:]])
#     mod_names = '_'.join([mod_names, model_name_file])
#     model_name = '-'.join([da_run[4:],
#                            icon_emvorado_run.split('/')[0][5:],
#                            icon_emvorado_run.split('/')[1][5:],
#                            spin_up_mm + 'min'])
#     # ------------------------------------ #
#     h_syn = syn_nc['height']
#     qvp_zh_syn = syn_nc['zrsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_zdr_syn = syn_nc['zdrsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_rho_syn = syn_nc['rhvsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_kdp_syn = syn_nc['kdpsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_d0q_syn = syn_nc['D0_r'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')  # * 1000
#     qvp_temp_syn = syn_nc['temp'].sel(
#         time=slice(date_start, date_end)).transpose(..., 'time') - 273.15
#     # ------------------------------------ #
#     current_row = current_row +1
#     add_colorbar = False
#     xlabel = None
#     if current_row == n_rows - 1:
#         add_colorbar = True
#         xlabel='UTC [mm-dd hh]'
#
#     # ------------------------------------ #
#     current_col = 0
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_zh_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_zh,
#         levels=header.levels_zh,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$Z_{H}$ [dBZ]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_zdr_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_zdr,
#         levels=header.levels_zdr,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$Z_{DR}$ [dB]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_kdp_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_kdp,
#         levels=header.levels_kdp,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$K_{DP}$ [°/km]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_rho_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_rhohv,
#         levels=header.levels_rhohv,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$\u03C1_{HV}$ [1]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     syn_nc.close()
#
# # --------------------------------------------------------------------------- #
# # QVPs SAVE                                                                   #
# # --------------------------------------------------------------------------- #
#
# hh_at=[2,4,6,8,10,12,14,16,18]
# hh_25=np.linspace(np.round(axs[-1,-1].get_xticks()[0]),
#                    np.round(axs[-1,-1].get_xticks()[0])+1,25,endpoint=True)
# str_hh_at=[str(z).zfill(2) for z in hh_at]
# axs[-1,-1].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-1].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-2].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-2].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-3].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-3].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-4].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-4].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
#
# if not os.path.exists(folder_plot):
#     os.makedirs(folder_plot)
#
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x4polmoms_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location + '_' +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x4polmoms_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.pdf', format='pdf', transparent=False)
# plt.close()
#
# # --------------------------------------------------------------------------- #
# # --------------------------------------------------------------------------- #
# # QVPs C: filter at 0.7                                                       #
# # --------------------------------------------------------------------------- #
# filter_entr = True
# filter_entr_at = 0.7
#
# # ------------------------------------ #
# # QVPs plot parameters                 #
# # ------------------------------------ #
# folder_plot = header.folder_qvp_plot
# mod_names = ''
# letters_i = 0
# n_rows = len(da_runs)+1
# n_cols = 4
# current_row = -1
# current_col = -1
# scale_font = 1.
# scale_numbers = 1.
# fig = plt.figure(figsize=(n_cols * 2.7, n_rows * 2.7), layout='constrained')
# gs = fig.add_gridspec(n_rows, n_cols, hspace=0.002,wspace=0.002)
# axs = gs.subplots(sharex=True, sharey=True)
# # ------------------------------------ #
#
# # --------------------------------------------------------------------------- #
# # QVPs OBS row 1                                                              #
# # --------------------------------------------------------------------------- #
# path_obs = "/".join([header.dir_obs_qvp + '*',
#                      year, year + '-' + mon,
#                      year + '-' + mon + '-' + day,
#                      location.lower(), mode + '*', sweep, ])
# path_obs = sorted(glob.glob(path_obs))
# if len(path_obs) == 1:
#     path_obs = path_obs[0]
#     path_obs = glob.glob(path_obs + '_'.join(
#         ['/ras*qvp*', '*_polmoms_nc_*' + date + '*' + location.lower() + '*']))
#     if len(path_obs) == 1:
#         path_obs = path_obs[0]
#         print(path_obs)
#     else:
#         print('no file or too many -> return')
# else:
#     print('no folder or too many-> return')
#
# obs_nc = xr.open_dataset(path_obs)
# if filter_entr:
#     obs_nc=obs_nc.where(obs_nc['min_entropy']>filter_entr_at)
#
# # ------------------------------------ #
# qvp_kdp_obs = obs_nc['KDP_NC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_zh_obs = obs_nc['ZH_AC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time')
# qvp_temp_obs = obs_nc['temp'].sel(
#     time=slice(date_start, date_end)).transpose(..., 'time') - 273.15
#
# # ------------------------------------ #
# current_row = 0
# # ------------------------------------ #
# current_col = 0
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_zh_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_zh,
#     levels=header.levels_zh,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$Z_{H}$ [dBZ]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_zdr_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_zdr,
#     levels=header.levels_zdr,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$Z_{DR}$ [dB]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_kdp_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_kdp,
#     levels=header.levels_kdp,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$K_{DP}$ [°/km]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# current_col = current_col + 1
# plot_qvp_of_polarimetric_variable(
#     mom=qvp_rho_obs,
#     cmap=header.cmap_radar,
#     norm=header.norm_rhohv,
#     levels=header.levels_rhohv,
#     # mom_cs=qvp_zh_obs,
#     levels_cs=np.arange(-50, 60, 5),
#     mom_cf=qvp_temp_obs,
#     levels_cf=np.arange(-50, 60, 5),
#     cbar_title='$\u03C1_{HV}$ [1]',
#     title='',
#     ax=axs[current_row, current_col],
#     mom_height_unit='km',
#     scale_font=scale_font,
#     scale_numbers=scale_numbers,
#     top_height=top_height,
#     xlabel=None,
#     ylabel='',
#     panel=letters[letters_i] + ') obs',
# )
# letters_i=letters_i+1
# # ------------------------------------ #
# obs_nc.close()
#
# # --------------------------------------------------------------------------- #
# # QVPs SYN row i                                                              #
# # --------------------------------------------------------------------------- #
# for da_run, icon_emvorado_run, spin_up_mm, short_name in zip(
#         da_runs, icon_emvorado_runs, spin_up_mms, short_names):
#     date_start = '-'.join([year, mon, day, hhmm_start_qvp])
#     date_end = '-'.join([year, mon, day, hhmm_end_qvp])
#     path_mod = '/'.join([header.dir_data_qvp + date, da_run, icon_emvorado_run,
#                          str(spin_up_mm) + 'min_spinup', 'QVP_' +
#                          str(elevation_deg) + '_Syn_' + location + '_' +
#                          date + '0000_' + date + '2355.nc'])
#     path_mod = sorted(glob.glob(path_mod))
#     if len(path_mod) == 1:
#         path_mod = path_mod[0]
#         syn_nc = xr.open_dataset(path_mod)
#         print(path_mod)
#     else:
#         path_mod = '/'.join([header.dir_data_qvp + date, da_run,
#                              icon_emvorado_run,
#                              str(spin_up_mm) + 'min_spinup', 'QVP_' +
#                              str(elevation_deg) + '_Syn_' + location + '_' +
#                              date + '*_' + date + '*.nc'])
#         path_mod = sorted(glob.glob(path_mod))
#         if (len(path_mod) < 4) and (len(path_mod) > 0):
#             syn_nc = xr.merge([xr.open_dataset(syn_nc_i) for
#                                syn_nc_i in path_mod])
#             print(path_mod)
#         else:
#             print('no file or too many -> return')
#             continue
#
#     if filter_entr:
#         syn_nc = syn_nc.where(syn_nc['min_entropy'] > filter_entr_at)
#
#     model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
#                                 icon_emvorado_run.split('/')[1][5:]])
#     mod_names = '_'.join([mod_names, model_name_file])
#     model_name = '-'.join([da_run[4:],
#                            icon_emvorado_run.split('/')[0][5:],
#                            icon_emvorado_run.split('/')[1][5:],
#                            spin_up_mm + 'min'])
#     # ------------------------------------ #
#     h_syn = syn_nc['height']
#     qvp_zh_syn = syn_nc['zrsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_zdr_syn = syn_nc['zdrsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_rho_syn = syn_nc['rhvsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_kdp_syn = syn_nc['kdpsim'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')
#     qvp_d0q_syn = syn_nc['D0_r'].sel(
#         time=slice(date_start, date_end)) \
#         .transpose(..., 'time')  # * 1000
#     qvp_temp_syn = syn_nc['temp'].sel(
#         time=slice(date_start, date_end)).transpose(..., 'time') - 273.15
#     # ------------------------------------ #
#     current_row = current_row +1
#     add_colorbar = False
#     xlabel = None
#     if current_row == n_rows - 1:
#         add_colorbar = True
#         xlabel='UTC [mm-dd hh]'
#
#     # ------------------------------------ #
#     current_col = 0
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_zh_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_zh,
#         levels=header.levels_zh,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$Z_{H}$ [dBZ]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_zdr_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_zdr,
#         levels=header.levels_zdr,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$Z_{DR}$ [dB]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_kdp_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_kdp,
#         levels=header.levels_kdp,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$K_{DP}$ [°/km]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     current_col = current_col +1
#     plot_qvp_of_polarimetric_variable(
#         mom=qvp_rho_syn,
#         cmap=header.cmap_radar,
#         norm=header.norm_rhohv,
#         levels=header.levels_rhohv,
#         # mom_cs=qvp_zh_syn,
#         levels_cs=np.arange(-50, 60, 5),
#         mom_cf=qvp_temp_syn,
#         levels_cf=np.arange(-50, 60, 5),
#         cbar_title='$\u03C1_{HV}$ [1]',
#         ax=axs[current_row, current_col],
#         scale_font=scale_font,
#         scale_numbers=scale_numbers,
#         top_height=top_height,
#         mom_height_unit='km',
#         add_colorbar=add_colorbar,
#         xlabel=xlabel,
#         ylabel='',
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     syn_nc.close()
#
# # --------------------------------------------------------------------------- #
# # QVPs SAVE                                                                   #
# # --------------------------------------------------------------------------- #
#
# hh_at=[2,4,6,8,10,12,14,16,18]
# hh_25=np.linspace(np.round(axs[-1,-1].get_xticks()[0]),
#                    np.round(axs[-1,-1].get_xticks()[0])+1,25,endpoint=True)
# str_hh_at=[str(z).zfill(2) for z in hh_at]
# axs[-1,-1].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-1].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-2].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-2].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-3].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-3].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
# axs[-1,-4].set_xticks(hh_25[hh_at],str_hh_at)
# axs[-1,-4].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
#
# if not os.path.exists(folder_plot):
#     os.makedirs(folder_plot)
#
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x4polmoms_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location + '_' +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x4polmoms_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.pdf', format='pdf', transparent=False)
# plt.close()
#
