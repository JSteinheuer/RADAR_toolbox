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
hhmm_start_qvp = '00:00'
hhmm_end_qvp = '20:00'
year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start_qvp])
date_end = '-'.join([year, mon, day, hhmm_end_qvp])
top_height = 8
mode = 'vol'
elevation_deg = 12
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])

# # ------------------------------------ #
# # CFTDs                                #
# # ------------------------------------ #
# hhmm_start_cftds = '00:00'
# hhmm_end_cftds = '23:59'
# # ------------------------------------ #
# # full: ------------------------------ #
# elevation_degs = [8,12,17]
# locations = list(rad_dict().keys())
# dates = ['20210714', '20210713']
# data_max = 125000
# data_max = None
# testing = False
# # testing: --------------------------- #
# # elevation_degs = [12,]
# # locations = ['ESS']  # TODO: remove
# # dates = ['20210714']  # TODO: remove
# # data_max = 2000  # TODO: remove
# # testing = True
# # ------------------------------------ #
#
# # CFADs ? ---------------------------- #
# vert_temp = False
# height_min = 0  # in km
# height_max = 10  # in km
# bins_height = 20
# # or CFTDs ! ------------------------- #
# vert_temp = True
# temp_min = -20
# # temp_max = 16
# # bins_temp = 18
# temp_max = 0
# bins_temp = 10
# # ------------------------------------ #

# ------------------------------------ #
# QVPs & CFTDs                         #
# ------------------------------------ #
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
filter_entr = False
filter_entr_at = 0
filter_moms = False
folder_plot = header.folder_plot + 'Paper_IV/'

# ------------------------------------ #
# MODELS                               #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00010000.2')
spin_up_mms.append('120')
short_names.append('R0E1')
colors.append('red')
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
# spin_up_mms.append('120')
# short_names.append('R0E2')
# colors.append('orange')
# ------------------------------------ #
# SYN data row 3                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R0E3')
# colors.append('green')
# # ------------------------------------ #
# # SYN data row 4                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.1/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R1E3')
# colors.append('magenta')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R2E3')
colors.append('cyan')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# QVPs                                                                        #
# --------------------------------------------------------------------------- #
# for hm in ['g','h','i', 'r', 's','c']:
# # for hm in ['g','h',]:
#     # ------------------------------------ #
#     # QVPs plot parameters                 #
#     # ------------------------------------ #

mod_names = ''
letters_i = 0
n_rows = len(da_runs)
n_cols = 3
current_row = -1
current_col = -1
scale_font = 1.
scale_numbers = 1.
fig = plt.figure(figsize=(n_cols * 2.7, n_rows * 2.7), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.002,wspace=0.002)
axs = gs.subplots(sharex=True, sharey=True)
# ------------------------------------ #

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

    #TODO
    # filter for low ice:
    # syn_nc = syn_nc.where(syn_nc['qg'] + syn_nc['qh'] + syn_nc['qi'] + syn_nc['qs'] > 1E-6)
    #TODO

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
    hm='s'
    qvp_1_syn = syn_nc['D0_' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_1_syn=qvp_1_syn.where(syn_nc['vol_q' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time'),np.nan)
    qvp_1_syn=qvp_1_syn.where(syn_nc['vol_qn' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')>0,np.nan)

    hm='g'
    qvp_2_syn = syn_nc['vol_q' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')*1000
    qvp_2_syn=qvp_2_syn.where(syn_nc['D0_' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')>0,np.nan)
    qvp_2_syn=qvp_2_syn.where(syn_nc['vol_qn' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')>0,np.nan)

    hm='i'
    qvp_3_syn = np.log10(syn_nc['vol_qn' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')/1000)
    qvp_3_syn=qvp_3_syn.where(syn_nc['D0_' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')>0,np.nan)
    qvp_3_syn=qvp_3_syn.where(syn_nc['vol_q' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time'),np.nan)

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
    hm='s'
    plot_qvp_of_polarimetric_variable(
        mom=qvp_1_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_d0_ice,
        levels=header.levels_d0_ice,
        mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 10),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$D_{0,\,hm}\,[mm]$'.replace('hm', hm),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        panel= letters[letters_i] +') ' + short_name + ' (snow)',
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    hm='g'
    plot_qvp_of_polarimetric_variable(
        mom=qvp_2_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_iwc,
        levels=header.levels_iwc,
        mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 10),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$vol_{q,\,hm}\,[g\,m^{-3}]$'.replace('hm', hm),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + short_name + ' (graupel)',
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    hm='i'
    plot_qvp_of_polarimetric_variable(
        mom=qvp_3_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_nt_iwc,
        levels=header.levels_nt_iwc,
        mom_cs=qvp_zh_syn,
        levels_cs=np.arange(-50, 60, 10),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$vol_{qn,\,hm}\,[log_{10}(L^{-1})]$'.replace('hm', hm),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + short_name + ' (ice)',
    )
    letters_i=letters_i+1
    syn_nc.close()

# --------------------------------------------------------------------------- #
# QVPs SAVE                                                                   #
# --------------------------------------------------------------------------- #
# hh_at=[0,2,4,6,8,10,12,14,16,18,20,22]
hh_at=np.arange(int(hhmm_start_qvp[:2])+1, int(hhmm_end_qvp[:2])+1, 2)
hh_25=np.linspace(np.round(axs[-1,-1].get_xticks()[0]),
                   np.round(axs[-1,-1].get_xticks()[0])+1,25,endpoint=True)
str_hh_at=[str(z).zfill(2) for z in hh_at]
axs[-1,-1].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-1].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-2].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-2].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-3].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-3].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot +
    '/QVPs_' + str(n_rows) + 'x3_isg_QVP_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location + '_' +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/QVPs_' + str(n_rows) + 'x3_isg_QVP_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location + '_' +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.pdf', format='pdf', transparent=True)
plt.close()

