#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 29.10.25                                                 #
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

# ------------------------------------ #
# MODEL (only 1 SYN here!)             #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R2E3')
colors.append('cyan')
da_run, icon_emvorado_run, spin_up_mm, short_name, color = (
    da_runs[0], icon_emvorado_runs[0], spin_up_mms[0], short_names[0], colors[0])

# ------------------------------------ #
# QVPs plot parameters                 #
# ------------------------------------ #
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2\u03B3\u03B4'
filter_entr = False
filter_entr_at = 0
filter_moms = False
folder_plot = header.folder_plot + 'Paper_IV/'
mod_names = ''
letters_i = 0
n_rows = 6
n_cols = 5
current_row = -1
current_col = -1
scale_font = 1.
scale_numbers = 1.
fig = plt.figure(figsize=(n_cols * 2.7, n_rows * 2.7), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.002,wspace=0.002)
axs = gs.subplots(sharex=True, sharey=True)

# ------------------------------------ #
# QVP plots                            #
# ------------------------------------ #
for hm in ['i','s','h', 'g', 'c','r']:
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
    qvp_zdr_syn = syn_nc['zdrsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_zh_syn = syn_nc['zrsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_1_syn = syn_nc['D0_' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_2_syn = syn_nc['vol_q' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')*1000
    qvp_3_syn = np.log10(syn_nc['vol_qn' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')/1000)
    qvp_4_syn = syn_nc['q' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')*1000
    qvp_5_syn = np.log10(syn_nc['qn' + hm].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')/1000)
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
    current_col = 4
    letters_i = current_col + current_row*(n_cols)
    plot_qvp_of_polarimetric_variable(
        mom=qvp_1_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_d0_ice,
        levels=header.levels_d0_ice,
        # mom_cs=qvp_zh_syn,
        # levels_cs=np.arange(-50, 60, 10),
        mom_cs=qvp_zdr_syn,
        levels_cs=np.arange(0, 6, 1),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$D_{0,\,hm}\,[mm]$'.replace('hm', 'x'),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + '$D_{0,\,hm}$'.replace('hm', hm),
    )
    # ------------------------------------ #
    current_col = 1
    letters_i = current_col + current_row*(n_cols)
    plot_qvp_of_polarimetric_variable(
        mom=qvp_2_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_iwc,
        levels=header.levels_iwc,
        # mom_cs=qvp_zh_syn,
        # levels_cs=np.arange(-50, 60, 10),
        mom_cs=qvp_zdr_syn,
        levels_cs=np.arange(0, 6, 1),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$vol_{q,\,hm}\,[g\,m^{-3}]$'.replace('hm', 'x'),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + '$vol_{q,\,hm}]$'.replace('hm', hm),
    )
    # ------------------------------------ #
    current_col = 3
    letters_i = current_col + current_row*(n_cols)
    plot_qvp_of_polarimetric_variable(
        mom=qvp_3_syn,
        cmap=header.cmap_radar,
        norm=header.norm_nt_iwc,
        levels=header.levels_nt_iwc,
        # mom_cs=qvp_zh_syn,
        # levels_cs=np.arange(-50, 60, 10),
        mom_cs=qvp_zdr_syn,
        levels_cs=np.arange(0, 6, 1),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$vol_{qn,\,hm}\,[log_{10}(L^{-1})]$'.replace('hm', 'x'),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + '$vol_{qn,\,hm}\,$'.replace('hm', hm),
    )
    # ------------------------------------ #
    current_col = 0
    letters_i = current_col + current_row*(n_cols)
    plot_qvp_of_polarimetric_variable(
        mom=qvp_4_syn,
        cmap=header.cmap_radar_white,
        norm=header.norm_iwc,
        levels=header.levels_iwc,
        # mom_cs=qvp_zh_syn,
        # levels_cs=np.arange(-50, 60, 10),
        mom_cs=qvp_zdr_syn,
        levels_cs=np.arange(0, 6, 1),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$mass_{q,\,hm}\,[g\,kg^{-1}]$'.replace('hm', 'x'),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        panel= letters[letters_i] +') ' + '$mass_{q,\,hm}$'.replace('hm', hm),
    )
    # ------------------------------------ #
    current_col = 2
    letters_i = current_col + current_row*(n_cols)
    plot_qvp_of_polarimetric_variable(
        mom=qvp_5_syn,
        cmap=header.cmap_radar,
        norm=header.norm_nt_iwc,
        levels=header.levels_nt_iwc,
        # mom_cs=qvp_zh_syn,
        # levels_cs=np.arange(-50, 60, 10),
        mom_cs=qvp_zdr_syn,
        levels_cs=np.arange(0, 6, 1),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$mass_{qn,\,hm}\,[log_{10}(g^{-1})]$'.replace('hm', 'x'),
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        panel= letters[letters_i] +') ' + '$mass_{qn,\,hm}$'.replace('hm', hm),
    )
    syn_nc.close()

# --------------------------------------------------------------------------- #
# QVPs SAVE                                                                   #
# --------------------------------------------------------------------------- #
# hh_at=[0,2,4,6,8,10,12,14,16,18,20,22]
hh_at=np.arange(int(hhmm_start_qvp[:2])+2, int(hhmm_end_qvp[:2])+2, 2)
hh_25=np.linspace(np.round(axs[-1,-1].get_xticks()[0]),
                   np.round(axs[-1,-1].get_xticks()[0])+1,25,endpoint=True)
str_hh_at=[str(z).zfill(2) for z in hh_at]
axs[-1,-1].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-1].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-2].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-2].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-3].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-3].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-4].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-4].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1,-5].set_xticks(hh_25[hh_at],str_hh_at)
axs[-1,-5].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot +
    '/QVPs_' + str(n_rows) + 'x' + str(n_cols) +'_ZDR_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location + '_' +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/QVPs_' + str(n_rows) + 'x' + str(n_cols) +'_ZDR_' +
#     str(elevation_deg) + '°_' +
#     date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
#     location +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at) + '_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     '.pdf', format='pdf', transparent=True)
plt.close()
