#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 13.02.25                                                 #
# plot_syn_RADAR_QVP_of_4_polmoms.py                                          #
#                                                                             #
# Run the QVP functions in PLOT_SYN_RADAR.py for generating specific QVP plot.#
# --------------------------------------------------------------------------- #

import sys
for entry in sys.path.copy():
    if '/RADAR_toolbox/' in entry:
        entry_folders = entry.split('/')
        index_mother = entry_folders.index('RADAR_toolbox') + 1
        sys.path.extend(['/'.join(entry_folders[:index_mother])])

import HEADER_RADAR_toolbox as header
import ColorBlindFriendlyRadarColorMaps as radar_colors
import os
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PROCESS_RADAR import d0_bringi
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
top_height = 6
mode = 'vol'
elevation_deg = 12
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])

# ------------------------------------ #
# CFTDs                                #
# ------------------------------------ #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
# ------------------------------------ #
# full: ------------------------------ #
elevation_degs = [8,12,17]
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
data_max = 125000
data_max = None
testing = False
# testing: --------------------------- #
# elevation_degs = [12,]
# locations = ['ESS']  # TODO: remove
# dates = ['20210714']  # TODO: remove
# data_max = 2000  # TODO: remove
# testing = True
# ------------------------------------ #

# CFADs ? ---------------------------- #
vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ! ------------------------- #
vert_temp = True
temp_min = 0.01
temp_max = 16.01
# bins_temp = 18
# temp_max = 0
bins_temp = 8
# ------------------------------------ #

# ------------------------------------ #
# QVPs & CFTDs                         #
# ------------------------------------ #
letters=('abcdefghijklmnopqrstuvwxyz'
         '\u03B1\u03B2\u03B3\u03B4\u03B5\u03B6\u03B7\u03B8\u03B9'
         '\u03BA\u03BB\u03BC\u03BD\u03BE\u03BF'
         '\u03C1\u03C2\u03C3\u03C4\u03C5\u03C6\u03C7\u03C8\u03C9')
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
colors2 = []
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_20010000.2')
# spin_up_mms.append('120')
# short_names.append('I1E1')
# colors.append('cyan')
# colors2.append('blue')
# ------------------------------------ #
# # SYN data row 2                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_20410000.2')
# spin_up_mms.append('120')
# short_names.append('I1E2')
# colors.append('green')
# colors2.append('darkgreen')
# # ------------------------------------ #
# # SYN data row 3                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_20510000.2')
# spin_up_mms.append('120')
# short_names.append('I1E3')
# # colors.append('yellow')
# colors2.append('orange')
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_20510000.2')
# spin_up_mms.append('120')
# short_names.append('I2E3')
# colorsX = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 5))
# colorsX[3] = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 7))[-3]
# colors.append(colorsX[np.array([ -2])])
# colors2.append('orange')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_20510840.2')
spin_up_mms.append('120')
short_names.append('I2E4')
colorsX = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 5))
colorsX[3] = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 7))[-3]
colors.append(colorsX[np.array([ -1])])
colorsX2 = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 15))
# colors2.append(colorsX2[np.array([ -1])])
colors2.append('magenta')
# ------------------------------------ #
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_20510840.2qnx')
# spin_up_mms.append('120')
# short_names.append('I2E4')
# colors.append('red')
# colors = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 5))
# colors[3] = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 7))[-3]
# colors = colors[np.array([-2, -1])]

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# QVPs                                                                        #
# --------------------------------------------------------------------------- #

# ------------------------------------ #
# QVPs plot parameters                 #
# ------------------------------------ #
mod_names = ''
letters_i = 0
n_rows = 2*len(da_runs)+1
n_cols = 3
current_row = -1
current_col = -1
scale_font = 1.
scale_numbers = 1.
factor=.75
fig = plt.figure(figsize=(factor*n_cols * 2.7, factor*n_rows * 2.8), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.0,wspace=0.05)
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

print('Nt_rain_qvp')
mom = (-2.37 + 0.1 * qvp_zh_obs -
                       2.89 * qvp_zdr_obs +
                       1.28 * qvp_zdr_obs ** 2 -
                       0.213 * qvp_zdr_obs ** 3)

mom=xr.where(qvp_temp_obs<0,np.nan,mom)
qvp_3_obs = mom

print('LWC_qvp')
mom = 10**(0.058 * qvp_zh_obs - 0.118 * qvp_zdr_obs - 2.36)
mom=xr.where(mom<0,np.nan,mom)
mom=xr.where(qvp_kdp_obs<0,np.nan,mom)
mom=xr.where(qvp_temp_obs<0,np.nan,mom)
qvp_2_obs = mom

# print('D0_bringi')
# mom = d0_bringi(qvp_zdr_obs)['d0']
print('Dm_bringi')
mom = d0_bringi(qvp_zdr_obs)['dm']
mom=xr.where(qvp_temp_obs<0,np.nan,mom)
qvp_1_obs = mom

extend='max'
extend='both'

# ------------------------------------ #
current_row = 0
# ------------------------------------ #
current_col = 0
plot_qvp_of_polarimetric_variable(
    mom=qvp_1_obs,
    cmap=radar_colors.cmap_radar_dm,
    norm=radar_colors.norm_dm_r,
    levels=radar_colors.levels_dm_r,
    mom_cs=qvp_zdr_obs,
    levels_cs=np.array([1]),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$D_{m,\,rain}\,[mm]$',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    extend=extend,
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
current_col = current_col + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_2_obs,
    cmap=radar_colors.cmap_radar_cont,
    norm=radar_colors.norm_cont,
    levels=radar_colors.levels_cont,
    mom_cs=qvp_zdr_obs,
    levels_cs=np.array([1]),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$RWC\,[g\,m^{-3}]$',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    ylabel='',
    extend=extend,
    panel=letters[letters_i] + ') obs',
)
letters_i=letters_i+1
# ------------------------------------ #
current_col = current_col + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_3_obs,
    cmap=radar_colors.cmap_radar_nt,
    norm=radar_colors.norm_nt,
    levels=radar_colors.levels_nt,
    mom_cs=qvp_zdr_obs,
    levels_cs=np.array([1]),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title='$N_{t,\,totice}\,[log_{10}(L^{-1})]$',
    title='',
    ax=axs[current_row, current_col],
    mom_height_unit='km',
    scale_font=scale_font,
    scale_numbers=scale_numbers,
    top_height=top_height,
    xlabel=None,
    ylabel='',
    extend=extend,
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
    if icon_emvorado_run[-3:] == 'qnx':
        path_mod = '/'.join(
            [header.dir_data_qvp + date, da_run, icon_emvorado_run[:-3],
             str(spin_up_mm) + 'min_spinup', 'QVPqnx_' +
             str(elevation_deg) + '_Syn_' + location + '_' +
             date + '0000_' + date + '2355.nc'])

    else:
        path_mod = '/'.join(
            [header.dir_data_qvp + date, da_run, icon_emvorado_run,
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
    qvp_kdp_syn = syn_nc['kdpsim'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_1_syn = syn_nc['D0_r'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')
    qvp_2_syn = syn_nc['vol_qr'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')*1000
    qvp_3_syn = np.log10(syn_nc['vol_qnr'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')/1000)

    qvp_3_syn=qvp_3_syn.where(qvp_2_syn>0,np.nan)
    qvp_3_syn=qvp_3_syn.where(qvp_1_syn>0,np.nan)

    qvp_2_syn=qvp_2_syn.where(qvp_3_syn,np.nan)
    qvp_2_syn=qvp_2_syn.where(qvp_1_syn>0,np.nan)

    qvp_1_syn=qvp_1_syn.where(qvp_3_syn,np.nan)
    qvp_1_syn=qvp_1_syn.where(qvp_2_syn>0,np.nan)

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
        mom=qvp_1_syn,
        cmap=radar_colors.cmap_radar_dm,
        norm=radar_colors.norm_dm_r,
        levels=radar_colors.levels_dm_r,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$D_{m,\,r}\,[mm]$',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        extend=extend,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_2_syn,
        cmap=radar_colors.cmap_radar_cont,
        norm=radar_colors.norm_cont,
        levels=radar_colors.levels_cont,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$RWC\,[g\,m^{-3}]$',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        extend=extend,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    current_col = current_col +1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_3_syn,
        cmap=radar_colors.cmap_radar_nt,
        norm=radar_colors.norm_nt,
        levels=radar_colors.levels_nt,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$N_{t,\,r}\,[log_{10}(L^{-1})]$',
        ax=axs[current_row, current_col],
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        xlabel=xlabel,
        ylabel='',
        extend=extend,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1

    # ----------------------------------------------------------------------- #
    # Syn Retrievals                                                          #
    # ----------------------------------------------------------------------- #

    print('Nt_rain_qvp')
    mom = (-2.37 + 0.1 * qvp_zh_syn -
           2.89 * qvp_zdr_syn +
           1.28 * qvp_zdr_syn ** 2 -
           0.213 * qvp_zdr_syn ** 3)

    mom = xr.where(qvp_temp_syn < 0, np.nan, mom)
    qvp_3_syn = mom

    print('LWC_qvp')
    mom = 10 ** (0.058 * qvp_zh_syn - 0.118 * qvp_zdr_syn - 2.36)
    mom = xr.where(mom < 0, np.nan, mom)
    mom = xr.where(qvp_kdp_syn < 0, np.nan, mom)
    mom = xr.where(qvp_temp_syn < 0, np.nan, mom)
    qvp_2_syn = mom

    print('Dm_bringi')
    mom = d0_bringi(qvp_zdr_syn)['dm']
    mom = xr.where(qvp_temp_syn < 0, np.nan, mom)
    qvp_1_syn = mom

    extend = 'max'
    extend = 'both'

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
        mom=qvp_1_syn,
        cmap=radar_colors.cmap_radar_dm,
        norm=radar_colors.norm_dm_r,
        levels=radar_colors.levels_dm_r,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$D_{m,\,rain}\,[mm]$',
        title='',
        ax=axs[current_row, current_col],
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        xlabel=xlabel,
        extend=extend,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #
    current_col = current_col + 1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_2_syn,
        cmap=radar_colors.cmap_radar_cont,
        norm=radar_colors.norm_cont,
        levels=radar_colors.levels_cont,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$RWC\,[g\,m^{-3}]$',
        title='',
        ax=axs[current_row, current_col],
        mom_height_unit='km',
        add_colorbar=add_colorbar,
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        xlabel=xlabel,
        ylabel='',
        extend=extend,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #
    current_col = current_col + 1
    plot_qvp_of_polarimetric_variable(
        mom=qvp_3_syn,
        cmap=radar_colors.cmap_radar_nt,
        norm=radar_colors.norm_nt,
        levels=radar_colors.levels_nt,
        mom_cs=qvp_zdr_syn,
        levels_cs=np.array([1]),
        mom_cf=qvp_temp_syn,
        levels_cf=np.arange(-50, 60, 5),
        cbar_title='$N_{t,\,totice}\,[log_{10}(L^{-1})]$',
        title='',
        ax=axs[current_row, current_col],
        add_colorbar=add_colorbar,
        mom_height_unit='km',
        scale_font=scale_font,
        scale_numbers=scale_numbers,
        top_height=top_height,
        xlabel=xlabel,
        ylabel='',
        extend=extend,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #

    syn_nc.close()

# --------------------------------------------------------------------------- #
# QVPs SAVE                                                                   #
# --------------------------------------------------------------------------- #
hh_at = np.arange(int(hhmm_start_qvp[:2]) + 2, int(hhmm_end_qvp[:2]) + 2, 2)
hh_25 = np.linspace(np.round(axs[-1, -1].get_xticks()[0]),
                    np.round(axs[-1, -1].get_xticks()[0]) + 1, 25,
                    endpoint=True)
str_hh_at = [str(z).zfill(2) for z in hh_at]
axs[-1, -1].set_xticks(hh_25[hh_at], str_hh_at, rotation=60)
axs[-1, -1].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1, -2].set_xticks(hh_25[hh_at], str_hh_at, rotation=60)
axs[-1, -2].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
axs[-1, -3].set_xticks(hh_25[hh_at], str_hh_at, rotation=60)
axs[-1, -3].set_xlabel('UTC [hh]', fontsize=12 * scale_font)
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot +
    '/QVPs_Fig_9_syn_ret_' + str(n_rows) + 'x3rainretrievals_QVP_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location + '_' +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/QVPs_Fig_9_syn_ret_' + str(n_rows) + 'x3rainretrievals_QVP_' +
    str(elevation_deg) + '°_' +
    date + '_' + hhmm_start_qvp + '-' + hhmm_end_qvp + '_' +
    location +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    '.pdf', format='pdf', transparent=True)
plt.close()


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# CFTDS plot parameters                #
# ------------------------------------ #
mod_names = ''
letters_i=0
n_rows = 2*len(da_runs) + 1 + 1  # add one once more for mean of all
n_cols = 3
factor=0.7
fig = plt.figure(figsize=(factor*n_cols * 2.8, factor*n_rows * 2.4),)
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.03,wspace=0.03)
axs = gs.subplots()
# ------------------------------------ #
ax_mean1=axs[-1,0]
ax_mean1.set_ylabel('temperature [°C]')
ax_mean1.set_xlabel('$D_{m,\,r}\,[mm]$')
ax_mean1.set_xlim([mom_plot_dict('Dm_r')['mom_min'],
                   mom_plot_dict('Dm_r')['mom_max']])
ax_mean1.set_ylim([temp_min,
                   temp_max])

ax_mean2 = axs[-1, 1]
ax_mean2.set_ylabel('temperature [°C]')
ax_mean2.set_xlabel('$RWC\,[g\,m^{-3}]$')
ax_mean2.set_xlim([mom_plot_dict('LWC')['mom_min'],
                   mom_plot_dict('LWC')['mom_max']])
ax_mean2.set_ylim([temp_min,
                   temp_max])

ax_mean3 = axs[-1, 2]
ax_mean3.set_ylabel('temperature [°C]')
ax_mean3.set_xlabel('$N_{t,\,r}\,[log_{10}(L^{-1})]$')
ax_mean3.set_xlim([mom_plot_dict('Nt_r')['mom_min'],
                   mom_plot_dict('Nt_r')['mom_max']])
ax_mean3.set_ylim([temp_min,
                   temp_max])
# --------------------------------------------------------------------------- #
# CFTDs OBS row 1                                                             #
# --------------------------------------------------------------------------- #
color='black'
# ------------------------------------ #
current_row = 0
current_col = 0
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_deg=elevation_degs,
    da_icon_emvorado_run=None,
    moment='Dm_bringi',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    filter_entr=filter_entr,
    filter_entr_at=filter_entr_at,
    filter_moms=filter_moms,
    ax=ax,
    save=False,
    color=color,
    plot_legend=False,
    plot_data=False,
    data_max=data_max,
    panel=letters[letters_i] + ') obs',
)
# ------------------------------------ #
letters_i=letters_i+1
y_bins = np.linspace(temp_min,temp_max,bins_temp+1)
y_step=y_bins[1]-y_bins[0]
y_mid = np.linspace(temp_min+1,temp_max-1,bins_temp)
quant_prof_obs1 = np.zeros([3, len(y_mid)])
quant_prof_obs1[:] = np.nan
mean_prof = np.repeat(np.nan, len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    if x_layer.size>3:
        wq = DescrStatsW(data=x_layer)
        quant_prof_obs1[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

# ax_mean1.plot(quant_prof_obs1[0, ], y_mid, color=color, ls='dashed', alpha=0.8,
#          linewidth=1, label='_nolegend_')
ax_mean1.plot(quant_prof_obs1[1, ], y_mid, color=color, ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4,label='obs')
# ax_mean1.plot(quant_prof_obs1[2, ], y_mid, color=color, ls='dashed',alpha=0.8,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',alpha=0.8,
#          linewidth=2, label='_nolegend_')
# --------------------------------------------------------------------------- #
current_col = current_col + 1
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_deg=elevation_degs,
    da_icon_emvorado_run=None,
    moment='LWC_qvp',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    filter_entr=filter_entr,
    filter_entr_at=filter_entr_at,
    filter_moms=filter_moms,
    ax=ax,
    save=False,
    color=color,
    plot_legend=True,
    plot_data=False,
    data_max=data_max,
    panel=letters[letters_i] + ') obs',
)
# ------------------------------------ #
letters_i=letters_i+1
y_bins = np.linspace(temp_min,temp_max,bins_temp+1)
y_step=y_bins[1]-y_bins[0]
y_mid = np.linspace(temp_min+1,temp_max-1,bins_temp)
quant_prof_obs2 = np.zeros([3, len(y_mid)])
quant_prof_obs2[:] = np.nan
mean_prof = np.repeat(np.nan, len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    if x_layer.size>3:
        wq = DescrStatsW(data=x_layer)
        quant_prof_obs2[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

# ax_mean2.plot(quant_prof_obs2[0, ], y_mid, color=color, ls='dashed',alpha=0.8,
#          linewidth=1, label='_nolegend_')
ax_mean2.plot(quant_prof_obs2[1, ], y_mid, color=color, ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4, label='obs')
# ax_mean2.plot(quant_prof_obs2[2, ], y_mid, color=color, ls='dashed',alpha=0.8,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',alpha=0.8,
#          linewidth=2, label='_nolegend_')
# --------------------------------------------------------------------------- #
current_col = current_col + 1
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_deg=elevation_degs,
    da_icon_emvorado_run=None,
    moment='Nt_rain_qvp',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    filter_entr=filter_entr,
    filter_entr_at=filter_entr_at,
    filter_moms=filter_moms,
    ax=ax,
    save=False,
    color=color,
    plot_legend=False,
    plot_data=False,
    data_max=data_max,
    data_label=True,
    panel=letters[letters_i] + ') obs',
)
# ------------------------------------ #
letters_i=letters_i+1
y_bins = np.linspace(temp_min,temp_max,bins_temp+1)
y_step=y_bins[1]-y_bins[0]
y_mid = np.linspace(temp_min+1,temp_max-1,bins_temp)
quant_prof_obs3 = np.zeros([3, len(y_mid)])
quant_prof_obs3[:] = np.nan
mean_prof = np.repeat(np.nan, len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    if x_layer.size>3:
        wq = DescrStatsW(data=x_layer)
        quant_prof_obs3[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

# ax_mean3.plot(quant_prof_obs3[0, ], y_mid, color=color, ls='dashed', alpha=0.8,
#          linewidth=1, label='_nolegend_')
ax_mean3.plot(quant_prof_obs3[1, ], y_mid, color=color, ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4, label='obs')
# ax_mean3.plot(quant_prof_obs3[2, ], y_mid, color=color, ls='dashed', alpha=0.8,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
#          linewidth=2, label='_nolegend_')
# --------------------------------------------------------------------------- #
# CFTDs CBAND SYN row i                                                       #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm, color, color2,short_name in zip(
        da_runs, icon_emvorado_runs, spin_up_mms, colors, colors2, short_names):
    da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
    model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
                                icon_emvorado_run.split('/')[1][5:]])
    mod_names = '_'.join([mod_names, model_name_file])
    model_name = '-'.join([da_run[4:],
                           icon_emvorado_run.split('/')[0][5:],
                           icon_emvorado_run.split('/')[1][5:],
                           spin_up_mm + 'min'])
    # ------------------------------------ #
    current_row = current_row + 1
    current_col = 0
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='D0_r',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color,
        plot_legend=False,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean1.plot(quant_prof[1,], y_mid, color=color, ls='solid',# TODO: mean and median swapped
                  linewidth=2+(n_rows-current_row-3)/4, label=short_name)
    ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')
    # ----------------------------------------------------------------------- #
    current_col = current_col + 1
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='vol_qr',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color,
        plot_legend=True,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean2.plot(quant_prof[1,], y_mid, color=color, ls='solid',# TODO: mean and median swapped
                  linewidth=2+(n_rows-current_row-3)/4, label=short_name)
    ax_mean2.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')
    # ----------------------------------------------------------------------- #
    current_col = current_col + 1
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='vol_qnr',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color,
        plot_legend=False,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean3.plot(quant_prof[1,], y_mid, color=color, ls='solid',# TODO: mean and median swapped
                  linewidth=2+(n_rows-current_row-3)/4, label=short_name)
    ax_mean3.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')



    # ----------------------------------------------------------------------- #
    # syn Retrievals                                                          #
    # ----------------------------------------------------------------------- #

    current_row = current_row + 1
    current_col = 0
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x, y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='Dm_bringi',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color2,
        plot_legend=False,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean1.plot(quant_prof[0,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean1.plot(quant_prof[1,], y_mid, color=color2, ls='solid',
                  # TODO: mean and median swapped
                  linewidth=2 + (n_rows - current_row - 3) / 4,
                  label='R('+short_name+')')
    ax_mean1.plot(quant_prof[2,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean1.plot(mean_prof, y_mid, color=color2, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')
    # ----------------------------------------------------------------------- #
    current_col = current_col + 1
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x, y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='LWC_qvp',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color2,
        plot_legend=True,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean2.plot(quant_prof[0,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean2.plot(quant_prof[1,], y_mid, color=color2, ls='solid',
                  # TODO: mean and median swapped
                  linewidth=2 + (n_rows - current_row - 3) / 4,
                  label='R('+short_name+')')
    ax_mean2.plot(quant_prof[2,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean2.plot(mean_prof, y_mid, color=color2, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')
    # ----------------------------------------------------------------------- #
    current_col = current_col + 1
    print(current_row)
    print(current_col)
    ax = axs[current_row, current_col]
    x, y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_degs,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='Nt_rain_qvp',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color2,
        plot_legend=False,
        plot_data=False,
        data_max=data_max,
        panel= letters[letters_i] +') R(' + short_name + ')' ,
    )
    letters_i = letters_i + 1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    quant_prof[:] = np.nan
    mean_prof = np.repeat(np.nan, len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        if x_layer.size > 3:
            wq = DescrStatsW(data=x_layer)
            quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                             return_pandas=False)
            mean_prof[t_i] = wq.mean

    ax_mean3.plot(quant_prof[0,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean3.plot(quant_prof[1,], y_mid, color=color2, ls='solid',
                  # TODO: mean and median swapped
                  linewidth=2 + (n_rows - current_row - 3) / 4,
                  label='R('+short_name+')')
    ax_mean3.plot(quant_prof[2,], y_mid, color=color2, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean3.plot(mean_prof, y_mid, color=color2, ls='solid', alpha=0.8,
    #               linewidth=2,label='_nolegend_')





# --------------------------------------------------------------------------- #
# CFTDs SAVE                                                                  #
# --------------------------------------------------------------------------- #

current_row=1

ax_mean1.plot(quant_prof_obs1[0, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')
ax_mean1.plot(quant_prof_obs1[1, ], y_mid, color='black', ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4, label='_nolegend_', alpha=.5)
ax_mean1.plot(quant_prof_obs1[2, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')

ax_mean2.plot(quant_prof_obs2[0, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')
ax_mean2.plot(quant_prof_obs2[1, ], y_mid, color='black', ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4, label='_nolegend_', alpha=.5)
ax_mean2.plot(quant_prof_obs2[2, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')

ax_mean3.plot(quant_prof_obs3[0, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')
ax_mean3.plot(quant_prof_obs3[1, ], y_mid, color='black', ls='solid',# TODO: mean and median swapped
         linewidth=2+(n_rows-current_row-3)/4, label='', alpha=.5)
ax_mean3.plot(quant_prof_obs3[2, ], y_mid, color='black', ls='dashed',alpha=0.8,
         linewidth=1, label='_nolegend_')

ax_mean1.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean1.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean1.legend()
letters_i=letters_i +1

ax_mean2.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
ax_mean2.legend()
letters_i=letters_i +1

ax_mean3.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean3.legend()
letters_i=letters_i +1


for i_r in range(n_rows):
    for i_c in range(n_cols):
        axs[i_r,i_c].set_title('')
        if i_r<n_rows-1:
            axs[i_r,i_c].set_xlabel('')
            axs[i_r, i_c].set_xticklabels('')

        if i_c>0:
            axs[i_r,i_c].set_ylabel('')
            axs[i_r, i_c].set_yticklabels('')

        # if i_c==1:
        #     axs[i_r,i_c].set_xlim([mom_plot_dict('ZDR')['mom_min'], 4.6])

cmap = mpl.cm.terrain_r  #TODO
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
norm = mpl.colors.Normalize(vmin=0, vmax=16)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
             extend='max')
axs[-1,-1].set_xlim([mom_plot_dict('Nt_r')['mom_min'],
                     mom_plot_dict('Nt_r')['mom_min']+
                     (mom_plot_dict('Nt_r')['mom_max']-
                      mom_plot_dict('Nt_r')['mom_min'])*.8])


# axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].set_xticks(
#     axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].get_xticks(),
#     [str(int(i/1000))+'k' for i in
#          axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].get_xticks()
#      ],color='gray')

gs.tight_layout(fig, rect=[0, 0, 0.5, 1.0])
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

if len(dates) == 1:
    dates_str = dates[0] + '_'
elif dates == ['20210714', '20210713'] or dates == ['20210713', '20210714']:
    dates_str = '20210713+14_'
else:
    dates_str = 'All_t_'

if len(locations) == 1:
    locations_str = locations[0] + '_'
else:
    locations_str = 'all_radars_'

plt.savefig(
    folder_plot +
    '/CFTDs_Fig_10_syn_ret_' + str(n_rows)+ 'x3rainetrievals_QVP_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_Fig_10_syn_ret_' + str(n_rows) + 'x3rainretrievals_QVP_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
plt.close()

# --------------------------------------------------------------------------- #
# CFTDs END                                                                   #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
