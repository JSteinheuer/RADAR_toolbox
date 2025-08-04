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

# ------------------------------------ #
# CFTDs                                #
# ------------------------------------ #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
elevation_degs = [8,12,17]
# elevation_degs = [12]

# ------------------------------------ #
# full: ------------------------------ #
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
data_max = 31000
data_max = 125000
testing = False
# testing: --------------------------- #
# locations = ['ESS']  # TODO: remove
# dates = ['20210714']  # TODO: remove
# data_max = 2200  # TODO: remove
# testing = True
# elevation_degs = [8,]
# ------------------------------------ #

# CFADs ? ---------------------------- #
vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ! ------------------------- #
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
# ------------------------------------ #

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
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
spin_up_mms.append('120')
short_names.append('R0E2')
colors.append('orange')
# ------------------------------------ #
# SYN data row 3                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R0E3')
colors.append('green')
# # ------------------------------------ #
# # SYN data row 4                       #
# # ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.1/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R1E3')
colors.append('magenta')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R2E3')
# colors.append('blue')
colors.append('cyan')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# QVPs                                                                        #
# --------------------------------------------------------------------------- #

# ------------------------------------ #
# QVPs plot parameters                 #
# ------------------------------------ #
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
    # title='Z$_{H}$ (C-band Obs. at ' + location + ')',
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
    # title='Z$_{DR}$ (C-band Obs. at ' + location + ')',
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
    # title='K$_{DP}$ (C-band Obs. at ' + location + ')',
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
    # title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
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
        # title='Z$_{H}$ (%s)' % model_name,
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
        # title='Z$_{DR}$ (%s)' % model_name,
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
        # title='K$_{DP}$ (%s)' % model_name,
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
        # title='$\u03C1_{hv}$ (%s)' % model_name,
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
hh_at=[2,4,6,8,10,12,14,16,18]
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
plt.savefig(
    folder_plot +
    '/QVPs_' + str(n_rows) + 'x4polmoms_' +
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
n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
n_cols = 4
# fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
# fig = plt.figure(figsize=(n_cols * 3.3, n_rows * 2.5), layout='constrained')
# fig = plt.figure(figsize=(n_cols * 2.8, n_rows * 2.5), layout='constrained')
fig = plt.figure(figsize=(n_cols * 2.8, n_rows * 2.2), layout='constrained')
gs = fig.add_gridspec(n_rows, n_cols, hspace=0.03,wspace=0.03)
axs = gs.subplots()
# ------------------------------------ #
ax_mean1=axs[-1,0]
ax_mean1.set_ylabel('temperature [°C]')
ax_mean1.set_xlabel('$Z_{H}$ [dBZ]')
ax_mean1.set_xlim([mom_plot_dict('ZH')['mom_min'],
                   mom_plot_dict('ZH')['mom_max']])
ax_mean2 = axs[-1, 1]
ax_mean2.set_ylabel('temperature [°C]')
ax_mean2.set_xlabel('$Z_{DR}$ [dB]')
ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],
                   mom_plot_dict('ZDR')['mom_max']])
ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],3.3])
ax_mean3 = axs[-1, 2]
ax_mean3.set_ylabel('temperature [°C]')
ax_mean3.set_xlabel('$K_{DP}$ [°/km]')
ax_mean3.set_xlim([mom_plot_dict('KDP')['mom_min'],
                   mom_plot_dict('KDP')['mom_max']])
ax_mean4 = axs[-1, 3]
ax_mean4.set_ylabel('temperature [°C]')
ax_mean4.set_ylabel('temperature [°C]')
ax_mean4.set_xlabel('$\u03C1_{HV}$ [1]')
ax_mean4.set_xlim([mom_plot_dict('RHOHV')['mom_min'],
                   mom_plot_dict('RHOHV')['mom_max']])
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
    moment='ZH_AC',
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    wq = DescrStatsW(data=x_layer)
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean1.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
# ax_mean1.plot(quant_prof[1, ], y_mid, color=color, ls='dashdot',
#          linewidth=2, label='_nolegend_')
ax_mean1.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
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
    moment='ZDR_AC_OC',
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    wq = DescrStatsW(data=x_layer)
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean2.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
# ax_mean2.plot(quant_prof[1, ], y_mid, color=color, ls='dashdot',
#          linewidth=2, label='_nolegend_')
ax_mean2.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
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
    moment='KDP_NC',
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
    plot_data=True,
    data_max=data_max,
    data_label=True,
    panel=letters[letters_i] + ') obs',
)
# ------------------------------------ #
letters_i=letters_i+1
y_bins = np.linspace(temp_min,temp_max,bins_temp+1)
y_step=y_bins[1]-y_bins[0]
y_mid = np.linspace(temp_min+1,temp_max-1,bins_temp)
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    wq = DescrStatsW(data=x_layer)
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean3.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
# ax_mean3.plot(quant_prof[1, ], y_mid, color=color, ls='dashdot',
#          linewidth=2, label='_nolegend_')
ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
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
    moment='RHOHV_NC2P',
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    x_layer=x[(y>y_mid[t_i]-y_step/2) * (y<=y_mid[t_i]+y_step/2)]
    wq = DescrStatsW(data=x_layer)
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean4.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
# ax_mean4.plot(quant_prof[1, ], y_mid, color=color, ls='dashdot',
#          linewidth=2, label='_nolegend_')
ax_mean4.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
# --------------------------------------------------------------------------- #
# CFTDs CBAND SYN row i                                                       #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm, color, short_name in zip(
        da_runs, icon_emvorado_runs, spin_up_mms, colors, short_names):
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
        moment='zrsim',
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
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        wq = DescrStatsW(data=x_layer)
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
    # ax_mean1.plot(quant_prof[1,], y_mid, color=color, ls='dashdot',
    #               linewidth=2, label='_nolegend_')
    ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
    ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
                  linewidth=2,
                  label=short_name)
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
        moment='zdrsim',
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
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        wq = DescrStatsW(data=x_layer)
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.3,
                  linewidth=1, label='_nolegend_')
    # ax_mean2.plot(quant_prof[1,], y_mid, color=color, ls='dashdot',
    #               linewidth=2, label='_nolegend_')
    ax_mean2.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
                  linewidth=1, label='_nolegend_')
    ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
                  linewidth=2,
                  label=short_name)
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
        moment='kdpsim',
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
        plot_data=True,
        data_max=data_max,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    y_bins = np.linspace(temp_min, temp_max, bins_temp + 1)
    y_step = y_bins[1] - y_bins[0]
    y_mid = np.linspace(temp_min + 1, temp_max - 1, bins_temp)
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        wq = DescrStatsW(data=x_layer)
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
    # ax_mean3.plot(quant_prof[1,], y_mid, color=color, ls='dashdot',
    #               linewidth=2, label='_nolegend_')
    ax_mean3.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
                  linewidth=1, label='_nolegend_')
    ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
                  linewidth=2,
                  label=short_name)
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
        moment='rhvsim',
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
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        x_layer = x[
            (y > y_mid[t_i] - y_step / 2) * (y <= y_mid[t_i] + y_step / 2)]
        wq = DescrStatsW(data=x_layer)
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean4.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
    # ax_mean4.plot(quant_prof[1,], y_mid, color=color, ls='dashdot',
    #               linewidth=2, label='_nolegend_')
    ax_mean4.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
    ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
                  linewidth=2,
                  label=short_name)

# --------------------------------------------------------------------------- #
# CFTDs SAVE                                                                  #
# --------------------------------------------------------------------------- #

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

ax_mean4.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean4.legend()
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

        if i_c==1:
            axs[i_r,i_c].set_xlim([mom_plot_dict('ZDR')['mom_min'], 4.6])

# cmap = mpl.cm.YlGnBu  #TODO
cmap = mpl.cm.terrain_r  #TODO
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
norm = mpl.colors.Normalize(vmin=0, vmax=16)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
             extend='max')
axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])


axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].set_xticks(
    axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].get_xticks(),
    [str(int(i/1000))+'k' for i in
         axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].get_xticks()
     ],color='gray')
    # ['0','25k','50k', '75k', '100k'],color='gray')

# ax2 = axs[0,2].twiny()
# ax2.set_xticks(axs[0,2].get_shared_y_axes().get_siblings(axs[0,2])[0].get_xticks(),
#                ['0','25k','50k', '75k', '100k'],color='gray')

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
    '/CFTDs_' + str(n_rows)+ 'x4polmoms_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=False, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_' + str(n_rows) + 'x4polmoms_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at) + '_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.png', format='png', transparent=False, dpi=300, bbox_inches='tight')
plt.close()

# --------------------------------------------------------------------------- #
# CFTDs END                                                                   #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
