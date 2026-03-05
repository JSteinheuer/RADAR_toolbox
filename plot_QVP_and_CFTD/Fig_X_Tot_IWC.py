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
        entry_folders=entry.split('/')
        index_mother=entry_folders.index('RADAR_toolbox')+1
        sys.path.extend(['/'.join(entry_folders[:index_mother])])

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import ColorBlindFriendlyRadarColorMaps as radar_colors
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
# CFTDs                                #
# ------------------------------------ #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
# full: ------------------------------ #
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
data_max = 125000
elevation_degs = [8,12,17]
data_max = None
testing = False
# testing: --------------------------- #
locations = ['ESS']  # TODO: remove
dates = ['20210714']  # TODO: remove
elevation_degs = [12,]
testing = True
# ------------------------------------ #

# CFADs ? ---------------------------- #
vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ! ------------------------- #
vert_temp = True
temp_min = -20.01
# temp_max = 16
# bins_temp = 18
temp_max = -.01
bins_temp = 10
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
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_20010000.2')
spin_up_mms.append('120')
short_names.append('I1E1')
# colors.append('cyan')
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_20510000.2')
spin_up_mms.append('120')
short_names.append('I2E3')
# colors.append('orange')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_20510840.2qnx')
spin_up_mms.append('120')
short_names.append('I2E4')
# colors.append('red')
# ------------------------------------ #
colors = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 5))
colors[3] = mpl.colormaps._cmaps['HomeyerRainbow'](np.linspace(0, 1, 7))[-3]
colors = colors[np.array([0,-2,-1])]

# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# CFTDS plot parameters                #
# ------------------------------------ #
mod_names = ''
letters_i=0
n_rows = 1
n_cols = 1
factor=0.7
fig = plt.figure(figsize=(factor*n_cols * 2.8, factor*n_rows * 2.8),)
gs = fig.add_gridspec(n_rows, n_cols)
axs = gs.subplots()
# ------------------------------------ #
ax_mean=axs
ax_mean.set_ylabel('temperature [°C]')
ax_mean.set_xlabel('$total IWC\,[g\,m^{-3}]$')
# ax_mean.set_xlim([mom_plot_dict('IWC')['mom_min'],
#                    mom_plot_dict('IWC')['mom_max']])
ax_mean.set_ylim([temp_min, temp_max])

# --------------------------------------------------------------------------- #
# CFTDs OBS row 1                                                             #
# --------------------------------------------------------------------------- #
color='black'
# ------------------------------------ #
x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_deg=elevation_degs,
    da_icon_emvorado_run=None,
    moment='IWC_qvp',
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
)
# ------------------------------------ #
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
    moment='IWC_qvp',
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
    moment='Nt_totice_qvp',
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
ax.set_xticks([-1,0,1], [-1,0,1])
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
# ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.8,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
#          linewidth=2, label='_nolegend_')
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
        moment='D0_totice',
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

    ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    ax_mean1.plot(quant_prof[1,], y_mid, color=color, ls='solid',# TODO: mean and median swapped
                  linewidth=2+(n_rows-current_row-3)/4, label=short_name)
    ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.8,
                  linewidth=1, label='_nolegend_')
    # ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid', alpha=0.8,
    #               linewidth=2, label='_nolegend_')
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
        moment='vol_qtotice',
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
        moment='vol_qntotice',
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
    ax.set_xticks([-1, 0, 1], [-1, 0, 1])
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
axs[-1,-1].set_xlim([mom_plot_dict('Nt_totice')['mom_min'],
                     mom_plot_dict('Nt_totice')['mom_min']+
                     (mom_plot_dict('Nt_totice')['mom_max']-
                      mom_plot_dict('Nt_totice')['mom_min'])*.8])


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
    '/CFTDs_' + str(n_rows)+ 'x3iceretrievals_QVP_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_' + str(n_rows) + 'x3iceretrievals_QVP_' +
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
