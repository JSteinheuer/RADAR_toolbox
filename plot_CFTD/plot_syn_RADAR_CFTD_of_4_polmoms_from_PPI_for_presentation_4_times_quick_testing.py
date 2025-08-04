#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.05.25                                                 #
# --------------------------------------------------------------------------- #

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable
from PLOT_SYN_RADAR import (plot_CFAD_or_CFTD_from_PPI_with_list,
                            mom_plot_dict,
                            plot_CFAD_or_CFTD_from_QVP_with_list,
                            plot_CFAD_or_CFTD_from_PPI_with_list_quick)
from SET_SYN_RADAR import rad_dict
from statsmodels.stats.weightstats import DescrStatsW

# --------------------------------------------------------------------------- #
#                              PLOT 0 entr                                    #
# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
# ------------------------------------ #
locations = ['ESS']
dates = ['20210714']
data_max = None
testing = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
elevation_degs = [12]
elevation_deg = 12
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
filter_moms = False
# with entropy filter: --------------- #
filter_entr = True
filter_entr_at = 0.7
# ------------------------------------ #
# MODELS                               #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
spin_up_mms.append('120')
short_names.append('R0E2')
colors.append('orange')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R2E3')
colors.append('blue')

# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# CFTDS plot parameters                #
# ------------------------------------ #
folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
letters_i=0
n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
n_cols = 4
fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
letters_i=letters_i +1

ax_mean2.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
ax_mean2.legend()
letters_i=letters_i +1

ax_mean3.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
letters_i=letters_i +1

ax_mean4.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
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

cmap = mpl.cm.YlGnBu
norm = mpl.colors.Normalize(vmin=0, vmax=15)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
             extend='max')
axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])

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
    '/CFTDs_' + str(n_rows)+ 'x4polmoms_QVP_' +
    str(elevation_deg) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_' + str(n_rows) + 'x4polmoms_QVP_' +
    str(elevation_deg) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
plt.close()

# --------------------------------------------------------------------------- #
#                              PLOT 0 no entr                                 #
# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
# ------------------------------------ #
locations = ['ESS']
dates = ['20210714']
data_max = None
testing = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
elevation_degs = [12]
elevation_deg = 12
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
filter_moms = False
# with entropy filter: --------------- #
filter_entr = False
filter_entr_at = 0
# ------------------------------------ #
# MODELS                               #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
spin_up_mms.append('120')
short_names.append('R0E2')
colors.append('orange')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R2E3')
colors.append('blue')

# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# CFTDS plot parameters                #
# ------------------------------------ #
folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
letters_i=0
n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
n_cols = 4
fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
    elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
        elevation_deg=elevation_deg,
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
letters_i=letters_i +1

ax_mean2.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
ax_mean2.legend()
letters_i=letters_i +1

ax_mean3.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
letters_i=letters_i +1

ax_mean4.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
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

cmap = mpl.cm.YlGnBu
norm = mpl.colors.Normalize(vmin=0, vmax=15)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
             extend='max')
axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])

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
    '/CFTDs_' + str(n_rows)+ 'x4polmoms_QVP_' +
    str(elevation_deg) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_' + str(n_rows) + 'x4polmoms_QVP_' +
    str(elevation_deg) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
plt.close()


# # --------------------------------------------------------------------------- #
# #                              PLOT A                                         #
# # --------------------------------------------------------------------------- #
# # INITIALIZATION                                                              #
# # --------------------------------------------------------------------------- #
# hhmm_start_cftds = '00:00'
# hhmm_end_cftds = '23:59'
# # ------------------------------------ #
# locations = ['ESS']
# dates = ['20210714']
# data_max = None
# testing = False
# height_min = 0  # in km
# height_max = 10  # in km
# bins_height = 20
# vert_temp = True
# temp_min = -20
# temp_max = 16
# bins_temp = 18
# elevation_degs = [12]
# letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
# filter_moms = False
# # ------------------------------------ #
# # MODELS                               #
# # ------------------------------------ #
# da_runs = []
# icon_emvorado_runs = []
# spin_up_mms = []
# short_names = []
# colors = []
# # ------------------------------------ #
# # SYN data row 2                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
# spin_up_mms.append('120')
# short_names.append('R0E2')
# colors.append('orange')
# # ------------------------------------ #
# # SYN data row 5                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R2E3')
# colors.append('blue')
#
# # --------------------------------------------------------------------------- #
# # CFTDs                                                                       #
# # --------------------------------------------------------------------------- #
# # ------------------------------------ #
# # CFTDS plot parameters                #
# # ------------------------------------ #
# folder_plot = header.folder_plot + 'CFADs'
# mod_names = ''
# letters_i=0
# n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
# n_cols = 4
# fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
# gs = fig.add_gridspec(n_rows, n_cols, hspace=0.03,wspace=0.03)
# axs = gs.subplots()
# # ------------------------------------ #
# ax_mean1=axs[-1,0]
# ax_mean1.set_ylabel('temperature [°C]')
# ax_mean1.set_xlabel('$Z_{H}$ [dBZ]')
# ax_mean1.set_xlim([mom_plot_dict('ZH')['mom_min'],
#                    mom_plot_dict('ZH')['mom_max']])
# ax_mean2 = axs[-1, 1]
# ax_mean2.set_ylabel('temperature [°C]')
# ax_mean2.set_xlabel('$Z_{DR}$ [dB]')
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],
#                    mom_plot_dict('ZDR')['mom_max']])
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],3.3])
# ax_mean3 = axs[-1, 2]
# ax_mean3.set_ylabel('temperature [°C]')
# ax_mean3.set_xlabel('$K_{DP}$ [°/km]')
# ax_mean3.set_xlim([mom_plot_dict('KDP')['mom_min'],
#                    mom_plot_dict('KDP')['mom_max']])
# ax_mean4 = axs[-1, 3]
# ax_mean4.set_ylabel('temperature [°C]')
# ax_mean4.set_xlabel('$\u03C1_{HV}$ [1]')
# ax_mean4.set_xlim([mom_plot_dict('RHOHV')['mom_min'],
#                    mom_plot_dict('RHOHV')['mom_max']])
# # --------------------------------------------------------------------------- #
# # CFTDs OBS row 1                                                             #
# # --------------------------------------------------------------------------- #
# color='black'
# # ------------------------------------ #
# current_row = 0
# current_col = 0
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZH_AC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean1.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZDR_AC_OC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=True,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean2.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='KDP_NC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=True,
#     data_max=data_max,
#     data_label=True,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean3.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='RHOHV_NC2P',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean4.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# # CFTDs CBAND SYN row i                                                       #
# # --------------------------------------------------------------------------- #
# for da_run, icon_emvorado_run, spin_up_mm, color, short_name in zip(
#         da_runs, icon_emvorado_runs, spin_up_mms, colors, short_names):
#     da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
#     model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
#                                 icon_emvorado_run.split('/')[1][5:]])
#     mod_names = '_'.join([mod_names, model_name_file])
#     model_name = '-'.join([da_run[4:],
#                            icon_emvorado_run.split('/')[0][5:],
#                            icon_emvorado_run.split('/')[1][5:],
#                            spin_up_mm + 'min'])
#     # ------------------------------------ #
#     current_row = current_row + 1
#     current_col = 0
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zdrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=True,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='kdpsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=True,
#         data_max=data_max,
#         data_label=True,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='rhvsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean4.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#
# # --------------------------------------------------------------------------- #
# # CFTDs SAVE                                                                  #
# # --------------------------------------------------------------------------- #
# ax_mean1.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean1.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean2.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean2.legend()
# letters_i=letters_i +1
#
# ax_mean3.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean4.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# for i_r in range(n_rows):
#     for i_c in range(n_cols):
#         axs[i_r,i_c].set_title('')
#         if i_r<n_rows-1:
#             axs[i_r,i_c].set_xlabel('')
#             axs[i_r, i_c].set_xticklabels('')
#
#         if i_c>0:
#             axs[i_r,i_c].set_ylabel('')
#             axs[i_r, i_c].set_yticklabels('')
#
#         if i_c==1:
#             axs[i_r,i_c].set_xlim([mom_plot_dict('ZDR')['mom_min'], 4.6])
#
# cmap = mpl.cm.YlGnBu
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
#              extend='max')
# axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])
#
# gs.tight_layout(fig, rect=[0, 0, 0.5, 1.0])
# if not os.path.exists(folder_plot):
#     os.makedirs(folder_plot)
#
# if len(dates) == 1:
#     dates_str = dates[0] + '_'
# elif dates == ['20210714', '20210713'] or dates == ['20210713', '20210714']:
#     dates_str = '20210713+14_'
# else:
#     dates_str = 'All_t_'
#
# if len(locations) == 1:
#     locations_str = locations[0] + '_'
# else:
#     locations_str = 'all_radars_'
#
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows)+ 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.pdf', format='pdf', transparent=True, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows) + 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
# plt.close()
#
# # --------------------------------------------------------------------------- #
# #                              PLOT B                                         #
# # --------------------------------------------------------------------------- #
# # INITIALIZATION                                                              #
# # --------------------------------------------------------------------------- #
# hhmm_start_cftds = '00:00'
# hhmm_end_cftds = '23:59'
# # ------------------------------------ #
# locations = ['ESS']
# dates = ['20210714']
# data_max = None
# testing = False
# height_min = 0  # in km
# height_max = 10  # in km
# bins_height = 20
# vert_temp = True
# temp_min = -20
# temp_max = 16
# bins_temp = 18
# elevation_degs = [5.5,4.5,3.5,2.5,1.5,0.5,8,12,17,25]
# letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
# filter_moms = False
# # ------------------------------------ #
# # MODELS                               #
# # ------------------------------------ #
# da_runs = []
# icon_emvorado_runs = []
# spin_up_mms = []
# short_names = []
# colors = []
# # ------------------------------------ #
# # SYN data row 2                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
# spin_up_mms.append('120')
# short_names.append('R0E2')
# colors.append('orange')
# # ------------------------------------ #
# # SYN data row 5                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R2E3')
# colors.append('blue')
#
# # --------------------------------------------------------------------------- #
# # CFTDs                                                                       #
# # --------------------------------------------------------------------------- #
# # ------------------------------------ #
# # CFTDS plot parameters                #
# # ------------------------------------ #
# folder_plot = header.folder_plot + 'CFADs'
# mod_names = ''
# letters_i=0
# n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
# n_cols = 4
# fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
# gs = fig.add_gridspec(n_rows, n_cols, hspace=0.03,wspace=0.03)
# axs = gs.subplots()
# # ------------------------------------ #
# ax_mean1=axs[-1,0]
# ax_mean1.set_ylabel('temperature [°C]')
# ax_mean1.set_xlabel('$Z_{H}$ [dBZ]')
# ax_mean1.set_xlim([mom_plot_dict('ZH')['mom_min'],
#                    mom_plot_dict('ZH')['mom_max']])
# ax_mean2 = axs[-1, 1]
# ax_mean2.set_ylabel('temperature [°C]')
# ax_mean2.set_xlabel('$Z_{DR}$ [dB]')
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],
#                    mom_plot_dict('ZDR')['mom_max']])
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],3.3])
# ax_mean3 = axs[-1, 2]
# ax_mean3.set_ylabel('temperature [°C]')
# ax_mean3.set_xlabel('$K_{DP}$ [°/km]')
# ax_mean3.set_xlim([mom_plot_dict('KDP')['mom_min'],
#                    mom_plot_dict('KDP')['mom_max']])
# ax_mean4 = axs[-1, 3]
# ax_mean4.set_ylabel('temperature [°C]')
# ax_mean4.set_xlabel('$\u03C1_{HV}$ [1]')
# ax_mean4.set_xlim([mom_plot_dict('RHOHV')['mom_min'],
#                    mom_plot_dict('RHOHV')['mom_max']])
# # --------------------------------------------------------------------------- #
# # CFTDs OBS row 1                                                             #
# # --------------------------------------------------------------------------- #
# color='black'
# # ------------------------------------ #
# current_row = 0
# current_col = 0
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZH_AC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean1.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZDR_AC_OC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=True,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean2.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='KDP_NC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=True,
#     data_max=data_max,
#     data_label=True,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean3.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='RHOHV_NC2P',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean4.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# # CFTDs CBAND SYN row i                                                       #
# # --------------------------------------------------------------------------- #
# for da_run, icon_emvorado_run, spin_up_mm, color, short_name in zip(
#         da_runs, icon_emvorado_runs, spin_up_mms, colors, short_names):
#     da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
#     model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
#                                 icon_emvorado_run.split('/')[1][5:]])
#     mod_names = '_'.join([mod_names, model_name_file])
#     model_name = '-'.join([da_run[4:],
#                            icon_emvorado_run.split('/')[0][5:],
#                            icon_emvorado_run.split('/')[1][5:],
#                            spin_up_mm + 'min'])
#     # ------------------------------------ #
#     current_row = current_row + 1
#     current_col = 0
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zdrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=True,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='kdpsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=True,
#         data_max=data_max,
#         data_label=True,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='rhvsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean4.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#
# # --------------------------------------------------------------------------- #
# # CFTDs SAVE                                                                  #
# # --------------------------------------------------------------------------- #
# ax_mean1.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean1.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean2.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean2.legend()
# letters_i=letters_i +1
#
# ax_mean3.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean4.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# for i_r in range(n_rows):
#     for i_c in range(n_cols):
#         axs[i_r,i_c].set_title('')
#         if i_r<n_rows-1:
#             axs[i_r,i_c].set_xlabel('')
#             axs[i_r, i_c].set_xticklabels('')
#
#         if i_c>0:
#             axs[i_r,i_c].set_ylabel('')
#             axs[i_r, i_c].set_yticklabels('')
#
#         if i_c==1:
#             axs[i_r,i_c].set_xlim([mom_plot_dict('ZDR')['mom_min'], 4.6])
#
# cmap = mpl.cm.YlGnBu
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
#              extend='max')
# axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])
#
# gs.tight_layout(fig, rect=[0, 0, 0.5, 1.0])
# if not os.path.exists(folder_plot):
#     os.makedirs(folder_plot)
#
# if len(dates) == 1:
#     dates_str = dates[0] + '_'
# elif dates == ['20210714', '20210713'] or dates == ['20210713', '20210714']:
#     dates_str = '20210713+14_'
# else:
#     dates_str = 'All_t_'
#
# if len(locations) == 1:
#     locations_str = locations[0] + '_'
# else:
#     locations_str = 'all_radars_'
#
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows)+ 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.pdf', format='pdf', transparent=True, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows) + 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
# plt.close()
#
# # --------------------------------------------------------------------------- #
# #                              PLOT C                                         #
# # --------------------------------------------------------------------------- #
# # INITIALIZATION                                                              #
# # --------------------------------------------------------------------------- #
# hhmm_start_cftds = '00:00'
# hhmm_end_cftds = '23:59'
# # ------------------------------------ #
# locations = ['ESS']
# dates = ['20210714']
# data_max = None
# testing = False
# height_min = 0  # in km
# height_max = 10  # in km
# bins_height = 20
# vert_temp = True
# temp_min = -20
# temp_max = 16
# bins_temp = 18
# elevation_degs = [12]
# letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
# filter_moms = False
# # ------------------------------------ #
# # MODELS                               #
# # ------------------------------------ #
# da_runs = []
# icon_emvorado_runs = []
# spin_up_mms = []
# short_names = []
# colors = []
# # ------------------------------------ #
# # SYN data row 1                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00010000.2')
# spin_up_mms.append('120')
# short_names.append('R0E1')
# colors.append('red')
# # ------------------------------------ #
# # SYN data row 2                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
# spin_up_mms.append('120')
# short_names.append('R0E2')
# colors.append('orange')
# # ------------------------------------ #
# # SYN data row 3                       #
# # ------------------------------------ #
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
# # ------------------------------------ #
# # SYN data row 5                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R2E3')
# colors.append('blue')
#
# # --------------------------------------------------------------------------- #
# # CFTDs                                                                       #
# # --------------------------------------------------------------------------- #
# # ------------------------------------ #
# # CFTDS plot parameters                #
# # ------------------------------------ #
# folder_plot = header.folder_plot + 'CFADs'
# mod_names = ''
# letters_i=0
# n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
# n_cols = 4
# fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
# gs = fig.add_gridspec(n_rows, n_cols, hspace=0.03,wspace=0.03)
# axs = gs.subplots()
# # ------------------------------------ #
# ax_mean1=axs[-1,0]
# ax_mean1.set_ylabel('temperature [°C]')
# ax_mean1.set_xlabel('$Z_{H}$ [dBZ]')
# ax_mean1.set_xlim([mom_plot_dict('ZH')['mom_min'],
#                    mom_plot_dict('ZH')['mom_max']])
# ax_mean2 = axs[-1, 1]
# ax_mean2.set_ylabel('temperature [°C]')
# ax_mean2.set_xlabel('$Z_{DR}$ [dB]')
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],
#                    mom_plot_dict('ZDR')['mom_max']])
# ax_mean2.set_xlim([mom_plot_dict('ZDR')['mom_min'],3.3])
# ax_mean3 = axs[-1, 2]
# ax_mean3.set_ylabel('temperature [°C]')
# ax_mean3.set_xlabel('$K_{DP}$ [°/km]')
# ax_mean3.set_xlim([mom_plot_dict('KDP')['mom_min'],
#                    mom_plot_dict('KDP')['mom_max']])
# ax_mean4 = axs[-1, 3]
# ax_mean4.set_ylabel('temperature [°C]')
# ax_mean4.set_xlabel('$\u03C1_{HV}$ [1]')
# ax_mean4.set_xlim([mom_plot_dict('RHOHV')['mom_min'],
#                    mom_plot_dict('RHOHV')['mom_max']])
# # --------------------------------------------------------------------------- #
# # CFTDs OBS row 1                                                             #
# # --------------------------------------------------------------------------- #
# color='black'
# # ------------------------------------ #
# current_row = 0
# current_col = 0
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZH_AC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean1.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='ZDR_AC_OC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=True,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean2.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='KDP_NC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=True,
#     data_max=data_max,
#     data_label=True,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean3.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# current_col = current_col + 1
# print(current_row)
# print(current_col)
# ax = axs[current_row, current_col]
# mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start_cftds,
#     hhmm_end=hhmm_end_cftds,
#     elevation_degs=elevation_degs,
#     da_icon_emvorado_run=None,
#     moment='RHOHV_NC2P',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
#     color=color,
#     plot_legend=False,
#     plot_data=False,
#     data_max=data_max,
#     panel=letters[letters_i] + ') obs',
# )
# # ------------------------------------ #
# letters_i=letters_i+1
# quant_prof = np.zeros([3, len(y_mid)])
# mean_prof = np.zeros(len(y_mid))
# for t_i in range(len(y_mid)):
#     wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#     quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                      return_pandas=False)
#     mean_prof[t_i] = wq.mean
#
# ax_mean4.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
#          linewidth=1, label='_nolegend_')
# ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#          linewidth=2, label='obs')
# # --------------------------------------------------------------------------- #
# # CFTDs CBAND SYN row i                                                       #
# # --------------------------------------------------------------------------- #
# for da_run, icon_emvorado_run, spin_up_mm, color, short_name in zip(
#         da_runs, icon_emvorado_runs, spin_up_mms, colors, short_names):
#     da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
#     model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
#                                 icon_emvorado_run.split('/')[1][5:]])
#     mod_names = '_'.join([mod_names, model_name_file])
#     model_name = '-'.join([da_run[4:],
#                            icon_emvorado_run.split('/')[0][5:],
#                            icon_emvorado_run.split('/')[1][5:],
#                            spin_up_mm + 'min'])
#     # ------------------------------------ #
#     current_row = current_row + 1
#     current_col = 0
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='zdrsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=True,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='kdpsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=True,
#         data_max=data_max,
#         data_label=True,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(quant_prof[2,], y_mid, color=color, ls='dashed', alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#     # ----------------------------------------------------------------------- #
#     current_col = current_col + 1
#     print(current_row)
#     print(current_col)
#     ax = axs[current_row, current_col]
#     mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
#         locations=locations,
#         dates=dates,
#         hhmm_start=hhmm_start_cftds,
#         hhmm_end=hhmm_end_cftds,
#         elevation_degs=elevation_degs,
#         da_icon_emvorado_run=da_icon_emvorado_run,
#         spin_up_mm=spin_up_mm,
#         moment='rhvsim',
#         vert_temp=vert_temp,
#         temp_min=temp_min,
#         temp_max=temp_max,
#         bins_temp=bins_temp,
#         height_min=height_min,  # in km
#         height_max=height_max,  # in km
#         bins_height=bins_height,
#         filter_moms=filter_moms,
#         ax=ax,
#         save=False,
#         color=color,
#         plot_legend=False,
#         plot_data=False,
#         data_max=data_max,
#         panel= letters[letters_i] +') ' + short_name,
#     )
#     letters_i=letters_i+1
#     # ------------------------------------ #
#     quant_prof = np.zeros([3, len(y_mid)])
#     mean_prof = np.zeros(len(y_mid))
#     for t_i in range(len(y_mid)):
#         wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
#         quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
#                                          return_pandas=False)
#         mean_prof[t_i] = wq.mean
#
#     ax_mean4.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(quant_prof[2,], y_mid, color=color, ls='dashed',alpha=0.3,
#                   linewidth=1, label='_nolegend_')
#     ax_mean4.plot(mean_prof, y_mid, color=color, ls='solid',
#                   linewidth=2,
#                   label=short_name)
#
# # --------------------------------------------------------------------------- #
# # CFTDs SAVE                                                                  #
# # --------------------------------------------------------------------------- #
# ax_mean1.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean1.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean2.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# ax_mean2.legend()
# letters_i=letters_i +1
#
# ax_mean3.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# ax_mean4.invert_yaxis()
# p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
# p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
# letters_i=letters_i +1
#
# for i_r in range(n_rows):
#     for i_c in range(n_cols):
#         axs[i_r,i_c].set_title('')
#         if i_r<n_rows-1:
#             axs[i_r,i_c].set_xlabel('')
#             axs[i_r, i_c].set_xticklabels('')
#
#         if i_c>0:
#             axs[i_r,i_c].set_ylabel('')
#             axs[i_r, i_c].set_yticklabels('')
#
#         if i_c==1:
#             axs[i_r,i_c].set_xlim([mom_plot_dict('ZDR')['mom_min'], 4.6])
#
# cmap = mpl.cm.YlGnBu
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
#              extend='max')
# axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])
#
# gs.tight_layout(fig, rect=[0, 0, 0.5, 1.0])
# if not os.path.exists(folder_plot):
#     os.makedirs(folder_plot)
#
# if len(dates) == 1:
#     dates_str = dates[0] + '_'
# elif dates == ['20210714', '20210713'] or dates == ['20210713', '20210714']:
#     dates_str = '20210713+14_'
# else:
#     dates_str = 'All_t_'
#
# if len(locations) == 1:
#     locations_str = locations[0] + '_'
# else:
#     locations_str = 'all_radars_'
#
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows)+ 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.pdf', format='pdf', transparent=True, bbox_inches='tight')
# plt.savefig(
#     folder_plot +
#     '/CFTDs_' + str(n_rows) + 'x4polmoms_PPI_' +
#     str(elevation_degs) + '°_' + dates_str + locations_str +
#     ['', 'mom_'][filter_moms] + mod_names[1:] +
#     ['', '_2e14'][testing] +
#     '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
# plt.close()

# --------------------------------------------------------------------------- #
#                              PLOT D                                         #
# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'
# ------------------------------------ #
locations = ['ESS']
dates = ['20210714']
data_max = None
testing = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
# elevation_degs = [12]
elevation_degs = [5.5,4.5,3.5,2.5,1.5,0.5,8,12,17,25]
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
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
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
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
colors.append('blue')

# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# CFTDS plot parameters                #
# ------------------------------------ #
folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
letters_i=0
n_rows = len(da_runs) + 1 + 1  # add one once more for mean of all
n_cols = 4
fig = plt.figure(figsize=(n_cols * 3, n_rows * 3))
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
mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_degs=elevation_degs,
    da_icon_emvorado_run=None,
    moment='ZH_AC',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean1.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean1.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean1.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
# --------------------------------------------------------------------------- #
current_col = current_col + 1
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_degs=elevation_degs,
    da_icon_emvorado_run=None,
    moment='ZDR_AC_OC',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean2.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean2.plot(quant_prof[2, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean2.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
# --------------------------------------------------------------------------- #
current_col = current_col + 1
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_degs=elevation_degs,
    da_icon_emvorado_run=None,
    moment='KDP_NC',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean

ax_mean3.plot(quant_prof[0, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean3.plot(quant_prof[2, ], y_mid, color=color, ls='dashed', alpha=0.3,
         linewidth=1, label='_nolegend_')
ax_mean3.plot(mean_prof, y_mid, color=color, ls='solid',
         linewidth=2, label='obs')
# --------------------------------------------------------------------------- #
current_col = current_col + 1
print(current_row)
print(current_col)
ax = axs[current_row, current_col]
mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_degs=elevation_degs,
    da_icon_emvorado_run=None,
    moment='RHOHV_NC2P',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
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
quant_prof = np.zeros([3, len(y_mid)])
mean_prof = np.zeros(len(y_mid))
for t_i in range(len(y_mid)):
    wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
    quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                     return_pandas=False)
    mean_prof[t_i] = wq.mean


ax_mean4.plot(quant_prof[0, ], y_mid, color=color, ls='dashed',alpha=0.3,
         linewidth=1, label='_nolegend_')
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
    mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_degs=elevation_degs,
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
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean1.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
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
    mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_degs=elevation_degs,
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
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean2.plot(quant_prof[0,], y_mid, color=color, ls='dashed', alpha=0.3,
                  linewidth=1, label='_nolegend_')
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
    mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_degs=elevation_degs,
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
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color,
        plot_legend=False,
        plot_data=True,
        data_max=data_max,
        data_label=True,
        panel= letters[letters_i] +') ' + short_name,
    )
    letters_i=letters_i+1
    # ------------------------------------ #
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean3.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
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
    mom_mid, y_mid, freq = plot_CFAD_or_CFTD_from_PPI_with_list_quick(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_degs=elevation_degs,
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
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid ,weights=freq[:,t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    ax_mean4.plot(quant_prof[0,], y_mid, color=color, ls='dashed',alpha=0.3,
                  linewidth=1, label='_nolegend_')
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
letters_i=letters_i +1

ax_mean2.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean2.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
ax_mean2.legend()
letters_i=letters_i +1

ax_mean3.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean3.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
letters_i=letters_i +1

ax_mean4.invert_yaxis()
p = plt.text(.04, .9, letters[letters_i] +')', transform=ax_mean4.transAxes)
p.set_bbox(dict(facecolor='white', alpha=0.8, linewidth=0.1))
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

cmap = mpl.cm.YlGnBu
norm = mpl.colors.Normalize(vmin=0, vmax=15)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs[-1,-1], orientation='vertical', label='frequency [%]',
             extend='max')
axs[-1,-1].set_xlim([mom_plot_dict('RHOHV')['mom_min'], 1.005])

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
    '/CFTDs_' + str(n_rows)+ 'x4polmoms_PPI_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/CFTDs_' + str(n_rows) + 'x4polmoms_PPI_' +
    str(elevation_degs) + '°_' + dates_str + locations_str +
    ['', 'mom_'][filter_moms] + mod_names[1:] +
    ['', '_2e14'][testing] +
    '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
plt.close()
