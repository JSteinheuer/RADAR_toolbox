#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 17.02.25                                                 #
# plot_syn_RADAR_CFTDs_of_4_polmoms.py                                        #
#                                                                             #
# Run the function in PLOT_SYN_RADAR.py for generating specific CFTD.         #
# --------------------------------------------------------------------------- #

import os
import HEADER_RADAR_toolbox as header
import matplotlib.pyplot as plt
import numpy as np
import warnings
from pathlib import Path
import glob

warnings.simplefilter('ignore')
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
# locations = ['PRO']
locations = list(rad_dict().keys())
dates = ['20170725']
# dates = [
#     '20170725',  # start this day
#     '20170810',
#     '20180809',
#     '20170719',
#     '20170720',
#     '20170724',
#     '20170726',
#     '20170727',
#     '20180728',
#     '20180923',
#     '20181202',
# ]
hhmm_start = '00:00'
hhmm_end = '23:55'
# hhmm_end = '09:55'
elevation_deg = 12
# filter_entr = True  # TODO
filter_entr = False  # TODO
# filter_moms = True  # TODO
filter_moms = False  # TODO
# ------------------------------------ #
# CFADs ?
# vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ?
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
da_runs = ['ASS_2211']
icon_emvorado_runs = ['MAIN_2411.0/EMVO_00000000.2']
spin_up_mms = ['60']
da_icon_emvorado_runs = [da_runs[-1] + '/' + icon_emvorado_runs[-1]]
if hhmm_end == '23:55':
    da_runs = []
    icon_emvorado_runs = []
    spin_up_mms = []
    da_icon_emvorado_runs = []

# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
da_runs.append('ASS_2211')
icon_emvorado_runs.append('MAIN_2211.0/EMVO_00000000.2')
spin_up_mms.append('60')
da_icon_emvorado_runs.append(da_runs[-1] + '/' + icon_emvorado_runs[-1])
# ------------------------------------ #
# SYN data row 3                       #
# ------------------------------------ #
da_runs.append('ASS_2211')
icon_emvorado_runs.append('MAIN_2308.1/EMVO_00400000.2')
spin_up_mms.append('60')
da_icon_emvorado_runs.append(da_runs[-1] + '/' + icon_emvorado_runs[-1])
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
da_runs.append('ASS_2211')
icon_emvorado_runs.append('MAIN_2401.1/EMVO_00510000.2')
spin_up_mms.append('60')
da_icon_emvorado_runs.append(da_runs[-1] + '/' + icon_emvorado_runs[-1])
# ------------------------------------ #
# SYN data row i                       #
# ------------------------------------ #
# da_runs.append('ASS_2211')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00000000.2')
# spin_up_mms.append('60')
# da_icon_emvorado_runs.append(da_runs[-1] + '/' + icon_emvorado_runs[-1])
# [...]
# ------------------------------------ #
# year = date[0:4]
# mon = date[4:6]
# day = date[6:8]
# date_start = '-'.join([year, mon, day, hhmm_start])
# date_end = '-'.join([year, mon, day, hhmm_end])
mode = 'vol'
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
n_rows = len(da_runs) + 1
n_cols = 4
# plt.figure(figsize=(n_cols * 8, n_rows * 6))
plt.figure(figsize=(n_cols * 5, n_rows * 4))  # TODO
n_i = 0
current_row = 0
current_col = 0

# --------------------------------------------------------------------------- #
# CBAND OBS row 1                                                             #
# --------------------------------------------------------------------------- #
current_row = current_row + 1
current_col = current_col + 1
n_i = (current_row - 1) * n_cols + current_col
print(n_i)
ax = plt.subplot(n_rows, n_cols, n_i)
plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
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
    filter_moms=filter_moms,
    ax=ax,
    save=False,
)
# ------------------------------------ #
current_col = current_col + 1
n_i = (current_row - 1) * n_cols + current_col
print(n_i)
ax = plt.subplot(n_rows, n_cols, n_i)
plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
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
    filter_moms=filter_moms,
    ax=ax,
    save=False,
)
# ------------------------------------ #
current_col = current_col + 1
n_i = (current_row - 1) * n_cols + current_col
print(n_i)
ax = plt.subplot(n_rows, n_cols, n_i)
plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
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
    filter_moms=filter_moms,
    ax=ax,
    save=False,
)
# ------------------------------------ #
current_col = current_col + 1
n_i = (current_row - 1) * n_cols + current_col
print(n_i)
ax = plt.subplot(n_rows, n_cols, n_i)
plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
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
    filter_moms=filter_moms,
    ax=ax,
    save=False,
)

# --------------------------------------------------------------------------- #
# CBAND SYN row i                                                             #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm in zip(
        da_runs, icon_emvorado_runs, spin_up_mms):
    da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
    current_col = 0
    model_name = '-'.join([da_run[4:],
                           icon_emvorado_run.split('/')[0][5:],
                           icon_emvorado_run.split('/')[1][5:],
                           spin_up_mm + 'min'])
    mod_names = '-'.join([mod_names, model_name])
    # ------------------------------------ #
    current_row = current_row + 1
    current_col = current_col + 1
    n_i = (current_row - 1) * n_cols + current_col
    print(n_i)
    ax = plt.subplot(n_rows, n_cols, n_i)
    plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start,
        hhmm_end=hhmm_end,
        elevation_deg=elevation_deg,
        da_icon_emvorado_run=da_icon_emvorado_run,
        moment='zrsim',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
    )
    # ------------------------------------ #
    current_col = current_col + 1
    n_i = (current_row - 1) * n_cols + current_col
    print(n_i)
    ax = plt.subplot(n_rows, n_cols, n_i)
    plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start,
        hhmm_end=hhmm_end,
        elevation_deg=elevation_deg,
        da_icon_emvorado_run=da_icon_emvorado_run,
        moment='zdrsim',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
    )
    # ------------------------------------ #
    current_col = current_col + 1
    n_i = (current_row - 1) * n_cols + current_col
    print(n_i)
    ax = plt.subplot(n_rows, n_cols, n_i)
    plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start,
        hhmm_end=hhmm_end,
        elevation_deg=elevation_deg,
        da_icon_emvorado_run=da_icon_emvorado_run,
        moment='kdpsim',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
    )
    # ------------------------------------ #
    current_col = current_col + 1
    n_i = (current_row - 1) * n_cols + current_col
    print(n_i)
    ax = plt.subplot(n_rows, n_cols, n_i)
    plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start,
        hhmm_end=hhmm_end,
        elevation_deg=elevation_deg,
        da_icon_emvorado_run=da_icon_emvorado_run,
        moment='rhvsim',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
    )

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

if len(dates) == 1:
    dates_str = dates[0]
else:
    dates_str = 'All_t_'

if len(locations) == 1:
    locations_str = locations[0]
else:
    locations_str = 'All_loc_'

plt.savefig(
    folder_plot + '/CFTD_4polmoms_' + str(elevation_deg) + '_' +
    dates_str + '_' + hhmm_start + '-' + hhmm_end + '_' +
    locations_str + ['', 'entr_'][filter_entr] + ['', 'mom_'][filter_moms] +
    mod_names +
    '.pdf', format='pdf', transparent=True)
plt.close()
