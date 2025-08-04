#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 02.04.25                                                 #
# plot_syn_RADAR_CFTDs_of_4_ICE.py                                            #
#                                                                             #
# --------------------------------------------------------------------------- #

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
# dates = ['20170725']
hhmm_start = '00:00'
hhmm_end = '23:55'
elevation_deg = 12
top_height = 8
# ------------------------------------ #
da_run = 'ASS_2411'
icon_emvorado_run ='MAIN_2411.3/EMVO_00510000.2'
spin_up_mm = '120'
# ------------------------------------ #
# filter_entr = True  # TODO
filter_entr = False  # TODO
filter_moms = False
# ------------------------------------ #
# CFADs ?
# vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ?
vert_temp = True
temp_min = -20
temp_max = 0
bins_temp = 10
# -----------------------
mode = 'vol'
elevation_deg = 12
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
folder_plot = header.folder_qvp_plot
mod_names = ''
n_rows = 2
n_cols = 3
plt.figure(figsize=(n_cols * 5, n_rows * 4))  # TODO
n_i = 0
current_row = 0
current_col = 0

# --------------------------------------------------------------------------- #
# CBAND OBS row 1 PPI                                                         #
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
    moment='Dm_totice_qvp',
    title= 'Dm (Obs QVP-based)',
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
    moment='IWC_qvp',
    title= 'IWC (Obs QVP-based)',
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
    moment='Nt_totice_qvp',
    title= 'Nt (Obs QVO-based)',
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
# CBAND OBS row 2 QVO                                                         #
# --------------------------------------------------------------------------- #
current_row = current_row + 1
current_col = 1
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
    moment='Dm_totice',
    title= 'Dm (Obs PPI-based)',
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
    moment='IWC',
    title= 'IWC (Obs PPI-based)',
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
    moment='Nt_totice',
    title= 'Nt (Obs PPI-based)',
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
# # ------------------------------------ #
# # --------------------------------------------------------------------------- #
# # CBAND SYN row 3 PPI                                                         #
# # --------------------------------------------------------------------------- #
# current_row = current_row + 1
# current_col = 1
# n_i = (current_row - 1) * n_cols + current_col
# print(n_i)
# ax = plt.subplot(n_rows, n_cols, n_i)
# plot_CFAD_or_CFTD_from_QVP_with_list(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     elevation_deg=elevation_deg,
#     da_icon_emvorado_run= da_run + '/' + icon_emvorado_run,
#     moment='D0_totice',
#     # title= 'Dm (Obs PPI-based)',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_entr=filter_entr,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
# )
# # ------------------------------------ #
# current_col = current_col + 1
# n_i = (current_row - 1) * n_cols + current_col
# print(n_i)
# ax = plt.subplot(n_rows, n_cols, n_i)
# plot_CFAD_or_CFTD_from_QVP_with_list(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     elevation_deg=elevation_deg,
#     da_icon_emvorado_run=da_run + '/' + icon_emvorado_run,
#     moment='vol_qtotice',
#     # title= 'IWC (Obs PPI-based)',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_entr=filter_entr,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
# )
# # ------------------------------------ #
# current_col = current_col + 1
# n_i = (current_row - 1) * n_cols + current_col
# print(n_i)
# ax = plt.subplot(n_rows, n_cols, n_i)
# plot_CFAD_or_CFTD_from_QVP_with_list(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     elevation_deg=elevation_deg,
#     da_icon_emvorado_run=da_run + '/' + icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment='vol_qntotice',
#     # title= 'Nt (Obs PPI-based)',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_entr=filter_entr,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
# )
# # ------------------------------------ #

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
    folder_plot + '/newCFTD_3ice_' + str(elevation_deg) + '_' +
    dates_str + '_' + hhmm_start + '-' + hhmm_end + '_' +
    locations_str + ['', 'entr_'][filter_entr] + ['', 'mom_'][filter_moms] +
    mod_names +
    '.pdf', format='pdf', transparent=True)
plt.savefig(
    folder_plot + '/newCFTD_3ice_' + str(elevation_deg) + '_' +
    dates_str + '_' + hhmm_start + '-' + hhmm_end + '_' +
    locations_str + ['', 'entr_'][filter_entr] + ['', 'mom_'][filter_moms] +
    mod_names +
    '.png', format='png', transparent=True)
plt.close()
