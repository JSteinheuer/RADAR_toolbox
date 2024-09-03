#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_syn_RADAR_CFTDs.py                                                     #
#                                                                             #
# Run the function in PLOT_SYN_RADAR.py for generating specific CFTD.         #
# --------------------------------------------------------------------------- #

# TODO: change plot_CFAD_or_CFTD_from_QVP towards:
# TODO:     - include entropy filter

import HEADER_RADAR_toolbox as header
import matplotlib.pyplot as plt
import numpy as np
import warnings
from pathlib import Path

warnings.simplefilter('ignore')
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# PARAMS

n_rows = 3
n_cols = 4
# plt.figure(figsize=(n_cols * 5, n_rows * 4))
# plt.figure(figsize=(n_cols * 4, n_rows * 3))
plt.figure(figsize=(n_cols * 4.5, n_rows * 3.6))
n_i = 0
current_row = 0
current_col = 0
mod_names = ''
filter_entr_ML = True # as email 15.2.24
# filter_entr_ML = False
filter_ML = True
filter2 = False
# filter2 = True

# All
date = '20170725'
dates = [date]
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
# vert_temp = False
vert_temp = True
temp_min = -2
temp_max = 16
bins_temp = 9
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# locations = ['PRO']
locations = list(rad_dict().keys())
locations.remove('FLD')
elevation_deg = 12
vmax = 20
vmax = 15
# vmax=1

# MOMENTS
moment_1s = 'zrsim'
moment_1o = 'zh'
mom_min_1 = 0
mom_max_1 = 40
bins_mom_1 = 40

moment_2s = 'zdrsim'
moment_2o = 'zdr'
mom_min_2 = -0.5
mom_max_2 = 3
bins_mom_2 = 40

moment_3s = 'kdpsim_ml_corrected'
moment_3o = 'KDP_ML_corrected'
mom_min_3 = 0
mom_max_3 = 0.3
bins_mom_3 = 40

moment_4s = 'rhvsim'
moment_4o = 'rho'
mom_min_4 = 0.951
mom_max_4 = 1.051
bins_mom_4 = 40

moment_5s = 'D0_r_Bringi_syn'
moment_5o = 'D0_r_Bringi_obs'
moment_5s2 = 'D0_r'
mom_min_5 = 0.2
mom_max_5 = 2.6
bins_mom_5 = 40

moment_7s = 'vol_qntotice'
moment_7o = 'nt_obs'
# mom_min_7 = 0.
# mom_max_7 = 1000
mom_min_7 = -2
mom_max_7 = 2
bins_mom_7 = 50

# moment_8s = 'vol_qtotice'
moment_8s = 'qtotice'
moment_8o = 'iwc_obs'
mom_min_8 = 0.
mom_max_8 = 0.0006
mom_max_8 = 0.6
bins_mom_8 = 50

moment_9s = 'D0_totice'
moment_9o = 'dm_obs'
mom_min_9 = 0.
mom_max_9 = 10
bins_mom_9 = 50

# --------------------------------------------------------------------------- #
# Obs row 1
folder_obs = '/automount/realpep/upload/s6toscha/Statistik/' + \
             'CBAND_OBS_FERTIG_BRINGI_BUCH/NEU_PHI_NEU_TIME/'
paths_in = [folder_obs + 'fin_qvp_' + location.lower() + date + '.nc'
            for location in locations for date in dates]

title = 'C-band observation'
current_row = current_row + 1

# Obs 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_1o,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_1,
    mom_min=mom_min_1,
    bins_mom=bins_mom_1,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_2o,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_2,
    mom_min=mom_min_2,
    bins_mom=bins_mom_2,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# # Obs 3
# current_col = 3
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     paths_in=paths_in,
#     title=title,
#     moment=moment_3o,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_3,
#     mom_min=mom_min_3,
#     bins_mom=bins_mom_3,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Obs 4
# current_col = 4
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     paths_in=paths_in,
#     title=title,
#     moment=moment_4o,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_4,
#     mom_min=mom_min_4,
#     bins_mom=bins_mom_4,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )

# Obs 5
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_5o,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_5,
    mom_min=mom_min_5,
    bins_mom=bins_mom_5,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# # Syn 7
# current_col = 7
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     paths_in=paths_in,
#     title=title,
#     moment=moment_7o,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_7,
#     mom_min=mom_min_7,
#     bins_mom=bins_mom_7,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 8
# current_col = 8
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     paths_in=paths_in,
#     title=title,
#     moment=moment_8o,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_8,
#     mom_min=mom_min_8,
#     bins_mom=bins_mom_8,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 9
# current_col = 9
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     paths_in=paths_in,
#     title=title,
#     moment=moment_9o,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_9,
#     mom_min=mom_min_9,
#     bins_mom=bins_mom_9,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
# --------------------------------------------------------------------------- #
# CBAND SYN 1 row 2

da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2211.0/EMVO_00000000.2'
spin_up_mm = '60'

folder_syn = header.dir_data_qvp
current_row = current_row + 1
model_name = '-'.join([da_run[4:],
                        icon_emvorado_run.split('/')[0][5:],
                        icon_emvorado_run.split('/')[1][5:],
                        spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])

# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_1s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_1,
    mom_min=mom_min_1,
    bins_mom=bins_mom_1,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_2s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_2,
    mom_min=mom_min_2,
    bins_mom=bins_mom_2,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# # Syn 3
# current_col = 3
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_3s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_3,
#     mom_min=mom_min_3,
#     bins_mom=bins_mom_3,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 4
# current_col = 4
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_4s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_4,
#     mom_min=mom_min_4,
#     bins_mom=bins_mom_4,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )

# Syn 5
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_5s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_5,
    mom_min=mom_min_5,
    bins_mom=bins_mom_5,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
    filter2=filter2
)

# Syn 5 2
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_5s2,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_5,
    mom_min=mom_min_5,
    bins_mom=bins_mom_5,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# # Syn 7
# current_col = 7
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_7s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_7,
#     mom_min=mom_min_7,
#     bins_mom=bins_mom_7,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 8
# current_col = 8
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_8s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_8,
#     mom_min=mom_min_8,
#     bins_mom=bins_mom_8,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 9
# current_col = 9
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_9s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_9,
#     mom_min=mom_min_9,
#     bins_mom=bins_mom_9,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
# --------------------------------------------------------------------------- #
# CBAND SYN 3 row 4

da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2401.1/EMVO_00510000.2'
spin_up_mm = '60'

folder_syn = header.dir_data_qvp
current_row = current_row + 1
model_name = '-'.join([da_run[4:],
                        icon_emvorado_run.split('/')[0][5:],
                        icon_emvorado_run.split('/')[1][5:],
                        spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])

# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_1s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_1,
    mom_min=mom_min_1,
    bins_mom=bins_mom_1,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_2s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_2,
    mom_min=mom_min_2,
    bins_mom=bins_mom_2,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# # Syn 3
# current_col = 3
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_3s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_3,
#     mom_min=mom_min_3,
#     bins_mom=bins_mom_3,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 4
# current_col = 4
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_4s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_4,
#     mom_min=mom_min_4,
#     bins_mom=bins_mom_4,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )

# Syn 5
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_5s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_5,
    mom_min=mom_min_5,
    bins_mom=bins_mom_5,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
    filter2=filter2
)

# Syn 5 2
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_5s2,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_5,
    mom_min=mom_min_5,
    bins_mom=bins_mom_5,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)
# # Syn 7
# current_col = 7
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_7s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_7,
#     mom_min=mom_min_7,
#     bins_mom=bins_mom_7,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 8
# current_col = 8
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_8s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_8,
#     mom_min=mom_min_8,
#     bins_mom=bins_mom_8,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
#
# # Syn 9
# current_col = 9
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_9s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_9,
#     mom_min=mom_min_9,
#     bins_mom=bins_mom_9,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

save_path = header.folder_plot + 'CFADs/'
Path(save_path).mkdir(parents=True, exist_ok=True)

if vert_temp:
    save_name = 'CFTD_'
else:
    save_name = 'CFAD_'
plt.savefig(
    save_path + save_name + str(elevation_deg) + '_' + date +
    '_' + hhmm_start + '-' + hhmm_end + '_all' +
    mod_names +
    ['', '_filtered'][filter_entr_ML] +
    ['', '_g1'][filter2]
    + '_D0_0°C.pdf', format='pdf', transparent=True)
plt.close()