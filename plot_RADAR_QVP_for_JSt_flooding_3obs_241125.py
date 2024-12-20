#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 31.01.24                                                 #
# plot_syn_RADAR_QVP.py                                                       #
#                                                                             #
# Run the QVP functions in PLOT_SYN_RADAR.py for generating specific QVP plot.#
# --------------------------------------------------------------------------- #

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable

# --------------------------------------------------------------------------- #

location = 'ESS'
date = '20210714'
hhmm_start = '00:00'
hhmm_end = '23:55'
elevation_deg = 12

year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])

folder_plot = header.folder_qvp_plot
mod_names = 'bb_comb_das_ras07-qvp5minng01_sweeph5onem_polmoms_nc_' \
            '07-202107140003-202107142358-ess-10410'
n_rows = 4
n_cols = 4
# plt.figure(figsize=(n_cols * 9, n_rows * 7))
plt.figure(figsize=(n_cols * 12, n_rows * 9))
n_i = 0
current_row = 0

# --------------------------------------------------------------------------- #
# CBAND OBS 1                                                                 #
# --------------------------------------------------------------------------- #

folder_obs = '/automount/data02/agradar/operation_hydrometeors/' \
             'data/obs_qvp/OpHymet2-case09-20210714/' \
             '2021/2021-07/2021-07-14/ess/vol5minng01/07/'

path = folder_obs + 'ras07-qvp5minng01_sweeph5onem_polmoms_bb_nc_' \
                    '07-202107140003-202107142358-ess-10410.hd5'
obs_nc = xr.open_dataset(path)

top_height = 10
h_obs = obs_nc['height']

qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_zh_obs = obs_nc['ZH_AC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_temp_obs = obs_nc['temp'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ BB (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_obs,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_obs,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
obs_nc.close()

# --------------------------------------------------------------------------- #
# CBAND OBS 2                                                                 #
# --------------------------------------------------------------------------- #

folder_obs = '/automount/data02/agradar/operation_hydrometeors/' \
             'data/obs_qvp/OpHymet2-case09-20210714/' \
             '2021/2021-07/2021-07-14/ess/vol5minng01/07/'

path = folder_obs + 'ras07-qvp5minng01_sweeph5onem_polmoms_comb_nc_' \
                    '07-202107140003-202107142358-ess-10410.hd5'
obs_nc = xr.open_dataset(path)

top_height = 10
h_obs = obs_nc['height']

qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_zh_obs = obs_nc['ZH_AC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_temp_obs = obs_nc['temp'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ SD/LR (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_obs,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_obs,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
obs_nc.close()

# --------------------------------------------------------------------------- #
# CBAND OBS 3                                                                 #
# --------------------------------------------------------------------------- #

folder_obs = '/automount/data02/agradar/operation_hydrometeors/' \
             'data/obs_qvp/OpHymet2-case09-20210714/' \
             '2021/2021-07/2021-07-14/ess/vol5minng01/07/'

path = folder_obs + 'ras07-qvp5minng01_sweeph5onem_polmoms_das_nc_' \
                    '07-202107140003-202107142358-ess-10410.hd5'
obs_nc = xr.open_dataset(path)

top_height = 10
h_obs = obs_nc['height']

qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_zh_obs = obs_nc['ZH_AC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_temp_obs = obs_nc['temp'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ DAS (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_obs,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_obs,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
obs_nc.close()

# --------------------------------------------------------------------------- #
# CBAND OBS 4                                                                 #
# --------------------------------------------------------------------------- #

folder_obs = '/automount/data02/agradar/operation_hydrometeors/' \
             'data/obs_qvp/OpHymet2-case09-20210714/' \
             '2021/2021-07/2021-07-14/ess/vol5minng01/07/'

path = folder_obs + 'ras07-qvp5minng01_sweeph5onem_polmoms_das2_nc_' \
                    '07-202107140003-202107142358-ess-10410.hd5'
obs_nc = xr.open_dataset(path)

top_height = 10
h_obs = obs_nc['height']

qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_kdp_obs = obs_nc['KDP_NC'].sel(
    time=slice(date_start, date_end)) \
    .transpose(..., 'time')  # raw from radar
qvp_zh_obs = obs_nc['ZH_AC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['ZDR_AC_OC'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['RHOHV_NC2P'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_temp_obs = obs_nc['temp'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_obs,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ DAS2 (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_obs,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_obs,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_obs,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_obs,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
obs_nc.close()

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot + 'QVP_' + str(n_rows) + '_' + str(elevation_deg) + '_' +
    date + '_' + hhmm_start + '-' + hhmm_end + '_' +
    location + mod_names +
    '.pdf', format='pdf', transparent=True)
plt.close()
