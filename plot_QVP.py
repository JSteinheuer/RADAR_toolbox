#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# plot_QVP.py                                                                 #
#                                                                             #
# Run the functions in PLOT_QVP.py for generating specific QVP plot.          #
# --------------------------------------------------------------------------- #

import os
import xarray as xr
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PLOT_QVP import plot_qvp_of_polarimetric_variable

# --------------------------------------------------------------------------- #

location = 'PRO'
date = '20170725'
filter = False
hhmm_start = '00:00'
# hhmm_end = '23:55'
hhmm_end = '10:00'

year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])

folder_plot = header.folder_qvp_plot
mod_names = ''
n_rows = 4
n_cols = 4
plt.figure(figsize=(n_cols * 9, n_rows * 7))
n_i = 0
current_row = 0

# --------------------------------------------------------------------------- #
# CBAND OBS 0                                                                 #
# --------------------------------------------------------------------------- #

folder_obs = '/automount/realpep/upload/s6toscha/Statistik/' + \
             'CBAND_OBS_FERTIG_BRINGI_BUCH/NEU_PHI_NEU_TIME/'
# TODO: TS style data?

path = folder_obs + 'fin_qvp_' + location.lower() + date + '.nc'
obs_nc = xr.open_dataset(path)

top_height = 10000
top_height = 8000

h_obs = obs_nc['height']
if filter:
    min_entropy_obs = obs_nc['min_entropy']
    qvp_zh_obs = obs_nc['zh']
    qvp_rho_obs = obs_nc['rho']
    qvp_kdp_obs = obs_nc['KDP_ML_corrected']
    height_ml_bottom_obs = obs_nc['height_ml_bottom_new_gia']
    obs_nc = xr.where(min_entropy_obs > 0.8, obs_nc, np.nan)
    obs_nc = xr.where(qvp_zh_obs > 0, obs_nc, np.nan)
    obs_nc = xr.where(qvp_rho_obs > 0.7, obs_nc, np.nan)
    obs_nc = xr.where(qvp_kdp_obs > 0.01, obs_nc, np.nan)
    obs_nc = xr.where(h_obs < height_ml_bottom_obs, obs_nc,
                      np.nan)
    qvp_kdp_obs = obs_nc['KDP_ML_corrected'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')  # ML corrected
else:
    qvp_kdp_obs = obs_nc['kdp'].sel(
        time=slice(date_start, date_end)) \
        .transpose(..., 'time')  # raw from radar

qvp_zh_obs = obs_nc['zh'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_zdr_obs = obs_nc['zdr'].sel(
    time=slice(date_start, date_end)).transpose(..., 'time')
qvp_rho_obs = obs_nc['rho'].sel(
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
    title='Z$_{DR}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
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
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
n_i = n_i + 1
obs_nc.close()

# --------------------------------------------------------------------------- #
# CBAND SYN 1                                                                 #
# --------------------------------------------------------------------------- #
folder_syn = header.dir_data_qvp
da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2308.1/EMVO_00500000.2'
spin_up_mm = '30'

date_start = '-'.join([year, mon, day, hhmm_start])
path = '/'.join([folder_syn, date, da_run, icon_emvorado_run,
                 str(spin_up_mm) + 'min_spinup', 'QVP_Syn_' +
                 location + '_' + date + '0000_' + date +
                 '2355.nc'])
syn_nc = xr.open_dataset(path)

top_height = 10
top_height = 8
h_syn = syn_nc['height']
if filter:
    min_entropy_syn = syn_nc['min_entropy']
    qvp_zh_syn = syn_nc['zrsim']
    qvp_rho_syn = syn_nc['rhvsim']
    qvp_kdp_syn = syn_nc['kdpsim']
    height_ml_bottom_syn = syn_nc['mlh_bottom']
    syn_nc = xr.where(min_entropy_syn > 0.8, syn_nc, np.nan)
    syn_nc = xr.where(qvp_zh_syn > 0, syn_nc, np.nan)
    syn_nc = xr.where(qvp_rho_syn > 0.7, syn_nc, np.nan)
    syn_nc = xr.where(qvp_kdp_syn > 0.01, syn_nc, np.nan)
    syn_nc = xr.where(h_syn < height_ml_bottom_syn, syn_nc,
                      np.nan)

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
    time=slice(date_start, date_end)) \
                   .transpose(..., 'time') - 273.15
# --------------------------------------------------------------------------- #
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name + ['', 'ML'][filter]])
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_syn,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_syn,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
syn_nc.close()

# --------------------------------------------------------------------------- #
# CBAND SYN 2                                                                 #
# --------------------------------------------------------------------------- #
folder_syn = header.dir_data_qvp
da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2308.1/EMVO_00600000.2'
spin_up_mm = '30'

date_start = '-'.join([year, mon, day, hhmm_start])
path = '/'.join([folder_syn, date, da_run, icon_emvorado_run,
                 str(spin_up_mm) + 'min_spinup', 'QVP_Syn_' +
                 location + '_' + date + '0000_' + date +
                 '2355.nc'])
syn_nc = xr.open_dataset(path)

top_height = 10
top_height = 8
h_syn = syn_nc['height']
if filter:
    min_entropy_syn = syn_nc['min_entropy']
    qvp_zh_syn = syn_nc['zrsim']
    qvp_rho_syn = syn_nc['rhvsim']
    qvp_kdp_syn = syn_nc['kdpsim']
    height_ml_bottom_syn = syn_nc['mlh_bottom']
    syn_nc = xr.where(min_entropy_syn > 0.8, syn_nc, np.nan)
    syn_nc = xr.where(qvp_zh_syn > 0, syn_nc, np.nan)
    syn_nc = xr.where(qvp_rho_syn > 0.7, syn_nc, np.nan)
    syn_nc = xr.where(qvp_kdp_syn > 0.01, syn_nc, np.nan)
    syn_nc = xr.where(h_syn < height_ml_bottom_syn, syn_nc,
                      np.nan)

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
    time=slice(date_start, date_end)) \
                   .transpose(..., 'time') - 273.15
# --------------------------------------------------------------------------- #
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name + ['', 'ML'][filter]])
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_syn,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_syn,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
syn_nc.close()
# --------------------------------------------------------------------------- #
# CBAND SYN 3                                                                 #
# --------------------------------------------------------------------------- #
folder_syn = header.dir_data_qvp
da_run = 'ASS_2109'
icon_emvorado_run = 'MAIN_2109.0/EMVO_00000000.2'
spin_up_mm = '30'

date_start = '-'.join([year, mon, day, hhmm_start])
path = '/'.join([folder_syn, date, da_run, icon_emvorado_run,
                 str(spin_up_mm) + 'min_spinup', 'QVP_Syn_' +
                 location + '_' + date + '0000_' + date +
                 '2355.nc'])
syn_nc = xr.open_dataset(path)

top_height = 10
top_height = 8
h_syn = syn_nc['height']
if filter:
    min_entropy_syn = syn_nc['min_entropy']
    qvp_zh_syn = syn_nc['zrsim']
    qvp_rho_syn = syn_nc['rhvsim']
    qvp_kdp_syn = syn_nc['kdpsim']
    height_ml_bottom_syn = syn_nc['mlh_bottom']
    syn_nc = xr.where(min_entropy_syn > 0.8, syn_nc, np.nan)
    syn_nc = xr.where(qvp_zh_syn > 0, syn_nc, np.nan)
    syn_nc = xr.where(qvp_rho_syn > 0.7, syn_nc, np.nan)
    syn_nc = xr.where(qvp_kdp_syn > 0.01, syn_nc, np.nan)
    syn_nc = xr.where(h_syn < height_ml_bottom_syn, syn_nc,
                      np.nan)

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
    time=slice(date_start, date_end)) \
                   .transpose(..., 'time') - 273.15
# --------------------------------------------------------------------------- #
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name + ['', 'ML'][filter]])
# --------------------------------------------------------------------------- #
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# --------------------------------------------------------------------------- #
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zh_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zh,
    levels=header.levels_zh,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Reflectivity [dBZ]',
    title='Z$_{H}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_zdr_syn,
    cmap=header.cmap_radar,
    norm=header.norm_zdr,
    levels=header.levels_zdr,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Differential Reflectivity [dB]',
    title='Z$_{DR}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_kdp_syn,
    cmap=header.cmap_radar,
    norm=header.norm_kdp,
    levels=header.levels_kdp,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Specific differential Phase [°/km]',
    title='K$_{DP}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
n_i = n_i + 1
plot_qvp_of_polarimetric_variable(
    mom=qvp_rho_syn,
    cmap=header.cmap_radar,
    norm=header.norm_rhohv,
    levels=header.levels_rhohv,
    mom_cs=qvp_zh_syn,
    levels_cs=np.arange(-50, 60, 5),
    mom_cf=qvp_temp_syn,
    levels_cf=np.arange(-50, 60, 5),
    cbar_title=r'Crosscorrelation Coefficient [1]',
    title='$\u03C1_{hv}$ (%s)' % model_name,
    ax=plt.subplot(n_rows, n_cols, n_i),
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
    mom_height_unit='km',
)
syn_nc.close()
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(
    folder_plot + date + '_' + hhmm_start + '-' + hhmm_end + '_' +
    location + '_' + mod_names +
    '.pdf', format='pdf', transparent=True)
plt.close()
