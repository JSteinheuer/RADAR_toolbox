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
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
location = 'ESS'
date = '20210714'
locations = [location]
dates = [date]
hhmm_start = '00:00'
hhmm_end = '21:00'
elevation_deg = 12
top_height = 8
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
# ------------------------------------ #
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
da_runs.append('ASS_2411')  # ASS_newererer
# MAIN_newererererRH8_MP-RUC1.0  /
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00010000.2')
spin_up_mms.append('120')
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
da_runs.append('ASS_2411')  # ASS_newererer
# MAIN_newererererRH8_MP-RUC1.0  / # EMVO_dynwetgrow-BBold_sig-ice-25
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510000.2')
spin_up_mms.append('120')
# # ------------------------------------ #
# # SYN data row 3                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')  # ASS_newererer
# # MAIN_newererererRH8_MP-RUC1.0  / # EMVO_
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510200.2')
# spin_up_mms.append('120')

# # ------------------------------------ #
# # SYN data row 4                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')  # ASS_newererer
# # MAIN_  / # EMVO_dynwetgrow-BBold_sig-ice-25
# icon_emvorado_runs.append('MAIN_2411.03/EMVO_00510000.2')
# spin_up_mms.append('120')
# # ------------------------------------ #
# # SYN data row 5                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')  # ASS_newererer
# # MAIN_  / #
# icon_emvorado_runs.append('MAIN_2411.03/EMVO_00510100.2')
# spin_up_mms.append('120')

# # ------------------------------------ #
# # SYN data row i                       #
# # ------------------------------------ #
# da_runs.append('ASS_2411')  # ASS_newererer
# # # MAIN_newererererRH8_MP-RUC1.0  / # EMVO-neu
# icon_emvorado_runs.append('MAIN_2411.1/EMVO_00410000.2')
# spin_up_mms.append('120')
# ------------------------------------ #
# SYN data row i                       #
# ------------------------------------ #
da_runs.append('ASS_2411')  # ASS_newererer
# # MAIN_newererererRH8_MP-RUC1.0  / # EMVO_dynwetgrow-BBold_sig-ice-25
icon_emvorado_runs.append('MAIN_2411.1/EMVO_00510000.2')
spin_up_mms.append('120')

# ------------------------------------ #
# SYN data row i                       #
# ------------------------------------ #
da_runs.append('ASS_2411')  # ASS_newererer
# MAIN_  / # EMVO_dynwetgrow-BBold_sig-ice-25
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')

# ------------------------------------ #
year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
mode = 'vol'
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
folder_plot = header.folder_qvp_plot
mod_names = ''
n_rows = len(da_runs)+1
n_cols = 4
# plt.figure(figsize=(n_cols * 9, n_rows * 7))
# plt.figure(figsize=(n_cols * 8, n_rows * 6))
plt.figure(figsize=(n_cols * 10, n_rows * 6))
n_i = 0
current_row = 0

# --------------------------------------------------------------------------- #
# CBAND OBS row 1                                                             #
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
current_row = current_row + 1
n_i = n_cols * (current_row - 1)
# ------------------------------------ #
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
    cbar_title=r'reflectivity [dBZ]',
    title='Z$_{H}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
# ------------------------------------ #
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
    cbar_title=r'diff. reflectivity [dB]',
    title='Z$_{DR}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
# ------------------------------------ #
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
    cbar_title=r'spec. diff. phase [°/km]',
    title='K$_{DP}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
# ------------------------------------ #
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
    cbar_title=r'.crosscorr. coeff. [1]',
    title='$\u03C1_{hv}$ (C-band Obs. at ' + location + ')',
    ax=plt.subplot(n_rows, n_cols, n_i),
    mom_height_unit='km',
    scale_font=2.6,
    scale_numbers=2,
    top_height=top_height,
)
# ------------------------------------ #
obs_nc.close()

# --------------------------------------------------------------------------- #
# CBAND SYN row i                                                             #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm in zip(
        da_runs, icon_emvorado_runs, spin_up_mms):
    current_row = current_row + 1
    n_i = n_cols * (current_row - 1)
    date_start = '-'.join([year, mon, day, hhmm_start])
    date_end = '-'.join([year, mon, day, hhmm_end])
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
    print(top_height)
    # ------------------------------------ #
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
        cbar_title=r'reflectivity [dBZ]',
        title='Z$_{H}$ (%s)' % model_name,
        ax=plt.subplot(n_rows, n_cols, n_i),
        scale_font=2.6,
        scale_numbers=2,
        top_height=top_height,
        mom_height_unit='km',
    )
    # ------------------------------------ #
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
        cbar_title=r'diff. reflectivity [dB]',
        title='Z$_{DR}$ (%s)' % model_name,
        ax=plt.subplot(n_rows, n_cols, n_i),
        scale_font=2.6,
        scale_numbers=2,
        top_height=top_height,
        mom_height_unit='km',
    )
    # ------------------------------------ #
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
        cbar_title=r'spec. diff. phase [°/km]',
        title='K$_{DP}$ (%s)' % model_name,
        ax=plt.subplot(n_rows, n_cols, n_i),
        scale_font=2.6,
        scale_numbers=2,
        top_height=top_height,
        mom_height_unit='km',
    )
    # ------------------------------------ #
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
        cbar_title=r'crosscorr. coeff. [1]',
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
    folder_plot + 'QVP_4polmoms_' + str(elevation_deg) + '_' +
    date + '_' + hhmm_start + '-' + hhmm_end + '_' +
    location + mod_names +
    '.png', format='png', transparent=True)
plt.savefig(
    folder_plot + 'QVP_4polmoms_' + str(elevation_deg) + '_' +
    date + '_' + hhmm_start + '-' + hhmm_end + '_' +
    location + mod_names +
    '.pdf', format='pdf', transparent=True)
plt.close()


# --------------------------------------------------------------------------- #

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
    model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
                                icon_emvorado_run.split('/')[1][5:]])
    mod_names = '_'.join([mod_names, model_name_file])
    model_name = '-'.join([da_run[4:],
                           icon_emvorado_run.split('/')[0][5:],
                           icon_emvorado_run.split('/')[1][5:],
                           spin_up_mm + 'min'])
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
plt.savefig(
    folder_plot + '/CFTD_4polmoms_' + str(elevation_deg) + '_' +
    dates_str + '_' + hhmm_start + '-' + hhmm_end + '_' +
    locations_str + ['', 'entr_'][filter_entr] + ['', 'mom_'][filter_moms] +
    mod_names +
    '.png', format='png', transparent=True)
plt.close()
