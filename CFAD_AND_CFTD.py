#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.02.24                                                 #
# CFAD_AND_CFTD.py                                                            #
#                                                                             #
# plotting routine for Continuous-frequency-altitude-diagram (CFAD) and       #
# Continuous-frequency-temperatur-diagram (CFTD).                             #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import pandas as pd
import numpy as np
import glob as glob
import wradlib as wrl
import random
import gc
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
# import bz2
# import pickle
# import _pickle as cPickle
# import hmcp
import xarray as xr
import warnings
from statsmodels.stats.weightstats import DescrStatsW
from pathlib import Path

warnings.simplefilter('ignore')


def plot_CFAD_or_CFTD_from_QVP(
        date='20170725',
        hhmm_start='00:00',
        hhmm_end='23:55',
        path_in=None,
        title=None,
        folder_syn=header.dir_data_qvp,
        location='PRO',
        elevation_deg=12,
        da_run='ASS_2211',
        icon_emvorado_run='MAIN_2401.3/EMVO_00500000.2',
        spin_up_mm='60',
        moment='zrsim',
        mom_min=0,
        mom_max=40,
        bins_mom=40,
        vert_temp=True,  # CFTD
        temp_min=-20,
        temp_max=16,
        bins_temp=18,
        height_min=0,  # in km
        height_max=10,  # in km
        bins_height=20,
        ax=None,
        save=False,
        save_path=header.folder_plot + 'CFADs/',
        save_name='test_CFAD'
):
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    date_start = '-'.join([year, mon, day, hhmm_start])
    date_end = '-'.join([year, mon, day, hhmm_end])
    if path_in:
        if not title:
            title = path_in.split('/')[-1]
    else:
        path_in = '/'.join([folder_syn, date, da_run, icon_emvorado_run,
                            str(spin_up_mm) + 'min_spinup', 'QVP_' +
                            str(elevation_deg) + '_Syn_' + location + '_' +
                            date + '0000_' + date + '2355.nc'])
        if not title:
            title = '-'.join([da_run[4:],
                              icon_emvorado_run.split('/')[0][5:],
                              icon_emvorado_run.split('/')[1][5:],
                              spin_up_mm + 'min'])

    # OPEN
    syn_nc = xr.open_dataset(path_in)
    syn_nc = syn_nc.sel(time=slice(date_start, date_end))
    mom = syn_nc[moment].transpose('time', ...).values
    if vert_temp:
        y = syn_nc.temp.transpose('time', ...).values
        y_min, y_max, bins_y = temp_min, temp_max, bins_temp
        if 'units' in syn_nc.temp.attrs:  # else: TS in °C
            if syn_nc.temp.units == 'K':
                y = y - 273.15
    else:
        y = np.array([syn_nc.height.values] * syn_nc.time.size) / 1000
        y_min, y_max, bins_y = height_min, height_max, bins_height
        if 'units' in syn_nc.height.attrs:  # else: TS in m
            if syn_nc.height.units == 'km':
                y = y * 1000

    mask = (y >= y_min) & (y <= y_max) & \
           (mom >= mom_min) & (mom <= mom_max)  # TODO: realistic values might
    # TODO: ... be set nan
    mom = mom[mask]
    y = y[mask]
    i_sort = np.argsort(y)
    mom = mom[i_sort]
    y = y[i_sort]

    # PLOT
    if ax is None:
        plt.figure(figsize=(6, 5))

    a = plt.hist(y, bins=bins_y, )
    weights = np.repeat(100 / a[0], np.int16(a[0]))
    h2d, mom2d, y2d, fg = plt.hist2d(mom, y, bins=[bins_mom, bins_y],
                                     range=[[mom_min, mom_max],
                                            [y_min, y_max]],
                                     weights=weights, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]')
    if vert_temp:
        plt.gca().invert_yaxis()
        plt.ylabel('temperature (°C)')
    else:
        plt.ylabel(r'height (km)')

    if 'standard_name' in syn_nc[moment].attrs:  # else: TS without any
        plt.xlabel(
            syn_nc[moment].standard_name + ' (' + syn_nc[moment].units + ')')
    else:
        plt.xlabel(moment)

    plt.title(title)
    y_mid = y2d[1:] / 2 + y2d[:-1] / 2
    mom_mid = mom2d[1:] / 2 + mom2d[:-1] / 2
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid, weights=h2d[:, t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    plt.plot(quant_prof[0,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.2}$')
    plt.plot(quant_prof[1,], y_mid, color='red', ls='solid',
             linewidth=2, label='$Q_{0.5}$')
    plt.plot(quant_prof[2,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.8}$')
    plt.plot(mean_prof, y_mid, color='orange', ls='solid',
             linewidth=2, label='$\mu$')
    plt.legend()
    syn_nc.close()
    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)


# --------------------------------------------------------------------------- #
# PARAMS

n_rows = 4
n_cols = 4
plt.figure(figsize=(n_cols * 5, n_rows * 4))
n_i = 0
current_row = 1
current_col = 1
mod_names = ''

# All
date = '20170725'
hhmm_start = '00:00'
hhmm_end = '10:00'
# vert_temp = False
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
location = 'PRO'
elevation_deg = 12

# MOMENTS
moment_1s = 'zrsim'
moment_1o = 'zh'
mom_min_1 = 0
mom_max_1 = 40
bins_mom_1 = 40

moment_2s = 'zdrsim'
moment_2o = 'zdr'
mom_min_2 = -0.5
mom_max_2 = 1.5
bins_mom_2 = 40

moment_3s = 'kdpsim_ml_corrected'
moment_3o = 'KDP_ML_corrected'
mom_min_3 = 0
mom_max_3 = 0.3
bins_mom_3 = 40

moment_4s = 'rhvsim'
moment_4o = 'rho'
mom_min_4 = 0.9601
mom_max_4 = 1.001
bins_mom_4 = 40

# Obs row 1
folder_obs = '/automount/realpep/upload/s6toscha/Statistik/' + \
             'CBAND_OBS_FERTIG_BRINGI_BUCH/NEU_PHI_NEU_TIME/'
path_in = folder_obs + 'fin_qvp_' + location.lower() + date + '.nc'
title = 'C-band observation'
current_row = 1

# Obs 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    path_in=path_in,
    title=title,
    moment=moment_1o,
    # moment=moment_1s,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    path_in=path_in,
    title=title,
    moment=moment_2o,
    # moment=moment_2s,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 3
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    path_in=path_in,
    title=title,
    moment=moment_3o,
    # moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    path_in=path_in,
    title=title,
    moment=moment_4o,
    # moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# --------------------------------------------------------------------------- #
# SYN row 2
da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2308.0/EMVO_00400000.2'
spin_up_mm = '60'
folder_syn = header.dir_data_qvp
current_row = 2
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])
# --------------------------------------------------------------------------- #
# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_1o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_2o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 3
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_3o,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_4o,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# --------------------------------------------------------------------------- #
# SYN row 3
da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2308.1/EMVO_00500000.2'
spin_up_mm = '60'
folder_syn = header.dir_data_qvp
current_row = 3
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])
# --------------------------------------------------------------------------- #

# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_1o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_2o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 3
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_3o,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_4o,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# --------------------------------------------------------------------------- #
# SYN row 4
da_run = 'ASS_2211'
icon_emvorado_run = 'MAIN_2401.3/EMVO_00500000.2'
spin_up_mm = '60'
folder_syn = header.dir_data_qvp
current_row = 4
model_name = '-'.join([da_run[4:],
                       icon_emvorado_run.split('/')[0][5:],
                       icon_emvorado_run.split('/')[1][5:],
                       spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])
# --------------------------------------------------------------------------- #

# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_1o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_2o,
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
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 3
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_3o,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    date=date,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    location=location,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    # path_in=path_in,
    # title=title,
    # moment=moment_4o,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

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
    save_path + save_name + str(elevation_deg) + '_' +
    date + '_' + hhmm_start + '-' + hhmm_end + '_' +
    location + mod_names +
    '.pdf', format='pdf', transparent=True)
plt.close()
