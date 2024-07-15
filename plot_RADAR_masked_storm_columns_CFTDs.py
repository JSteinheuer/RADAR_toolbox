#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 19.06.24                                                 #
# plot_RADAR_masked_storm_columns_CFTDs.py                                    #
#                                                                             #
# masking for every sweep weather in the column is a storm bin (max(zh)>40),  #
# it is storm surrounding 'envelope' (40>max(zh)>25), or neither.             #
# --------------------------------------------------------------------------- #

import datatree as dttree
from statsmodels.stats.weightstats import DescrStatsW
import numpy as np
import pandas as pd
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr
from scipy.ndimage import uniform_filter, gaussian_filter
import time
import warnings
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy
import sklearn
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore")
from pathlib import Path
from PLOT_RADAR import plot_CFAD_or_CFTD_from_2_arrays


# J. Steinheuer
def mask_storm_columns(date="20220520", location='fld',
                       zh_thresholds=[0, 80],
                       zdr_thresholds=[-1, np.nan],
                       kdp_thresholds=[-2, np.nan],
                       rho_thresholds=[.8, np.nan],
                       dist_thresholds=[5, 100],
                       dir_data_obs=header.dir_data_obs,
                       core_zh=45,
                       envelope_zh=25,
                       ):
    """
    .

    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
    dir_data_obs : directory to search for input cases
                  (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
    """
    # initialization
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    mode = 'vol'
    mask_times = np.array([])
    mask_lat = np.array([])
    mask_lon = np.array([])
    mask_alt = np.array([])
    mask_temp = np.array([])
    mask_kdp = np.array([])
    mask_zdr = np.array([])
    mask_zh = np.array([])
    mask_rho = np.array([])
    init = 0
    for elevation_deg in header.ELEVATIONS_ALL:
        sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                                   float(elevation_deg))[0][0])
        path_in = "/".join([dir_data_obs + '*',
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep,
                            'ras*_polmoms_*'])
        files = sorted(glob.glob(path_in))
        if not files:
            print('no input data *_allmoms_*')
            print('return')
        else:
            path_in = files[0]

        # load
        data = dttree.open_datatree(path_in)[
            'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
        data = data.transpose('time', 'azimuth', 'range')
        data = data.assign({'h_dist': data.range * np.cos(
            data.elevation.values * 2 * np.pi / 360)})

        # cut donut around radar
        data = data.where(data.h_dist >= dist_thresholds[0] * 1000)
        data = data.where(data.h_dist <= dist_thresholds[1] * 1000)

        # all coordinates
        all_lat = data.lat.expand_dims(dim={'time': data.time.values},
                                       axis=0).values.ravel()
        all_lon = data.lon.expand_dims(dim={'time': data.time.values},
                                       axis=0).values.ravel()
        all_alt = data.alt.expand_dims(dim={'time': data.time.values,
                                            'azimuth': data.azimuth.values},
                                       axis=[0, 1]).values.ravel()
        all_times = data.time.expand_dims(dim={'range': data.range.values,
                                               'azimuth': data.azimuth.values},
                                          axis=[2, 1]).values.ravel()

        # all_variables
        all_kdp = data.KDP_NC.values.ravel()
        all_rho = data.RHOHV_NC2P.values.ravel()
        all_zdr = data.ZDR_AC_OC.values.ravel()
        all_zh = data.ZH_AC.values.ravel()
        all_temp = data.temp.values.ravel()

        # mask due to threshold and what is already nan
        mask = (all_zh >= zh_thresholds[0]) & \
               (all_zh <= zh_thresholds[1]) & \
               (all_kdp >= kdp_thresholds[0]) & \
               (all_kdp <= kdp_thresholds[1]) & \
               (all_zdr >= zdr_thresholds[0]) & \
               (all_zdr <= zdr_thresholds[1]) & \
               (all_rho >= rho_thresholds[0]) & \
               (all_rho <= rho_thresholds[1])

        # apply masking
        mask_times = np.append(mask_times, np.round(
            (all_times[mask] - np.datetime64('2021-01-01')
             ).astype(int) * 1e-9, 0).astype(int))
        mask_lat = np.append(mask_lat, all_lat[mask])
        mask_lon = np.append(mask_lon, all_lon[mask])
        mask_alt = np.append(mask_alt, all_alt[mask])
        mask_temp = np.append(mask_temp, all_temp[mask])
        mask_kdp = np.append(mask_kdp, all_kdp[mask])
        mask_zdr = np.append(mask_zdr, all_zdr[mask])
        mask_zh = np.append(mask_zh, all_zh[mask])
        mask_rho = np.append(mask_rho, all_rho[mask])

        if init == 0:
            init = 1
            # src grid in lat x lon x time
            # ~converting km to deg lon/lat
            lon_w = data.longitude.values - dist_thresholds[1] / \
                    (110 * np.cos(data.latitude.values * 2 * np.pi / 360))
            lon_e = data.longitude.values + dist_thresholds[1] / \
                    (110 * np.cos(data.latitude.values * 2 * np.pi / 360))
            lat_n = data.latitude.values + dist_thresholds[1] / (110)
            lat_s = data.latitude.values - dist_thresholds[1] / (110)

            grid_lat = np.linspace(lat_s, lat_n, dist_thresholds[1])
            grid_lon = np.linspace(lon_w, lon_e, dist_thresholds[1])
            # times always in middle of 5min, i.e. 12:37:30
            grid_times = np.round(
                (data.time.values - np.datetime64('2021-01-01')
                 ).astype(int) * 1e-9, 0).astype(int) - \
                         np.round(
                             (data.time.values - np.datetime64('2021-01-01')
                              ).astype(int) * 1e-9, 0).astype(int) % 300 + 150

        data.close()
        print(path_in.split('/')[-1] + ' loaded')

    id_lat = np.digitize(mask_lat, grid_lat)
    id_lon = np.digitize(mask_lon, grid_lon)
    id_times = np.digitize(mask_times, grid_times)
    id_llt = id_lat + id_lon * 100 + id_times * 10000

    sort_llt = np.argsort(id_llt)
    id_llt = id_llt[sort_llt]
    mask_times = mask_times[sort_llt]
    mask_lat = mask_lat[sort_llt]
    mask_lon = mask_lon[sort_llt]
    mask_alt = mask_alt[sort_llt]
    mask_temp = mask_temp[sort_llt]
    mask_kdp = mask_kdp[sort_llt]
    mask_zdr = mask_zdr[sort_llt]
    mask_zh = mask_zh[sort_llt]
    mask_rho = mask_rho[sort_llt]

    max_zh = np.repeat(0.0, mask_zh.size)
    id_llt_diff = np.append(1, id_llt[1:] - id_llt[:-1])
    i_starters = np.where(id_llt_diff > 0)[0]
    i_enders = np.append(i_starters[1:], id_llt_diff.size)
    for i_st, i_en in zip(i_starters, i_enders):
        max_zh[i_st:i_en] = np.max(mask_zh[i_st:i_en])

    column_type = np.repeat(0, mask_zh.size)
    column_type[max_zh >= core_zh] = 3
    column_type[max_zh < core_zh] = 2
    column_type[max_zh < envelope_zh] = 1
    return dict(mask_times=mask_times,
                mask_lat=mask_lat,
                mask_lon=mask_lon,
                mask_alt=mask_alt,
                mask_temp=mask_temp,
                mask_kdp=mask_kdp,
                mask_zdr=mask_zdr,
                mask_zh=mask_zh,
                mask_rho=mask_rho,
                max_zh=max_zh,
                column_type=column_type,
                )


# --------------------------------------------------------------------------- #
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# # SET PARAMS:
# DATES = [
#     "20210604",  # case01
#     "20210620", "20210621",  # case02
#     "20210628", "20210629",  # case03
#     "20220519", "20220520",  # case04
#     "20220623", "20220624", "20220625",  # case05
#     "20220626", "20220627", "20220628",  # case06+07
#     "20220630", "20220701",  # case08
#     "20210714",  # case09
#     "20221222",  # case10
# ]
# LOCATIONS = [
#     'asb', 'boo', 'drs', 'eis', 'ess', 'fbg',
#     'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
#     'oft', 'pro', 'ros', 'tur', 'umd',
# ]
# ELEVATIONS = np.array([
#     5.5,
#     4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0,
#     25.0,
# ])
# overwrite = False
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
# for date in DATES:
#     for location in LOCATIONS:
date = "20220520"
location = 'fld'
zh_thresholds = [0, 80]
zdr_thresholds = [-1, 10]
kdp_thresholds = [-2, 10]
rho_thresholds = [.8, 2]
dist_thresholds = [5, 100]
dir_data_obs = header.dir_data_obs
core_zh = 45
envelope_zh = 25
# --------------------------------------------------------------------------- #
masked_storm = mask_storm_columns(date=date, location=location,
                                  zh_thresholds=zh_thresholds,
                                  zdr_thresholds=zdr_thresholds,
                                  kdp_thresholds=kdp_thresholds,
                                  rho_thresholds=rho_thresholds,
                                  dist_thresholds=dist_thresholds,
                                  dir_data_obs=dir_data_obs,
                                  core_zh=core_zh,
                                  envelope_zh=envelope_zh)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
vert_temp = True
temp_min = -20
temp_max = 18
bins_temp = 19
height_min = 0
height_max = 10
bins_height = 20
vmax = None
vmax = 10
# ax = None
save = False
save_name = 'dummy'
save_path = header.folder_plot + 'CFADs/'
# --------------------------------------------------------------------------- #
n_rows = 3
n_cols = 4
n_i = 0
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
# --------------------------------------------------------------------------- #
mom_min = 0
mom_max = 60
bins_mom = 60
# --------------------------------------------------------------------------- #
n_i = 2 * n_cols + 1
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm core ZH [dBZ]'
# save_name = 'storm_core_zh'
x = masked_storm['mask_zh'][masked_storm['column_type'] == 3]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 3]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 0 * n_cols + 1
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'all ZH [dBZ]'
# save_name = 'all_data_zh'
x = masked_storm['mask_zh']
y = masked_storm['mask_temp']
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 1 * n_cols + 1
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm envelope ZH [dBZ]'
# save_name = 'storm_envelope_zh'
x = masked_storm['mask_zh'][masked_storm['column_type'] == 2]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 2]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
mom_min = -0.5
mom_max = 3
bins_mom = 70
# --------------------------------------------------------------------------- #
n_i = 2 * n_cols + 2
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm core ZDR [dB]'
# save_name = 'storm_core_zdr'
x = masked_storm['mask_zdr'][masked_storm['column_type'] == 3]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 3]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 0 * n_cols + 2
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'all ZDR [dB]'
# save_name = 'all_data_zdr'
x = masked_storm['mask_zdr']
y = masked_storm['mask_temp']
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 1 * n_cols + 2
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm envelope ZDR [db]'
# save_name = 'storm_envelope_zdr'
x = masked_storm['mask_zdr'][masked_storm['column_type'] == 2]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 2]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
# mom_min = -.2
# mom_max = .6
# bins_mom = 80
mom_min = -.1
mom_max = .5
bins_mom = 60
# --------------------------------------------------------------------------- #
n_i = 2 * n_cols + 3
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm core KDP [°]'
# save_name = 'storm_core_kdp'
x = masked_storm['mask_kdp'][masked_storm['column_type'] == 3]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 3]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 0 * n_cols + 3
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'all KDP [°]'
# save_name = 'all_data_kdp'
x = masked_storm['mask_kdp']
y = masked_storm['mask_temp']
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 1 * n_cols + 3
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm envelope KDP [°]'
# save_name = 'storm_envelope_kdp'
x = masked_storm['mask_kdp'][masked_storm['column_type'] == 2]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 2]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
# mom_min = .8
# mom_max = 1.05
# bins_mom = 82
# mom_min = .9
# mom_max = 1.01
# bins_mom = 55
mom_min = .931
mom_max = 1.01
bins_mom = 70
# --------------------------------------------------------------------------- #
n_i = 2 * n_cols + 4
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm core RHO [1]'
# save_name = 'storm_core_rho'
x = masked_storm['mask_rho'][masked_storm['column_type'] == 3]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 3]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
n_i = 0 * n_cols + 4
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'all RHO [1]'
# save_name = 'all_data_rho'
x = masked_storm['mask_rho']
y = masked_storm['mask_temp']
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
save = True
save_name = 'Vol_CFTDs_' + location.upper() + '_' + date + '_storm'
# --------------------------------------------------------------------------- #
n_i = 1 * n_cols + 4
ax = plt.subplot(n_rows, n_cols, n_i)
title = 'storm envelope RHO [1]'
# save_name = 'storm_envelope_rho'
x = masked_storm['mask_rho'][masked_storm['column_type'] == 2]
y = masked_storm['mask_temp'][masked_storm['column_type'] == 2]
plot_CFAD_or_CFTD_from_2_arrays(x, y, title=title, mom_min=mom_min,
                                mom_max=mom_max, bins_mom=bins_mom,
                                vert_temp=vert_temp, temp_min=temp_min,
                                temp_max=temp_max, bins_temp=bins_temp,
                                height_min=height_min, height_max=height_max,
                                bins_height=bins_height, vmax=vmax, ax=ax,
                                save=save, save_path=save_path,
                                save_name=save_name)
# --------------------------------------------------------------------------- #
