#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_syn_RADAR_pseudoRHI.py                                                 #
#                                                                             #
# plot pseudoRHI scans of moments                                             #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import xarray as xr
import datatree as dttree
import numpy as np
import matplotlib.pyplot as plt
import wradlib as wrl
import glob
import pandas as pd
import os
import matplotlib as mpl
from PLOT_SYN_RADAR import plot_syn_pseudoRHI

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20220520"
location = 'fld'
# time_i = 180
time_utc = 1455
time_utc = 1500
# time_utc = 1505
azimuth_deg = 315
azimuth_deg = 305
azimuth_deg = 303.5  # c1
azimuth_deg = 230.5  # c2
pdf_or_png = 'png'
folder_plot = header.folder_rhi_plot
da_run = 'ASS_2405'
icon_emvorado_run = 'MAIN_2405.1/EMVO_00400000.2'
spin_up_mm = '60'
# --------------------------------------------------------------------------- #
# folder and file search
year = date[0:4]
mon = date[4:6]
day = date[6:8]
hhmm_start = str(time_utc - time_utc % 600).zfill(4)
hhmm_end = str(time_utc - time_utc % 600 + 555).zfill(4)
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
path = '/'.join([header.dir_data_vol[:-1], date, da_run,
                 icon_emvorado_run, str(spin_up_mm) +
                 'min_spinup/EMV_Vol_' + location.upper() + '_' +
                 date + hhmm_start + '_' + date + hhmm_end + '.nc'])
syn_nc = xr.open_dataset(path)
# --------------------------------------------------------------------------- #

elevation_degs = np.array([
    0.5, 1.5, 2.5, 3.5, 4.5, 5.5,
    8., 12., 17., 25.
])
range_max = 200
# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# plot parameters
n_rows = 1#elevation_degs.size
n_cols = 4
n_i_zh = 0
n_i_zdr = 1*1
n_i_rho = 1*2
n_i_kdp = 1*3
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
# --------------------------------------------------------------------------- #
# loop over all elevations:
# --------------------------------------------------------------------------- #

# ------------------------------------------------------------------------#
print(azimuth_deg)
# ------------------------------------------------------------------------#
# sweep
# nc_file_comb = syn_nc.sel(azimuth=azimuth_deg, method='nearest')
# time
dti = pd.date_range(date+hhmm_start, date+hhmm_end,
                    freq="5min", inclusive='both')
time_i = np.where(dti == pd.to_datetime(date+str(time_utc).zfill(4)))[0][0]
time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
nc_file_comb = syn_nc.isel(time=time_i)
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# plot moments
moment = 'zrsim'
n_i_zh = n_i_zh + 1
ax = plt.subplot(n_rows, n_cols, n_i_zh)
title = moment + ' at ' + str(azimuth_deg) + '° ' + \
        location.upper() + ' ' + time_UTC

plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
         range_max=range_max)
# ----------------------------------------------------------------------- #
moment = 'zdrsim'
n_i_zdr = n_i_zdr + 1
ax = plt.subplot(n_rows, n_cols, n_i_zdr)
title = moment + ' at ' + str(azimuth_deg) + '° ' + \
        location.upper() + ' ' + time_UTC

plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
         range_max=range_max)
# ----------------------------------------------------------------------- #
moment = 'rhvsim'
n_i_rho = n_i_rho + 1
ax = plt.subplot(n_rows, n_cols, n_i_rho)
title = moment + ' at ' + str(azimuth_deg) + '° ' + \
        location.upper() + ' ' + time_UTC

plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
         range_max=range_max)
# ------------------------------------------------------------------------#
moment = 'kdpsim'
n_i_kdp = n_i_kdp + 1
ax = plt.subplot(n_rows, n_cols, n_i_kdp)
title = moment + ' at ' + str(azimuth_deg) + '° ' + \
        location.upper() + ' ' + time_UTC

plot_syn_pseudoRHI(nc_file_comb, ax, azimuth_deg, moment, title=title,
         range_max=range_max)

###

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
str_mod = '_'.join(path.split('/')[-5:-2]) + '_' + str(spin_up_mm) + 'min.'
str_location = '_'.join([location.upper(), date+str(time_utc).zfill(4)])
# if sum(include_sweep) == 1:
file_out = folder_plot + 'SYN_pseudoRHI_' + str(azimuth_deg) + \
           '°_' + str_location + '_' + str_mod + pdf_or_png
# else:
#     file_out = folder_plot + 'SYN_VOL_' + str(n_rows) + 'x' + \
#                str(n_cols) + '_' + str_location + '_' + str_mod + pdf_or_png

plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(file_out, format=pdf_or_png, transparent=True)
plt.close()
