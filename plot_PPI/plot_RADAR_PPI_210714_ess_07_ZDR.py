#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 18.11.24                                                 #
# plot_RADAR_PPI.py                                                           #
#                                                                             #
# [...] Description here [...]                                                #
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
from PLOT_RADAR import plot_PPI

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20210714"
location = 'ess'
elevation_deg = 12
mode = 'vol'  # 'pcp'
time_i = 210
time_i = 109
time_i = 110
time_i2 = time_i+1
time_i3 = time_i+2
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
# PPI 1 row: ZDR BB 211, ZDR BB 212, ZDR BB 213
file_in_1 = 'ras07-vol5minng01_sweeph5onem_polmoms_nc_07-202107140003-202107142358-ess-10410_backup.hd5'
# PPI 2 row: ZDR new 211, ZDR new 212, ZDR new 213
file_in_2 = 'ras07-vol5minng01_sweeph5onem_polmoms_nc_07-202107140003-202107142358-ess-10410.hd5'
# PPI 2 row: ZDR raw 211, ZDR raw 212, ZDR raw 213
file_in_3 = 'ras07-vol5minng01_sweeph5onem_allmoms_07-202107140003-202107142358-ess-10410.hd5'
# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# a                                                                           #
# --------------------------------------------------------------------------- #
n_rows = 2
n_cols = 3
n_i = 0
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
time_start = date + file_in_1.split('-')[-4][-4:]
time_end = date + file_in_1.split('-')[-3][-4:]
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
time_UTC2 = dti[time_i2].strftime('%Y-%m-%d %H:%M')
time_UTC3 = dti[time_i3].strftime('%Y-%m-%d %H:%M')
folder_in = "/".join([header.dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
nc_files_1 = glob.glob(folder_in + '/' + file_in_1)
nc_file_1 = nc_files_1[0]
if len(nc_files_1) > 1:
    print('more files')
elif len(nc_files_1) == 0:
    print('none files')

nc_files_2 = glob.glob(folder_in + '/' + file_in_2)
nc_file_2 = nc_files_2[0]
if len(nc_files_2) > 1:
    print('more files')
elif len(nc_files_2) == 0:
    print('none files')
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_1, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_1, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_2.replace('.hd5', '_ZDR_' + time_UTC2.replace(' ', '_') +
                             '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# b                                                                           #
# --------------------------------------------------------------------------- #
n_i = 0
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
time_start = date + file_in_1.split('-')[-4][-4:]
time_end = date + file_in_1.split('-')[-3][-4:]
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
time_UTC2 = dti[time_i2].strftime('%Y-%m-%d %H:%M')
time_UTC3 = dti[time_i3].strftime('%Y-%m-%d %H:%M')
folder_in = "/".join([header.dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
nc_files_1 = glob.glob(folder_in + '/' + file_in_1)
nc_file_1 = nc_files_1[0]
if len(nc_files_1) > 1:
    print('more files')
elif len(nc_files_1) == 0:
    print('none files')

nc_files_3 = glob.glob(folder_in + '/' + file_in_3)
nc_file_3 = nc_files_3[0]
if len(nc_files_3) > 1:
    print('more files')
elif len(nc_files_3) == 0:
    print('none files')
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_1, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_1, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'ZDR'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_3, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_3, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_3, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_3.replace('.hd5', '_ZDR_b_' + time_UTC2.replace(' ', '_') +
                             '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# c                                                                           #
# --------------------------------------------------------------------------- #
n_rows = 3
n_cols = 5
n_i = 0
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(5 * n_cols, 4 * n_rows))
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
time_start = date + file_in_1.split('-')[-4][-4:]
time_end = date + file_in_1.split('-')[-3][-4:]
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
time_UTC = dti[time_i].strftime('%Y-%m-%d %H:%M')
time_UTC2 = dti[time_i2].strftime('%Y-%m-%d %H:%M')
time_UTC3 = dti[time_i3].strftime('%Y-%m-%d %H:%M')
folder_in = "/".join([header.dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
nc_files_2 = glob.glob(folder_in + '/' + file_in_2)
nc_file_2 = nc_files_2[0]
if len(nc_files_2) > 1:
    print('more files')
elif len(nc_files_2) == 0:
    print('none files')
# --------------------------------------------------------------------------- #
moment = 'ZH_AC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'KDP_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'PHI_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'ZH_AC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'KDP_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'PHI_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC2
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC2
plot_PPI(nc_file_2, ax, time_i2, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'ZH_AC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'ZDR_AC_OC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'KDP_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
moment = 'PHI_NC'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC3
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC3
plot_PPI(nc_file_2, ax, time_i3, moment, title=title, range_max=40)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_2.replace('.hd5', '_ZDR_c_' + time_UTC2.replace(' ', '_') +
                             '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #

