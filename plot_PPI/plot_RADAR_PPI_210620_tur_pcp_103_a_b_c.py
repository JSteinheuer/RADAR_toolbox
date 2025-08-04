#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
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
# case (adjust!)
date = "20210620"
location = 'tur'
elevation_deg = 5.5
mode = 'pcp'
time_i = 103
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
# PPI 1 row: DBZH, ZDR, RHOHV
file_in_1 = 'ras11-pcpng10_sweeph5allm_allmoms_00-202106201400-202106202355-tur-10832.hd5'
# PPI 2 row: phi_c, PHI_NC, KDP
file_in_2 = 'ras11-pcpng10_sweeph5allm_kdp_nc_00-202106201400-202106202355-tur-10832.hd5'
# RHOHV_NC2P
file_in_3 = 'ras11-pcpng10_sweeph5allm_rhohv_nc_00-202106201400-202106202355-tur-10832.hd5'
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

nc_files_3 = glob.glob(folder_in + '/' + file_in_3)
if len(nc_files_3) > 1:
    print('more files')
elif len(nc_files_3) == 0:
    print('none files')

nc_file_3 = nc_files_3[0]
# --------------------------------------------------------------------------- #
moment = 'DBZH'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_1, ax, time_i, moment, title=title)
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
plot_PPI(nc_file_1, ax, time_i, moment, title=title)
# --------------------------------------------------------------------------- #
moment = 'RHOHV'
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_3, ax, time_i, moment, title=title)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'UPHIDP'
moment = 'PHI_C'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_2.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
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
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20210620"
location = 'tur'
elevation_deg = 5.5
mode = 'pcp'
time_i = 103
range_max = 50
# pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
# PPI 1 row: DBZH, ZDR, RHOHV
file_in_1 = 'ras11-pcpng10_sweeph5allm_allmoms_00-202106201400-202106202355-tur-10832.hd5'
# PPI 2 row: phi_c, PHI_NC, KDP
file_in_2 = 'ras11-pcpng10_sweeph5allm_kdp_nc_00-202106201400-202106202355-tur-10832.hd5'
# RHOHV_NC2P
file_in_3 = 'ras11-pcpng10_sweeph5allm_rhohv_nc_00-202106201400-202106202355-tur-10832.hd5'
# case (adjust!):
# --------------------------------------------------------------------------- #
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
moment = 'DBZH'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=range_max)
# --------------------------------------------------------------------------- #
moment = 'RHOHV'
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_3, ax, time_i, moment, title=title, range_max=range_max)
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'UPHIDP'
moment = 'PHI_C'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_2.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
                             '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_50km_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# c                                                                          #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20210620"
location = 'tur'
elevation_deg = 5.5
mode = 'pcp'
time_i = 103
range_max = 50
# pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
# PPI 1 row: DBZH, ZDR, RHOHV
file_in_1 = 'ras11-pcpng10_sweeph5allm_allmoms_00-202106201400-202106202355-tur-10832.hd5'
# PPI 2 row: phi_c, PHI_NC, KDP
file_in_2 = 'ras11-pcpng10_sweeph5allm_kdp_nc_00-202106201400-202106202355-tur-10832.hd5'
# RHOHV_NC2P
file_in_3 = 'ras11-pcpng10_sweeph5allm_rhohv_nc_00-202106201400-202106202355-tur-10832.hd5'
# case (adjust!):
# --------------------------------------------------------------------------- #
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
moment = 'DBZH'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_1, ax, time_i, moment, title=title, range_max=range_max)
# --------------------------------------------------------------------------- #
moment = 'RHOHV'
moment = 'RHOHV_NC2P'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_3, ax, time_i, moment, title=title, range_max=range_max,
         levels=[.7, .95, .998],
         norm=mpl.colors.BoundaryNorm([.7, .95, .998], 2), cmap='jet'
         )
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
moment = 'UPHIDP'
moment = 'PHI_C'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + time_UTC
else:
    title = moment + ' at pcp ' + \
            location.upper() + ' ' + time_UTC
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
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
plot_PPI(nc_file_2, ax, time_i, moment, title=title, range_max=range_max)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = file_in_2.replace('.hd5', '_' + time_UTC.replace(' ', '_') +
                             '.' + pdf_or_png)
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + 'PPI_50km_rho_thr_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #



