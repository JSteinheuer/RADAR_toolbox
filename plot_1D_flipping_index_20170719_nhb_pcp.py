#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_PPI.py                                                            #
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
from PLOT_PPI import plot_PPI

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# case (adjust!):
date = "20170719"
location = 'nhb'
elevation_deg = 5.5
mode = 'pcp'
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot
file_in_1 = 'ras07-pcpng01_sweeph5onem_kdp_nc_00-201707190000-201707192355-nhb-10605.hd5'
# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
n_rows = 1
n_cols = 1
n_i = 0
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(6 * n_cols, 4 * n_rows))
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
time_start = date + file_in_1.split('-')[-4][-4:]
time_end = date + file_in_1.split('-')[-3][-4:]
folder_in = "/".join([header.dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
nc_files_1 = glob.glob(folder_in + '/' + file_in_1)
if len(nc_files_1) > 1:
    print('more files')
elif len(nc_files_1) == 0:
    print('none files')

nc_file_1 = nc_files_1[0]
# --------------------------------------------------------------------------- #
moment = 'PHI_FLIPPED_INDEX'
n_i = n_i + 1
ax = plt.subplot(n_rows, n_cols, n_i)
sweep = nc_file_1.split('/')[-2]
vol = dttree.open_datatree(nc_file_1)[
    'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
if ax is None:
    plt.figure(figsize=(5, 4))
    ax = plt.gca()

vol[moment].plot(ax=ax, ylim=[-100, 12000], color='blue')
ax.axhline(y=1000, color='orange')
ax.axhline(y=0, color='grey')
ax.set_ylabel("$\Delta \Phi$ [°]")
ax.set_xlabel("UTC [mm-dd hh]")
if mode == 'vol':
    title = moment + ' at ' + str(elevation_deg) + '° ' + \
            location.upper() + ' ' + date
if mode == 'pcp':
    title = moment + ' for precip. scan at ' + \
            location.upper() + ' ' + date

plt.title(title)
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #
file_out = (file_in_1.replace('.hd5', '_').replace(' ', '_')
            + '.pdf').replace('kdp_nc_', '')
plt.tight_layout()
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)

plt.savefig(folder_plot + moment + '_' + file_out, format=pdf_or_png,
            transparent=True)
plt.close()
# --------------------------------------------------------------------------- #
