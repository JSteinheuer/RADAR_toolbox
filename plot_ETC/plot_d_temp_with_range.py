#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 21.02.24                                                 #
# plot_d_temp_with_range.py                                                   #
#                                                                             #
# [...] Description here [...]                                                #
# --------------------------------------------------------------------------- #

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr
import wradlib as wrl
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# SET PARAMS:

ELEVATIONS = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0])
# MODE = ['pcp', 'vol']
date = '20210604'
hh = 12
azim = 172.5
location = 'ess'
mode = 'vol'
overwrite = True
# --------------------------------------------------------------------------- #
# plt.figure(figsize=(14, 9))
plt.figure(figsize=(12, 7))
# elevation_deg = 0.5
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
for elevation_deg in ELEVATIONS:

    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        print('return')

    path_in = "/".join([header.dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*ERA5_temp*'])
    files = sorted(glob.glob(path_in))
    if not files:
        print('No input: ' + path_in + ' -> continue')
        print('return')
    else:
        path_in = files[0]
        path_out = path_in.replace('_allmoms_', '_ERA5_temp_')

    if not overwrite and os.path.exists(path_out):
        print('exists: ' + path_out + ' -> continue')
        print('return')

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')

    dtemp = data.temp_beambottom - data.temp_beamtop

    r = data.range.values / 1000
    dt = dtemp.isel(time=hh - 1).sel(azimuth=azim).values
    alt = data.alt.values / 1000

    ax = plt.subplot(2, 1, 1)
    plt.plot(r, dt, label=str(elevation_deg) + '°')

    ax = plt.subplot(2, 1, 2),
    plt.plot(dt, alt, label=str(elevation_deg) + '°')

ax = plt.subplot(2, 1, 1)
plt.xlabel('range [km]')
plt.ylabel('$\Delta$T')
plt.legend()

ax = plt.subplot(2, 1, 2)
plt.ylabel('height [km]')
plt.xlabel('$\Delta$T')
plt.legend()

plt.tight_layout()
plt.savefig(header.folder_plot + 'd_temp_ ' + location + '_' + date + '-' +
            str(hh) + '_az-' + str(int(round(azim, 0))) + '.pdf',
            format='pdf', transparent=True)
