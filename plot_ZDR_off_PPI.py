#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 19.11.24                                                 #
# >script_name.py<                                                            #
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
# case (adjust!):
date = "20210714"
location = 'ess'
elevation_deg = 12
mode = 'vol'  # 'pcp'
overwrite = False  # TODO
dir_data_obs = header.dir_data_obs
other_zdr_off_day = ''  # TODO
n_zdr_lowest = 1000
pdf_or_png = 'png'
folder_plot = header.folder_ppi_plot

file_in = 'ras07-vol5minng01_sweeph5onem_zdr_off_07-202107140003-202107142358-ess-10410.hd5'
# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
n_rows = 1
n_cols = 1
n_i = 0
# --------------------------------------------------------------------------- #
year = date[0:4]
mon = date[4:6]
day = date[6:8]
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
if mode == 'pcp':
    if sweep != '00':
        print('return')

path_in = "/".join([dir_data_obs + '*',
                    year, year + '-' + mon,
                    year + '-' + mon + '-' + day,
                    location, mode + '*', sweep,
                    'ras*_zdr_off_*'])
files = sorted(glob.glob(path_in))
if not files:
    path_in = "/".join([dir_data_obs +
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_zdr_off_*'])
    files = sorted(glob.glob(path_in))
if not files:
    print('return')
else:
    path_in = files[0]

if other_zdr_off_day:
    year_of = other_zdr_off_day[0:4]
    mon_of = other_zdr_off_day[4:6]
    day_of = other_zdr_off_day[6:8]
else:
    year_of = date[0:4]
    mon_of = date[4:6]
    day_of = date[6:8]

# time_start = date + file_in.split('-')[-4][-4:]
# time_end = date + file_in.split('-')[-3][-4:]
# dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
# dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
# dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
folder_in = "/".join([dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])


nc_files = glob.glob(folder_in + '/' + file_in)
nc_file = nc_files[0]

sweep = nc_file.split('/')[-2]
vol = dttree.open_datatree(nc_file)[
    'sweep_' + str(int(sweep))].to_dataset().chunk('-1')

path_in_bb = "/".join([dir_data_obs + '*',
                       year_of, year_of + '-' + mon_of,
                       year_of + '-' + mon_of + '-' + day_of,
                       location, '90grad*', '00',
                       'ras*_zdr_off_*'])
files = sorted(glob.glob(path_in_bb))
if not files:
    path_in_bb = "/".join([dir_data_obs +
                           year_of, year_of + '-' + mon_of,
                           year_of + '-' + mon_of + '-' + day_of,
                           location, '90grad*', '00',
                           'ras*_zdr_off_*'])
    files = sorted(glob.glob(path_in_bb))
if not files:
    print('no birdbath input data *90grad*_zdr_off_*')
else:
    path_in_bb = files[0]
    data_bb = dttree.open_datatree(path_in_bb)[
        'sweep_0'].to_dataset().chunk(-1)

# ppi = vol.isel(time=time_i)
# ppi = ppi.transpose('azimuth', 'range')

zdr_off_lr_ppi=vol.zdr_off_lr_ppi
zdr_off_lr_n_ppi=vol.zdr_off_lr_n_ppi
zdr_off_sd_ppi=vol.zdr_off_sd_ppi
zdr_off_sd_n_ppi=vol.zdr_off_sd_n_ppi


fig, ax1 = plt.subplots(figsize=(13, 6))
zdr_off_sd_n_ppi.plot(color='green', alpha=.2,#, ls='dotted',
                      label='count small drops')
zdr_off_lr_n_ppi.plot(color='blue', alpha=.2,#, ls='dotted',
                      label='count light rain')
vol.zdr_off_das_n_ppi.plot(color='orange', alpha=.2,#, ls='dotted',
                      label='count dry snow')
ax1.hlines(n_zdr_lowest, color='grey', ls='dotted', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1])

zdr_off_sd_ppi_filtered=zdr_off_sd_ppi.where(zdr_off_sd_n_ppi>n_zdr_lowest)
zdr_off_lr_ppi_filtered=zdr_off_lr_ppi.where(zdr_off_lr_n_ppi>n_zdr_lowest)
zdr_off_combined_ppi=zdr_off_sd_ppi_filtered.where(
    zdr_off_sd_n_ppi>zdr_off_lr_n_ppi, zdr_off_lr_ppi_filtered)
of_all = np.nansum(vol.zdr_off_sd_n_ppi * vol.zdr_off_sd_ppi +
                   vol.zdr_off_lr_n_ppi * vol.zdr_off_lr_ppi) / \
                     np.nansum(vol.zdr_off_sd_n_ppi + vol.zdr_off_lr_n_ppi)

zdr_off_combined_ppi=zdr_off_combined_ppi.where(
    np.isnan(zdr_off_combined_ppi)== False, of_all)
ax2 = ax1.twinx()
ax2.hlines(0, color='grey', ls='dotted', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1])
ax2.hlines(data_bb.zdr_off_bb.data, color='red', ls='dashed', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1], label = 'bird bath constant')
ax2.hlines(vol.zdr_off_sd.data, color='green', ls='dashed', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1], label = 'small drops constant')
ax2.hlines(vol.zdr_off_lr.data, color='blue', ls='dashed', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1], label = 'small drops constant')
ax2.hlines(vol.zdr_off_das.data, color='orange', ls='dashed', xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[1], label = 'dry snow constant')
zdr_off_sd_ppi_filtered.plot(color='green', linewidth=5, alpha=.4,
                             label='small drops offset')
zdr_off_lr_ppi_filtered.plot(color='blue', linewidth=5, alpha=.4,
                             label='light rain offset')


# das down
vol['zdr_off_das_ppi']= vol.zdr_off_das_ppi.where(vol.zdr_off_das_n_ppi>n_zdr_lowest, vol.zdr_off_das.data)
vol['zdr_off_das_ppi']=vol.zdr_off_das_ppi-vol.zdr_off_das+vol.zdr_off_lr

vol.zdr_off_das_ppi.plot(color='orange', linewidth=5, alpha=.4,
                             label='dry snow offset')
zdr_off_combined_ppi.plot(color='black', linewidth=2,
                          label='combined offset')
ax2.set_ylim(-.7, .7)
ax2.set_xlim(18822, 18823)
ax2.legend(loc='upper right')
ax1.legend(loc='upper left')

ax1.set_title('')
ax2.set_title('')
# ax2.set_title('ZDR offset')
ax1.set_xlabel("UTC [mm-dd hh]")
ax1.set_ylabel("counts [1]")
ax2.set_ylabel("offset [dB]")
plt.tight_layout()
plt.savefig(folder_plot + 'ZDR_offset_' + file_in.replace('.hd5', '.pdf'),
            format='pdf')
