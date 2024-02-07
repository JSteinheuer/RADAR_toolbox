#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 05.02.24                                                 #
# DWD_obs_to_MIUB_obs_3_ERA5_temp_to_RAD.py                                   #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 3: write hourly temperature values to each radar and its grid. ERA5 3D #
#         temperature fields from download_ERA5_data_DE.py                    #
# --------------------------------------------------------------------------- #

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr
from SYN_RADAR_1_CREATE_VOLUME_SCAN import get_lon_lat_alt, ipol_fc_to_radgrid

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils

# --------------------------------------------------------------------------- #
# SET PARAMS:

DATES = ["20210604",  # case01
         "20210620", "20210621",  # case02
         "20210628", "20210629",  # case03
         "20220519", "20220520",  # case04
         "20220623", "20220624", "20220625",  # case05
         "20220626", "20220627", "20220628",  # case06+07
         "20220630", "20220701",  # case08
         "20210714",  # case09
         "20221222",  # case10
         ]
LOCATIONS = ['asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld', 'hnr', 'isn',
             'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd', ]
ELEVATIONS = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0])
MODE = ['pcp', 'vol']
overwrite = False

# --------------------------------------------------------------------------- #

# [...]

# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# LOCATIONS = ['pro']
# ELEVATIONS = np.array([5.5])
# MODE = ['vol']
# overwrite = True

# date = '20210604'
date = '20210714'
location = 'pro'
# location = 'hnr'
location = 'ess'
elevation_deg = 5.5
elevation_deg = 4.5
elevation_deg = 8
# elevation_deg = 12
# elevation_deg = 25
mode = 'vol'
overwrite = True

# for date in DATES:
#     for location in LOCATIONS:
# for elevation_deg in ELEVATIONS:
# #             for mode in MODE:

dir_data_obs = header.dir_data_obs
dir_data_era5 = header.dir_data_era5
# # J. Steinheuer
# def correct_rho_hv(date, location , elevation_deg=5.5, mode='vol',
#                    overwrite=False, dir_data_obs=header.dir_data_obs):
#     """
#     Load *all_moms*-file of one day and correct for rho_hv.
#
#     Parameter
#     ---------
#     date : 'yyyymmdd' datestring.
#     location : 'rrr' 3-letter string for radar location.
#     elevation_deg : elevation in degrees, set to 5.5 for precipitation scan
#                     (as this is the sweep 0 for the volume).
#     mode : set 'vol' for volume and 'pcp' for precipitation.
#     overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
#     dir_data_obs : directory to search for input cases
#                   (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
#     """
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
# if mode == 'pcp' and sweep != '00':
# return

path_in = "/".join([dir_data_obs + '*',
                    year, year + '-' + mon,
                    year + '-' + mon + '-' + day,
                    location, mode + '*', sweep,
                    'ras*_allmoms_*'])
files = sorted(glob.glob(path_in))
if not files:
    print('No input: ' + path_in + ' -> continue')
    # return
else:
    path_in = files[0]
    path_out = path_in.replace('_allmoms_', '_ERA5_temp_')

if not overwrite and os.path.exists(path_out):
    print('exists: ' + path_out + ' -> continue')
    # return

path_in_t = dir_data_era5 + date + '-3D-T-q-ml.nc'
path_in_z = dir_data_era5 + date + '-3D-z.nc'
if not overwrite and os.path.exists(path_out):
    print('missing ERA5: ' + path_in_t + ' -> continue')
    # return

if not overwrite and os.path.exists(path_out):
    print('missing ERA5: ' + path_in_z + ' -> continue')
    # return

# new datasets
data = dttree.open_datatree(path_in)[
    'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
data = data.drop_vars(list(data.data_vars.keys()))
data_t_era5 = dttree.open_datatree(path_in_t).to_dataset()
data_z_era5 = dttree.open_datatree(path_in_z).to_dataset()

r = data.range  # TODO /2 leads to more egg artifacts ?! hmmmmm
az = data.azimuth
el = data.elevation
sitecoords = [data.longitude.values,
              data.latitude.values,
              data.altitude.values]

lon, lat, alt = get_lon_lat_alt(r, az, el, sitecoords)

# data_t_era5 = data_t_era5.sel(lat=slice(lat.max()+.25, lat.min()-.25),
#                               lon=slice(lon.min()-.25, lon.max()+.25))
# data_z_era5 = data_z_era5.sel(lat=slice(lat.max()+.25, lat.min()-.25),
#                               lon=slice(lon.min()-.25, lon.max()+.25))

# variables
data['lat'] = (['range', 'azimuth'],
               lat.transpose(),
               dict(standard_name='latitude',
                    units='degrees_north'))
data['lon'] = (['range', 'azimuth'],
               lon.transpose(),
               dict(standard_name='longitude',
                    units='degrees_east'))
data['alt'] = (['range'],
               alt[0, :],
               dict(standard_name='altitude',
                    comments='height above mean sea level',
                    units='m'))
data['time'] = data_z_era5['time']
# dummy_tar = np.empty(shape=[data.time.size, data.azimuth.size, # changed
#                             data.range.size, ])  # changed
dummy_tra = np.empty(shape=[data.time.size, data.range.size,  # changed
                            data.azimuth.size, ])  # changed
# dummy_tar[:] = np.nan  # changed
dummy_tra[:] = np.nan  # changed
# data['temp'] = (['time', 'azimuth', 'range'], dummy_tar.copy(), # changed
data['temp'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),  # changed
                dict(standard_name='air temperature', units='K'))
# shape_ar = (data.azimuth.size, data.range.size) # changed
shape_ra = (data.range.size, data.azimuth.size)  # changed

t_i=0
# for t_i in range(data['time'].size):
# data.lon  # (range: 720, azimuth: 360)
# data.lat  # (range: 720, azimuth: 360)
# data.alt  # (range: 720)

# grid for searching later the closest ERA% cells
rad_lon = data['lon'].data.flatten()
rad_lat = data['lat'].data.flatten()
rad_alt = np.repeat(data['alt'].data[:, np.newaxis],
                    data['azimuth'].shape, axis=1).flatten()

era5_lon, era5_lat = np.meshgrid(data_z_era5.lon, data_z_era5.lat)
era5_alt = data_z_era5.z.isel(time=t_i)/9.81
era5_lon = era5_lon.flatten()
era5_lat = era5_lat.flatten()
era5_z = era5_alt.data.reshape(era5_alt.shape[0],
                               era5_alt.shape[1] * era5_alt.shape[2])

func_ipol, mask = ipol_fc_to_radgrid(
    # lon with shape (65,258775): # old
    # lon with shape (137, 2013):
    np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
    # lat with shape (65,258775)# old
    # lon with shape (137, 2013):
    np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
    # alt with shape (137, 2013):
    era5_z,
    rad_lon, rad_lat, rad_alt,  # (259200,)
    method='Linear',
)  # mask.shape (137, 2013)

data['temp'][t_i, :, :] = func_ipol(
    data_t_era5['t'].data[t_i, :, :, :].reshape(
        era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
    [mask]).reshape(shape_ra)  # changed

mom_use = [x for x in list(data.keys())]
for mom in mom_use:
    data[mom].encoding["coordinates"] = "time range azimuth"

dtree = dttree.DataTree(name="root")
dttree.DataTree(data, name=f"sweep_{int(sweep)}", parent=dtree)
print('saving: ... ' + path_out + ' ...')
dtree.load().to_netcdf(path_out)
data.close()
print('saved:  ' + path_out + ' !')
