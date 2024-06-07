#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 05.02.24                                                 #
# DWD_obs_to_MIUB_obs_3_ERA5_temp.py                                          #
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
import wradlib as wrl

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
# from radar_processing_scripts import utils


# --------------------------------------------------------------------------- #
# FUNCTIONS: J. Steinheuer


# J. Steinheuer
def era5_temp(date, location, elevation_deg=5.5, mode='vol',
              overwrite=False, dir_data_obs=header.dir_data_obs,
              dir_data_era5=header.dir_data_era5):
    """
    Load one *all_moms*-file (of one day and sweep) and create on the
    radar grid all temperature values taken from ERA5.

    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    elevation_deg : elevation in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    mode : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
    dir_data_obs : directory to search for input cases
                  (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
    dir_data_era5 : directory to search for era5 files.
    """
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        return

    path_in = "/".join([dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_allmoms_*'])
    files = sorted(glob.glob(path_in))
    if not files:
        print('No input: ' + path_in + ' -> continue')
        return
    else:
        path_in = files[0]
        path_out = path_in.replace('_allmoms_', '_ERA5_temp_')

    if not overwrite and os.path.exists(path_out):
        print('exists: ' + path_out + ' -> continue')
        return

    path_in_any = '/'.join(path_in.split('/')[:-1]) + '/*any*'
    files_any = sorted(glob.glob(path_in_any))
    if not files_any:
        path_in_any = "/".join([header.dir_data_obs_realpep + '' +
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep, 'ras*'])
        files_any = sorted(glob.glob(path_in_any))
        if not files_any:
            print('nothing found -> bw=1')
            bw = 1
        else:
            path_in_any = files_any[0]
            bw = dttree.open_datatree(path_in_any)['how'].attrs['beamwidth']

    else:
        path_in_any = files_any[0]
        bw = dttree.open_datatree(path_in_any)['how'].attrs['beamwidth']

    bw = round(bw, 3)
    print(date + ' ' + location + ' bw=' + str(bw))
    path_in_t = dir_data_era5 + date + '-3D-T-q-ml.nc'
    path_in_z = dir_data_era5 + date + '-3D-z.nc'
    if not os.path.exists(path_in_t):
        print('missing ERA5: ' + path_in_t + ' -> continue')
        return

    if not os.path.exists(path_in_z):
        print('missing ERA5: ' + path_in_z + ' -> continue')
        return

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data = data.drop_vars(list(data.data_vars.keys()))
    data_t_era5 = dttree.open_datatree(path_in_t).to_dataset()
    data_z_era5 = dttree.open_datatree(path_in_z).to_dataset()
    r = data.range
    az = data.azimuth
    el = data.elevation
    sitecoords = [data.longitude.values,
                  data.latitude.values,
                  data.altitude.values]
    lon, lat, alt = get_lon_lat_alt(r, az, el, sitecoords)
    # variables
    data['lat'] = (['range', 'azimuth'],
                   lat.transpose(),
                   dict(standard_name='latitude',
                        units='degrees_north'))
    data['lon'] = (['range', 'azimuth'],
                   lon.transpose(),
                   dict(standard_name='longitude',
                        units='degrees_east'))
    if mode == 'pcp':
        data['alt'] = (['range', 'azimuth'],
                       alt.transpose(),
                       dict(standard_name='altitude',
                            comments='height above mean sea level',
                            units='m'))
    else:
        data['alt'] = (['range'],
                       alt[0, :],
                       dict(standard_name='altitude',
                            comments='height above mean sea level',
                            units='m'))

    data['time'] = data_t_era5['time']
    dummy_tra = np.empty(shape=[data.time.size, data.range.size,
                                data.azimuth.size, ])
    dummy_tra[:] = np.nan
    data['temp'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                    dict(standard_name='air temperature', units='K'))
    data['temp_beamtop'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                    dict(standard_name='air temperature', units='K',
                         comment='beambroadening induced temperature of '
                                 'bin volume top (approximated for pw=0)'))
    data['temp_beambottom'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                    dict(standard_name='air temperature', units='K',
                         comment='beambroadening induced temperature of '
                                 'bin volume bottom (approximated for pw=0)'))
    shape_ra = (data.range.size, data.azimuth.size)
    for t_i in range(data['time'].size):
        # grid for searching later the closest ERA5 cells
        rad_lon = data['lon'].data.flatten()
        rad_lat = data['lat'].data.flatten()
        if mode == 'pcp':
            rad_alt = data['alt'].data.flatten()
        else:
            rad_alt = np.repeat(data['alt'].data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        era5_lon, era5_lat = np.meshgrid(data_z_era5.lon, data_z_era5.lat)
        era5_alt = data_z_era5.z.isel(time=t_i)/9.81
        era5_lon = era5_lon.flatten()
        era5_lat = era5_lat.flatten()
        era5_z = era5_alt.data.reshape(era5_alt.shape[0],
                                       era5_alt.shape[1] * era5_alt.shape[2])
        # temp center
        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
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
            [mask]).reshape(shape_ra)

        # temp top
        beamradius = wrl.util.half_power_radius(data.range, bw)
        if mode == 'pcp':
            rad_alt = (data['alt']+beamradius).data.flatten()
        else:
            rad_alt = np.repeat((data['alt']+beamradius).data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
            # lon with shape (137, 2013):
            np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
            # alt with shape (137, 2013):
            era5_z,
            rad_lon, rad_lat, rad_alt,  # (259200,)
            method='Linear',
        )  # mask.shape (137, 2013)
        data['temp_beamtop'][t_i, :, :] = func_ipol(
            data_t_era5['t'].data[t_i, :, :, :].reshape(
                era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
            [mask]).reshape(shape_ra)

        # temp bottom
        if mode == 'pcp':
            rad_alt = (data['alt']+beamradius).data.flatten()
        else:
            rad_alt = np.repeat((data['alt']-beamradius).data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
            # lon with shape (137, 2013):
            np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
            # alt with shape (137, 2013):
            era5_z,
            rad_lon, rad_lat, rad_alt,  # (259200,)
            method='Linear',
        )  # mask.shape (137, 2013)
        data['temp_beambottom'][t_i, :, :] = func_ipol(
            data_t_era5['t'].data[t_i, :, :, :].reshape(
                era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
            [mask]).reshape(shape_ra)

    mom_use = [x for x in list(data.keys())]
    for mom in mom_use:
        data[mom].encoding["coordinates"] = "time range azimuth"

    dtree = dttree.DataTree(name="root")
    dttree.DataTree(data, name=f"sweep_{int(sweep)}", parent=dtree)
    print('saving: ... ' + path_out + ' ...')
    dtree.load().to_netcdf(path_out)
    data.close()
    print('saved:  ' + path_out + ' !')


# --------------------------------------------------------------------------- #
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20210604",  # case01
    "20210620", "20210621",  # case02
    "20210628", "20210629",  # case03
    "20220519", "20220520",  # case04
    "20220623", "20220624", "20220625",  # case05
    "20220626", "20220627", "20220628",  # case06+07
    "20220630", "20220701",  # case08
    "20210714",  # case09
    "20221222",  # case10
]
LOCATIONS = [
    'asb', 'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,
    # 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,  # TODO
])
MODE = [
    'pcp',
    # 'vol',  # TODO
    # '90grad',  # TODO
]
# overwrite = False  # TODO
overwrite = True  # TODO
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                era5_temp(date=date, location=location,
                          elevation_deg=elevation_deg,
                          mode=mode, overwrite=overwrite)

# --------------------------------------------------------------------------- #
# OLD CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20170719",
    "20170725",
]
LOCATIONS = [
    'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,
    # 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,  # TODO
])
MODE = [
    'pcp',
    # 'vol',  # TODO
    # '90grad',  # TODO
]
# overwrite = False  # TODO
overwrite = True  # TODO
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                era5_temp(date=date, location=location,
                          elevation_deg=elevation_deg,
                          mode=mode, overwrite=overwrite)

# --------------------------------------------------------------------------- #
# CONTINUE?
# import DWD_obs_to_MIUB_obs_4_correct_phi_kdp
import DWD_obs_to_MIUB_obs_5_calibrate_zdr  # TODO
