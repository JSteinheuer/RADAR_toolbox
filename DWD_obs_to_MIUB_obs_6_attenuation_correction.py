#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 27.05.24                                                 #
# DWD_obs_to_MIUB_obs_6_attenuation_correction.py                             #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 6: attenuation correction.                                             #
#         see Ryzhkov and Zrnic (2019): 6.4 (pp. 162, 167)                    #
# --------------------------------------------------------------------------- #

import datatree as dttree
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
warnings.filterwarnings("ignore")


# J. Steinheuer
def correct_zh_zdr(swp_cf, uh_tresh=0, rho_tresh=0.8, snr_tresh=15,
                   alpha=0.08, beta=0.02,
                   until_temp_beambottom=273.15+4):
    """
    ...
    """
    # 1: thresholding
    swp_mask = swp_cf.where(#(swp_cf.PHI_NC >= 0) &
                            (swp_cf.DBZH > uh_tresh) &
                            # (swp_cf.RHOHV > rho_tresh) &
                            # (swp_cf.SNRH > snr_tresh) &
                            np.isnan(swp_cf.CMAP)
                            )
    swp_mask['PHI_NC'] = xr.where(swp_mask.PHI_NC < 0, 0, swp_mask.PHI_NC)
    swp_last_l = swp_mask.where(
        (swp_mask.temp_beambottom >= until_temp_beambottom) &
        (swp_mask.temp_beambottom < until_temp_beambottom + 1))  # TODO
    phi_const = swp_last_l.PHI_NC.median(dim="range", skipna=True)
    phi_const = xr.where(np.isnan(phi_const), swp_mask.where(
        swp_mask.temp_beambottom >= until_temp_beambottom).PHI_NC.max(
        dim="range", skipna=True), phi_const)
    phi_const = xr.where(np.isnan(phi_const), swp_mask.where(
        swp_mask.temp_beambottom < until_temp_beambottom).PHI_NC.min(
        dim="range", skipna=True), phi_const)
    phi_4ac = xr.where(swp_mask.temp_beambottom >= until_temp_beambottom,
                       swp_mask.PHI_NC, phi_const)
    zh_ac = swp_mask.DBZH + phi_4ac * alpha
    zdr_ac = swp_mask.ZDR + phi_4ac * beta
    zdr_ac.attrs["long_name"] = 'reflectivity factor attenuation corrected'
    zdr_ac.attrs["short_name"] = 'ZH ac'
    zdr_ac.attrs["units"] = 'dBZ'
    swp_cf = swp_cf.assign(ZH_AC=zh_ac)
    zdr_ac.attrs["long_name"] = 'Log differential reflectivity ' + \
                                'attenuation corrected'
    zdr_ac.attrs["short_name"] = 'ZDR ac'
    zdr_ac.attrs["units"] = 'dB'
    swp_cf = swp_cf.assign(ZDR_AC=zdr_ac)
    swp_cf = swp_cf.assign(PHI_4AC=phi_4ac)
    return swp_cf


# J. Steinheuer
def attenuation_correction(date, location, elevation_deg=5.5, mode='vol',
                           overwrite=False, dir_data_obs=header.dir_data_obs):
    """
    .

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
    """
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp':
        if sweep != '00':
            return

        print('\nstart: ' + date + ' ' + location + ' pcp')
    else:
        print('\nstart: ' + date + ' ' + location + ' ' +
              str(elevation_deg))

    path_in = "/".join([dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_allmoms_*'])
    files = sorted(glob.glob(path_in))
    if not files:
        path_in = "/".join([dir_data_obs +
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep,
                            'ras*_allmoms_*'])
        files = sorted(glob.glob(path_in))
    if not files:
        print('no input data *_allmoms_*')
        return
    else:
        path_in = files[0]

    path_out = path_in.replace('_allmoms_', '_zh_zdr_ac_')
    if os.path.isfile(path_out) and not overwrite:
        print(path_out + ' exists;\n' + ' ... set: > ' +
              'overwrite = True < for recalculation')
        return

    path_kdp = path_in.replace('_allmoms_', '_kdp_nc_')
    if not os.path.exists(path_kdp):
        print('not exists: ' + path_kdp + ' -> continue')
        return

    path_rho = path_in.replace('_allmoms_', '_rhohv_nc_')
    if not os.path.exists(path_rho):
        print('not exists: ' + path_rho + ' -> continue')
        return

    path_temp = path_in.replace('_allmoms_', '_ERA5_temp_')
    if not os.path.exists(path_rho):
        print('not exists: ' + path_temp + ' -> continue')
        return

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_rho = dttree.open_datatree(path_rho)[
        'sweep_' + str(int(sweep))].to_dataset()
    data_kdp = dttree.open_datatree(path_kdp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_temp = dttree.open_datatree(path_temp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_temp2 = data_temp.interp(
        coords=data.drop(['longitude', 'latitude',
                          'altitude', 'elevation']).coords,
        method='nearest')

    data.RHOHV.values = data_rho.RHOHV_NC2P.values
    data = data.assign({'SNRH': data_rho.SNRH})
    data = data.assign({'PHI_NC': data_kdp.PHI_NC})
    data = data.assign({'temp_beambottom': data_temp2.temp_beambottom})
    data = data.assign({'temp_beamtop': data_temp.temp_beamtop})
    remo_var = list(data.data_vars.keys())
    data = data.transpose('time', 'azimuth', 'range')
    data = correct_zh_zdr(data)
    mom_use = [x for x in list(data.keys())]
    for mom in mom_use:
        data[mom].encoding["coordinates"] = \
            "time azimuth range"
    data = data.drop_vars(remo_var)
    dtree = dttree.DataTree(name="root")
    dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                    parent=dtree)
    print('saving: ... ' + path_out.split('/')[-1] + ' ...')
    dtree.load().to_netcdf(path_out)
    data.close()
    data_kdp.close()
    data_temp.close()


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
    4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,
])
MODE = [
    'pcp',
    'vol',
]
overwrite = False
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                attenuation_correction(date=date, location=location,
                                       elevation_deg=elevation_deg, mode=mode,
                                       overwrite=overwrite,
                                       dir_data_obs=header.dir_data_obs)

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
    4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,
])
MODE = [
    'pcp',
    'vol',
]
MODE = [
    'pcp',
    'vol',
]
overwrite = False
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                attenuation_correction(date=date, location=location,
                                       elevation_deg=elevation_deg,
                                       mode=mode,
                                       overwrite=overwrite,
                                       dir_data_obs=header.dir_data_obs)

# --------------------------------------------------------------------------- #
# CONTINUE?
# import DWD_obs_to_MIUB_obs_7.....
