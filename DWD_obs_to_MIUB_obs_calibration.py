#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 18.01.24                                                 #
# DWD_obs_to_MIUB_obs_calibration.py                                          #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 2: correct for rho_hv.                                                 #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/build_radar_database/correct_rhohv.py      #
# --------------------------------------------------------------------------- #
"""
@author: jgiles
Script for noise-correcting RHOHV.
# """

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils


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
LOCATIONS = ['boo', 'eis', 'fld', 'mem', 'neu', 'ros', 'tur', 'umd',
             'drs', 'ess', 'fbg', 'hnr', 'isn', 'nhb', 'oft', 'pro'
             ]
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['vol', 'pcp']
overwrite = False

import time
time.sleep(60*60*24)
# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# LOCATIONS = ['pro']
# ELEVATIONS = np.array([12])
# # MODE = ['pcp']

# date = '20210604'
# # date = '20210714'
# location = 'pro'
# elevation_deg = 12
# mode = ['vol']

for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:

                year = date[0:4]
                mon = date[4:6]
                day = date[6:8]
                sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                                           float(elevation_deg))[0][0])
                if mode == 'pcp':
                    sweep = '00'

                path_in = "/".join([header.dir_data_obs + '*' + date,
                                    year, year + '-' + mon,
                                    year + '-' + mon + '-' + day,
                                    location, mode + '*', sweep,
                                    'ras*_allmoms_*'])
                files = sorted(glob.glob(path_in))
                if not files:
                    path_in = "/".join([header.dir_data_obs +
                                        year, year + '-' + mon,
                                        year + '-' + mon + '-' + day,
                                        location, mode + '*', sweep,
                                        'ras*_allmoms_*'])
                    files = sorted(glob.glob(path_in))
                if not files:
                    print('no input data *_allmoms_*')
                    continue
                else:
                    path_in = files[0]
                    path_out = path_in.replace('_allmoms_', '_rhohv_nc_')

                if overwrite or os.path.exists(path_out):
                    print(path_out.split('/')[-1] + ' is already existing; so'
                                                    ' continue with the next'
                                                    ' file.')
                    continue

                dbzh_names = ["DBZH"]  # names to look for the DBZH variable
                rhohv_names = ["RHOHV"]  # same but for RHOHV
                data = dttree.open_datatree(path_in)[
                    'sweep_' + str(int(sweep))].to_dataset()
                # get RHOHV name
                for X_RHO in rhohv_names:
                    if X_RHO in data.data_vars:
                        break

                # get DBZH name
                for X_DBZH in dbzh_names:
                    if X_DBZH in data.data_vars:
                        break

                # check that the variables actually exist, otherwise continue
                if X_DBZH not in data.data_vars:
                    print("DBZH not found in data")
                    coninue
                if X_RHO not in data.data_vars:
                    print("RHOHV not found in data")
                    sys.exit("RHOHV not found in data.")
                    coninue

                rho_nc = utils.calculate_noise_level(
                    data[X_DBZH], data[X_RHO], noise=(-45, -15, 1))

                print('linear fit for nois levels for ' + path_out + ' ...')
                # lets do a linear fit for every noise level
                fits = []
                for nn, rhon in enumerate(rho_nc[0]):
                    merged = xr.merge(rhon)
                    rhonc_snrh = xr.DataArray(
                        merged.RHOHV_NC.values.flatten(),
                        coords={"SNRH": merged.SNRH.values.flatten()})
                    try:
                        fits.append(float(rhonc_snrh.where(
                            (0 < rhonc_snrh.SNRH) &
                            (rhonc_snrh.SNRH < 20) &
                            (rhonc_snrh > 0.7)
                        ).polyfit("SNRH", deg=1, skipna=True
                                  ).polyfit_coefficients[0].values))
                    except:
                        # if it does not work, just attach nan
                        fits.append(np.nan)

                # checking which fit has the slope closest to zero
                try:
                    bci = np.nanargmin(np.abs(np.array(fits)))
                except ValueError:
                    # if all slopes are nan, no correction good enough, abort
                    print("Could not calculate noise correction (possibly "
                          "due to not enough data points). Aborting...")
                    continue

                # get the best noise correction level according to bci
                ncl = np.arange(-45, -15, 1)[bci]
                # merge into a single array
                rho_nc_out = xr.merge(rho_nc[0][bci])
                # add noise correction level as attribute
                rho_nc_out.attrs["noise correction level"] = ncl
                # Just in case, calculate again for a NCL slightly lower (2%),
                # in case the automatically-selected one is too strong
                rho_nc2 = utils.noise_correction2(
                    data[X_DBZH], data[X_RHO], ncl * 1.02)
                rho_nc_out2 = xr.merge(rho_nc2)
                rho_nc_out2.attrs["noise correction level"] = ncl * 1.02
                rho_nc_out["RHOHV_NC"].encoding = data[X_RHO].encoding
                rho_nc_out2["RHOHV_NC"].encoding = data[X_RHO].encoding
                rho_nc_out.to_netcdf(path_out)
                rho_nc_out2.to_netcdf(path_out.replace(
                    'rhohv_nc', 'rhohv_nc_2percent'))

