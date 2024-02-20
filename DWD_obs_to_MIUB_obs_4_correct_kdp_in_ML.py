#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.01.24                                                 #
# DWD_obs_to_MIUB_obs_4_correct_kdp_in_ML.py                                  #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 4: correct KDP in ML.                                                  #
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
LOCATIONS = ['asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld',  'hnr', 'isn',
             'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd',
             ]
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']
# overwrite = True
overwrite = False

# import time
# time.sleep(60*60*24)

# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# LOCATIONS = ['pro']
# ELEVATIONS = np.array([12])
# # MODE = ['pcp']

date = '20210604'
# date = '20210714'
location = 'pro'
elevation_deg = 5.5
mode = 'vol'

# for date in DATES:
#     for location in LOCATIONS:
#         for mode in MODE:
#             for elevation_deg in ELEVATIONS:

year = date[0:4]
mon = date[4:6]
day = date[6:8]
sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
if mode == 'pcp' and sweep != '00':
    print('pcp only 00')
    # continue

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
    # continue
else:
    path_in = files[0]
    path_out = path_in.replace('_allmoms_', '_rhohv_nc_')

