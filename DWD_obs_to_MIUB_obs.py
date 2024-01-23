#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 10.01.24                                                 #
# DWD_obs_to_MIUB_obs.py                                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 1: Load all moments and all times (of day) in one file.                #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/build_radar_database/concat_dwd_data_to_d* #
# --------------------------------------------------------------------------- #
"""
@author: jgiles
This script takes all dwd radar files from a folder (for one elevation) and
merges them into a single file combining all moments along all timesteps.
Then saves the resulting dataset into a new file with the same naming
style but with "allmoms" instead of the moment name. Additionally, it saves
either a true.txt or false.txt file alongside, if the data fulfills certain
condition, as an attempt to check if there is actually something interesting
in that period of data.
"""

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
from pathlib import Path
import os

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
LOCATIONS = ['asb', 'boo', 'eis', 'fld', 'mem', 'neu', 'ros', 'tur', 'umd',
             'drs', 'ess', 'fbg', 'hnr', 'isn', 'nhb', 'oft', 'pro'
             ]
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']
# overwrite = True
overwrite = False

# START: Loop over cases, dates, and radars:

# DATES = ['20210604']
# # DATES = ['20210714']
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
                if mode == 'pcp' and sweep != '00':
                    continue

                path_in = "/".join([header.dir_data_obs + '*',
                                    year, year + '-' + mon,
                                    year + '-' + mon + '-' + day,
                                    location, mode + '*', sweep, 'ras*'])
                files = sorted(glob.glob(path_in))
                if not files:
                    path_in = "/".join([header.dir_data_obs_realpep + '' +
                                        year, year + '-' + mon,
                                        year + '-' + mon + '-' + day,
                                        location, mode + '*', sweep, 'ras*'])
                    files = sorted(glob.glob(path_in))
                    path_out = '/'.join((files[0].split('/'))[:-1])
                    path_out = path_out.replace(header.dir_data_obs_realpep,
                                                header.dir_data_obs)
                else:
                    path_out = '/'.join((files[0].split('/'))[:-1])
                files_temp = []
                for file in files:
                    if not 'allmoms' in file:
                        if not 'rhohv_nc' in file:
                            files_temp.append(file)
                    # else:
                    #     print(file)

                files = files_temp
                name = files[0].split("/")[-1].split("_")
                t_start = files[0].split("/")[-1].split("-")[2][:12]
                t_end = files[-1].split("/")[-1].split("-")[2][:12]
                name[-2] = "allmoms"
                name_m1 = name[-1]
                name_m1 = name_m1.replace('-hd5', '.hd5')
                name_m1 = name_m1.split("-")
                name_m1[1] = t_start + '-' + t_end
                name[-1] = '-'.join(name_m1)
                name_out = ("_".join(name))
                file_out = '/'.join([path_out, name_out])
                Path(path_out).mkdir(parents=True, exist_ok=True)
                if not overwrite and os.path.exists(file_out):
                    print(name_out + ' is already existing;'
                                     ' so continue with the next file.')
                    continue

                ds = utils.load_dwd_raw(files)
                dtree = dttree.DataTree(name="root")
                if "longitude" in ds:
                    # rename the variable
                    # ds = ds.rename({"longitude": "longitude_loc"})
                    ds['longitude_loc'] = ([], ds['longitude'].data,
                                           ds['longitude'].attrs)
                    ds = ds.drop_vars({'longitude'})
                if "latitude" in ds:
                    # rename the variable
                    ds['latitude_loc'] = ([], ds['latitude'].data,
                                          ds['latitude'].attrs)
                    ds = ds.drop_vars({'latitude'})
                if "altitude" in ds:
                    # rename the variable
                    ds['altitude_loc'] = ([], ds['altitude'].data,
                                          ds['altitude'].attrs)
                    ds = ds.drop_vars({'altitude'})
                if "fixed_angle" in ds:  # rename the variable
                    ds = ds.rename({"fixed_angle": "sweep_fixed_angle"})

                dttree.DataTree(ds, name=f"sweep_{int(sweep)}", parent=dtree)
                print('saving ' + file_out + ' ...')
                dtree.load().to_netcdf(file_out)

# import os
#
# folder = "/automount/agradar/operation_hydrometeors/data/obs/"
# files = sorted(glob.glob(folder + '*/*/*/*/*/*/*/*'))
# for file in files:
#     sweep_folder = file.split('/')[-2]
#     sweep_file = file.split('/')[-1].split('_')[-1][:2]
#     if 'allmoms' in file:
#         print(file)
#         file_new = '/'.join(file.split('/')[:-2]) + '/' + sweep_file + \
#                    '/' + file.split('/')[-1]
#         # print(file_new)
#         # os.system('rm ' + file )

# import DWD_obs_to_MIUB_obs_calibration

# sweep=4
# path_in="/automount/agradar/operation_hydrometeors/data/obs/" \
#         "OpHymet2-case10-20221222/2022/2022-12/2022-12-22/boo/" \
#         "vol5minng10/04/" \
#         "ras11-vol5minng10_sweeph5allm_allmoms_04-202212220002-202212221742-boo-10132.hd5"
# ddata = dttree.open_datatree(path_in)[
#                     'sweep_' + str(int(sweep))].to_dataset().chunk('auto')


# sweep=4
# path_in="/automount/agradar/operation_hydrometeors/data/obs/" \
#         "OpHymet2-case10-20221222/2022/2022-12/2022-12-22/boo/" \
#         "vol5minng10/04/" \
#         "ras11-vol5minng10_sweeph5allm_allmoms_04-202212220002-202212221742-boo-10132.hd5"
# ddata = dttree.open_datatree(path_in)[
#                     'sweep_' + str(int(sweep))].to_dataset().chunk('auto')

# folder = "/automount/agradar/operation_hydrometeors/data/obs/*"
# files = sorted(glob.glob(folder + '*/*/*/*/*/*/*/*'))
# for file in files:
#     if 'rhohv_nc' in file:
#         print(file)
#         os.system('rm ' + file )
