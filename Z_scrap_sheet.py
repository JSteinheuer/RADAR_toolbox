#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 09.01.25                                                 #
# SET_SYN_RADAR.py                                                            #
#                                                                             #
# Functions to calculate synthetic volume scans from EMVORADO and ICON.       #
# --------------------------------------------------------------------------- #

import multiprocessing as mp
import numpy as np
import glob
import pandas as pd
from pathlib import Path
import os
from osgeo import osr
import wradlib as wrl
import xarray as xr
import datetime as dt
import pyinterp
import time as time_p
from SET_SYN_RADAR import \
    create_8_vol_nc_of_day, \
    create_8_vol_nc_of_day_paralell, \
    create_vol_nc, \
    create_vol_nc_new, \
    rad_dict

import HEADER_RADAR_toolbox as header

# import HEADER_RADAR_toolbox as header
# --------------------------------------------------------------------------- #
# preamble necessary for wrl.georef.reproject: Tell the shell where to find   #
# the projection maps.                                                        #
import sys




# --------------------------------------------------------------------------- #
# specifications                                                              #
# --------------------------------------------------------------------------- #

# time
day = '20210714'
time_start_pd = pd.to_datetime(day, format="%Y%m%d") + pd.Timedelta('12h')
# time_end_pd = time_start_pd + pd.Timedelta('5min')
time_end_pd = time_start_pd + pd.Timedelta('6h')
# time_end_pd = time_start_pd + pd.Timedelta('3h')
time_start = time_start_pd.strftime('%Y%m%d%H')
time_end = time_end_pd.strftime('%Y%m%d%H')
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H")
dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='left')

# run specifications
da_run = 'ASS_2411'
icon_run = 'MAIN_2411.1'
emvorado_run = 'EMVO_00510000.2'
icon_emvorado_run = 'MAIN_2411.1/' + emvorado_run
spin_up_mm = 120
spin_up_mm = str(spin_up_mm)

# radar
radar_loc = 'ESS'
radar_id = rad_dict()[radar_loc]

# ICON, EMVORADO, new?
include_icon = True
include_emv = False

# include_icon = False
# include_emv = True
#
include_icon = True
include_emv = True

overwrite = False
overwrite='2025-01-28'

# folders
# dir_data_in = header.dir_data_mod
dir_data_in = '/automount/data02/agradar/operation_hydrometeors/data/mod/'
dir_data_out = '/automount/data02/agradar/operation_hydrometeors/data/Syn_vol/'
dir_out = dir_data_out + dti[0].strftime('%Y%m%d') + '/' + \
          da_run + '/' + icon_emvorado_run + '/' + \
          str(spin_up_mm) + 'min_spinup/'
if not include_icon:
    file_out = 'EMV_Vol_'
    if not include_emv:
        print('Nothing to do. Please include ICON or EMVORADO!')
        print('return')
elif not include_emv:
    file_out = 'ICON_Vol_'
    dir_out = dir_out.replace(icon_emvorado_run, icon_run + '/ICONdata')
else:
    file_out = 'Syn_Vol_'

file_out = file_out + radar_loc + '_' + dti[0].strftime('%Y%m%d%H%M') + \
           '_' + dti[-1].strftime('%Y%m%d%H%M') + '_' + \
           dt.datetime.now().strftime('%Y%m%d_%H%M%S') + '.nc'
if os.path.isfile(dir_out + file_out) and not overwrite:
    print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
    print(file_out + ' exists;\n' +
          ' ... set: > overwrite = True < for recalculation')
    print('___________________________')
    print('return')


# include_icon = True
# include_emv = False
# create_vol_nc(time_start=time_start, time_end=time_end,
#                    dir_data_in=header.dir_data_mod,
#                    dir_data_out=header.dir_data_vol,
#                    radar_loc=radar_loc, radar_id= rad_dict()[radar_loc], spin_up_mm=spin_up_mm,
#                    da_run=da_run, icon_run=icon_run,
#                    icon_emvorado_run=icon_emvorado_run,
#                    overwrite=overwrite, include_icon=include_icon,
#                   include_emv=include_emv)
include_icon = False
include_emv = True
create_vol_nc(time_start=time_start, time_end=time_end,
                   dir_data_in=header.dir_data_mod,
                   dir_data_out=header.dir_data_vol,
                   radar_loc=radar_loc, radar_id= rad_dict()[radar_loc], spin_up_mm=spin_up_mm,
                   da_run=da_run, icon_run=icon_run,
                   icon_emvorado_run=icon_emvorado_run,
                   overwrite=overwrite, include_icon=include_icon,
                  include_emv=include_emv)
