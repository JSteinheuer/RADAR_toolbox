#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.05.25                                                 #
# set_syn_RADAR_1_create_volume_scan_LIFT.py                                  #
#                                                                             #
# Run the functions from SET_SYN_RADAR_1_CREATE_VOLUME_SCAN.py:               #
# to calculate synthetic volume scans from EMVORADO and ICON                  #
# --------------------------------------------------------------------------- #

import os
import pandas as pd
import HEADER_RADAR_toolbox as header
from SET_SYN_RADAR import \
    create_8_vol_nc_of_day, \
    create_8_vol_nc_of_day_paralell, \
    create_vol_nc, \
    rad_dict, \
    create_8_vol_nc_of_day_cdo

# --------------------------------------------------------------------------- #
# 25.05.25   # TODO
overwrite_EMV = False
overwrite_ICON = False
# radar_locs = list(rad_dict().keys())
radar_locs = ['OFT', 'FBG', 'MEM', 'TUR']
spin_up_mm = 120
for day in [
    '20230814',
    '20230817',
]:
    da_run = 'ASS_RUCorig'
    icon_run = 'MAIN_2411.3'
    for emvorado_run in [
        'EMVO_00510000.2',
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            for i in [11, 12, 13]: # range(24):
                time_start = (pd.to_datetime(day, format="%Y%m%d") +
                              pd.Timedelta('1h')*i)
                time_end = time_start + pd.Timedelta('1h')
                print('________________________________________')
                print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
                      str(spin_up_mm) + '_spinup/')
                create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
                              time_end=time_end.strftime('%Y%m%d%H'),
                              spin_up_mm=spin_up_mm, da_run=da_run,
                              icon_run=icon_run,
                              icon_emvorado_run=icon_emvorado_run,
                              dir_data_in=header.dir_data_mod,
                              dir_data_out=header.dir_data_vol,
                              radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                              overwrite=overwrite_ICON,
                              include_icon=True, include_emv=False,
                              icon_folder='ICONdata.SW')
                print('________________________________________')
                print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
                      str(spin_up_mm) + '_spinup/')
                create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
                              time_end=time_end.strftime('%Y%m%d%H'),
                              spin_up_mm=spin_up_mm, da_run=da_run,
                              icon_run=icon_run,
                              icon_emvorado_run=icon_emvorado_run,
                              dir_data_in=header.dir_data_mod,
                              dir_data_out=header.dir_data_vol,
                              radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                              overwrite=overwrite_EMV,
                              include_icon=False, include_emv=True,
                              icon_folder='ICONdata.SW')


# --------------------------------------------------------------------------- #