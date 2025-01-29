#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# set_syn_RADAR_1_create_volume_scan.py                                       #
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

dir_data_in = header.dir_data_mod

# --------------------------------------------------------------------------- #
# start rerun 28.01.25 for all ICON volumes
overwrite_EMV = False
# overwrite_ICON = False
overwrite_ICON = '2025-01-28'
# radar_locs = list(rad_dict().keys())
radar_locs = ['ESS']
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# 11.11.24  # start rerun 28.01.25
# icon + emvorado mixtures
for day in [
    # '20210713',  # not there!
    '20210714'
]:
    da_run = 'ASS_2407'
    icon_run = 'MAIN_2405.3'
    emvorado_run = 'EMVO_00510000.2'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    spin_up_mm = 120
    # create_8_vol_nc_of_day(day=day, da_run=da_run,
    create_8_vol_nc_of_day(day=day, da_run=da_run,
                           icon_run=icon_run,
                           icon_emvorado_run=icon_emvorado_run,
                           spin_up_mm=spin_up_mm,
                           radar_locs=radar_locs,
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol,
                           overwrite_EMV=overwrite_EMV,
                           overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# 16.12.24  # start rerun 28.01.25
# icon + emvorado mixtures
for day in [
    '20210713',
    '20210714'
]:
    da_run = 'ASS_2411'
    icon_run = 'MAIN_2411.1'
    for emvorado_run in [
        'EMVO_00510000.2',
        'EMVO_00410000.2'
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        spin_up_mm = 120
        create_8_vol_nc_of_day(day=day, da_run=da_run,
                               icon_run=icon_run,
                               icon_emvorado_run=icon_emvorado_run,
                               spin_up_mm=spin_up_mm,
                               radar_locs=radar_locs,
                               dir_data_in=header.dir_data_mod,
                               dir_data_out=header.dir_data_vol,
                               overwrite_EMV=overwrite_EMV,
                               overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# 18.12.24 + 14.01.25  # start rerun 28.01.25
# icon + emvorado mixtures
for day in [
    '20210713',
    '20210714'
]:
    da_run = 'ASS_2411'
    for icon_run in [
        'MAIN_2411.6',
        'MAIN_2411.61',
    ]:
        for emvorado_run in [
            'EMVO_00510000.2',
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            spin_up_mm = 120
            create_8_vol_nc_of_day(day=day, da_run=da_run,
                                   icon_run=icon_run,
                                   icon_emvorado_run=icon_emvorado_run,
                                   spin_up_mm=spin_up_mm,
                                   radar_locs=radar_locs,
                                   dir_data_in=header.dir_data_mod,
                                   dir_data_out=header.dir_data_vol,
                                   overwrite_EMV=overwrite_EMV,
                                   overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# TODO: set up for BoXPol
# # 07.01.25
# # icon + emvorado mixtures
# for day in [
#     '20210713',
#     '20210714'
# ]:
#     da_run = 'ASS_2411'
#     icon_run = 'MAIN_2411.1'
#     icon_run = 'MAIN_2411.6'
#     radar_loc = 'BOX'
#     for emvorado_run in [
#         'EMVO_00510000.X',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         spin_up_mm = 120
#         time_start = pd.to_datetime(day, format="%Y%m%d")
#         for i in range(4):
#             time_end = time_start + pd.Timedelta('6h')
#             print('________________________________________')
#             print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
#                   str(spin_up_mm) + '_spinup/')
#             create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
#                           time_end=time_end.strftime('%Y%m%d%H'),
#                           spin_up_mm=spin_up_mm, da_run=da_run,
#                           icon_run=icon_run,
#                           icon_emvorado_run=icon_emvorado_run,
#                           dir_data_in=dir_data_in,
#                           dir_data_out=header.dir_data_vol,
#                           radar_loc=radar_loc,
#                           radar_id=rad_dict(xband_res=100)[radar_loc],
#                           include_icon=True, include_emv=False,
#                           method='Nearest')
#             print('________________________________________')
#             print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
#                   str(spin_up_mm) + '_spinup/')
#             create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
#                           time_end=time_end.strftime('%Y%m%d%H'),
#                           spin_up_mm=spin_up_mm, da_run=da_run,
#                           icon_run=icon_run,
#                           icon_emvorado_run=icon_emvorado_run,
#                           dir_data_in=dir_data_in,
#                           dir_data_out=header.dir_data_vol,
#                           radar_loc=radar_loc,
#                           radar_id=rad_dict(xband_res=125)[radar_loc],
#                           include_icon=False, include_emv=True,
#                           method='Nearest')
#             time_start = time_end


# --------------------------------------------------------------------------- #
# CONTINUE?
# import process_syn_RADAR_QVP_from_volume_scan_JSt_ahrflood_241111
