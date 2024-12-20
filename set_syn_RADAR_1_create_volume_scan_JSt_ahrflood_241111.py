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
import HEADER_RADAR_toolbox as header
from SET_SYN_RADAR import \
    create_8_vol_nc_of_day, \
    create_8_vol_nc_of_day_paralell, \
    create_vol_nc, \
    rad_dict

dir_data_in = header.dir_data_mod

# --------------------------------------------------------------------------- #
# 11.11.24
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
    radar_locs = list(rad_dict().keys())
    # radar_locs = ['ESS']
    # create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
    create_8_vol_nc_of_day(day=day, da_run=da_run,
                           icon_run=icon_run,
                           icon_emvorado_run=icon_emvorado_run,
                           spin_up_mm=spin_up_mm,
                           radar_locs=radar_locs,
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# 16.12.24
# icon + emvorado mixtures
for day in [
    '20210713',
    '20210714'
]:
    da_run = 'ASS_2411'
    icon_run = 'MAIN_2411.1'
    # emvorado_run = 'EMVO_00510000.2'
    # emvorado_run = 'EMVO_00410000.2'
    for emvorado_run in [
        'EMVO_00510000.2',
        'EMVO_00410000.2'
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        spin_up_mm = 120
        radar_locs = list(rad_dict().keys())
        # radar_locs = ['ESS']
        # create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
        create_8_vol_nc_of_day(day=day, da_run=da_run,
                               icon_run=icon_run,
                               icon_emvorado_run=icon_emvorado_run,
                               spin_up_mm=spin_up_mm,
                               radar_locs=radar_locs,
                               dir_data_in=header.dir_data_mod,
                               dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# 18.12.24
# icon + emvorado mixtures
for day in [
    '20210713',
    '20210714'
]:
    da_run = 'ASS_2411'
    icon_run = 'MAIN_2411.6'
    # emvorado_run = 'EMVO_00510000.2'
    for emvorado_run in [
        'EMVO_00510000.2',
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        spin_up_mm = 120
        radar_locs = list(rad_dict().keys())
        # radar_locs = ['ESS']
        # create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
        create_8_vol_nc_of_day(day=day, da_run=da_run,
                               icon_run=icon_run,
                               icon_emvorado_run=icon_emvorado_run,
                               spin_up_mm=spin_up_mm,
                               radar_locs=radar_locs,
                               dir_data_in=header.dir_data_mod,
                               dir_data_out=header.dir_data_vol)
