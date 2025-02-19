#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# set_syn_RADAR_1_create_volume_scan_bringibook.py                            #
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

# --------------------------------------------------------------------------- #
# Consider rerun for all ICON volumes with new interpolation method
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = list(rad_dict().keys())
spin_up_mm = 30

# --------------------------------------------------------------------------- #
# 10.07.24
for day in [
    '20170719',
    '20170720',
    '20170724',
    '20170726',
    '20170727',
    '20180728',
    '20180923',
    '20181202',
]:
    for da_run in [
        'ASS_2111',  # ASS_old
    ]:
        for icon_run in [
            'MAIN_2203.0',  # MAIN_old
        ]:
            for emvorado_run in [
                'EMVO_00000000.2',  # EMVO_BBold
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
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
# 10.07.24
for day in [
    '20170725',
    '20170810',
    '20180809',
]:
    for da_run in [
        'ASS_2211',  # ASS_new
    ]:
        for icon_run in [
            'MAIN_2211.0',  # MAIN_new
        ]:
            for emvorado_run in [
                'EMVO_00000000.2',  # EMVO_BBold
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
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
# Consider rerun for all ICON volumes with new interpolation method
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
radar_locs = ['PRO']
spin_up_mm = 30

# --------------------------------------------------------------------------- #
# 21.11.23
for day in [
    '20170725',
    '20180809',
]:
    for da_run in [
        'ASS_2109',  # ASS_vold
    ]:
        for icon_run in [
            'MAIN_2109.0',  # MAIN_vold
        ]:
            for emvorado_run in [
                'EMVO_00000000.2',  # EMVO_BBold
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
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

