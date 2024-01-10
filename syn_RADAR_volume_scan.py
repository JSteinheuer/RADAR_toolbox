#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# syn_RADAR_volume_scan.py                                                    #
#                                                                             #
# Run the functions from SYN_RADAR_VOLUME_SCAN.py:                            #
# to calculate synthetic volume scans from EMVORADO and ICON                  #
# --------------------------------------------------------------------------- #

import os
import HEADER_RADAR_toolbox as header
from SYN_RADAR_VOLUME_SCAN import \
    create_8_vol_nc_of_day, \
    create_8_vol_nc_of_day_paralell, \
    create_vol_nc, \
    rad_dict

dir_data_in = header.dir_data_mod

# --------------------------------------------------------------------------- #
# icon + emvorado mixtures

day = '20170725'
for emvorado_run in [
    'EMVO_00500000.2',
    'EMVO_00600000.2',
]:
    da_run = 'ASS_2211'
    icon_run = 'MAIN_2308.0'
    icon_emvorado_run = 'MAIN_2308.1/' + emvorado_run
    spin_up_mm = 30
    # radar_locs=list(rad_dict().keys())[1:2]
    radar_locs = ['PRO']
    # create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
    #                        icon_emvorado_run=icon_emvorado_run,
    #                        spin_up_mm=spin_up_mm,
    #                        radar_locs=radar_locs,
    #                        dir_data_in=header.dir_data_mod,
    #                        dir_data_out=header.dir_data_vol)
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# icon + emvorado mixtures

day = '20170725'
for emvorado_run in [
    'EMVO_00000000.2',
]:
    da_run = 'ASS_2109'
    icon_run = 'MAIN_2109.0'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    spin_up_mm = 30
    radar_locs = ['PRO']
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# all days, all ass, all main, all emv, only PRO
# TODO: following produces daily files if previous day had some data

for day in [
    '20170725',  # start this day
    '20170719',
    '20170720',
    '20170724',
    '20170726',
    '20170727',
    '20170810',
    '20180728',
    '20180809',
    '20180923',
    '20181202',
]:
    for da_run in [
        'ASS_2109',  # ASS_vold,
        'ASS_2111',  # ASS_old,
        'ASS_2211',  # ASS_new,
    ]:
        for icon_run in [
            'MAIN_2109.0',  # MAIN_vold,
            'MAIN_2203.0',  # MAIN_old,
            'MAIN_2211.0',  # MAIN_new,
            'MAIN_2308.0',  # MAIN_newer,
            'MAIN_2308.1',  # MAIN_newer_new2mom-MP,
        ]:
            for emvorado_run in [
                'EMVO_00000000.2',  # 'EMVO_BBold',
                'EMVO_00200000.2',  # 'EMVO_BB-ML',
                'EMVO_00300000.2',  # 'EMVO_no-melt',
                'EMVO_00401000.2',  # 'EMVO_BBnew2momMP_SSDB',
                'EMVO_00400000.2',  # 'EMVO_BBnew2momMP',
                'EMVO_00100000.2',  # 'EMVO_BBnew',

            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
                spin_up_mm = 30
                # radar_locs=list(rad_dict().keys())[1:2]
                radar_locs = ['PRO']
                if os.path.isdir('/'.join([dir_data_in, day,
                                           da_run, icon_emvorado_run])):
                    a = 1 + 1
                    # create_8_vol_nc_of_day(day=day, da_run=da_run,
                    #                        icon_run=icon_run,
                    #                        icon_emvorado_run=
                    #                        icon_emvorado_run,
                    #                        spin_up_mm=spin_up_mm,
                    #                        radar_locs=radar_locs,
                    #                        dir_data_in=header.dir_data_mod,
                    #                        dir_data_out=header.dir_data_vol)
                    create_8_vol_nc_of_day_paralell(
                        day=day,
                        da_run=da_run,
                        icon_run=icon_run,
                        icon_emvorado_run=icon_emvorado_run,
                        spin_up_mm=spin_up_mm,
                        radar_locs=radar_locs,
                        dir_data_in=header.dir_data_mod,
                        dir_data_out=header.dir_data_vol
                    )
                else:
                    if os.path.isdir('/'.join([dir_data_in, str(int(day) - 1),
                                               da_run, icon_emvorado_run])):
                        print('Not existing: ' + '/'.join([dir_data_in, day,
                                                           da_run,
                                                           icon_emvorado_run]))

# --------------------------------------------------------------------------- #
# icon + emvorado mixtures

day = '20170725'
for emvorado_run in [
    'EMVO_00500000.2',
    'EMVO_00600000.2',
]:
    da_run = 'ASS_2211'
    icon_run = 'MAIN_2308.0'
    icon_emvorado_run = 'MAIN_2308.1/' + emvorado_run
    spin_up_mm = 60
    radar_locs = ['PRO']
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# reference

day = '20170725'
for emvorado_run in [
    'EMVO_00000000.2',
]:
    da_run = 'ASS_2109'
    icon_run = 'MAIN_2109.0'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    spin_up_mm = 60
    radar_locs = ['PRO']
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
