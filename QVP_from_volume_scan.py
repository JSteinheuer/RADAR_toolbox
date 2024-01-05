#!/usr/bin/env python3.11
###############################################################################
# Julian Steinheuer; November 2023                                            #
# QVP_from_volume_scan.py                                                     #
#                                                                             #
# Run the functions from QVP_FROM_VOLUME_SCAN.py:                             #
# calculate QVPs from given synthetic (EMVORADO) volume scans                 #
###############################################################################

import HEADER_RADAR_toolbox as header
from QVP_FROM_VOLUME_SCAN import qvp_from_syn_vol


day = '20170725'
da_run = 'ASS_2211'
icon_run = 'MAIN_2308.0'
qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                 icon_emvorado_run='MAIN_2308.1/EMVO_00500000.2',
                 dir_data_in=header.dir_data_syn,
                 dir_data_out=header.dir_data_qvp)
qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                 icon_emvorado_run='MAIN_2308.1/EMVO_00600000.2',
                 dir_data_in=header.dir_data_syn,
                 dir_data_out=header.dir_data_qvp)

###############################################################################

# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                  spin_up_mm=60,
#                  icon_emvorado_run='MAIN_2308.1/EMVO_00500000.2',
#                  dir_data_in=header.dir_data_syn,
#                  dir_data_out=header.dir_data_qvp)
# qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                  spin_up_mm=60,
#                  icon_emvorado_run='MAIN_2308.1/EMVO_00600000.2',
#                  dir_data_in=header.dir_data_syn,
#                  dir_data_out=header.dir_data_qvp)

###############################################################################

# all days, all ass, all main, all emv, only PRO
for day in [
    '20170725',  # start this day
    '20170719',
    '20170720',
    '20170724',
    # '20170725',
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
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                    icon_emvorado_run=icon_run + '/' + emvorado_run,
                    dir_data_in=header.dir_data_syn,
                    dir_data_out=header.dir_data_qvp)
