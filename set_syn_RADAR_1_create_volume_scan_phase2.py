#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# set_syn_RADAR_1_create_volume_scan_phase2.py                                #
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
# 19.07.24
# PRISTINE case
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS', 'NHB']
spin_up_mm = 120
for day in [
    '20181223',
    '20181224'
]:
    da_run = 'ASS_2407'
    icon_run = 'MAIN_2405.3'
    for emvorado_run in [
        'EMVO_00510000.2',
        'EMVO_00512000.2',
        'EMVO_00513000.2',
        'EMVO_00513900.2',  # 29.10.25
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        create_8_vol_nc_of_day(day=day, da_run=da_run,
                               icon_run=icon_run,
                               icon_emvorado_run=icon_emvorado_run,
                               spin_up_mm=spin_up_mm,
                               radar_locs=radar_locs,
                               dir_data_in=header.dir_data_mod,
                               # dir_data_out=header.dir_data_vol,
                               dir_data_out='/automount/agradar/operation_hydrometeors/data/Syn_vol/',
                               overwrite_EMV=overwrite_EMV,
                               overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# 14.06.24
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = list(rad_dict().keys())
spin_up_mm = 60
for day in [
    '20220519',
    '20220520',
]:
    da_run = 'ASS_2405'
    icon_run = 'MAIN_2405.1'
    emvorado_run = 'EMVO_00400000.2'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
                           icon_emvorado_run=icon_emvorado_run,
                           spin_up_mm=spin_up_mm,
                           radar_locs=radar_locs,
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol,
                           overwrite_EMV=overwrite_EMV,
                           overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# 04.07.24
overwrite_EMV = False
# overwrite_ICON = False
overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = list(rad_dict().keys())
radar_locs = ['FLD', 'ESS', 'HNR']  # TODO: recalculate everything?
spin_up_mm = 120
for day in [
    '20220519',
    '20220520',
]:
    da_run = 'ASS_2405'
    icon_run = 'MAIN_2405.1'
    emvorado_run = 'EMVO_00400000.2'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
                           icon_emvorado_run=icon_emvorado_run,
                           spin_up_mm=spin_up_mm,
                           radar_locs=radar_locs,
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol,
                           overwrite_EMV=overwrite_EMV,
                           overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # TODO: 2025 not done yet
# overwrite_EMV = False
# # overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# radar_locs = ['FLD', 'ESS', 'HNR']
# spin_up_mm = 120
# for day in [
#     '20220519',
#     '20220520',
# ]:
#     da_run = 'ASS_2405'
#     icon_run = 'MAIN_2405.1'
#     emvorado_run = 'EMVO_00510000.2'
#     icon_emvorado_run = 'MAIN_2405.3' + '/' + emvorado_run  # mixture
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol,
#                            overwrite_EMV=overwrite_EMV,
#                            overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
