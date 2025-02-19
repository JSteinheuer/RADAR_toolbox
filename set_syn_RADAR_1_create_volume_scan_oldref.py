#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# set_syn_RADAR_1_create_volume_scan_oldref.py                                #
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
# 10.02.25
overwrite_EMV = False
overwrite_ICON = False
radar_locs = list(rad_dict().keys())
spin_up_mm = 60
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new,
]:
    for icon_run in [
        'MAIN_2411.0',  # ICON-Stand in 2411 im Setup 0 (bzw alte Mikrophysik)
    ]:
        for emvorado_run in [
            'EMVO_00000000.2',  # EMVO_BBold
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run  # mixture
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
# --------------------------------------------------------------------------- #
# 21.11.23
# Consider rerun for all ICON volumes with new interpolation method
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = ['PRO']
spin_up_mm = 30
day = '20170725'
for da_run in [
    'ASS_2109',  # ASS_vold,
]:
    for icon_run in [
        'MAIN_2109.0',  # MAIN_vold,
    ]:
        for emvorado_run in [
            'EMVO_00000000.2',  # 'EMVO_BBold',
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
# # 23.11.23
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = ['PRO']
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2109',  # ASS_vold,
# ]:
#     for icon_run in [
#         'MAIN_2109.0',  # MAIN_vold,
#     ]:
#         for emvorado_run in [
#             'EMVO_00000000.2',  # 'EMVO_BBold',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# 21.11.23
# Consider rerun for all ICON volumes with new interpolation method
overwrite_EMV = False
overwrite_ICON = False
# overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
radar_locs = list(rad_dict().keys())
spin_up_mm = 30
day = '20170725'
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
# # 21.11.23
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = ['PRO']
# spin_up_mm = 30
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2211.0',  # MAIN_new
#     ]:
#         for emvorado_run in [
#             'EMVO_00100000.2',  # EMVO_BBnew
#             'EMVO_00200000.2',  # EMVO_BB-ML
#             'EMVO_00300000.2',  # EMVO_no-melt
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 08.02.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2211.0',  # MAIN_new
#     ]:
#         for emvorado_run in [
#             'EMVO_00000000.2',  # EMVO_BBold
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 21.11.23
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = ['PRO']
# spin_up_mm = 30
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00400000.2',  # EMVO_BBnew2momMP
#             'EMVO_00401000.2',  # EMVO_BBnew2momMP_SSDB
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00600000.2',  # EMVO_dynwetgrow-wgEMAmai
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 08.02.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = ['PRO']
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00400000.2',  # EMVO_BBnew2momMP
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00600000.2',  # EMVO_dynwetgrow-wgEMAmai
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 08.02.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00400000.2',  # EMVO_BBnew2momMP
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 17.06.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00510000.2',  # EMVO_dynwetgrow-BBold_sig-ice-25
#         ]:
#             # MAIN_newererRH8_new2mom-MP
#             icon_emvorado_run = 'MAIN_2401.1/' + emvorado_run  # mixture
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 07.02.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new,
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
#             'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
#         ]:
#             # MAIN_newererRH8limrime_new2mom-MP_limrime
#             icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run  # mixture
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 26.02.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = ['PRO']
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new,
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00500004.2',  # EMVO_dynwetgrow-BBold no snow
#             'EMVO_00500006.2',  # EMVO_dynwetgrow-BBold no hail
#         ]:
#             # MAIN_newererRH8limrime_new2mom-MP_limrime
#             icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run  # mixture
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 21.03.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new,
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
#             'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
#         ]:
#             # MAIN_newererRH8limrime_new2mom-MP_limrime+alpha100
#             icon_emvorado_run = 'MAIN_2401.4' + '/' + emvorado_run  # mixture
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 21.03.24
# # Consider rerun for all ICON volumes with new interpolation method
# overwrite_EMV = False
# overwrite_ICON = False
# # overwrite_ICON = '2025-01-28'  # TODO: recalculate everything?
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# day = '20170725'
# for da_run in [
#     'ASS_2211',  # ASS_new,
# ]:
#     for icon_run in [
#         'MAIN_2308.1',  # MAIN_newer_new2mom-MP
#     ]:
#         for emvorado_run in [
#             'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
#             'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
#             'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
#         ]:
#             # MAIN_newererRH8limrime_new2mom-MP_limrime+thr-rho-sg1e0
#             icon_emvorado_run = 'MAIN_2401.5' + '/' + emvorado_run  # mixture
#             create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                    icon_run=icon_run,
#                                    icon_emvorado_run=icon_emvorado_run,
#                                    spin_up_mm=spin_up_mm,
#                                    radar_locs=radar_locs,
#                                    dir_data_in=header.dir_data_mod,
#                                    dir_data_out=header.dir_data_vol,
#                                    overwrite_EMV=overwrite_EMV,
#                                    overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #

