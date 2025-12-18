#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 24.11.11                                                 #
# set_syn_RADAR_1_create_volume_scan_ahrflood.py                              #
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
# # 20.02.25   # Done 24.02.25
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.0'  # MAIN_newererererRH8_MP-RUC1.0
#     for emvorado_run in [
#         'EMVO_00010000.2',
#         'EMVO_00510000.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 16.12.24  # start rerun 28.01.25  # Done 24.02.25
# overwrite_EMV = False
# overwrite_ICON = '2025-01-28'
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.1'  # MAIN_newererererRH8_MP-RUC1.0
#     for emvorado_run in [
#         'EMVO_00510000.2',  # EMVO-ganzneu
#         'EMVO_00410000.2',  # EMVO-neu
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 21.02.25   # Done 06.03.25
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.3'
#     for emvorado_run in [
#         'EMVO_00510000.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 24.02.25 2  # Done 07.03.25
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.0'
#     for emvorado_run in [
#         'EMVO_00510200.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 24.02.25 1  # Done
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())  # TODO: '20210713'
# spin_up_mm = 120
# for day in [
#     '20210714',  # Done 06.03.25
#     '20210713',  # TODO
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.03'
#     for emvorado_run in [
#         'EMVO_00510000.2',  # Done 14.03.25
#         'EMVO_00510100.2',  # ICON done
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# # 18.12.24 + 14.01.25  # start rerun 28.01.25
# # INP study
# overwrite_EMV = False
# overwrite_ICON = '2025-01-28'
# radar_locs = ['ESS']  # Done 24.02.25
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     for icon_run in [
#         'MAIN_2411.6',  # INP*5
#         'MAIN_2411.61',  # INP*100
#     ]:
#         for emvorado_run in [
#             'EMVO_00510000.2',  # EMVO-ganzneu
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
# # 11.11.24  # start rerun 28.01.25   # Done 24.02.25
# # old runs
# overwrite_EMV = False
# overwrite_ICON = '2025-01-28'
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714'
# ]:
#     da_run = 'ASS_2407'  # ASS_newerer (old reference)
#     icon_run = 'MAIN_2405.3'  # MAIN_newerererRH8_new2mom-MP_RUC2.0
#     emvorado_run = 'EMVO_00510000.2'  # EMVO_dynwetgrow-BBold_sig-ice-25
#     icon_emvorado_run = icon_run + '/' + emvorado_run
#     create_8_vol_nc_of_day(day=day, da_run=da_run,
#                            icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol,
#                            overwrite_EMV=overwrite_EMV,
#                            overwrite_ICON=overwrite_ICON)
#
# --------------------------------------------------------------------------- #
# # 08.05.25   # Done 12.05.25
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.0'  # MAIN_newererererRH8_MP-RUC1.0
#     for emvorado_run in [
#         'EMVO_00410000.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)
#
# # --------------------------------------------------------------------------- #
# # 13.05.25 # Kobras runs  # TODO
# overwrite_EMV = False
# overwrite_ICON = False
# radar_locs = list(rad_dict().keys())# TODO
# # radar_locs = ['ESS']
# spin_up_mm = 120
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.3'
#     for emvorado_run in [
#         'EMVO_00511300.2',
#         'EMVO_00521300.2',
#         'EMVO_01521300.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         create_8_vol_nc_of_day(day=day, da_run=da_run,
#                                icon_run=icon_run,
#                                icon_emvorado_run=icon_emvorado_run,
#                                spin_up_mm=spin_up_mm,
#                                radar_locs=radar_locs,
#                                dir_data_in=header.dir_data_mod,
#                                dir_data_out=header.dir_data_vol,
#                                overwrite_EMV=overwrite_EMV,
#                                overwrite_ICON=overwrite_ICON)

# --------------------------------------------------------------------------- #
# 17.12.25 # new thresholds to get rid of graupel oversizes
overwrite_EMV = False
overwrite_ICON = False
# radar_locs = list(rad_dict().keys())# TODO
radar_locs = ['ESS']
spin_up_mm = 120
for day in [
    '20210714',
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.3'
    for emvorado_run in [
        'EMVO_00510040.2',
        'EMVO_00510030.2',
        'EMVO_00510020.2',
        'EMVO_00510010.2',
        'EMVO_005100NN.2',
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

# CONTINUE?
# import process_syn_RADAR_QVP_from_volume_scan_ahrflood
