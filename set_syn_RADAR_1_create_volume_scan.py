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
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# for emvorado_run in [
#     'EMVO_00500000.2',
#     'EMVO_00600000.2',
# ]:
#     icon_emvorado_run = 'MAIN_2308.1/' + emvorado_run
#     spin_up_mm = 30
#     # radar_locs=list(rad_dict().keys())[1:2]
#     radar_locs = ['PRO']
#     # create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#     #                        icon_emvorado_run=icon_emvorado_run,
#     #                        spin_up_mm=spin_up_mm,
#     #                        radar_locs=radar_locs,
#     #                        dir_data_in=header.dir_data_mod,
#     #                        dir_data_out=header.dir_data_vol)
#     create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
#                                     icon_run=icon_run,
#                                     icon_emvorado_run=icon_emvorado_run,
#                                     spin_up_mm=spin_up_mm,
#                                     radar_locs=radar_locs,
#                                     dir_data_in=header.dir_data_mod,
#                                     dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2109'
# icon_run = 'MAIN_2109.0'
# for emvorado_run in [
#     'EMVO_00000000.2',
# ]:
#     icon_emvorado_run = icon_run + '/' + emvorado_run
#     spin_up_mm = 30
#     radar_locs = ['ESS']
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)
#     # create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
#     #                                 icon_run=icon_run,
#     #                                 icon_emvorado_run=icon_emvorado_run,
#     #                                 spin_up_mm=spin_up_mm,
#     #                                 radar_locs=radar_locs,
#     #                                 dir_data_in=header.dir_data_mod,
#     #                                 dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # all days, all ass, all main, all emv, only PRO
# # TODO: following produces daily files if previous day had some data
#
# for day in [
#     '20170725',  # start this day
#     '20170719',
#     '20170720',
#     '20170724',
#     '20170726',
#     '20170727',
#     '20170810',
#     '20180728',
#     '20180809',
#     '20180923',
#     '20181202',
# ]:
#     for da_run in [
#         'ASS_2109',  # ASS_vold,
#         'ASS_2111',  # ASS_old,
#         'ASS_2211',  # ASS_new,
#     ]:
#         for icon_run in [
#             'MAIN_2109.0',  # MAIN_vold,
#             'MAIN_2203.0',  # MAIN_old,
#             'MAIN_2211.0',  # MAIN_new,
#             'MAIN_2308.0',  # MAIN_newer,
#             'MAIN_2308.1',  # MAIN_newer_new2mom-MP,
#         ]:
#             for emvorado_run in [
#                 'EMVO_00000000.2',  # 'EMVO_BBold',
#                 'EMVO_00200000.2',  # 'EMVO_BB-ML',
#                 'EMVO_00300000.2',  # 'EMVO_no-melt',
#                 'EMVO_00401000.2',  # 'EMVO_BBnew2momMP_SSDB',
#                 'EMVO_00400000.2',  # 'EMVO_BBnew2momMP',
#                 'EMVO_00100000.2',  # 'EMVO_BBnew',
#
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 spin_up_mm = 30
#                 # radar_locs=list(rad_dict().keys())[1:2]
#                 radar_locs = ['PRO']
#                 if os.path.isdir('/'.join([dir_data_in, day,
#                                            da_run, icon_emvorado_run])):
#                     # create_8_vol_nc_of_day(day=day, da_run=da_run,
#                     #                        icon_run=icon_run,
#                     #                        icon_emvorado_run=
#                     #                        icon_emvorado_run,
#                     #                        spin_up_mm=spin_up_mm,
#                     #                        radar_locs=radar_locs,
#                     #                        dir_data_in=header.dir_data_mod,
#                     #                        dir_data_out=
#                     #                        header.dir_data_vol)
#                     create_8_vol_nc_of_day_paralell(
#                         day=day,
#                         da_run=da_run,
#                         icon_run=icon_run,
#                         icon_emvorado_run=icon_emvorado_run,
#                         spin_up_mm=spin_up_mm,
#                         radar_locs=radar_locs,
#                         dir_data_in=header.dir_data_mod,
#                         dir_data_out=header.dir_data_vol
#                     )
#                 else:
#                     if os.path.isdir('/'.join([dir_data_in,
#                                               str(int(day) - 1),
#                                                da_run, icon_emvorado_run])):
#                         print('Not existing: ' +
#                               '/'.join([dir_data_in, day, da_run,
#                                        icon_emvorado_run]))

# --------------------------------------------------------------------------- #
# # 10.7.24
# # all days, Bringi Book runs
#
# for day in [
#     '20170725',  # start this day
#     '20170719',
#     '20170720',
#     '20170724',
#     '20170726',
#     '20170727',
#     '20170810',
#     '20180728',
#     '20180809',
#     '20180923',
#     '20181202',
# ]:
#     for da_run in [
#         'ASS_2111',  # ASS_old,
#     ]:
#         for icon_run in [
#             'MAIN_2203.0',  # MAIN_old,
#         ]:
#             for emvorado_run in [
#                 'EMVO_00000000.2',  # 'EMVO_BBold',
#
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 spin_up_mm = 30
#                 radar_locs=list(rad_dict().keys())
#                 if os.path.isdir('/'.join([dir_data_in, day,
#                                            da_run, icon_emvorado_run])):
#                     create_8_vol_nc_of_day_paralell(
#                         day=day,
#                         da_run=da_run,
#                         icon_run=icon_run,
#                         icon_emvorado_run=icon_emvorado_run,
#                         spin_up_mm=spin_up_mm,
#                         radar_locs=radar_locs,
#                         dir_data_in=header.dir_data_mod,
#                         dir_data_out=header.dir_data_vol
#                     )
#                 else:
#                     if os.path.isdir('/'.join([dir_data_in,
#                                               str(int(day) - 1),
#                                                da_run, icon_emvorado_run])):
#                         print('Not existing: ' +
#                               '/'.join([dir_data_in, day, da_run,
#                                        icon_emvorado_run]))

# --------------------------------------------------------------------------- #
# # 19.7.24
# # all days, Bringi Book runs the seconds ......-.-
#
# for day in [
#     '20170725',
#     '20170810',
#     '20180809',
# ]:
#     for da_run in [
#         'ASS_2211',  # ASS_old,
#     ]:
#         for icon_run in [
#             'MAIN_2211.0',  # MAIN_old,
#         ]:
#             for emvorado_run in [
#                 'EMVO_00000000.2',  # 'EMVO_BBold',
#
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 spin_up_mm = 30
#                 radar_locs=list(rad_dict().keys())
#                 if os.path.isdir('/'.join([dir_data_in, day,
#                                            da_run, icon_emvorado_run])):
#                     create_8_vol_nc_of_day_paralell(
#                         day=day,
#                         da_run=da_run,
#                         icon_run=icon_run,
#                         icon_emvorado_run=icon_emvorado_run,
#                         spin_up_mm=spin_up_mm,
#                         radar_locs=radar_locs,
#                         dir_data_in=header.dir_data_mod,
#                         dir_data_out=header.dir_data_vol
#                     )
#                 else:
#                     if os.path.isdir('/'.join([dir_data_in,
#                                               str(int(day) - 1),
#                                                da_run, icon_emvorado_run])):
#                         print('Not existing: ' +
#                               '/'.join([dir_data_in, day, da_run,
#                                        icon_emvorado_run]))

# --------------------------------------------------------------------------- #
# 19.7.24
# icon + emvorado mixtures
for day in [
    '20181223',
    '20181224'
]:
    da_run = 'ASS_2407'
    icon_run = 'MAIN_2405.3'
    emvorado_run = 'EMVO_00510000.2'
    icon_emvorado_run = 'MAIN_2405.3/' + emvorado_run
    spin_up_mm = 120
    radar_locs = ['ESS', 'NHB']
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)
    radar_locs = list(rad_dict().keys())
    create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
                                    icon_run=icon_run,
                                    icon_emvorado_run=icon_emvorado_run,
                                    spin_up_mm=spin_up_mm,
                                    radar_locs=radar_locs,
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # reference
#
# day = '20170725'
# da_run = 'ASS_2109'
# icon_run = 'MAIN_2109.0'
# for emvorado_run in [
#     'EMVO_00000000.2',
# ]:
#     icon_emvorado_run = icon_run + '/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = ['PRO']
#     create_8_vol_nc_of_day_paralell(day=day, da_run=da_run,
#                                     icon_run=icon_run,
#                                     icon_emvorado_run=icon_emvorado_run,
#                                     spin_up_mm=spin_up_mm,
#                                     radar_locs=radar_locs,
#                                     dir_data_in=header.dir_data_mod,
#                                     dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 31.1.24 # line 4
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_run = 'EMVO_00500000.2'
# icon_emvorado_run = 'MAIN_2401.3/' + emvorado_run
# spin_up_mm = 60
# # radar_locs = ['PRO']
# radar_locs = list(rad_dict().keys())
# create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                        icon_emvorado_run=icon_emvorado_run,
#                        spin_up_mm=spin_up_mm,
#                        radar_locs=radar_locs,
#                        dir_data_in=header.dir_data_mod,
#                        dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 31.1.24 # line 2
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_run = 'EMVO_00400000.2'
# icon_emvorado_run = icon_run + '/' + emvorado_run
# spin_up_mm = 60
# # radar_locs = ['PRO']
# radar_locs = list(rad_dict().keys())
# create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                        icon_emvorado_run=icon_emvorado_run,
#                        spin_up_mm=spin_up_mm,
#                        radar_locs=radar_locs,
#                        dir_data_in=header.dir_data_mod,
#                        dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 7.2.24 # line 3
# # icon + emvorado mixtures
#
# day = '20170725'
# emvorado_run = 'EMVO_00500000.2'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# icon_emvorado_run = 'MAIN_2308.1/' + emvorado_run
# spin_up_mm = 60
# # radar_locs = ['PRO']
# radar_locs = list(rad_dict().keys())
# create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                        icon_emvorado_run=icon_emvorado_run,
#                        spin_up_mm=spin_up_mm,
#                        radar_locs=radar_locs,
#                        dir_data_in=header.dir_data_mod,
#                        dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 7.2.24 # line new 2?!
# # icon + emvorado mixtures
#
# day = '20170725'
# emvorado_run = 'EMVO_00000000.2'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2211.0'
# icon_emvorado_run = 'MAIN_2211.0/' + emvorado_run
# spin_up_mm = 60
# # radar_locs = ['PRO']
# radar_locs = list(rad_dict().keys())
# create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                        icon_emvorado_run=icon_emvorado_run,
#                        spin_up_mm=spin_up_mm,
#                        radar_locs=radar_locs,
#                        dir_data_in=header.dir_data_mod,
#                        dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 26.2.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500002.2',
#                  'EMVO_00500004.2',
#                  'EMVO_00500005.2',
#                  'EMVO_00500006.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.3/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = ['PRO']
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 21.3.24 - 1 of 3
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500002.2',
#                  'EMVO_00500000.2',
#                  'EMVO_00500004.2',
#                  'EMVO_00500005.2',
#                  'EMVO_00500006.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.3/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 21.3.24  - 2 of 3
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500002.2',
#                  'EMVO_00500005.2',
#                  'EMVO_00500000.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.4/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 21.3.24 - 3 of 3
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500002.2',
#                  'EMVO_00500005.2',
#                  'EMVO_00500000.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.5/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 14.6.24
# # icon + emvorado
#
# for day in [
#     '20220519',
#     '20220520',
# ]:
#     da_run = 'ASS_2405'
#     emvorado_run = 'EMVO_00400000.2'
#     icon_run = 'MAIN_2405.1'
#     icon_emvorado_run = 'MAIN_2405.1/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# processing: 28.08.24
# icon + emvorado mixtures

day = '20170725'
da_run = 'ASS_2211'
icon_run = 'MAIN_2308.0'
emvorado_runs = ['EMVO_00500000.2',
                 'EMVO_00510000.2']
for emvorado_run in emvorado_runs:
    icon_emvorado_run = 'MAIN_2401.1/' + emvorado_run
    spin_up_mm = 120
    radar_locs = list(rad_dict().keys())
    radar_locs = ['PRO']
    create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
                           icon_emvorado_run=icon_emvorado_run,
                           spin_up_mm=spin_up_mm,
                           radar_locs=radar_locs,
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 27.06.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500000.2',
#                  'EMVO_00510000.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.1/' + emvorado_run
#     spin_up_mm = 60
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# # processing: 04.7.24
# # icon + emvorado
#
# for day in [
#     '20220519',
#     '20220520',
# ]:
#     da_run = 'ASS_2405'
#     emvorado_run = 'EMVO_00400000.2'
#     icon_run = 'MAIN_2405.1'
#     icon_emvorado_run = 'MAIN_2405.1/' + emvorado_run
#     spin_up_mm = 120
#     radar_locs = list(rad_dict().keys())
#     create_8_vol_nc_of_day(day=day, da_run=da_run, icon_run=icon_run,
#                            icon_emvorado_run=icon_emvorado_run,
#                            spin_up_mm=spin_up_mm,
#                            radar_locs=radar_locs,
#                            dir_data_in=header.dir_data_mod,
#                            dir_data_out=header.dir_data_vol)

# --------------------------------------------------------------------------- #
# processing: 16.7.24 method='Linear'
# icon + emvorado

# import pandas as pd
# method = 'Linear'
#
# day = '20220520'
# da_run = 'ASS_2405'
# emvorado_run = 'EMVO_00400000.2'
# icon_run = 'MAIN_2405.1'
# icon_emvorado_run = 'MAIN_2405.1/' + emvorado_run
# spin_up_mm = 60
# radar_loc = 'FLD'
# time_start = pd.to_datetime(day, format="%Y%m%d")
# for i in range(4):
#     time_end = time_start + pd.Timedelta('6h')
#     print('________________________________________')
#     print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
#           str(spin_up_mm) + '_spinup/')
#     if i == 2:
#         create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
#                       time_end=time_end.strftime('%Y%m%d%H'),
#                       spin_up_mm=spin_up_mm, da_run=da_run,
#                       icon_run=icon_run,
#                       icon_emvorado_run=icon_emvorado_run,
#                       dir_data_in=dir_data_in,
#                       dir_data_out=header.dir_data_vol,
#                       radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
#                       include_icon=True, include_emv=False,
#                       method=method)
#         print('________________________________________')
#         print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
#               str(spin_up_mm) + '_spinup/')
#         create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
#                       time_end=time_end.strftime('%Y%m%d%H'),
#                       spin_up_mm=spin_up_mm, da_run=da_run,
#                       icon_run=icon_run,
#                       icon_emvorado_run=icon_emvorado_run,
#                       dir_data_in=dir_data_in,
#                       dir_data_out=header.dir_data_vol,
#                       radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
#                       include_icon=False, include_emv=True,
#                       method=method)
#
#     time_start = time_end

# --------------------------------------------------------------------------- #
