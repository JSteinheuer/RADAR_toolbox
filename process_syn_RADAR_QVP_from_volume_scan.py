#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# process_syn_RADAR_QVP_from_volume_scan.py                                   #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# # one day
#
# day = '20170725'
# da_run = 'ASS_2109'
# icon_run = 'MAIN_2109.0'
# emvorado_run = 'EMVO_00000000.2'
# icon_emvorado_run = icon_run + '/' + emvorado_run
# spin_up_mm = 30
# elevation_deg = 12
# radar_loc = 'ESS'
# qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                  icon_emvorado_run=icon_emvorado_run,
#                  spin_up_mm=spin_up_mm,
#                  elevation_deg=elevation_deg, radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # one day
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_run = 'EMVO_00500000.2'
# icon_emvorado_run = 'MAIN_2308.1' + '/' + emvorado_run
# spin_up_mm = 60
# elevation_deg = 12
# radar_loc = 'PRO'
# qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                  icon_emvorado_run=icon_emvorado_run,
#                  spin_up_mm=spin_up_mm,
#                  elevation_deg=elevation_deg, radar_loc=radar_loc)
#
# emvorado_run = 'EMVO_00600000.2'
# icon_emvorado_run = 'MAIN_2308.1' + '/' + emvorado_run
# qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                  icon_emvorado_run=icon_emvorado_run,
#                  spin_up_mm=spin_up_mm,
#                  elevation_deg=elevation_deg, radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # all days, all ass, all main, all emv, only PRO
#
# for day in [
#     '20170725',  # start this day
#     '20170719',
#     '20170720',
#     '20170724',
#     # '20170725',
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
#                 elevation_deg = 12
#                 radar_loc = 'PRO'
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  elevation_deg=elevation_deg,
#                                  spin_up_mm=spin_up_mm,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 10.7.24
# all days, Bringi Book runs

for day in [
    '20170725',  # start this day
    '20170719',
    '20170720',
    '20170724',
    '20170726',
    '20170727',  # TODO!
    '20170810',
    '20180728',
    '20180809',
    '20180923',
    '20181202',
]:
    for da_run in [
        'ASS_2111',  # ASS_old,
        'ASS_2211',
    ]:
        for icon_run in [
            'MAIN_2203.0',  # MAIN_old,
            'MAIN_2211.0',
        ]:
            for emvorado_run in [
                'EMVO_00000000.2',  # 'EMVO_BBold',

            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
                spin_up_mm = 30
                elevation_deg = 12
                radar_locs = list(rad_dict().keys())
                for radar_loc in radar_locs:
                    qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                     icon_emvorado_run=icon_emvorado_run,
                                     elevation_deg=elevation_deg,
                                     spin_up_mm=spin_up_mm,
                                     radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 31.1.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_run = 'EMVO_00400000.2'
# icon_emvorado_run = icon_run + '/' + emvorado_run
# spin_up_mm = 60
# elevation_deg = 12
# radar_locs = list(rad_dict().keys())
# for radar_loc in radar_locs:
#     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                      icon_emvorado_run=icon_emvorado_run,
#                      elevation_deg=elevation_deg,
#                      spin_up_mm=spin_up_mm,
#                      radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 31.1.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_run = 'EMVO_00500000.2'
# icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run
# spin_up_mm = 60
# elevation_deg = 12
# radar_locs = list(rad_dict().keys())
# for radar_loc in radar_locs:
#     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                      icon_emvorado_run=icon_emvorado_run,
#                      elevation_deg=elevation_deg,
#                      spin_up_mm=spin_up_mm,
#                      radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 7.2.24
# # icon + emvorado mixtures
#
# day = '20170725'
# emvorado_run = 'EMVO_00500000.2'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# icon_emvorado_run = 'MAIN_2308.1' + '/' + emvorado_run
# spin_up_mm = 60
# elevation_deg = 12
# radar_locs = list(rad_dict().keys())
# for radar_loc in radar_locs:
#     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                      icon_emvorado_run=icon_emvorado_run,
#                      elevation_deg=elevation_deg,
#                      spin_up_mm=spin_up_mm,
#                      radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 7.2.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2211.0'
# emvorado_run = 'EMVO_00000000.2'
# icon_emvorado_run = 'MAIN_2211.0' + '/' + emvorado_run
# spin_up_mm = 60
# elevation_deg = 12
# radar_locs = list(rad_dict().keys())
# for radar_loc in radar_locs:
#     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                      icon_emvorado_run=icon_emvorado_run,
#                      elevation_deg=elevation_deg,
#                      spin_up_mm=spin_up_mm,
#                      radar_loc=radar_loc)

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
#     icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run
#     spin_up_mm = 60
#     elevation_deg = 12
#     radar_loc = 'PRO'
#     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                      icon_emvorado_run=icon_emvorado_run,
#                      elevation_deg=elevation_deg,
#                      spin_up_mm=spin_up_mm,
#                      radar_loc=radar_loc)

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
#     icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run
#     spin_up_mm = 60
#     elevation_deg = 12
#     # radar_locs = ['PRO']
#     radar_locs = list(rad_dict().keys())
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 21.3.24 - 2 of 3
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500002.2',
#                  'EMVO_00500005.2',
#                  'EMVO_00500000.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.4' + '/' + emvorado_run
#     spin_up_mm = 60
#     elevation_deg = 12
#     # radar_locs = ['PRO']
#     radar_locs = list(rad_dict().keys())
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)

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
#     icon_emvorado_run = 'MAIN_2401.5' + '/' + emvorado_run
#     spin_up_mm = 60
#     elevation_deg = 12
#     # radar_locs = ['PRO']
#     radar_locs = list(rad_dict().keys())
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# # processing: 17.6.24
# # icon + emvorado mixtures
#
# day = '20170725'
# da_run = 'ASS_2211'
# icon_run = 'MAIN_2308.0'
# emvorado_runs = ['EMVO_00500000.2',
#                  'EMVO_00510000.2']
# for emvorado_run in emvorado_runs:
#     icon_emvorado_run = 'MAIN_2401.1' + '/' + emvorado_run
#     spin_up_mm = 60
#     elevation_deg = 12
#     radar_locs = list(rad_dict().keys())
#     # radar_locs = ['PRO']
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # processing: 17.06.24
# # one case fld
#
# da_run = 'ASS_2405'
# icon_run = 'MAIN_2405.1'
# emvorado_run = 'EMVO_00400000.2'
# icon_emvorado_run = 'MAIN_2405.1' + '/' + emvorado_run
# spin_up_mm = 60
# for day in [
#     '20220520',
#     '20220519',
# ]:
#     for elevation_deg in header.ELEVATIONS_ALL:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          elevation_deg=elevation_deg, radar_loc='FLD')
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          n_lowest=5,
#                          lowest_zh=20, highest_zh=40,
#                          elevation_deg=elevation_deg, radar_loc='FLD')
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          n_lowest=5,
#                          lowest_zh=0,
#                          elevation_deg=elevation_deg, radar_loc='FLD')
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          n_lowest=5,
#                          lowest_zh=0, highest_zh=20,
#                          elevation_deg=elevation_deg, radar_loc='FLD')
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          n_lowest=5,
#                          lowest_zh=40, highest_zh=80,
#                          elevation_deg=elevation_deg, radar_loc='FLD')
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          n_lowest=5,
#                          lowest_zh=0, highest_zh=40,
#                          elevation_deg=elevation_deg, radar_loc='FLD')

# --------------------------------------------------------------------------- #
