#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# process_syn_RADAR_QVP_from_volume_scan_bringibook.py                        #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# # 10.07.24
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 30
# elevation_deg = 12
# for day in [
#     '20170719',
#     '20170720',
#     '20170724',
#     '20170726',
#     '20170727',
#     '20180728',
#     '20180923',
#     '20181202',
# ]:
#     for da_run in [
#         'ASS_2111',  # ASS_old
#     ]:
#         for icon_run in [
#             'MAIN_2203.0',  # MAIN_old
#         ]:
#             for emvorado_run in [
#                 'EMVO_00000000.2',  # EMVO_BBold
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 for radar_loc in radar_locs:
#                     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                      icon_emvorado_run=icon_emvorado_run,
#                                      spin_up_mm=spin_up_mm,
#                                      elevation_deg=elevation_deg,
#                                      radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # 10.07.24
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# elevation_deg = 12
# for day in [
#     '20170725',
#     '20170810',
#     '20180809',
# ]:
#     for da_run in [
#         'ASS_2211',  # ASS_new
#     ]:
#         for icon_run in [
#             'MAIN_2211.0',  # MAIN_new
#         ]:
#             for emvorado_run in [
#                 'EMVO_00000000.2',  # EMVO_BBold
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 for radar_loc in radar_locs:
#                     qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                      icon_emvorado_run=icon_emvorado_run,
#                                      spin_up_mm=spin_up_mm,
#                                      elevation_deg=elevation_deg,
#                                      radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 21.11.23
radar_locs = ['PRO']
spin_up_mm = 30
elevation_deg = 12
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
                for radar_loc in radar_locs:
                    qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                     icon_emvorado_run=icon_emvorado_run,
                                     spin_up_mm=spin_up_mm,
                                     elevation_deg=elevation_deg,
                                     radar_loc=radar_loc)

# --------------------------------------------------------------------------- #

