#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# process_syn_RADAR_QVP_from_volume_scan_oldref.py                            #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# 11.02.25
# radar_locs = list(rad_dict().keys())
radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
overwrite = '2025-02-14'
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
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 overwrite=overwrite,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# 21.11.23
# skipped ASS_2109 30min
# --------------------------------------------------------------------------- #
# 23.11.23
# skipped ASS_2109 60min
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# 21.11.23
# # radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
# spin_up_mm = 30
# elevation_deg = 12
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
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # 21.11.23
# radar_locs = ['PRO']
# spin_up_mm = 30
# elevation_deg = 12
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
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 08.02.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
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
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 21.11.23
# radar_locs = ['PRO']
# spin_up_mm = 30
# elevation_deg = 12
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
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 08.02.24
radar_locs = ['PRO']
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00400000.2',  # EMVO_BBnew2momMP
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
            'EMVO_00600000.2',  # EMVO_dynwetgrow-wgEMAmai
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 08.02.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00400000.2',  # EMVO_BBnew2momMP
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# 17.06.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
            'EMVO_00510000.2',  # EMVO_dynwetgrow-BBold_sig-ice-25
        ]:
            # MAIN_newererRH8_new2mom-MP
            icon_emvorado_run = 'MAIN_2401.1/' + emvorado_run  # mixture
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 07.02.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new,
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
            'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
            'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
        ]:
            # MAIN_newererRH8limrime_new2mom-MP_limrime
            icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run  # mixture
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 26.02.24
radar_locs = ['PRO']
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new,
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00500004.2',  # EMVO_dynwetgrow-BBold no snow
            'EMVO_00500006.2',  # EMVO_dynwetgrow-BBold no hail
        ]:
            # MAIN_newererRH8limrime_new2mom-MP_limrime
            icon_emvorado_run = 'MAIN_2401.3' + '/' + emvorado_run  # mixture
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 21.03.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new,
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
            'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
            'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
        ]:
            # MAIN_newererRH8limrime_new2mom-MP_limrime+alpha100
            icon_emvorado_run = 'MAIN_2401.4' + '/' + emvorado_run  # mixture
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 21.03.24
radar_locs = list(rad_dict().keys())
# radar_locs = ['PRO']  # TODO
spin_up_mm = 60
elevation_deg = 12
day = '20170725'
for da_run in [
    'ASS_2211',  # ASS_new,
]:
    for icon_run in [
        'MAIN_2308.1',  # MAIN_newer_new2mom-MP
    ]:
        for emvorado_run in [
            'EMVO_00500000.2',  # EMVO_dynwetgrow-BBold
            'EMVO_00500002.2',  # EMVO_dynwetgrow-BBold no rain
            'EMVO_00500005.2',  # EMVO_dynwetgrow-BBold no graupel
        ]:
            # MAIN_newererRH8limrime_new2mom-MP_limrime+thr-rho-sg1e0
            icon_emvorado_run = 'MAIN_2401.5' + '/' + emvorado_run  # mixture
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #

