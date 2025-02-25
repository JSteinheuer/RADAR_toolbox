#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 24.11.24                                                 #
# process_syn_RADAR_QVP_from_volume_scan_ahrflood.py                          #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# 20.02.25    # TODO
# radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    # '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.0'  # MAIN_newererererRH8_MP-RUC1.0
    for emvorado_run in [
        'EMVO_00010000.2',
        'EMVO_00510000.2',  # EMVO_dynwetgrow-BBold_sig-ice-25
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                             icon_emvorado_run=icon_emvorado_run,
                             spin_up_mm=spin_up_mm,
                             elevation_deg=elevation_deg,
                             overwrite=overwrite,
                             radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 21.02.25  # TODO
# radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    # '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.3'
    for emvorado_run in [
        'EMVO_00510000.2',
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                             icon_emvorado_run=icon_emvorado_run,
                             spin_up_mm=spin_up_mm,
                             elevation_deg=elevation_deg,
                             overwrite=overwrite,
                             radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 24.02.25 1  # TODO
# radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    # '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.03'
    for emvorado_run in [
        'EMVO_00510000.2',
        'EMVO_00510100.2',
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                             icon_emvorado_run=icon_emvorado_run,
                             spin_up_mm=spin_up_mm,
                             elevation_deg=elevation_deg,
                             overwrite=overwrite,
                             radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 24.02.25 2  # TODO
# radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.0'
    for emvorado_run in [
        'EMVO_00510200.2',
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                             icon_emvorado_run=icon_emvorado_run,
                             spin_up_mm=spin_up_mm,
                             elevation_deg=elevation_deg,
                             overwrite=overwrite,
                             radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 16.12.24  # start rerun 28.01.25
# radar_locs = list(rad_dict().keys())  # TODO
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    icon_run = 'MAIN_2411.1'  # MAIN_newererererRH8_MP-RUC1.0
    for emvorado_run in [
        'EMVO_00510000.2',  # EMVO-ganzneu
        'EMVO_00410000.2',  # EMVO-neu
    ]:
        icon_emvorado_run = icon_run + '/' + emvorado_run
        for radar_loc in radar_locs:
            qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                             icon_emvorado_run=icon_emvorado_run,
                             spin_up_mm=spin_up_mm,
                             elevation_deg=elevation_deg,
                             overwrite=overwrite,
                             radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 18.12.24 + 14.01.25  # start rerun 28.01.25 # TODO: rerun ?!
# INP study
radar_locs = ['ESS']
spin_up_mm = 120
overwrite = '2025-02-13'
elevation_deg = 12
for day in [
    '20210714',
    '20210713',  # TODO
]:
    da_run = 'ASS_2411'  # ASS_newererer
    for icon_run in [
        'MAIN_2411.6',  # INP*5
        'MAIN_2411.61',  # INP*100
    ]:
        for emvorado_run in [
            'EMVO_00510000.2',  # EMVO-ganzneu
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 icon_emvorado_run=icon_emvorado_run,
                                 spin_up_mm=spin_up_mm,
                                 elevation_deg=elevation_deg,
                                 overwrite=overwrite,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# 11.11.24  # TODO: rerun ?!
# old runs
radar_locs = list(rad_dict().keys())
spin_up_mm = 120
overwrite = False
elevation_deg = 12
for day in [
    '20210714'
]:
    da_run = 'ASS_2407'  # ASS_newerer (old reference)
    icon_run = 'MAIN_2405.3'  # MAIN_newerererRH8_new2mom-MP_RUC2.0
    emvorado_run = 'EMVO_00510000.2'  # EMVO_dynwetgrow-BBold_sig-ice-25
    icon_emvorado_run = icon_run + '/' + emvorado_run
    for radar_loc in radar_locs:
        qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                         icon_emvorado_run=icon_emvorado_run,
                         spin_up_mm=spin_up_mm,
                         elevation_deg=elevation_deg,
                         overwrite=overwrite,
                         radar_loc=radar_loc)

# --------------------------------------------------------------------------- #

