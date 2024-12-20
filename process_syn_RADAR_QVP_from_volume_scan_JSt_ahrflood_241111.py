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
# 11.1.24
# Ahr flooding processing

for da_run in [
    'ASS_2407',
]:
    for icon_run in [
        'MAIN_2405.3',
    ]:
        for emvorado_run in [
            'EMVO_00510000.2',
        ]:
            for day in [
                '20210714',
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
                spin_up_mm = 120
                elevation_deg = 12
                radar_locs = list(rad_dict().keys())
                # radar_locs = ['ESS']
                for radar_loc in radar_locs:
                    print(radar_loc)
                    qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                     icon_emvorado_run=icon_emvorado_run,
                                     elevation_deg=elevation_deg,
                                     spin_up_mm=spin_up_mm,
                                     radar_loc=radar_loc)


# --------------------------------------------------------------------------- #
# 16.12.24
# Ahr flooding processing

for da_run in [
    'ASS_2411',
]:
    for icon_run in [
        'MAIN_2411.1',
    ]:
        for emvorado_run in [
            'EMVO_00510000.2',
            'EMVO_00410000.2'
        ]:
            for day in [
                '20210714',
                # '20210713',
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
                spin_up_mm = 120
                elevation_deg = 12
                radar_locs = list(rad_dict().keys())
                # radar_locs = ['ESS']
                for radar_loc in radar_locs:
                    print(radar_loc)
                    qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                     icon_emvorado_run=icon_emvorado_run,
                                     elevation_deg=elevation_deg,
                                     spin_up_mm=spin_up_mm,
                                     radar_loc=radar_loc)


# --------------------------------------------------------------------------- #

# 18.12.24
# Ahr flooding processing

for da_run in [
    'ASS_2411',
]:
    for icon_run in [
        'MAIN_2411.6',
    ]:
        for emvorado_run in [
            'EMVO_00510000.2',
        ]:
            for day in [
                '20210714',
                # '20210713',
            ]:
                icon_emvorado_run = icon_run + '/' + emvorado_run
                spin_up_mm = 120
                elevation_deg = 12
                radar_locs = list(rad_dict().keys())
                # radar_locs = ['ESS']
                for radar_loc in radar_locs:
                    print(radar_loc)
                    qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                     icon_emvorado_run=icon_emvorado_run,
                                     elevation_deg=elevation_deg,
                                     spin_up_mm=spin_up_mm,
                                     radar_loc=radar_loc)


# --------------------------------------------------------------------------- #
