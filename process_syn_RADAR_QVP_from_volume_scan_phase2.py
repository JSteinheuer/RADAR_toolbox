#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 22.07.24                                                 #
# process_syn_RADAR_QVP_from_volume_scan_phase2.py                            #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict


# --------------------------------------------------------------------------- #
# # 22.07.24
# # PRISTINE runs
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# elevation_deg = 12
# for day in [
#     '20181223',
#     '20181224'
# ]:
#     da_run = 'ASS_2407'
#     icon_run = 'MAIN_2405.3'
#     for emvorado_run in [
#         'EMVO_00510000.2',
#         'EMVO_00510200.2',
#         'EMVO_00510300.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         for radar_loc in radar_locs:
#             qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                              icon_emvorado_run=icon_emvorado_run,
#                              elevation_deg=elevation_deg,
#                              spin_up_mm=spin_up_mm,
#                              radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# processing: 17.06.24
# only FLD but all sweeps
radar_loc = 'FLD'
spin_up_mm = 60
for day in [
    '20220520',
    '20220519',
]:
    da_run = 'ASS_2405'
    icon_run = 'MAIN_2405.1'
    emvorado_run = 'EMVO_00400000.2'
    icon_emvorado_run = icon_run + '/' + emvorado_run
    for elevation_deg in header.ELEVATIONS_ALL:
        qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                         icon_emvorado_run=icon_emvorado_run,
                         elevation_deg=elevation_deg,
                         spin_up_mm=spin_up_mm,
                         radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # TODO: 2025 not done yet
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 60
# elevation_deg = 12
# for day in [
#     '20220519',
#     '20220520',
# ]:
#     da_run = 'ASS_2405'
#     icon_run = 'MAIN_2405.1'
#     emvorado_run = 'EMVO_00400000.2'
#     icon_emvorado_run = icon_run + '/' + emvorado_run
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # TODO: 2025 not done yet
# radar_locs = list(rad_dict().keys())
# radar_locs = ['FLD', 'ESS', 'HNR']
# spin_up_mm = 120
# elevation_deg = 12
# for day in [
#     '20220519',
#     '20220520',
# ]:
#     da_run = 'ASS_2405'
#     icon_run = 'MAIN_2405.1'
#     emvorado_run = 'EMVO_00510000.2'
#     icon_emvorado_run = 'MAIN_2405.3' + '/' + emvorado_run  # mixture
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          elevation_deg=elevation_deg,
#                          spin_up_mm=spin_up_mm,
#                          radar_loc=radar_loc)

# --------------------------------------------------------------------------- #

