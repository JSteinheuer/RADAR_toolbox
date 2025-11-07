#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.10.25                                                 #
# process_syn_RADAR_QVP_from_volume_scan_PRISTINE25.py                        #
#                                                                             #
# Run the functions from SET_SYN_RADAR.py to calculate QVPs from given        #
# synthetic (EMVORADO) volume scans                                           #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_SYN_RADAR import qvp_from_syn_vol
from SET_SYN_RADAR import rad_dict


# --------------------------------------------------------------------------- #
# 15.10.25
# PRISTINE runs
# radar_locs = list(rad_dict().keys())
radar_locs = ['ESS', 'NHB']
spin_up_mm = 120
elevation_degs = [12, 8, 17]
for elevation_deg in elevation_degs:
    for day in [
        '20181223',
        # '20181224'
    ]:
        da_run = 'ASS_2407'
        icon_run = 'MAIN_2405.3'
        for emvorado_run in [
            'EMVO_00510000.2',
            # 'EMVO_00512000.2',
            'EMVO_00513000.2',
            'EMVO_00513900.2',
        ]:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            for radar_loc in radar_locs:
                qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
                                 dir_data_in='/automount/agradar/operation_hydrometeors/data/Syn_vol/',
                                 icon_emvorado_run=icon_emvorado_run,
                                 elevation_deg=elevation_deg,
                                 spin_up_mm=spin_up_mm,
                                 lowest_rhv=0,
                                 lowest_zh=-10, lowest_kdp=None,
                                 overwrite=False,
                                 radar_loc=radar_loc)

# --------------------------------------------------------------------------- #

