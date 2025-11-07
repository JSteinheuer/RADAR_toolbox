#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 14.02.25                                                 #
# process_RADAR_QVP_from_volume_scan.py                                       #
#                                                                             #
# Functions to calculate QVPs from given volume scans.                        #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PROCESS_RADAR import qvp_from_radar_PPIs, sector_qvp_from_radar_PPIs
import numpy as np
# --------------------------------------------------------------------------- #

elevation_deg = 12
overwrite = '2025-02-14'
lowest_rhv = 0.7
lowest_kdp = 0
lowest_kdp = None
highest_kdp = 10
n_lowest = 30
lowest_zh = 0
lowest_zh = -10
highest_zh = None
lowest_snr = 10
lowest_snr = 5

# --------------------------------------------------------------------------- #
# 14.02.25
DATES = [
    # "20210714",  # case09
    # "20210713",  # case09
    # "20170725",  # caseX -> old OP HM 1 case
    "20181223",
    ]
LOCATIONS = [
    'ess','nhb',
    # 'pro',
    # 'asb', 'boo', 'drs', 'eis', 'fbg',
    # 'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    # 'oft', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    # 5.5,  # for pcp
    # 5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8.0,
    12.0,
    17.0,
    # 25.0,
    # 5.5,  # for 90grad
])
overwrite = '2025-02-24'
overwrite = '2025-10-01'
for location in LOCATIONS:
    for date in DATES:
        for elevation_deg in ELEVATIONS:
            print(location + ' ' + date + ' ' + str(elevation_deg))
            qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
                                location=location,
                                overwrite=overwrite, lowest_rhv=lowest_rhv,
                                lowest_snr=lowest_snr,
                                lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
                                lowest_zh=lowest_zh, highest_zh=highest_zh,
                                n_lowest=n_lowest)

# min_azimuth=210
# max_azimuth=240
# n_lowest_sector=5
# for location in ['ess']:
#     for date in ["20181223"]:
#         print(location + ' ' + date + ' ' + str(elevation_deg))
#         sector_qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
#                                    location=location,
#                                    overwrite=overwrite, lowest_rhv=lowest_rhv,
#                                    lowest_snr=lowest_snr,
#                                    lowest_kdp=lowest_kdp,
#                                    highest_kdp=highest_kdp,
#                                    lowest_zh=lowest_zh,
#                                    highest_zh=highest_zh,
#                                    min_azimuth=min_azimuth,
#                                    max_azimuth=max_azimuth,
#                                    n_lowest=n_lowest_sector)
#
# min_azimuth=360-15
# max_azimuth=15
# n_lowest_sector=5
# for location in ['nhb']:
#     for date in ["20181223"]:
#         print(location + ' ' + date + ' ' + str(elevation_deg))
#         sector_qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
#                                    location=location,
#                                    overwrite=overwrite, lowest_rhv=lowest_rhv,
#                                    lowest_snr=lowest_snr,
#                                    lowest_kdp=lowest_kdp,
#                                    highest_kdp=highest_kdp,
#                                    lowest_zh=lowest_zh,
#                                    highest_zh=highest_zh,
#                                    min_azimuth=min_azimuth,
#                                    max_azimuth=max_azimuth,
#                                    n_lowest=n_lowest_sector)
# --------------------------------------------------------------------------- #
# # TODO but SNRs is missing
# # 14.02.25
# DATES = [
#     "20181223",  "20181224",  # case -> PRISTINE  # TODO: missing SNRH
# ]
# LOCATIONS = [
#     'ess', 'pro',
#     'asb', 'boo', 'drs', 'eis', 'fbg',
#     'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
#     'oft', 'ros', 'tur', 'umd',
# ]
# for date in DATES:
#     for location in LOCATIONS:
#         print(location + ' ' + date + '' + str(elevation_deg))
#         qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
#                             location=location,
#                             overwrite=overwrite, lowest_rhv=lowest_rhv,
#                             lowest_snr=lowest_snr,
#                             lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
#                             lowest_zh=lowest_zh, highest_zh=highest_zh,
#                             n_lowest=n_lowest)

# --------------------------------------------------------------------------- #
# # 2025 TODO ALL; but SNRs might be missing and do you net it?
# # 14.02.25
# DATES = [
#     "20210604",  # case01
#     "20210620", "20210621",  # case02
#     "20210628", "20210629",  # case03
#     "20220519", "20220520",  # case04
#     "20220623", "20220624", "20220625",  # case05
#     "20220626", "20220627", "20220628",  # case06+07
#     "20220630", "20220701",  # case08
#     "20210713",  # case09
#     "20210714",  # case09
#     "20221222",  # case10
#     "20170719",  # caseX -> old OP HM 1 case
#     "20170725",  # caseX -> old OP HM 1 case
#     "20181223",  "20181224",  # case -> PRISTINE  # TODO: missing SNRH
# ]
# LOCATIONS = [
#     'ess', 'pro',
#     'asb', 'boo', 'drs', 'eis', 'fbg',
#     'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
#     'oft', 'ros', 'tur', 'umd',
# ]
# for date in DATES:
#     for location in LOCATIONS:
#         print(location + ' ' + date + '' + str(elevation_deg))
#         qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
#                             location=location,
#                             overwrite=overwrite, lowest_rhv=lowest_rhv,
#                             lowest_snr=lowest_snr,
#                             lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
#                             lowest_zh=lowest_zh, highest_zh=highest_zh,
#                             n_lowest=n_lowest)

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #

# overwrite = '2025-02-14'
# lowest_rhv = 0.7
# lowest_kdp = 0
# lowest_kdp = None
# highest_kdp = 10
# n_lowest = 30
# lowest_zh = 0
# highest_zh = None
# lowest_snr = 10
#
# # --------------------------------------------------------------------------- #
# # 14.02.25
# DATES = [
#     "20210714",  # case09
#     # "20210713",  # case09
# ]
# LOCATIONS = [
#     'ess',
#     # 'pro',
#     # 'asb', 'boo', 'drs', 'eis', 'fbg',
#     # 'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
#     # 'oft', 'ros', 'tur', 'umd',
# ]
# overwrite = '2025-02-24'
# for elevation_deg in [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]:
#     for location in LOCATIONS:
#         for date in DATES:
#             print(location + ' ' + date + ' ' + str(elevation_deg))
#             qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
#                                 location=location,
#                                 overwrite=overwrite, lowest_rhv=lowest_rhv,
#                                 lowest_snr=lowest_snr,
#                                 lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
#                                 lowest_zh=lowest_zh, highest_zh=highest_zh,
#                                 n_lowest=n_lowest)

# --------------------------------------------------------------------------- #