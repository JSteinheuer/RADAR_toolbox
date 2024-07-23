#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 22.07.24                                                 #
# DWD_obs_to_MIUB_obs_1234567_full_processing_JM_240722.py                    #
#                                                                             #
# run DWD_OBS_TO_MIUB.py                                                      #
# --------------------------------------------------------------------------- #


import HEADER_RADAR_toolbox as header
import numpy as np

from DWD_OBS_TO_MIUB_OBS import \
    load_all_moms, \
    correct_rho_hv, \
    download_ERA5_temp, era5_temp, \
    correct_phi_kdp, \
    calibrate_zdr, \
    attenuation_correction, \
    combine_pol_mom_nc

# --------------------------------------------------------------------------- #
#  ! SET THE PARAMETERS !                                                     #
# --------------------------------------------------------------------------- #
# Parameters for ALL STEPS:
DATES = [
    # "20210604",  # case01
    # "20210620", "20210621",  # case02
    # "20210628", "20210629",  # case03
    # "20220519", "20220520",  # case04
    # "20220623", "20220624", "20220625",  # case05
    # "20220626", "20220627", "20220628",  # case06+07
    # "20220630", "20220701",  # case08
    # "20210714",  # case09
    # "20221222",  # case10
    # "20170719",  # caseX -> old OP HM 1 case
    # "20170725",  # caseX -> old OP HM 1 case
    "20181223",  "20181224",  # case -> PRISTINE
]
LOCATIONS = [
    # 'asb', 'boo', 'drs', 'eis',
    # 'ess',
    # 'fbg',
    # 'fld', 'hnr', 'isn', 'mem', 'neu',
    'nhb',
    # 'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,  # for pcp
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8.0, 12.0, 17.0, 25.0,
    5.5,  # for 90grad
])
MODE = [
    'pcp',
    'vol', 'vol', 'vol', 'vol', 'vol', 'vol',
    'vol', 'vol', 'vol', 'vol',
    '90grad',
]
ELEVATIONS = [12]
MODE = ['vol']
overwrite = False
# --------------------------------------------------------------------------- #
# Parameters for STEP 1
moments = ['CMAP', 'DBSNRH', 'DBZH', 'RHOHV', 'UPHIDP', 'ZDR', 'SNRHC',
           'VRADH', ]
# --------------------------------------------------------------------------- #
# Parameters for STEP 4
uh_tresh = 0
rho_tresh = 0.8
snr_tresh = 15
win_r = 25  # VP or PARK?!
win_azi = None
rng = 3000
wkdp_light = 9  # PARK: light filtering
wkdp_heavy = 25  # PARK: heavy filtering
parts = 6
merge = True
remove_parts = True
# --------------------------------------------------------------------------- #

# OSError: [Errno -51] NetCDF: Unknown file format: b'/automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-caseX-20181223/2018/2018-12/2018-12-23/nhb/vol5minng01/04/ras11-vol5minng01_sweeph5allm_rhohv_nc_04-201812230102-201812232357-nhb-10605.hd5'
# s6justei@sv47:~$ /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_1234567_full_processing_JM_240722_nhb.py

# --------------------------------------------------------------------------- #
# Parameters for STEP 7
other_zdr_off_day = ''
n_zdr_lowest = 1000
std_zdr_highest = 2
method_zdr_priorities = ['BB', 'SD_V', 'LR_V' 'SD_I', 'LR_I']
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:                                  #
# --------------------------------------------------------------------------- #
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg, mode in zip(ELEVATIONS, MODE):
            # STEP 1
            load_all_moms(date=date, location=location,
                          elevation_deg=elevation_deg,
                          mode=mode, moments=moments,
                          overwrite=overwrite,
                          dir_data_obs=header.dir_data_obs,
                          dir_data_obs_realpep=header.dir_data_obs_realpep)
            # STEP 2
            if mode != '90grad':
                correct_rho_hv(date=date, location=location,
                               elevation_deg=elevation_deg,
                               mode=mode, overwrite=overwrite,
                               dir_data_obs=header.dir_data_obs)

            # STEP 3
            download_ERA5_temp(date=date, overwrite=overwrite,
                               dir_out=header.dir_data_era5)
            era5_temp(date=date, location=location,
                      elevation_deg=elevation_deg,
                      mode=mode, overwrite=overwrite,
                      dir_data_obs=header.dir_data_obs,
                      dir_data_era5=header.dir_data_era5)
            # STEP 4
            if mode != '90grad':
                correct_phi_kdp(date=date, location=location,
                                elevation_deg=elevation_deg, mode=mode,
                                overwrite=overwrite,
                                dir_data_obs=header.dir_data_obs,
                                parts=parts, merge=merge,
                                remove_parts=remove_parts,
                                uh_tresh=uh_tresh,
                                rho_tresh=rho_tresh,
                                snr_tresh=snr_tresh,
                                win_r=win_r, win_azi=win_azi, rng=rng,
                                wkdp_light=wkdp_light,
                                wkdp_heavy=wkdp_heavy)

            # STEP 5
            calibrate_zdr(date=date, location=location,
                          elevation_deg=elevation_deg, mode=mode,
                          overwrite=overwrite)
            # STEP 6
            if mode != '90grad':
                attenuation_correction(date=date, location=location,
                                       elevation_deg=elevation_deg,
                                       mode=mode,
                                       overwrite=overwrite,
                                       dir_data_obs=header.dir_data_obs)

        for elevation_deg, mode in zip(ELEVATIONS, MODE):  # for Vol zdr calib.
            # STEP 7
            if mode != '90grad':
                other_zdr_off_day = ''
                n_zdr_lowest = 1000
                std_zdr_highest = 2
                method_zdr_priorities = ['BB', 'SD_V', 'LR_V' 'SD_I', 'LR_I']
                combine_pol_mom_nc(date=date, location=location,
                                   elevation_deg=elevation_deg, mode=mode,
                                   overwrite=overwrite,
                                   n_zdr_lowest=n_zdr_lowest,
                                   std_zdr_highest=std_zdr_highest,
                                   other_zdr_off_day=other_zdr_off_day,
                                   method_zdr_priorities=method_zdr_priorities,
                                   dir_data_obs=header.dir_data_obs)

# --------------------------------------------------------------------------- #
