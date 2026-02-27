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
from PROCESS_SYN_RADAR import qvp_from_syn_vol, qvp_from_syn_vol_with_qn_qnx
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 20.02.25  # R0E1  # R0E3  # 10.03.25 Done
# radar_locs = list(rad_dict().keys())
# # radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12 # DONE!l
# for elevation_deg in [12, 8, 17]: # DONE!
#     for day in [
#         '20210714',
#         '20210713',
#     ]:
#         da_run = 'ASS_2411'  # ASS_newererer
#         icon_run = 'MAIN_2411.0'  # MAIN_newererererRH8_MP-RUC1.0
#         for emvorado_run in [
#             'EMVO_00010000.2',
#             'EMVO_00510000.2',  # EMVO_dynwetgrow-BBold_sig-ice-25
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 16.12.24  # R1E3  # start rerun 28.01.25  # 10.03.25 Done
# radar_locs = list(rad_dict().keys())
# # radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12 # DONE!
# for elevation_deg in [12, 8, 17]: # DONE!
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.1'  # MAIN_newererererRH8_MP-RUC1.0
#     for emvorado_run in [
#         'EMVO_00510000.2',  # EMVO-ganzneu
#         # 'EMVO_00410000.2',  # EMVO-neu  # not in path of matrix for 8° + 17°
#     ]:
#         for day in [
#             '20210714',
#             '20210713',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 24.02.25 2  # 11.03.25 Done
# radar_locs = list(rad_dict().keys())
# # radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.0'
#     for emvorado_run in [
#         'EMVO_00510200.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         for radar_loc in radar_locs:
#             qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                              icon_emvorado_run=icon_emvorado_run,
#                              spin_up_mm=spin_up_mm,
#                              elevation_deg=elevation_deg,
#                              overwrite=overwrite,
#                              radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 24.02.25 1  # TODO
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12
# for day in [
#     '20210714',  # 12.03.25 Done
#     '20210713',  # TODO
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.03'
#     for emvorado_run in [
#         'EMVO_00510000.2',
#         'EMVO_00510100.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         for radar_loc in radar_locs:
#             qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                              icon_emvorado_run=icon_emvorado_run,
#                              spin_up_mm=spin_up_mm,
#                              elevation_deg=elevation_deg,
#                              overwrite=overwrite,
#                              radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 18.12.24 + 14.01.25  # start rerun 28.01.25 # TODO: rerun ?!
# # INP study
# radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12
# for day in [
#     '20210714',
#     '20210713',  # TODO
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     for icon_run in [
#         'MAIN_2411.6',  # INP*5
#         'MAIN_2411.61',  # INP*100
#     ]:
#         for emvorado_run in [
#             'EMVO_00510000.2',  # EMVO-ganzneu
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 11.11.24  # TODO: rerun ?!
# # old runs
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# overwrite = False
# elevation_deg = 12
# for day in [
#     '20210714'
# ]:
#     da_run = 'ASS_2407'  # ASS_newerer (old reference)
#     icon_run = 'MAIN_2405.3'  # MAIN_newerererRH8_new2mom-MP_RUC2.0
#     emvorado_run = 'EMVO_00510000.2'  # EMVO_dynwetgrow-BBold_sig-ice-25
#     icon_emvorado_run = icon_run + '/' + emvorado_run
#     for radar_loc in radar_locs:
#         qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                          icon_emvorado_run=icon_emvorado_run,
#                          spin_up_mm=spin_up_mm,
#                          elevation_deg=elevation_deg,
#                          overwrite=overwrite,
#                          radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# #  08.05.25  # R0E2  # Done 12.05.25
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12 # DONE!
# for elevation_deg in [12, 8, 17]:
#     for day in [
#         '20210714',
#         '20210713',
#     ]:
#         da_run = 'ASS_2411'  # ASS_newererer
#         icon_run = 'MAIN_2411.0'
#         for emvorado_run in [
#             'EMVO_00410000.2',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 21.02.25  # R2E3  # 10.03.25 Done
# # 22.01.26  # R2E3  new #
# radar_locs = list(rad_dict().keys())
# # radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12 # DONE!
# for elevation_deg in [12, 8, 17]:
#     for day in [
#         '20210714',
#         '20210713',
#     ]:
#         da_run = 'ASS_2411'  # ASS_newererer
#         icon_run = 'MAIN_2411.3'
#         for emvorado_run in [
#             'EMVO_00510000.2',
#             'EMVO_20510000.2',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 13.05.25  # Kobras runs  # Done!
# radar_locs = list(rad_dict().keys())
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12
# for day in [
#     '20210714',
#     '20210713',
# ]:
#     da_run = 'ASS_2411'  # ASS_newererer
#     icon_run = 'MAIN_2411.3'
#     for emvorado_run in [
#         'EMVO_00511300.2',
#         'EMVO_00521300.2',
#         'EMVO_01521300.2',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         for radar_loc in radar_locs:
#             qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                              icon_emvorado_run=icon_emvorado_run,
#                              spin_up_mm=spin_up_mm,
#                              elevation_deg=elevation_deg,
#                              overwrite=overwrite,
#                              radar_loc=radar_loc)
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 10.06.25  # R0E2  # 11.06.25 Done
# # radar_locs = list(rad_dict().keys())
# radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# # elevation_deg = 12 # DONE!
# for elevation_deg in [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]:
#     for day in [
#         '20210714',
#         # '20210713',
#     ]:
#         da_run = 'ASS_2411'  # ASS_newererer
#         icon_run = 'MAIN_2411.0'
#         for emvorado_run in [
#             'EMVO_00410000.2',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 10.06.25  # R2E3  # 11.06.25 Done
# # radar_locs = list(rad_dict().keys())
# radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# # elevation_deg = 12 # DONE!
# for elevation_deg in [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]:
#     for day in [
#         '20210714',
#         # '20210713',
#     ]:
#         da_run = 'ASS_2411'  # ASS_newererer
#         icon_run = 'MAIN_2411.3'
#         for emvorado_run in [
#             'EMVO_00510000.2',
#         ]:
#             icon_emvorado_run = icon_run + '/' + emvorado_run
#             for radar_loc in radar_locs:
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# # 19.12.25  # R2E3
# radar_locs = ['ESS']
# spin_up_mm = 120
# overwrite = '2025-02-13'
# elevation_deg = 12 # DONE!
# # for elevation_deg in [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]:
# for day in [
#     '20210714',
#     # '20210713',
# ]:
#     da_run = 'ASS_2411'
#     icon_run = 'MAIN_2411.3'
#     for emvorado_run in [
#         'EMVO_00510010.2',
#         'EMVO_00510020.2',
#         'EMVO_00510030.2',
#         'EMVO_00510040.2',
#         'EMVO_00510050.2',
#         'EMVO_00510060.2',
#         'EMVO_00510070.2',
#         'EMVO_00510000.2GF',
#     ]:
#         icon_emvorado_run = icon_run + '/' + emvorado_run
#         for radar_loc in radar_locs:
#             qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                              icon_emvorado_run=icon_emvorado_run,
#                              spin_up_mm=spin_up_mm,
#                              elevation_deg=elevation_deg,
#                              overwrite=overwrite,
#                              radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # 27.01.26  # R2E3  all new
# spin_up_mm = 120
# overwrite = '2025-02-13'
# radar_locs=['ESS', 'NHB', 'FLD', 'FBG'] + list(rad_dict().keys())
# for elevation_deg in [12, 8, 17]:
#     for radar_loc in radar_locs:
#         for day in [
#             '20210714',
#             '20210713',
#         ]:
#             da_run = 'ASS_2411'
#             icon_run = 'MAIN_2411.3'
#             for emvorado_run in [
#                 'EMVO_20510810.2',
#                 'EMVO_20510820.2',
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)



# --------------------------------------------------------------------------- #
# # 27.01.26  # R2E3  all new
# spin_up_mm = 120
# overwrite = '2025-02-13'
# radar_locs=['ESS', 'NHB', 'FLD', 'FBG'] + list(rad_dict().keys())
# for elevation_deg in [12, 8, 17]:
#     for radar_loc in radar_locs:
#         for day in [
#             '20210714',
#             '20210713',
#         ]:
#             da_run = 'ASS_2411'
#             icon_run = 'MAIN_2411.3'
#             for emvorado_run in [
#                 'EMVO_20510000.2',
#                 'EMVO_20510810.2',
#                 'EMVO_20510820.2',
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)
#
#             icon_run = 'MAIN_2411.0'
#             for emvorado_run in [
#                 'EMVO_20510000.2',
#                 'EMVO_20010000.2',
#                 'EMVO_20410000.2',
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)
#
#             icon_run = 'MAIN_2411.1'
#             for emvorado_run in [
#                 'EMVO_20510000.2',
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # 18.02.26  # R2E4,5
# spin_up_mm = 120
# overwrite = '2025-02-13'
# radar_locs=['ESS', 'FLD', 'NHB', 'FBG']
# for elevation_deg in [12]:
#     for radar_loc in radar_locs:
#         for day in [
#             '20210714',
#         ]:
#             da_run = 'ASS_2411'
#             icon_run = 'MAIN_2411.3'
#             for emvorado_run in [
#                 'EMVO_20510830.2',
#                 'EMVO_20510840.2',
#             ]:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)

# --------------------------------------------------------------------------- #
# # 18.02.26  # R2E4,5
# spin_up_mm = 120
# overwrite = '2025-02-13'
# radar_locs=['ESS', 'NHB', 'FLD', 'FBG'] + list(rad_dict().keys())
# for elevation_deg in [12, 8, 17]:
#     for day in [
#         '20210714',
#         '20210713',
#     ]:
#         da_run = 'ASS_2411'
#         icon_run = 'MAIN_2411.3'
#         for emvorado_run in [
#             'EMVO_20510830.2',
#             'EMVO_20510840.2',
#         ]:
#             for radar_loc in radar_locs:
#                 icon_emvorado_run = icon_run + '/' + emvorado_run
#                 qvp_from_syn_vol(day=day, da_run=da_run, icon_run=icon_run,
#                                  icon_emvorado_run=icon_emvorado_run,
#                                  spin_up_mm=spin_up_mm,
#                                  elevation_deg=elevation_deg,
#                                  overwrite=overwrite,
#                                  radar_loc=radar_loc)
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# 25.02.26  # R2E4,5 for ICON
spin_up_mm = 120
overwrite = '2025-02-13'
radar_locs=['ESS', 'NHB', 'FLD', 'FBG'] + list(rad_dict().keys())
for elevation_deg in [12, 8, 17]:
    for day in [
        '20210714',
        '20210713',
    ]:
        da_run = 'ASS_2411'
        icon_run = 'MAIN_2411.3'
        emvorado_run = 'EMVO_20510840.2'
        for radar_loc in radar_locs:
            icon_emvorado_run = icon_run + '/' + emvorado_run
            qvp_from_syn_vol_with_qn_qnx(
                day=day, da_run=da_run, icon_run=icon_run,
                icon_emvorado_run=icon_emvorado_run,
                spin_up_mm=spin_up_mm,
                elevation_deg=elevation_deg,
                overwrite=overwrite,
                radar_loc=radar_loc,
                qn_i=1e-0,
                q_i=1e-7,
                qn_c=1e-0,
                q_c=1e-7,
                qn_r=1e-0,
                q_r=1e-7,
                qn_s=1e-0,
                q_s=1e-7,
                qn_g=1e-0,
                q_g=1e-7,
                qn_h=1e-3,
                q_h=1e-5,
            )
# --------------------------------------------------------------------------- #