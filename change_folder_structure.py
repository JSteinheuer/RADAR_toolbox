import os
import glob

###############################################################################

# # day = '20170719'
# # day = '20170720'
# # day = '20170724'
# # day = '20170726'
# # day = '20170727'
# # day = '20180728'
# day = '20181202'
# # put 3h chuncks in .. /ASS_old/MAIN_old/ICONdata ...
#
# # cp and mv to emvorado
# dir_data = '/automount/agradar/operation_hydrometeors/data/mod/' + day + \
#            '/ASS_old/MAIN_old/'
# emvorado_run = 'EMVO_BBold'
# # os.system('mkdir ' + dir_data + emvorado_run)
# os.system(f"mkdir {dir_data}{emvorado_run}")
# # files = glob.glob(dir_data + 'ICONdata/**/*', recursive=True)
# files = sorted(glob.glob(f"{dir_data}ICONdata/**/*", recursive=True))
# for file in files:
#     if os.path.isfile(file):
#         if 'end_LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'NAMELIST_' in file or \
#                 'nml.atmo.log_main' in file or \
#                 'nml.emvorado.log_main' in file or \
#                 'qsubw_log_LL.icon_main' in file:
#             os.system('cp ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#         elif '/det/' in file:
#             os.system('mv ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#     else:  # directory
#         os.system('mkdir ' + file.replace('ICONdata', emvorado_run))
#
# # check!
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         # os.system('rm -r ' + file)
#
# # delete carfully:
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         os.system('rm -r ' + file)

###############################################################################

# # day = '20170725'
# day = '20180809'
# # put 3h chuncks in .. /ASS_vold/MAIN_vold/ICONdata ...
#
# # cp and mv to emvorado
# dir_data = '/automount/agradar/operation_hydrometeors/data/mod/' + day + \
#            '/ASS_vold/MAIN_vold/'
# emvorado_run = 'EMVO_BBold'
# # os.system('mkdir ' + dir_data + emvorado_run)
# os.system(f"mkdir {dir_data}{emvorado_run}")
# # files = glob.glob(dir_data + 'ICONdata/**/*', recursive=True)
# files = sorted(glob.glob(f"{dir_data}ICONdata/**/*", recursive=True))
# for file in files:
#     if os.path.isfile(file):
#         if 'end_LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'NAMELIST_' in file or \
#                 'nml.atmo.log_main' in file or \
#                 'nml.emvorado.log_main' in file or \
#                 'qsubw_log_LL.icon_main' in file:
#             os.system('cp ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run
#                                    ).replace('/radarout_DOM01', ''))
#         elif '/det/' in file:
#             os.system('mv ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run
#                                    ).replace('/radarout_DOM01', ''))
#     else:  # directory
#         os.system('mkdir ' + file.replace('ICONdata', emvorado_run
#                                           ).replace('/radarout_DOM01', ''))
#
# # check!
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         # os.system('rm -r ' + file)
#
# # delete carfully:
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         os.system('rm -r ' + file)

###############################################################################

# day = '20170810'
# # put 3h chuncks in .. /ASS_new/MAIN_new/ICONdata ...
#
# # cp and mv to emvorado
# dir_data = '/automount/agradar/operation_hydrometeors/data/mod/' + day + \
#            '/ASS_new/MAIN_new/'
# emvorado_run = 'EMVO_BBold'
# # os.system('mkdir ' + dir_data + emvorado_run)
# os.system(f"mkdir {dir_data}{emvorado_run}")
# # files = glob.glob(dir_data + 'ICONdata/**/*', recursive=True)
# files = sorted(glob.glob(f"{dir_data}ICONdata/**/*", recursive=True))
# for file in files:
#     if os.path.isfile(file):
#         if 'end_LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'NAMELIST_' in file or \
#                 'nml.atmo.log_main' in file or \
#                 'nml.emvorado.log_main' in file or \
#                 'qsubw_log_LL.icon_main' in file:
#             os.system('cp ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#         elif '/det/' in file:
#             os.system('mv ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#     else:  # directory
#         os.system('mkdir ' + file.replace('ICONdata', emvorado_run))
#
# # check!
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         # os.system('rm -r ' + file)
#
# # delete carfully:
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         os.system('rm -r ' + file)

###############################################################################

# day = '20170725'
# # put 3h chuncks in .. /ASS_new/MAIN_new/ICONdata ...
#
# # cp and mv to emvorado
# dir_data = '/automount/agradar/operation_hydrometeors/data/mod/' + day + \
#            '/ASS_new/MAIN_new/'
# emvorado_run = 'EMVO_BBold'
# # os.system('mkdir ' + dir_data + emvorado_run)
# os.system(f"mkdir {dir_data}{emvorado_run}")
# # files = glob.glob(dir_data + 'ICONdata/**/*', recursive=True)
# files = sorted(glob.glob(f"{dir_data}ICONdata/**/*", recursive=True))
# for file in files:
#     if os.path.isfile(file):
#         if 'end_LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'LL.icon_main' in file or \
#                 'icon_master.namelist_main' in file or \
#                 'NAMELIST_' in file or \
#                 'nml.atmo.log_main' in file or \
#                 'nml.emvorado.log_main' in file or \
#                 'qsubw_log_LL.icon_main' in file:
#             os.system('cp ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#         elif '/det/' in file:
#             os.system('mv ' + file + ' ' +
#                       file.replace('ICONdata', emvorado_run))
#     else:  # directory
#         os.system('mkdir ' + file.replace('ICONdata', emvorado_run))
#
# # check!
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         # os.system('rm -r ' + file)
#
# # delete carfully:
# for file in files:
#     if file[-4:] == '/det':
#         print(file)
#         print(glob.glob(file + '/**', recursive=True))
#         os.system('rm -r ' + file)

###############################################################################

# # add QVP elevation degree info
#
# dir_data = '/automount/agradar/operation_hydrometeors/data/QVP/**'
# files = glob.glob(dir_data, recursive=True)
# for file in files:
#     if 'QVP_Syn' in file:
#         print(file)
#         print(file.replace('QVP_Syn', 'QVP_12_Syn'))
#         os.system('mv ' + file + ' ' +
#                   file.replace('QVP_Syn', 'QVP_12_Syn'))

###############################################################################
