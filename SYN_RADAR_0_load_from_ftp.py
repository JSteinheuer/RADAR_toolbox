#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 13.06.24                                                 #
# SYN_RADAR_0_load_from_ftp.py                                                #
#                                                                             #
# Untar files from JM from ftp to                                             #
# /automount/data02/agradar/operation_hydrometeors/data/mod/                  #
# --------------------------------------------------------------------------- #


import os
import glob
import HEADER_RADAR_toolbox as header

folder_ftp = '/automount/ftp/wwwgast/spp-prom/'
folder_mod = header.dir_data_mod

# /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/SYN_RADAR_0_load_from_ftp.py

# --------------------------------------------------------------------------- #
folder_20240613 = folder_ftp + 'BACYdata/OH2/'
files = sorted(glob.glob(folder_20240613 + '*'))
for file in files:
    if file[-4:] == '.tgz':
        print('tar zxvf ' + file + ' -C ' + folder_mod)
        os.system('tar zxvf ' + file + ' -C ' + folder_mod)

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# folder_20240613_b = folder_ftp + 'BACYdata/OH1/'
# files = sorted(glob.glob(folder_20240613_b + '*'))
# for file in files:
#     if file[-4:] == '.tgz':
#         print('tar zxvf ' + file + ' -C ' + folder_mod)
#         os.system('tar zxvf ' + file + ' -C ' + folder_mod)

# --------------------------------------------------------------------------- #

