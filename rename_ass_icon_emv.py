#!/usr/bin/env python3.11
###############################################################################
# Julian Steinheuer; December 2023                                            #
# rename_ass_icon_emv.py                                                      #
#                                                                             #
# Following the new convention from JM, meaningfull names for the old runs    #
# are given.                                                                  #
###############################################################################

import glob
import os

dir_data_in = '/automount/agradar/operation_hydrometeors/data/mod/'
dir_ass = glob.glob(dir_data_in + '*/*')
for folder in dir_ass:
    if 'ASS_vold' in folder:
        print(folder)
        os.rename(folder, folder.replace('ASS_vold', 'ASS_2109'))
    elif 'ASS_old' in folder:
        print(folder)
        os.rename(folder, folder.replace('ASS_old', 'ASS_2111'))
    elif 'ASS_new' in folder:
        print(folder)
        os.rename(folder, folder.replace('ASS_new', 'ASS_2211'))
    else:
        print(folder + ' not changed')

dir_data_in = '/automount/agradar/operation_hydrometeors/data/mod/'
dir_icon = glob.glob(dir_data_in + '*/*/*')
for folder in dir_icon:
    if 'MAIN_vold' in folder:
        print(folder)
        os.rename(folder, folder.replace('MAIN_vold', 'MAIN_2109.0'))
    elif 'MAIN_old' in folder:
        print(folder)
        os.rename(folder, folder.replace('MAIN_old', 'MAIN_2203.0'))
    elif 'MAIN_newer_new2mom-MP' in folder:  # TODO
        print(folder)
        os.rename(folder, folder.replace('MAIN_newer_new2mom-MP', 'MAIN_2308.1'))
    elif 'MAIN_newer' in folder:  # TODO
        print(folder)
        os.rename(folder, folder.replace('MAIN_newer', 'MAIN_2308.0'))
    elif 'MAIN_new' in folder:  # TODO
        print(folder)
        os.rename(folder, folder.replace('MAIN_new', 'MAIN_2211.0'))
    else:
        print(folder + ' not changed')

dir_data_in = '/automount/agradar/operation_hydrometeors/data/mod/'
dir_emvorado = glob.glob(dir_data_in + '*/*/*/*')
for folder in dir_emvorado:
    if 'EMVO_att-no_BBold.Xband' in folder:
        print(folder, folder.replace('EMVO_att-no_BBold.Xband', 'EMVO_00000000.X'))
        os.rename(folder, folder.replace('EMVO_att-no_BBold.Xband', 'EMVO_00000000.X'))
    elif 'EMVO_att-yes_BBold.Xband' in folder:
        print(folder, folder.replace('EMVO_att-yes_BBold.Xband', 'EMVO_10000000.X'))
        os.rename(folder, folder.replace('EMVO_att-yes_BBold.Xband', 'EMVO_10000000.X'))
    elif 'EMVO_BBold' in folder:
        print(folder, folder.replace('EMVO_BBold', 'EMVO_00000000.2'))
        os.rename(folder, folder.replace('EMVO_BBold', 'EMVO_00000000.2'))
    elif 'EMVO_att-no_BB-ML.Xband' in folder:
        print(folder, folder.replace('EMVO_att-no_BB-ML.Xband', 'EMVO_00200000.X'))
        os.rename(folder, folder.replace('EMVO_att-no_BB-ML.Xband', 'EMVO_00200000.X'))
    elif 'EMVO_att-yes_BB-ML.Xband' in folder:
        print(folder, folder.replace('MVO_att-yes_BB-ML.Xband', 'EMVO_10200000.X'))
        os.rename(folder, folder.replace('MVO_att-yes_BB-ML.Xband', 'EMVO_10200000.X'))
    elif 'EMVO_BB-ML' in folder:
        print(folder, folder.replace('EMVO_BB-ML', 'EMVO_00200000.2'))
        os.rename(folder, folder.replace('EMVO_BB-ML', 'EMVO_00200000.2'))
    elif 'EMVO_att-no_BBnew.Xband' in folder:
        print(folder, folder.replace('EMVO_att-no_BBnew.Xband', 'EMVO_00100000.X'))
        os.rename(folder, folder.replace('EMVO_att-no_BBnew.Xband', 'EMVO_00100000.X'))
    elif 'EMVO_att-yes_BBnew.Xband' in folder:
        print(folder, folder.replace('MVO_att-yes_BBnew.Xband', 'EMVO_10100000.X'))
        os.rename(folder, folder.replace('MVO_att-yes_BBnew.Xband', 'EMVO_10100000.X'))
    elif 'EMVO_no-melt.Xband' in folder:
        print(folder, folder.replace('EMVO_no-melt.Xband', 'EMVO_00300000.X'))
        os.rename(folder, folder.replace('EMVO_no-melt.Xband', 'EMVO_00300000.X'))
    elif 'EMVO_no-melt' in folder:
        print(folder, folder.replace('EMVO_no-melt', 'EMVO_00300000.2'))
        os.rename(folder, folder.replace('EMVO_no-melt', 'EMVO_00300000.2'))
    elif 'EMVO_BBnew2momMP_SSDB' in folder:  # TODO
        print(folder, folder.replace('EMVO_BBnew2momMP_SSDB', 'EMVO_00401000.2'))
        os.rename(folder, folder.replace('EMVO_BBnew2momMP_SSDB', 'EMVO_00401000.2'))
    elif 'EMVO_BBnew2momMP' in folder:  # TODO
        print(folder, folder.replace('EMVO_BBnew2momMP', 'EMVO_00400000.2'))
        os.rename(folder, folder.replace('EMVO_BBnew2momMP', 'EMVO_00400000.2'))
    elif 'EMVO_BBnew' in folder:  # TODO
        print(folder, folder.replace('EMVO_BBnew', 'EMVO_00100000.2'))
        os.rename(folder, folder.replace('EMVO_BBnew', 'EMVO_00100000.2'))
    else:
        print(folder + ' not changed')

dir_data_in = '/automount/agradar/operation_hydrometeors/data/mod/'
dir_RADOLAN = glob.glob(dir_data_in + '*/*/*/*/*/*/*')
for file in dir_RADOLAN:
    if file[-7:] == 'RADOLAN':
        print(file)
        os.rename(file, file + '.nc')


# TODO: maybe adjust main0*00******* folder to main0*00
