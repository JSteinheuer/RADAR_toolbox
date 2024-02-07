#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.01.24                                                 #
# Z_scrap_sheet.py                                                            #
#                                                                             #
# Make Tests, scrapping code pieces, etc ...                                  #
# --------------------------------------------------------------------------- #

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
from pathlib import Path
import os

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils
# import os


# folder = "/automount/agradar/operation_hydrometeors/data/obs/"
# files = sorted(glob.glob(folder + '*/*/*/*/*/*/*/*'))
# for file in files:
#     sweep_folder = file.split('/')[-2]
#     sweep_file = file.split('/')[-1].split('_')[-1][:2]
#     if 'rhohv_nc' in file:
#         print(file)
# #         file_new = '/'.join(file.split('/')[:-2]) + '/' + sweep_file + \
# #                    '/' + file.split('/')[-1]
# #         print(file_new)
# #         os.system('rm ' + file)


# sweep=4
# path_in="/automount/agradar/operation_hydrometeors/data/obs/" \
#         "OpHymet2-case10-20221222/2022/2022-12/2022-12-22/boo/" \
#         "vol5minng10/04/" \
#         "ras11-vol5minng10_sweeph5allm_allmoms_" \
#         "04-202212220002-202212221742-boo-10132.hd5"
# ddata = dttree.open_datatree(path_in)[
#                     'sweep_' + str(int(sweep))].to_dataset().chunk('auto')


# sweep=4
# path_in="/automount/agradar/operation_hydrometeors/data/obs/" \
#         "OpHymet2-case10-20221222/2022/2022-12/2022-12-22/boo/" \
#         "vol5minng10/04/" \
#         "ras11-vol5minng10_sweeph5allm_allmoms_" \
#         "04-202212220002-202212221742-boo-10132.hd5"
# ddata = dttree.open_datatree(path_in)[
#                     'sweep_' + str(int(sweep))].to_dataset().chunk('auto')


# folder = "/automount/agradar/operation_hydrometeors/data/obs/*"
# files = sorted(glob.glob(folder + '*/*/*/*/*/*/*/*'))
# for file in files:
#     if 'rhohv_nc' in file:
#         print(file)
#         os.system('rm ' + file)

# Re=6371000
# h=500
# alt=Re*h/(Re-h)

# import numpy as np
# import matplotlib.pyplot as plt
# import wradlib as wrl
# import warnings
# import cmweather
# import xarray as xr
# import datatree as dttree
#
# # warnings.filterwarnings("ignore")
# # try:
# #     get_ipython().run_line_magic("matplotlib inline")
# # except:
# #     plt.ion()

# import wradlib as wrl
# # import warnings
# # import cmweather
# import xarray as xr
# import datatree as dttree
# filename = '/automount/agradar/operation_hydrometeors/data/obs/2021/2021-07/2021-07-14/pro/vol5minng01/07/' \
#            'ras07-vol5minng01_sweeph5onem_allmoms_07-202107140003-202107142358-pro-10392.hd5'
#
# # data = dttree.open_datatree(filename)['sweep_7'].to_dataset().chunk('auto')
# vol = dttree.open_datatree(filename)['sweep_7'].to_dataset().chunk('auto')
# # vol = vol.rename(station_longitude='longitude')
# # vol = vol.rename(station_latitude='latitude')
# # vol = vol.rename(station_height='altitude')
# # ppi = vol.sel(elevation=12).isel(time=60)
# ppi = vol.isel(time=60)
# ppi = ppi.transpose('azimuth', 'range')
# img0 = wrl.georef.create_xarray_dataarray(ppi.ZDR)
# img = img0.wrl.georef.georeference()
# pm = img0.wrl.vis.plot()
#
# filename = '/automount/agradar/operation_hydrometeors/data/obs/' \
#            'OpHymet2-case03-20210628/2021/2021-06/2021-06-28/nhb/' \
#            'vol5minng10/09/ras11-vol5minng10_sweeph5allm_allmoms_09' \
#            '-202106281204-202106282359-nhb-10605.hd5'
#
# data = dttree.open_datatree(filename)[
#     'sweep_9'].to_dataset().chunk('auto')
#
#
# # ppi_ts = ppi_ts.rename(station_longitude='longitude')
# # ppi_ts = ppi_ts.rename(station_latitude='latitude')
# # ppi_ts = ppi_ts.rename(station_height='altitude')
# ppi = ppi_ts.isel(time=60)
# ppi = ppi.transpose('azimuth', 'range')
# img0 = wrl.georef.create_xarray_dataarray(ppi.zrsim)
# img = img0.wrl.georef.georeference()
# pm = img0.wrl.vis.plot()
#
# # s6toscha: 32,2 TiB (35.378.542.263.321)
# # op_hydro: 12,9 TiB (14.164.418.575.977)

# /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220519 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y/%Y-%m/%Y-%m-%d" -v

# /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220623 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case05-20220623 -d "%Y/%Y-%m/%Y-%m-%d" -v
# /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220626 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case06+07-20220626 -d "%Y/%Y-%m/%Y-%m-%d" -v
# /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220630 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case08-20220630 -d "%Y/%Y-%m/%Y-%m-%d" -v

# os.system('/home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220519 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y/%Y-%m/%Y-%m-%d" -v')
# os.system('/home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220623 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case05-20220623 -d "%Y/%Y-%m/%Y-%m-%d" -v')
# os.system('/home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220626 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case06+07-20220626 -d "%Y/%Y-%m/%Y-%m-%d" -v')
# os.system('/home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220630 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case08-20220630 -d "%Y/%Y-%m/%Y-%m-%d" -v')