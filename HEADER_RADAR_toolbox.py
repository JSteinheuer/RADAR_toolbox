#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 08.01.24                                                 #
# HEADER_RADAR_toolbox.py                                                     #
#                                                                             #
# header for hard coded paths, colors, etc, ...                               #
# --------------------------------------------------------------------------- #

import getpass
import os
import sys
import matplotlib as mpl
import numpy as np

# --------------------------------------------------------------------------- #
# preamble necessary for wrl.georef.reproject: Tell the shell where to find   #
# the projection maps.                                                        #
path = sys.executable.split("/")[:-2]
path.extend(["share", "proj"])
path = "/".join(path)
os.environ["PROJ_LIB"] = path
os.environ["PROJ_NETWORK"] = 'ON'
# preamble ends.                                                              #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Where are you Julian?                                                       #
# --------------------------------------------------------------------------- #

if getpass.getuser() == 's6justei':  # Bonner tower/network
    dir_data_vol = '/automount/agradar/operation_hydrometeors/data/Syn_vol/'
    dir_data_qvp = '/automount/agradar/operation_hydrometeors/data/QVP/'
    dir_data_mod = '/automount/agradar/operation_hydrometeors/data/mod/'
    dir_data_era5 = '/automount/agradar/operation_hydrometeors/data/ERA5/'
    dir_projects = '/automount/user/s6justei/PyCharm/PyCharmProjects/'
    dir_data_obs = '/automount/agradar/operation_hydrometeors/data/obs/'
    dir_data_obs_realpep = '/automount/realpep/upload/RealPEP-SPP/DWD-CBand/'
    folder_plot = '/automount/agradar/operation_hydrometeors/plots/'
    folder_qvp_plot = '/automount/agradar/operation_hydrometeors/plots/QVPs/'
    folder_ppi_plot = '/automount/agradar/operation_hydrometeors/plots/PPIs/'


# elif getpass.getuser() == 'julian':  # personal notebook
#   [...]


# --------------------------------------------------------------------------- #
# Colors                                                                      #
# --------------------------------------------------------------------------- #

colors_radar = np.array(
    [[0.00, 1.00, 1.00], [0.00, 0.70, 0.93], [0.00, 0.00, 1.00],  # blues
     [0.50, 1.00, 0.00], [0.40, 0.80, 0.00], [0.27, 0.55, 0.00],  # greens
     [1.00, 1.00, 0.00], [0.80, 0.80, 0.00], [1.00, 0.65, 0.00],  # yellows
     [1.00, 0.27, 0.00], [0.80, 0.22, 0.00], [0.55, 0.15, 0.00],  # reds
     [1.00, 0.00, 1.00], [0.58, 0.44, 0.86]])  # pinks
cmap_radar = mpl.colors.ListedColormap(colors_radar)

# Zh
levels_zh = [-10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
norm_zh = mpl.colors.BoundaryNorm(levels_zh, len(levels_zh) - 1)

# ZDR
levels_zdr = [-1, -.1, 0, .1, .2, .3, .4, .5, .6, .8, 1, 2, 3]
norm_zdr = mpl.colors.BoundaryNorm(levels_zdr, len(levels_zdr) - 1)

# KDP
# levels_kdp = [-2, -1, 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9,
#               1., 1.2, 1.4, 1.7]
levels_kdp = [-.5, -.1, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 2, 3]
norm_kdp = mpl.colors.BoundaryNorm(levels_kdp, len(levels_kdp) - 1)

# RHOHV
levels_rhohv = [.7, .8, .85, .9, .92, .94, .95, .96, .97, .98, .99, .995, .998]
norm_rhohv = mpl.colors.BoundaryNorm(levels_rhohv, len(levels_rhohv) - 1)

# D0
levels_d0 = np.arange(0.2, 2.8, .2)
norm_d0 = mpl.colors.BoundaryNorm(levels_d0, len(levels_d0) - 1)

# --------------------------------------------------------------------------- #
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
# --------------------------------------------------------------------------- #




