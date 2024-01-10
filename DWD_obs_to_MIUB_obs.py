#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 10.01.24                                                 #
# DWD_obs_to_MIUB_obs.py                                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header

import sys
sys.path.insert(0, header.dir_projects + 'RADAR_toolbox/radar_processing_scripts')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
import radar_processing_scripts.utils as utils


