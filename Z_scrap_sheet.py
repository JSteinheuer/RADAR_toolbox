#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.01.24                                                 #
# Z_scrap_sheet.py                                                            #
#                                                                             #
# Make Tests, scrapping code pieces, etc ...                                  #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import datatree as dttree
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os
from osgeo import osr
import pandas as pd
from pathlib import Path
from scipy.ndimage import uniform_filter, gaussian_filter
import sys
import time
import warnings
warnings.filterwarnings("ignore")
import wradlib as wrl
from xhistogram.xarray import histogram
import xarray as xr
import xradar as xd


