#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# PLOT_PPI.py                                                            #
#                                                                             #
# [...] Description here [...]                                                #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import xarray as xr
import datatree as dttree
import numpy as np
import matplotlib.pyplot as plt
import wradlib as wrl
import glob
import pandas as pd
import os
import matplotlib as mpl


# --------------------------------------------------------------------------- #

def plot_PPI(nc_file,
             ax=None,
             time_i=0,
             moment='KDP_NC',
             cmap=None,  # 'jet',
             levels=None,
             norm=None,
             title=None,
             ):
    # get colour conventions from header:
    if moment in ['KDP_NC', 'kdp']:
        if levels is None:
            levels = header.levels_kdp
        if norm is None:
            norm = header.norm_kdp
        if cmap is None:
            cmap = header.cmap_radar
    elif moment in ['rho', 'RHOHV']:
        if levels is None:
            levels = header.levels_rhohv
        if norm is None:
            norm = header.norm_rhohv
        if cmap is None:
            cmap = header.cmap_radar
    elif moment in ['DBZH', 'zh']:
        if levels is None:
            levels = header.levels_zh
        if norm is None:
            norm = header.norm_zh
        if cmap is None:
            cmap = header.cmap_radar
    elif moment in ['ZDR']:
        if levels is None:
            levels = header.levels_zdr
        if norm is None:
            norm = header.norm_zdr
        if cmap is None:
            cmap = header.cmap_radar
    elif moment in ['PHI_NC', 'phi_c', 'UPHIDP', 'PHI_C']:
        if levels is None:
            n_color = 14

            # step=1
            step = 2
            # step=3

            levels = np.arange(-step * 2, (n_color - 4) * step + 1, step)
            # levels = np.arange(-5, 22, 1)
            # levels = np.arange(-1, 26, 1)
            # levels = np.arange(-2, 24, 2)
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar
    else:
        if cmap is None:
            cmap = 'jet'

    sweep = nc_file.split('/')[-2]
    vol = dttree.open_datatree(nc_file)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    ppi = vol.isel(time=time_i)
    ppi = ppi.transpose('azimuth', 'range')
    img = wrl.georef.create_xarray_dataarray(data=ppi[moment],
                                             r=ppi.range.values / 1000,
                                             phi=ppi.azimuth.values,
                                             theta=ppi.elevation.values.item(),
                                             site=[ppi.longitude.values,
                                                   ppi.latitude.values,
                                                   ppi.altitude.values])
    img = img.wrl.georef.georeference()
    if ax is None:
        plt.figure(figsize=(5, 4))
        ax = plt.gca()

    img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm, extend='both')
    if title:
        plt.title(title)

    ax.set_xlabel("easting [km]")
    ax.set_ylabel("northing [km]")
