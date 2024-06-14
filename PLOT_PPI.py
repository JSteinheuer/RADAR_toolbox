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
             range_max=None,
             ):
    sweep = nc_file.split('/')[-2]
    vol = dttree.open_datatree(nc_file)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    ppi = vol.isel(time=time_i)
    ppi = ppi.transpose('azimuth', 'range')

    if type(moment) == list:
        for moment_i in moment:
            if moment_i in ppi.variables:
                moment = moment_i
                break

    # get colour conventions from header:
    # if moment in ['KDP_NC', 'kdp']:
    if 'kdp' in moment.lower():
        if levels is None:
            levels = header.levels_kdp
        if norm is None:
            norm = header.norm_kdp
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['rho', 'RHOHV', 'RHOHV_NC', 'RHOHV_NC2P']:
    elif 'rho' in moment.lower():
        if levels is None:
            levels = header.levels_rhohv
        if norm is None:
            norm = header.norm_rhohv
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['DBZH', 'zh']:
    elif 'zh' in moment.lower():
        if levels is None:
            levels = header.levels_zh
        if norm is None:
            norm = header.norm_zh
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['ZDR']:
    elif 'zdr' in moment.lower():
        if levels is None:
            levels = header.levels_zdr
        if norm is None:
            norm = header.norm_zdr
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['PHI_NC', 'phi_c', 'UPHIDP', 'PHI_C']:
    elif 'phi' in moment.lower():
        if levels is None:
            levels = header.levels_phi
        if norm is None:
            norm = header.norm_phi
        if cmap is None:
            cmap = header.cmap_radar_smooth
        # if levels is None:
        #     n_color = 14
        #
        #     # step=1
        #     step = 2
        #     # step=3
        #
        #     levels = np.arange(-step * 2, (n_color - 4) * step + 1, step)
        #     # levels = np.arange(-5, 22, 1)
        #     # levels = np.arange(-1, 26, 1)
        #     # levels = np.arange(-2, 24, 2)
        # if norm is None:
        #     norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        # if cmap is None:
        #     cmap = header.cmap_radar
    elif 'vrad' in moment.lower():
        if levels is None:
            levels = np.arange(-20, 21, 0.1)
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar_smooth
    else:
        if cmap is None:
            cmap = 'jet'

    theta = ppi.elevation.values
    if theta.size == 1:
        theta = theta.item()

    img = wrl.georef.create_xarray_dataarray(data=ppi[moment],
                                             r=ppi.range.values / 1000,
                                             phi=ppi.azimuth.values,
                                             theta=theta,
                                             # theta=ppi.elevation.values.flatten(),#[:],#.item(),
                                             site=[ppi.longitude.values,
                                                   ppi.latitude.values,
                                                   ppi.altitude.values]
                                             )
    img = img.wrl.georef.georeference()
    if ax is None:
        plt.figure(figsize=(5, 4))
        ax = plt.gca()

    if range_max:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         extend='both',
                         ylim=[-range_max, range_max],
                         xlim=[-range_max, range_max])
    else:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         extend='both')

    if title:
        plt.title(title)

    ax.set_xlabel("easting [km]")
    ax.set_ylabel("northing [km]")


def plot_PPI_temp_ring(nc_file,
                       ax=None,
                       time_i=0,
                       temp=275.15,
                       temp_thickness=0.2,
                       moment='temp',
                       title=None,
                       range_max=None,
                       ):
    sweep = nc_file.split('/')[-2]
    vol = dttree.open_datatree(nc_file)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    ppi = vol.isel(time=time_i)
    ppi = ppi.transpose('azimuth', 'range')
    ppi_ring = ppi.where(abs(ppi[moment] - temp) < temp_thickness / 2)
    theta = ppi.elevation.values
    if theta.size == 1:
        theta = theta.item()

    img = wrl.georef.create_xarray_dataarray(data=ppi_ring[moment],
                                             r=ppi.range.values / 1000,
                                             phi=ppi.azimuth.values,
                                             theta=theta,
                                             # theta=ppi.elevation.values.flatten(),#.item(),
                                             site=[ppi.longitude.values,
                                                   ppi.latitude.values,
                                                   ppi.altitude.values]
                                             )
    img = img.wrl.georef.georeference()
    if ax is None:
        plt.figure(figsize=(5, 4))
        ax = plt.gca()

    colors = np.array(['white', 'black', 'black', 'white'])
    cmap = mpl.colors.ListedColormap(colors)
    levels = [-1000, -999, 9999, 10000]
    norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
    if range_max:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         ylim=[-range_max, range_max],
                         xlim=[-range_max, range_max],
                         add_colorbar=False)
    else:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         add_colorbar=False)
    if title:
        plt.title(title)

    ax.set_xlabel("easting [km]")
    ax.set_ylabel("northing [km]")
