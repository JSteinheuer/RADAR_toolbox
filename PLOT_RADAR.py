#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# PLOT_RADAR.py                                                               #
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
from pathlib import Path
import os
import matplotlib as mpl
from statsmodels.stats.weightstats import DescrStatsW


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
    if 'kdp' in moment.lower():
        if levels is None:
            levels = header.levels_kdp
        if norm is None:
            norm = header.norm_kdp
        if cmap is None:
            cmap = header.cmap_radar
    elif 'rho' in moment.lower():
        if levels is None:
            levels = header.levels_rhohv
        if norm is None:
            norm = header.norm_rhohv
        if cmap is None:
            cmap = header.cmap_radar
    elif 'zh' in moment.lower():
        if levels is None:
            levels = header.levels_zh
        if norm is None:
            norm = header.norm_zh
        if cmap is None:
            cmap = header.cmap_radar
    elif 'zdr' in moment.lower():
        if levels is None:
            levels = header.levels_zdr
        if norm is None:
            norm = header.norm_zdr
        if cmap is None:
            cmap = header.cmap_radar
    elif 'phi' in moment.lower():
        if levels is None:
            levels = header.levels_phi
        if norm is None:
            norm = header.norm_phi
        if cmap is None:
            cmap = header.cmap_radar_smooth
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
        #     cmap = header.cmap_radar
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

    theta = ppi.elevation.values
    if theta.size == 1:
        theta = theta.item()

    img = wrl.georef.create_xarray_dataarray(data=ppi[moment],
                                             r=ppi.range.values / 1000,
                                             phi=ppi.azimuth.values,
                                             theta=theta,
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


def plot_CFAD_or_CFTD_from_2_arrays(x, y,
                                    title=None,
                                    mom_min=0,
                                    mom_max=40,
                                    bins_mom=40,
                                    vert_temp=True,  # CFTD
                                    temp_min=-20,
                                    temp_max=16,
                                    bins_temp=18,
                                    height_min=0,  # in km
                                    height_max=10,  # in km
                                    bins_height=20,
                                    vmax=None,
                                    ax=None,
                                    save=False,
                                    save_path=header.folder_plot + 'CFADs/',
                                    save_name='test_CFAD'
                                    ):
    if vert_temp:
        y = y - 273.15
        y_min, y_max, bins_y = temp_min, temp_max, bins_temp
    else:
        y_min, y_max, bins_y = height_min, height_max, bins_height
        y = y * 1000

    x = x[y >= y_min]
    y = y[y >= y_min]
    x = x[y <= y_max]
    y = y[y <= y_max]
    if ax is None:
        plt.figure(figsize=(6, 5))

    i_sort = np.argsort(y)
    y = y[i_sort]
    x = x[i_sort]
    a = plt.hist(y, bins=bins_y, range=(y_min, y_max))
    weights = np.repeat(100 / a[0], np.int32(a[0]))

    # PLOT
    if vmax:
        extend = 'max'
    else:
        extend = 'neither'

    h2d, mom2d, y2d, fg = plt.hist2d(x, y, bins=[bins_mom, bins_y],
                                     range=[[mom_min, mom_max],
                                            [y_min, y_max]], vmax=vmax,
                                     weights=weights, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]', extend=extend)
    if vert_temp:
        plt.gca().invert_yaxis()
        plt.ylabel('temperature (°C)')
    else:
        plt.ylabel(r'height (km)')

    plt.title(title + '; n = ' +
              str(int(len(x[~np.isnan(x)])/1000)) + 'k')
    y_mid = y2d[1:] / 2 + y2d[:-1] / 2
    mom_mid = mom2d[1:] / 2 + mom2d[:-1] / 2
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid, weights=h2d[:, t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    plt.plot(quant_prof[0, ], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.2}$')
    plt.plot(quant_prof[1, ], y_mid, color='red', ls='solid',
             linewidth=2, label='$Q_{0.5}$')
    plt.plot(quant_prof[2, ], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.8}$')
    plt.plot(mean_prof, y_mid, color='orange', ls='solid',
             linewidth=2, label='$\mu$')
    plt.legend()
    if vmax and np.max(h2d) > vmax:
        if np.max(h2d) > 100:
            print('above 100% :' + str(np.max(h2d)))
        if ax:
            plt.text(1.04, 1.03, min(np.round(np.max(h2d), 1), 100.0),
                     transform=ax.transAxes)
        else:  # TODO: check!
            plt.text((1.04*mom_max-.04*mom_min), (-1.08*y_max+.08*y_min),
                     min(np.round(np.max(h2d), 1), 100.0))

    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)
        # plt.close()