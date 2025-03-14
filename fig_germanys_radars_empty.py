#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 20.02.25                                                 #
# fig_germanys_radars.py                                                      #
#                                                                             #
# [...] Description here [...]                                                #
# --------------------------------------------------------------------------- #

import cartopy.crs as ccrs
import cartopy
import gc
import glob
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from osgeo import osr
import pandas as pd
import warnings
import wradlib as wrl
import xarray as xr

warnings.filterwarnings('ignore')
xr.set_options(keep_attrs=True)


# --------------------------------------------------------------------------- #


def get_radar_locations():
    """
    Return DWD radars information

     Args:

    Returns:
        Dictionary of with 16 dictionaries with attributes to each DWD radar
        (due to 24.4.18).
    """

    radars = {}

    radar = {'name': 'ASR Borkum',
             'wmo': 10103,
             'lon': 6.748292,
             'lat': 53.564011,
             'alt': 36}
    radars['ASB'] = radar

    radar = {'name': 'Boostedt',
             'wmo': 10132,
             'lon': 10.046899,
             'lat': 54.004381,
             'alt': 125}
    radars['BOO'] = radar

    radar = {'name': 'Dresden',
             'wmo': 10488,
             'lon': 13.768639,
             'lat': 51.124639,
             'alt': 263}
    radars['DRS'] = radar

    radar = {'name': 'Eisberg',
             'wmo': 10780,
             'lon': 12.402788,
             'lat': 49.540667,
             'alt': 799}
    radars['EIS'] = radar

    radar = {'name': 'Essen',
             'wmo': 10410,
             'lon': 6.967111,
             'lat': 51.405649,
             'alt': 185}
    radars['ESS'] = radar

    radar = {'name': 'Feldberg',
             'wmo': 10908,
             'lon': 8.003611,
             'lat': 47.873611,
             'alt': 1516}
    radars['FBG'] = radar

    radar = {'name': 'Flechtdorf',
             'wmo': 10440,
             'lon': 8.801998,
             'lat': 51.311197,
             'alt': 628}
    radars['FLD'] = radar

    radar = {'name': 'Hannover',
             'wmo': 10339,
             'lon': 9.694533,
             'lat': 52.460083,
             'alt': 98}
    radars['HNR'] = radar

    radar = {'name': 'Isen',
             'wmo': 10873,
             'lon': 12.101779,
             'lat': 48.174705,
             'alt': 678}
    radars['ISN'] = radar

    radar = {'name': 'Memmingen',
             'wmo': 10950,
             'lon': 10.219222,
             'lat': 48.042145,
             'alt': 724}
    radars['MEM'] = radar

    radar = {'name': 'Neuhaus',
             'wmo': 10557,
             'lon': 11.135034,
             'lat': 50.500114,
             'alt': 880}
    radars['NEU'] = radar

    radar = {'name': 'Neuheilenbach',
             'wmo': 10605,
             'lon': 6.548328,
             'lat': 50.109656,
             'alt': 586}
    radars['NHB'] = radar

    radar = {'name': 'Offenthal',
             'wmo': 10629,
             'lon': 8.712933,
             'lat': 49.984745,
             'alt': 246}
    radars['OFT'] = radar

    radar = {'name': 'Proetzel',
             'wmo': 10392,
             'lon': 13.858212,
             'lat': 52.648667,
             'alt': 194}
    radars['PRO'] = radar

    radar = {'name': 'Rostock',
             'wmo': 10169,
             'lon': 12.058076,
             'lat': 54.175660,
             'alt': 37}
    radars['ROS'] = radar

    radar = {'name': 'Tuerkheim',
             'wmo': 10832,
             'lon': 9.782675,
             'lat': 48.585379,
             'alt': 768}
    radars['TUR'] = radar

    radar = {'name': 'Ummendorf',
             'wmo': 10356,
             'lon': 11.176091,
             'lat': 52.160096,
             'alt': 185}
    radars['UMM'] = radar

    return radars


# ------------------------------------ #
# Load RW of two flooding days
folder_rw = header.folder_radolan + 'RW202107/'
rw_files = glob.glob(folder_rw + '*')
for rw_file in rw_files:
    if (int(rw_file.split('/')[-1].split('-')[2]) > 2107150000) or \
            (int(rw_file.split('/')[-1].split('-')[2]) < 2107130000):
        print(rw_file)
        rw_files.remove(rw_file)

ds_m = wrl.io.open_radolan_mfdataset(rw_files)
ds_sum = ds_m.sum(dim='time')
ds_count = ds_m.count(dim='time')
ds_sum = xr.where(ds_count > 99999990, ds_sum, np.nan)

# ------------------------------------ #
# map with catropy
# KM: RADOLAN Stereographische Projection see:
# https://docs.wradlib.org/en/2.2.0/notebooks/fileio/radolan/radolan_grid.html#Polar-Stereographic-Projection
if ds_m.formatversion < 5:
    map_proj = ccrs.Stereographic(
        true_scale_latitude=60.0,
        central_latitude=90.0,
        central_longitude=10.0,
        globe=ccrs.Globe(ellipse="sphere")  # for RW formatversion < 5

    )
else:
    map_proj = ccrs.Stereographic(
        true_scale_latitude=60.0,
        central_latitude=90.0,
        central_longitude=10.0,
        globe=ccrs.Globe(ellipse="WGS84")  # for RW formatversion >= 5
    )

# # ------------------------------------ #
# # create radolan projection object
# proj_stereo = wrl.georef.create_osr("dwd-radolan")  # TODO: not used
# # create wgs84 projection object
# proj_wgs = osr.SpatialReference()  # TODO: used
# proj_wgs.ImportFromEPSG(4326)  # TODO

# --------------------------------------------------------------------------- #
# plot RW
cmap = mpl.cm.plasma
bounds = np.concatenate(
    (np.arange(0, 50, 10),
     np.arange(50, 225, 25),))
norm = mpl.colors.BoundaryNorm(bounds, len(bounds), extend='max')
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection=map_proj)
plot = ds_sum.RW.plot(ax=ax,
                      cmap=cmap, norm=norm)
plot.colorbar.set_ticks(bounds)
plot.colorbar.set_label('rain [mm]')
plt.title('polarimetric C-band radar network of DWD')

# plot borders etc
#ax = plt.gca()
ax.add_feature(cartopy.feature.OCEAN, facecolor="lightblue")
ax.add_feature(cartopy.feature.LAND)
ax.coastlines(color='cyan')
ax.add_feature(cartopy.feature.BORDERS, color='orange')
gl = ax.gridlines(draw_labels=True, y_inline=False)
gl.xlabels_top = False
gl.ylabels_right = False
xlim = ax.set_xlim()
ylim = ax.set_ylim()

# plot radar rings
radars = get_radar_locations()
for RRR in radars:
    radar = radars[RRR]
    site = (radar['lon'], radar['lat'], radar['alt'])
    map_trans = ccrs.AzimuthalEquidistant(
        central_latitude=radar['lat'],
        central_longitude=radar['lon'],
    )
    r = np.arange(1, 151) * 1000
    az = np.linspace(0, 359, 360)
    polygons, rad_proj = wrl.georef.spherical_to_polyvert(
        r, az, 0, site)
    polygons = polygons[..., 0:2]  # not heights
    polygons.shape = (len(az), len(r), 5, 2)  # az x ra x 5 x lo/la
    polygons = polygons[:, -1, :, :]  # only last range
    polycoll = mpl.collections.PolyCollection(
        polygons, closed=True, edgecolors='lightgray',
        transform=map_trans
    )
    ax.add_collection(polycoll, autolim=True)

# plot radars
for RRR in radars:
    radar = radars[RRR]
    site = (radar['lon'], radar['lat'], radar['alt'])
    plt.scatter(radar['lon'], radar['lat'],
                marker='H', s=40, color='grey',
                transform=ccrs.PlateCarree())
    plt.annotate(RRR, (radar['lon'], radar['lat']),
                 color='white', fontsize=12,
                 bbox=dict(facecolor='black', alpha=.25),
                 transform=ccrs.PlateCarree())

ax.set_xlim(xlim)
ax.set_ylim(ylim)

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

plt.tight_layout()
plt.savefig(
    header.folder_plot + 'Radar_network_Germany_empty.pdf', format='pdf',
    transparent=True)
plt.savefig(
    header.folder_plot + 'Radar_network_Germany_empty.png', format='png',
    transparent=True)
# plt.close()
