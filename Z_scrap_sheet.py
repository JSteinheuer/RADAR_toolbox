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

# ###############################################################################
# # arround isel z ~5000
# index=29
# index=28
# nc_file = '/automount/data02/agradar/operation_hydrometeors/data/mod/20220519/ASS_2405/MAIN_2405.1/ICONdata/20220520120000/main0200/fc_R19B07.20220520140000_00010000.nc'
# data = xr.open_dataset(nc_file)
# nc_vgrid ='/automount/data02/agradar/operation_hydrometeors/data/mod/grid/vgrd_R19B07.RADOLAN.nc'
# data_v=xr.open_dataset(nc_vgrid)
# nc_hgrid ='/automount/data02/agradar/operation_hydrometeors/data/mod/grid/hgrd_R19B07.RADOLAN.nc'
# data_h=xr.open_dataset(nc_hgrid)
# data=data.assign({'clon':(('ncells'),np.rad2deg(data_h.clon.values[:]))},)
# data=data.assign({'clat':(('ncells'),np.rad2deg(data_h.clat.values[:]))},)
# data=data.assign({'z_mc':(('height','ncells'),data_v.z_mc.values[0,:])},)
# data=data.assign({'z_ifc':(('height_2','ncells'),data_v.z_ifc.values[0,:])},)
#
# data=data.where(data.clon<11.4)
# data=data.where(data.clon>6.2)
# data=data.where(data.clat>49.7)
# data=data.where(data.clat<52.9)
#
# # z_mc_level=5000
# # z_mc_thick_half=100
# # # data_plot=data.sel(z_mc=5000,method='nearest') # not woriking no coord
#
# # # is it height?
# # data_plot = data.where(data.z_mc < z_mc_level+z_mc_thick_half)
# # data_plot = data_plot.where(data_plot.z_mc > z_mc_level-z_mc_thick_half)
#
# # # is it height2?
# # data_plot = data.where(data.z_ifc < z_mc_level+z_mc_thick_half)
# # data_plot = data_plot.where(data_plot.z_ifc > z_mc_level-z_mc_thick_half)
#
#
# x=data.clon.values
# y=data.clat.values
# w=data.w.values[0,index,:]
# z=data.z_ifc.values[index,:]
#
# mask=~np.isnan(x) & ~np.isnan(y) & ~np.isnan(w)\
#      # & (z < z_mc_level+z_mc_thick_half) &\
#      # (z > z_mc_level-z_mc_thick_half)
# x=x[mask]
# y=y[mask]
# w=w[mask]
# # xr.plot.scatter(data_plot, x='clon', y='clat', hue='w' )
# # plt.scatter(x,y,c=w, s=1, marker ='v',vmin=-2,vmax=2)
# plt.scatter(x,y,c=w, s=1, marker ='v',vmin=-5,vmax=5,cmap='jet')
# plt.colorbar()
# plt.tight_layout
# plt.savefig(header.folder_plot + 'fc_w_at_zi='+ str(index)+'.pdf', format='pdf',
#             transparent=True)
# plt.close()


###############################################################################
# arround z=5000
nc_file = '/automount/data02/agradar/operation_hydrometeors/data/mod/20220519/ASS_2405/MAIN_2405.1/ICONdata/20220520120000/main0200/fc_R19B07.20220520140000_00010000.nc'
data = xr.open_dataset(nc_file)
nc_vgrid = '/automount/data02/agradar/operation_hydrometeors/data/mod/grid/vgrd_R19B07.RADOLAN.nc'
data_v = xr.open_dataset(nc_vgrid)
nc_hgrid = '/automount/data02/agradar/operation_hydrometeors/data/mod/grid/hgrd_R19B07.RADOLAN.nc'
data_h = xr.open_dataset(nc_hgrid)
data = data.assign({'clon': (('ncells'), np.rad2deg(data_h.clon.values[:]))}, )
data = data.assign({'clat': (('ncells'), np.rad2deg(data_h.clat.values[:]))}, )
data = data.assign(
    {'z_mc': (('height', 'ncells'), data_v.z_mc.values[0, :])}, )
data = data.assign(
    {'z_ifc': (('height_2', 'ncells'), data_v.z_ifc.values[0, :])}, )

data = data.where(data.clon < 11.4)
data = data.where(data.clon > 6.2)
data = data.where(data.clat > 49.7)
data = data.where(data.clat < 52.9)

z_mc_level = 5000
# z_mc_level=4000
# z_mc_level=3750
z_mc_thick_half = 100
# # data_plot=data.sel(z_mc=5000,method='nearest') # not woriking no coord

# # is it height?
# data_plot = data.where(data.z_mc < z_mc_level+z_mc_thick_half)
# data_plot = data_plot.where(data_plot.z_mc > z_mc_level-z_mc_thick_half)

# # is it height2?
# data_plot = data.where(data.z_ifc < z_mc_level+z_mc_thick_half)
# data_plot = data_plot.where(data_plot.z_ifc > z_mc_level-z_mc_thick_half)


x = np.repeat(data.clon.values, data.height_2.size)  # .ravel()
y = np.repeat(data.clat.values, data.height_2.size)  # .ravel()
w = data.w.values.transpose().ravel()
z = data.z_ifc.values.transpose().ravel()

mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(w) & \
       (z < z_mc_level + z_mc_thick_half) & \
       (z > z_mc_level - z_mc_thick_half)
x = x[mask]
y = y[mask]
w = w[mask]
# xr.plot.scatter(data_plot, x='clon', y='clat', hue='w' )
fig = plt.figure(figsize=(10, 8))
plt.scatter(x, y, c=w, s=3, marker='.', vmin=-5, vmax=5, cmap='jet')

# plt.scatter([8.25, 7.25], [51.5, 50.3], s=80, facecolors='none',)
plt.plot([8.25, 7.], [51.5, 50.4], 'o', markersize=40,  mfc='none', color='black')
plt.colorbar(label='w [m/s]', extend='both')
# plt.tight_layout()
plt.title('w at [' + str(z_mc_level - z_mc_thick_half) + ', ' +
          str(z_mc_level + z_mc_thick_half) + '] m ')
plt.ylabel('lat [°]')
plt.xlabel('lon [°]')
plt.savefig(header.folder_plot + 'fc_w_at' + str(z_mc_level) + 'm.pdf',
            format='pdf',
            transparent=True)
plt.savefig(header.folder_plot + 'fc_w_at' + str(z_mc_level) + 'm.png',
            format='png',
            transparent=True)
plt.close()

###############################################################################
