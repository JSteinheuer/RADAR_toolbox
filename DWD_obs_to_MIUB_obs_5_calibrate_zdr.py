#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.01.24                                                 #
# DWD_obs_to_MIUB_obs_4_correct_kdp_in_ML.py                                  #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 4: correct KDP in ML.                                                  #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/build_radar_database/correct_rhohv.py      #
# --------------------------------------------------------------------------- #

import datatree as dttree
import numpy as np
import pandas as pd
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr
from scipy.ndimage import uniform_filter, gaussian_filter
import time
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# --------------------------------------------------------------------------- #

# V. Pejcic
def my_cmap(colmap='jet', csteps=10, setunder=False, ucolor='magenta'):
    """
    Get a colormap with discrete steps
    colmap ::: python colormap
    csteps ::: number of discrete colors
    """
    import matplotlib.cm as cm
    my_cmap = cm.get_cmap(colmap, csteps)
    if setunder:
        my_cmap.set_under(ucolor)

    return my_cmap


# V. Pejcic
def hist_2d(A,B, bins1=35, bins2=35, mini=1, maxi=None,
            cmap='jet', colsteps=30, alpha=1, fsize=15, colbar=True):
    """
    # Histogram 2d Quicklooks
    # ------------------------
    Plotting 2d Histogramm of two varibles
    # Input
    # -----
    A,B  ::: Variables
    bins1, bins2 ::: x, y bins
    mini, maxi  ::: min and max
    cmap  ::: colormap
    colsteps  ::: number of cmap steps
    alpha  ::: transperency
    fsize  ::: fontsize
    # Output
    # ------
    2D Histogramm Plot
    """
    from matplotlib.colors import LogNorm
    # discret cmap
    cmap = plt.cm.get_cmap(cmap, colsteps)
    # mask array
    m = ~np.isnan(A) & ~np.isnan(B)
    plt.hist2d(A[m], B[m], bins=(bins1, bins2), cmap=cmap,
               norm=LogNorm(vmin=mini, vmax=maxi), alpha=alpha)
    if colbar:
        cb = plt.colorbar(shrink=1, pad=0.01)
        # cb.set_label('number of samples', fontsize=fsize)
        cb.ax.set_title('#', fontsize=fsize)
        cb.ax.tick_params(labelsize=fsize)
        # cb.formatter = LogFormatterExponent(base=10) # 10 is the default
        # cb.update_ticks()

    plt.xticks(fontsize=fsize)
    plt.yticks(fontsize=fsize)


# V. Pejcic
def cal_zhzdr_lightrain(ZH, ZDR, plot=True):
    """
    ZH-ZDR Consistency in light rain
    AR p.155-156
    """
    zdr_zh_20 = np.nanmedian(ZDR[(ZH >= 19) & (ZH < 21)])
    zdr_zh_22 = np.nanmedian(ZDR[(ZH >= 21) & (ZH < 23)])
    zdr_zh_24 = np.nanmedian(ZDR[(ZH >= 23) & (ZH < 25)])
    zdr_zh_26 = np.nanmedian(ZDR[(ZH >= 25) & (ZH < 27)])
    zdr_zh_28 = np.nanmedian(ZDR[(ZH >= 27) & (ZH < 29)])
    zdr_zh_30 = np.nanmedian(ZDR[(ZH >= 29) & (ZH < 31)])
    zdroffset = np.nansum([zdr_zh_20-.23, zdr_zh_22-.27, zdr_zh_24-.32,
                           zdr_zh_26-.38, zdr_zh_28-.46, zdr_zh_30-.55])/6.
    nm = (ZH >= 19) & (ZH < 31) & (~np.isnan(ZH))
    if plot:
        plt.figure(figsize=(8, 3))

        plt.subplot(1, 2, 1)
        hist_2d(ZH, ZDR,
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        plt.plot([20, 22, 24, 26, 28, 30],
                 [.23, .27, .33, .40, .48, .56], color='black')
        plt.title('Non-calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)

        plt.subplot(1, 2, 2)
        hist_2d(ZH, ZDR-zdroffset,
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        plt.plot([20, 22, 24, 26, 28, 30],
                 [.23, .27, .33, .40, .48, .56], color='black')
        plt.title('Calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)
        plt.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                         'dB\n' + r'$N$: '+str(np.sum(nm)))
        plt.tight_layout()
        plt.show()

    return zdroffset, np.sum(nm)


# V. Pejcic
def cal_zhzdr_smalldrops(ZH, ZDR, band='S', plot=True):
    """
    Daniel zhzdr_for_small_drops ...
    """
    zdr_bar = {'X': 0.165, 'C': 0.183, 'S': 0.176}
    zdr_zh_1 = np.nanmedian(ZDR[(ZH >= 0) & (ZH < 2)])
    zdr_zh_3 = np.nanmedian(ZDR[(ZH >= 2) & (ZH < 4)])
    zdr_zh_5 = np.nanmedian(ZDR[(ZH >= 4) & (ZH < 6)])
    zdr_zh_7 = np.nanmedian(ZDR[(ZH >= 6) & (ZH < 8)])
    zdr_zh_9 = np.nanmedian(ZDR[(ZH >= 8) & (ZH < 10)])
    zdr_zh_11 = np.nanmedian(ZDR[(ZH >= 10) & (ZH < 12)])
    zdr_zh_13 = np.nanmedian(ZDR[(ZH >= 12) & (ZH < 14)])
    zdr_zh_15 = np.nanmedian(ZDR[(ZH >= 14) & (ZH < 16)])
    zdr_zh_17 = np.nanmedian(ZDR[(ZH >= 16) & (ZH < 18)])
    zdr_zh_19 = np.nanmedian(ZDR[(ZH >= 18) & (ZH < 20)])
    zdroffset = np.nansum([zdr_zh_1-zdr_bar[band],
    zdr_zh_3-zdr_bar[band],
    zdr_zh_5-zdr_bar[band],
    zdr_zh_7-zdr_bar[band],
    zdr_zh_9-zdr_bar[band],
    zdr_zh_11-zdr_bar[band],
    zdr_zh_13-zdr_bar[band],
    zdr_zh_15-zdr_bar[band],
    zdr_zh_17-zdr_bar[band],
    zdr_zh_19-zdr_bar[band], ]) / 10.
    nm = (ZH >= 0) & (ZH < 20)&(~np.isnan(ZH))
    if plot:
        plt.figure(figsize=(8, 3))

        plt.subplot(1, 2, 1)
        hist_2d(ZH.flatten(), ZDR.flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Non-calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)

        plt.subplot(1, 2, 2)
        hist_2d(ZH.flatten(), (ZDR-zdroffset).flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)
        plt.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset,3)) +
                         'dB\n' + r'$N$: '+str(np.sum(nm)))
        plt.tight_layout()
        plt.show()

    return zdroffset, np.sum(nm)


# V. Pejcic -> T. Scharbach
def cal_zhzdr_smalldrops2(ZH, ZDR, band='S', plot=True):
    """
    Daniel zhzdr_for_small_drops ...
    """
    zdr_bar = {'X': 0.165, 'C': 0.183, 'S': 0.176}
    zdr_zh_1 = np.nanmedian(ZDR[(ZH >= 0) & (ZH < 20)])
    zdroffset = np.nansum(zdr_zh_1-zdr_bar[band])
    nm = (ZH >= 0) & (ZH < 20) & (~np.isnan(ZH))
    if plot:
        plt.figure(figsize=(8, 3))

        plt.subplot(1, 2, 1)
        hist_2d(ZH.flatten(), ZDR.flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Non-calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)

        plt.subplot(1, 2, 2)
        hist_2d(ZH.flatten(), (ZDR-zdroffset).flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)
        plt.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                         'dB\n' + r'$N$: '+str(np.sum(nm)))
        plt.tight_layout()
        plt.show()
    return zdroffset, np.sum(nm)


# J. Steinheuer
def cal_zdr_lightrain(swp_cf, band='C', plot=True):
    """
    ZH-ZDR Consistency in light rain
    AR p.155-156
    """
    zdr_bar = {'S': [0.23, 0.27, 0.32, 0.38, 0.46, 0.55],
               'C': [0.23, 0.27, 0.33, 0.40, 0.48, 0.56],
               'X': [0.23, 0.28, 0.33, 0.41, 0.49, 0.58]}
    swp_mask = swp_cf.where((swp_cf.temp_beamtop > 4 + 273.15) &
                            (swp_cf.RHOHV > 0.98) &
                            np.isnan(swp_cf.CMAP))
    zdr_zh_20 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 19) &
                                            (swp_mask.DBZH < 21)).ZDR)
    zdr_zh_22 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 21) &
                                          (swp_mask.DBZH < 23)).ZDR)
    zdr_zh_24 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 23) &
                                          (swp_mask.DBZH < 25)).ZDR)
    zdr_zh_26 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 25) &
                                          (swp_mask.DBZH < 27)).ZDR)
    zdr_zh_28 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 27) &
                                          (swp_mask.DBZH < 29)).ZDR)
    zdr_zh_30 = np.nanmedian(swp_mask.where((swp_mask.DBZH >= 29) &
                                          (swp_mask.DBZH < 31)).ZDR)
    zdroffset = np.nansum(
        [zdr_zh_20 - zdr_bar[band][0], zdr_zh_22 - zdr_bar[band][1],
         zdr_zh_24 - zdr_bar[band][2], zdr_zh_26 - zdr_bar[band][3],
         zdr_zh_28 - zdr_bar[band][4], zdr_zh_30 - zdr_bar[band][5]]) / 6.
    nm = np.sum(((swp_mask.DBZH >= 19) & (swp_mask.DBZH < 31) &
                 (~np.isnan(swp_mask.DBZH)))).values.item()
    if plot:
        plt.figure(figsize=(8, 3))

        plt.subplot(1, 2, 1)
        hist_2d(swp_mask.where((swp_mask.DBZH >= 19) & (swp_mask.DBZH < 31)
                               ).DBZH.values.flatten(),
                swp_mask.where((swp_mask.DBZH >= 19) & (swp_mask.DBZH < 31)
                               ).ZDR.values.flatten(),
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        plt.plot([20, 22, 24, 26, 28, 30], zdr_bar[band], color='black')
        plt.title('Non-calibrated $Z_{DR}$ (used)')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)

        plt.subplot(1, 2, 2)
        hist_2d(swp_cf.DBZH.values.flatten(),
                swp_cf.ZDR.values.flatten() - zdroffset,
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        plt.plot([20, 22, 24, 26, 28, 30], zdr_bar[band], color='black')
        plt.title('Calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)
        plt.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                         'dB\n' + r'$N$: ' + str(nm))
        plt.tight_layout()
        plt.show()

    return zdroffset, nm


# J. Steinheuer
def cal_zdr_smalldrops(swp_cf, band='C', plot=True):
    """
    Daniel zhzdr_for_small_drops ...
    """
    zdr_bar = {'X': 0.165, 'C': 0.183, 'S': 0.176}
    swp_mask = swp_cf.where((swp_cf.DBZH > 0) &
                            (swp_cf.DBZH < 20) &
                            (swp_cf.RHOHV > 0.98) &
                            (swp_cf.temp_beamtop > 4 + 273.15) &
                            np.isnan(swp_cf.CMAP))
    zdr_zh_1 = np.nanmedian(swp_mask.ZDR)
    zdroffset = np.nansum(zdr_zh_1-zdr_bar[band])
    nm = np.sum(~np.isnan(swp_mask.ZDR.values))
    if plot:
        plt.figure(figsize=(8, 3))

        plt.subplot(1, 2, 1)
        hist_2d(swp_mask.DBZH.values.flatten(), swp_mask.ZDR.values.flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Non-calibrated $Z_{DR}$ (used)')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)

        plt.subplot(1, 2, 2)
        hist_2d(swp_cf.DBZH.values.flatten(),
                swp_cf.ZDR.values.flatten()-zdroffset,
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        plt.axhline(zdr_bar[band], color='black')
        plt.title('Calibrated $Z_{DR}$')
        plt.xlabel(r'$Z_H$', fontsize=15)
        plt.ylabel(r'$Z_{DR}$', fontsize=15)
        plt.grid(which='both', color='black', linestyle=':', alpha=0.5)
        plt.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                         'dB\n' + r'$N$: '+str(nm))

        plt.tight_layout()
        plt.show()

        return zdroffset, nm


# swp_cf=data.copy()
# # J. Steinheuer
# def zdr_at_zero_elev(swp_cf):
zdr = swp_cf.ZDR
el = swp_cf.elevation
if np.unique(el.values).size == 1:
    el = np.unique(el.values)
else:
    el = el  # TODO



# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
# --------------------------------------------------------------------------- #
DATES = [
    "20221222",  # case10"
    # "20210604",  # case01
    # "20210714",  # case09
    # "20210620", "20210621",  # case02
    # "20210628", "20210629",  # case03
    # "20220519", "20220520",  # case04
    # "20220623", "20220624", "20220625",  # case05
    # "20220626", "20220627", "20220628",  # case06+07
    # "20220630", "20220701",  # case08
    # "20170719",  # case_MA
]
LOCATIONS = [
    'oft',
    # 'umd',
    # 'ess',
    # 'pro',
    # 'tur',
    # 'asb', 'boo', 'drs', 'eis', 'fbg',
    # 'fld', 'hnr', 'isn', 'mem', 'neu',
    # 'nhb',
    # 'ros',
]
ELEVATIONS = ELEVATIONS_ALL.copy()
ELEVATIONS = [5.5]
MODE = [
        'pcp',
        # 'vol'
        ]
# --------------------------------------------------------------------------- #
merge = True
remove_parts = True
overwrite = False
# overwrite = True
parts = 6
print('Departing into: ' + str(parts))
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
# for date in DATES:
#     for location in LOCATIONS:
#         for mode in MODE:
#             for elevation_deg in ELEVATIONS:

date = DATES[0]
location = LOCATIONS[0]
mode = MODE[0]
elevation_deg = ELEVATIONS[0]
year = date[0:4]
mon = date[4:6]
day = date[6:8]
sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
if mode == 'pcp':
    if sweep != '00':
        print('break')
        # break

    print('\nstart: ' + date + ' ' + location + ' pcp')
else:
    print('\nstart: ' + date + ' ' + location + ' ' +
          str(elevation_deg))

# ------------------------------------------------------------------------#
# folder and file search
folder_plot = header.folder_ppi_plot
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
year = date[0:4]
mon = date[4:6]
day = date[6:8]
folder_in = "/".join([header.dir_data_obs + '*',
                      year, year + '-' + mon,
                      year + '-' + mon + '-' + day,
                      location, mode + '*', sweep])
nc_file_mom = glob.glob(folder_in + '/*allmoms*')
if len(nc_file_mom) > 1:
    print('mom: too many files')
elif len(nc_file_mom) == 0:
    print('mom: no files')
else:
    nc_file_mom = nc_file_mom[0]
nc_file_temp = glob.glob(folder_in + '/*temp*')
if len(nc_file_temp) > 1:
    print('temp: too many files')
elif len(nc_file_temp) == 0:
    print('temp: no files')
else:
    nc_file_temp = nc_file_temp[0]
nc_file_rho = glob.glob(folder_in + '/*rhohv_nc*')
if len(nc_file_rho) > 1:
    print('rho: too many files')
elif len(nc_file_rho) == 0:
    print('rho: no files')
    # nc_file_rho = nc_file_mom
else:
    nc_file_rho = nc_file_rho[0]
nc_file_kdp = glob.glob(folder_in + '/*kdp_nc*')
if len(nc_file_kdp) > 1:
    print('kdp: too many files')
    nc_file_kdp = nc_file_mom
elif len(nc_file_kdp) == 0:
    print('kdp: no files')
    # nc_file_kdp = nc_file_mom
else:
    nc_file_kdp = nc_file_kdp[0]
# ----------------------------------------------------------------------- #
# time
time_start = date + nc_file_mom.split('-')[-4][-4:]
time_end = date + nc_file_mom.split('-')[-3][-4:]
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
# ----------------------------------------------------------------------- #


data = dttree.open_datatree(nc_file_mom)[
    'sweep_' + str(int(sweep))].to_dataset()#.chunk('auto')
data_rho = dttree.open_datatree(nc_file_rho)[
    'sweep_' + str(int(sweep))].to_dataset()#.chunk('auto')
data_temp = dttree.open_datatree(nc_file_temp)[
    'sweep_' + str(int(sweep))].to_dataset()#.chunk('auto')
#
data_temp2 = data_temp.interp(
    coords=data.drop(['longitude','latitude',
                      'altitude','elevation']).coords, method='nearest')
#
data.RHOHV.values = data_rho.RHOHV_NC2P.values
data = data.assign({'SNRH': data_rho.SNRH})
data = data.assign({'temp_beamtop': data_temp2.temp_beamtop})
remo_var = list(data.data_vars.keys())
data = data.transpose('time', 'azimuth', 'range')
# data = data.sel(time=slice(dti_start, dti_end))

#TODO

# cal_zhzdr_lightrain(data.DBZH.values[0, :, :], data.ZDR.values[0, :, :],
#                     plot=True)
# cal_zhzdr_lightrain(data.DBZH.values[:, :, :], data.ZDR.values[:, :, :],
#                     plot=True)
#
#
# cal_zhzdr_smalldrops(data.DBZH.values[0, :, :], data.ZDR.values[0, :, :],
#                      band='C', plot=True)
# cal_zhzdr_smalldrops(data.DBZH.values[:, :, :], data.ZDR.values[:, :, :],
#                      band='C', plot=True)
#
# cal_zhzdr_smalldrops2(data.DBZH.values[0, :, :], data.ZDR.values[0, :, :],
#                      band='C', plot=True)
# cal_zhzdr_smalldrops2(data.DBZH.values[:, :, :], data.ZDR.values[:, :, :],
#                      band='C', plot=True)
#
cal_zdr_smalldrops(data, band='C', plot=True)
# cal_zhzdr_lightrain(data.DBZH.values[:, :, :], data.ZDR.values[:, :, :],
#                     plot=True)
cal_zdr_lightrain(data, band='C', plot=True)


# path_out = nc_file_mom.replace('_allmoms_', '_zdr_c_')
# mom_use = [x for x in list(data.keys())]
# for mom in mom_use:
#     data[mom].encoding["coordinates"] = \
#         "time azimuth range"
# data = data.drop_vars(remo_var)
# dtree = dttree.DataTree(name="root")
# dttree.DataTree(data, name=f"sweep_{int(sweep)}",
#                 parent=dtree)
# print('saving: ... ' + path_out.split('/')[-1] + ' ...')
# dtree.load().to_netcdf(path_out)
# data.close()
# print('saved: ' + path_out + ' !')


