# !/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 22.05.24                                                 #
# DWD_obs_to_MIUB_obs_5_calibrate_zdr.py                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 5: calibrate ZDR.                                                      #
#         Adapted from Velibor Pejcic                                         #
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
from matplotlib.colors import LogNorm
warnings.filterwarnings("ignore")


# V. Pejcic
def cal_zhzdr_lightrain(ZH, ZDR, plot=[True, True], axes=None):
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
    # valid for S-Band !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zdroffset = np.nansum([zdr_zh_20 - .23, zdr_zh_22 - .27, zdr_zh_24 - .32,
                           zdr_zh_26 - .38, zdr_zh_28 - .46,
                           zdr_zh_30 - .55]) / 6.
    # valid for S-Band !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nm = (ZH >= 19) & (ZH < 31) & (~np.isnan(ZH))
    if plot[0]:
        if axes is None:
            plt.figure(figsize=(8, 3))
            ax = plt.subplot(1, 2, 1)
        else:
            ax = axes[0]

        hist_2d(ZH, ZDR, ax=ax,
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        ax.plot([20, 22, 24, 26, 28, 30],
                [.23, .27, .33, .40, .48, .56], color='black')
        ax.set_title('Non-calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)

    if plot[1]:
        if axes is None:
            ax = plt.subplot(1, 2, 2)
        else:
            ax = axes[1]

        hist_2d(ZH, ZDR - zdroffset, ax=ax,
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        ax.plot([20, 22, 24, 26, 28, 30],
                [.23, .27, .33, .40, .48, .56], color='black')
        ax.set_title('Calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(np.sum(nm)))

    if sum(plot) > 0:
        plt.tight_layout()
        plt.show()

    return zdroffset, np.sum(nm)


# V. Pejcic
def cal_zhzdr_smalldrops(ZH, ZDR, band='S', plot=[True, True], axes=None):
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
    zdroffset = np.nansum([zdr_zh_1 - zdr_bar[band],
                           zdr_zh_3 - zdr_bar[band],
                           zdr_zh_5 - zdr_bar[band],
                           zdr_zh_7 - zdr_bar[band],
                           zdr_zh_9 - zdr_bar[band],
                           zdr_zh_11 - zdr_bar[band],
                           zdr_zh_13 - zdr_bar[band],
                           zdr_zh_15 - zdr_bar[band],
                           zdr_zh_17 - zdr_bar[band],
                           zdr_zh_19 - zdr_bar[band], ]) / 10.
    nm = (ZH >= 0) & (ZH < 20) & (~np.isnan(ZH))
    if plot[0]:
        if axes is None:
            plt.figure(figsize=(8, 3))
            ax = plt.subplot(1, 2, 1)
        else:
            ax = axes[0]

        hist_2d(ZH.flatten(), ZDR.flatten(), ax=ax,
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        ax.set_title('Non-calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)

    if plot[1]:
        if axes is None:
            ax = plt.subplot(1, 2, 2)
        else:
            ax = axes[1]

        hist_2d(ZH.flatten(), (ZDR - zdroffset).flatten(), ax=ax,
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        ax.set_title('Calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(np.sum(nm)))

    if sum(plot) > 0:
        plt.tight_layout()
        plt.show()

    return zdroffset, np.sum(nm)


# V. Pejcic -> T. Scharbach
def cal_zhzdr_smalldrops2(ZH, ZDR, band='S', plot=[True, True], axes=None):
    """
    Daniel zhzdr_for_small_drops ...
    """
    zdr_bar = {'X': 0.165, 'C': 0.183, 'S': 0.176}
    zdr_zh_1 = np.nanmedian(ZDR[(ZH >= 0) & (ZH < 20)])
    zdroffset = np.nansum(zdr_zh_1 - zdr_bar[band])
    nm = (ZH >= 0) & (ZH < 20) & (~np.isnan(ZH))
    if plot[0]:
        if axes is None:
            plt.figure(figsize=(8, 3))
            ax = plt.subplot(1, 2, 1)
        else:
            ax = axes[0]

        hist_2d(ZH.flatten(), ZDR.flatten(), ax=ax,
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        ax.set_title('Non-calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)

    if plot[1]:
        if axes is None:
            ax = plt.subplot(1, 2, 2)
        else:
            ax = axes[1]

        hist_2d(ZH.flatten(), (ZDR - zdroffset).flatten(),
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        ax.set_title('Calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(np.sum(nm)))

    if sum(plot) > 0:
        plt.tight_layout()
        plt.show()

    return zdroffset, np.sum(nm)


# V. Pejcic
def hist_2d(A, B, bins1=35, bins2=35, mini=1, maxi=None, ax=None, cmap='jet',
            alpha=1, fsize=15, colorbar=True, density=True):
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
    m = ~np.isnan(A) & ~np.isnan(B)
    if density:
        norm = LogNorm(vmin=0.000001, vmax=1.01)
    else:
        norm = LogNorm(vmin=mini, vmax=maxi)

    if ax:
        h = ax.hist2d(A[m], B[m], bins=(bins1, bins2), cmap=cmap,
                      density=density, norm=norm, alpha=alpha)
    else:
        h = plt.hist2d(A[m], B[m], bins=(bins1, bins2), cmap=cmap,
                       density=density, norm=norm, alpha=alpha)

    if colorbar:
        cb = plt.colorbar(h[3], shrink=1, pad=0.01, ax=ax)
        if not density:
            cb.ax.set_title('#', fontsize=fsize)

        cb.ax.tick_params(labelsize=fsize)


# J. Steinheuer
def cal_zdr_lightrain(swp_cf, band='C', plot=[True, True],
                      colorbar=[True, True], axes=None):
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
    if plot[0]:
        if axes is None:
            plt.figure(figsize=(8, 3))
            ax = plt.subplot(1, 2, 1)
        else:
            ax = axes[0]

        hist_2d(swp_mask.where((swp_mask.DBZH >= 19) & (swp_mask.DBZH < 31)
                               ).DBZH.values.flatten(),
                swp_mask.where((swp_mask.DBZH >= 19) & (swp_mask.DBZH < 31)
                               ).ZDR.values.flatten() - zdroffset,
                ax=ax,
                colorbar=colorbar[0],
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        ax.plot([20, 22, 24, 26, 28, 30], zdr_bar[band], color='black')
        # ax.set_title('Calibrated $Z_{DR}$ (used)')
        ax.set_xlabel(r'$Z_H [dBZ]$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR} [dB]$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm), loc='upper left')

    if plot[1]:
        if axes is None:
            ax = plt.subplot(1, 2, 2)
        else:
            ax = axes[1]

        hist_2d(swp_cf.DBZH.values.flatten(),
                swp_cf.ZDR.values.flatten() - zdroffset,
                ax=ax,
                colorbar=colorbar[1],
                bins1=np.arange(0, 40, 1),
                bins2=np.arange(-1, 3, .1))
        ax.plot([20, 22, 24, 26, 28, 30], zdr_bar[band], color='black')
        # ax.set_title('Calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm), loc='upper left')

    if sum(plot) > 0:
        plt.tight_layout()
        # plt.show()

    return zdroffset, nm


# J. Steinheuer
def cal_zdr_smalldrops(swp_cf, band='C', plot=[True, True],
                       colorbar=[True, True], axes=None):
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
    zdroffset = np.nansum(zdr_zh_1 - zdr_bar[band])
    nm = np.sum(~np.isnan(swp_mask.ZDR.values))
    if plot[0]:
        if axes is None:
            plt.figure(figsize=(8, 3))
            ax = plt.subplot(1, 2, 1)
        else:
            ax = axes[0]

        hist_2d(swp_mask.DBZH.values.flatten(),
                swp_mask.ZDR.values.flatten() - zdroffset,
                ax=ax,
                colorbar=colorbar[0],
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        # ax.set_title('Calibrated $Z_{DR}$ (used)')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm))

    if plot[1]:
        if axes is None:
            ax = plt.subplot(1, 2, 2)
        else:
            ax = axes[1]

        hist_2d(swp_cf.DBZH.values.flatten(),
                swp_cf.ZDR.values.flatten() - zdroffset,
                ax=ax,
                colorbar=colorbar[1],
                bins1=np.arange(0, 40, 1), bins2=np.arange(-1, 3, .1))
        ax.axhline(zdr_bar[band], color='black')
        # ax.set_title('Calibrated $Z_{DR}$')
        ax.set_xlabel(r'$Z_H$', fontsize=15)
        ax.set_ylabel(r'$Z_{DR}$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm))

    if sum(plot) > 0:
        plt.tight_layout()
        # plt.show()

    return zdroffset, nm


# J. Steinheuer
def zdr_at_zero_elev(swp_cf):
    """
    rearange Eq. (8) from Sanchez-Rivas, Rico-Ramirez 2022
    (resp. from Chandrasekar(2001)).
    """
    zdr = swp_cf.ZDR
    el = swp_cf.elevation
    zdr_lin = 10 ** (zdr / 10)
    el_rad = el * np.pi / 180
    zdr_lin_0 = (zdr_lin * np.cos(el_rad) ** 4) / \
                (zdr_lin * np.sin(el_rad) ** 4 -
                 2 * zdr_lin ** (1 / 2) * np.sin(el_rad) ** 2 + 1)
    # zdr_lin_02 = (zdr_lin * np.cos(el_rad)**4) / \
    #              (zdr_lin * np.sin(el_rad)**4 +
    #               2 * zdr_lin**(1 / 2) * np.sin(el_rad)**2 + 1)
    zdr_0 = 10 * np.log10(zdr_lin_0)
    zdr_0.attrs["long_name"] = 'equivalent log differential reflectivity at 0°'
    zdr_0.attrs["short_name"] = 'equivalent ZDR 0°'
    zdr_0.attrs["units"] = 'dB'
    swp_cf = swp_cf.assign(ZDR0=zdr_0)
    return swp_cf


# J. Steinheuer
def cal_zdr_birdbath(swp_cf, plot=True, ax=None):
    """
    ...
    """
    print(1)
    swp_cf = swp_cf.chunk(chunks=-1)
    print(1)
    zh=swp_cf.DBZH.values.flatten()
    print(1)
    zdr=swp_cf.ZDR.values.flatten()
    print(2)
    rho=swp_cf.RHOHV.values.flatten()
    print(3)
    cmap=swp_cf.CMAP.values.flatten()
    print(4)
    mask=~np.isnan(zdr)
    print(5)
    zdr=zdr[mask]
    print(6)
    zh=zh[mask]
    print(7)
    rho=rho[mask]
    print(8)
    cmap=cmap[mask]
    print(9)
    zdr=zdr[(zh > 0) & (zh < 40) & (rho > 0.98) & (np.isnan(cmap))]
    print(10)
    zdroffset = np.nanmedian(zdr)
    print(11)
    zdroffset_sd = np.nanstd(zdr)
    print(12)
    nm = np.sum(~np.isnan(zdr))


    # print(0)
    # # swp_mask=swp_cf.chunk(chunks=1).load()
    # print(-3)
    # # swp_mask.load()
    # # swp_mask = swp_cf.where((swp_cf.DBZH > 0) &
    # #                         (swp_cf.DBZH < 40) &
    # #                         (swp_cf.RHOHV > 0.98) &
    # #                         np.isnan(swp_cf.CMAP))
    # print(1)
    # # swp_mask = swp_cf.where(swp_mask.DBZH > 0)
    # swp_mask = swp_cf.where(swp_cf.DBZH > 0)
    # print(2)
    # swp_mask = swp_mask.where(swp_mask.DBZH < 40)
    # print(3)
    # swp_mask = swp_mask.where(swp_mask.RHOHV > 0.98)
    # print(4)
    # swp_mask = swp_mask.where(np.isnan(swp_mask.CMAP))
    # print(5)
    # # zdr=swp_mask.ZDR.data.flatten()
    # # print(-1)
    # # zdr=zdr[~np.isnan(zdr)]
    # # print(0)
    # # zdroffset = np.nanmedian(zdr)
    # # print(1)
    # # zdroffset_sd = np.nanstd(zdr)
    # # print(2)
    # # nm = len(zdr)
    # # pr1nt(00)
    # zdroffset = np.nanmedian(swp_mask.ZDR)
    # print(11)
    # zdroffset_sd = np.nanstd(swp_mask.ZDR)
    # print(22)
    # nm = np.sum(~np.isnan(swp_mask.ZDR.values))
    if plot:
        if ax is None:
            plt.figure(figsize=(4, 3))
            ax = plt.subplot(1, 1, 1)
        # ax.hist(swp_cf.ZDR.values.flatten(),
        #          bins=np.arange(-1, 3, .025), density=True, alpha=.5)
        # ax.hist(swp_mask.ZDR.values.flatten(),
        #         bins=np.arange(-1, 3, .025))
        ax.hist(zdr,
                bins=np.arange(-1, 3, .025))
        ax.axvline(zdroffset, color='black')
        ax.set_title('Birdbath $90°$')
        ax.set_xlabel(r'$Z_{DR}$', fontsize=15)
        ax.set_ylabel(r'$counts$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm) + '\n sd: ' +
                        str(np.round(zdroffset_sd, 3)) + 'dB')
        plt.tight_layout()
        # plt.show()

    return zdroffset, nm, zdroffset_sd


# J. Steinheuer
def calibrate_zdr(date, location, elevation_deg=5.5, mode='pcp',
                  overwrite=False):
    """
    calibrate zdr.
    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    elevation_degs : elevations in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    modes : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool;, if *zdr_off*-output exists, it can be
                       overwritten
    """
    print(mode + ' ' +
          (str(elevation_deg) if (mode == 'vol') else '') +
          ' ' + location + ' ' + date)
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    folder_in = "/".join([header.dir_data_obs + '*', year,
                          year + '-' + mon,
                          year + '-' + mon + '-' + day,
                          location, mode + '*', sweep])
    nc_file_mom = glob.glob(folder_in + '/*allmoms*')
    if len(nc_file_mom) > 1:
        print('mom: too many files')
        return
    elif len(nc_file_mom) == 0:
        print('mom: no files')
        return
    else:
        nc_file_mom = nc_file_mom[0]

    path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
    if overwrite or not os.path.exists(path_out_nc) or plot:
        nc_file_temp = glob.glob(folder_in + '/*temp*')
        if len(nc_file_temp) > 1:
            print('temp: too many files')
            if mode != '90grad':
                return

        elif len(nc_file_temp) == 0:
            print('temp: no files')
            if mode != '90grad':
                return

        else:
            nc_file_temp = nc_file_temp[0]

        nc_file_rho = glob.glob(folder_in + '/*rhohv_nc*')
        if len(nc_file_rho) > 1:
            print('rho: too many files')
            if mode != '90grad':
                return

        elif len(nc_file_rho) == 0:
            print('rho: no files')
            if mode != '90grad':
                return

        else:
            nc_file_rho = nc_file_rho[0]

        # --------------------------------------------------------------- #
        # vol/pcp or bird bad                                             #
        # --------------------------------------------------------------- #
        if mode == '90grad':
            data = dttree.open_datatree(nc_file_mom)[
                'sweep_' + str(int(sweep))].to_dataset()
            print('go')
            bb_off, bb_nm, bb_sd = cal_zdr_birdbath(data, plot=[False, False],
                                                    ax=None)
            path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
            if not overwrite and os.path.exists(path_out_nc):
                print('exists: ' + path_out_nc + ' -> continue')
            else:
                remo_var = list(data.data_vars.keys())
                remo_var.remove('ZDR')
                data = data.drop_vars(remo_var)
                data['zdr_off_bb'] = bb_off
                data['zdr_off_bb'].attrs["long_name"] = 'ZDR offset ' + \
                                                        'from bird bath'
                data['zdr_off_bb'].attrs["short_name"] = 'ZDR off BB'
                data['zdr_off_bb'].attrs["units"] = 'dB'
                data['zdr_off_bb'].attrs["comment"] = 'to subtract ' + \
                                                      'from ZDR'
                data['zdr_off_sd_bb'] = bb_sd
                data['zdr_off_sd_bb'].attrs["long_name"] = \
                    'ZDR offset standard deviation from bird bath'
                data['zdr_off_sd_bb'].attrs["short_name"] = 'ZDR off sd BB'
                data['zdr_off_sd_bb'].attrs["units"] = 'dB'
                data['zdr_off_bb_n'] = bb_nm
                data['zdr_off_bb_n'].attrs["long_name"] = 'number ' + \
                                                          'values for BB'
                data['zdr_off_bb_n'].attrs["short_name"] = 'nm for BB'
                data['zdr_off_bb_n'].attrs["units"] = '1'
                dtree = dttree.DataTree(name="root")
                dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                parent=dtree)
                print('saving: ... ' + path_out_nc.split('/')[-1] + ' ...')
                dtree.load().to_netcdf(path_out_nc)

            data.close()
        else:
            data = dttree.open_datatree(nc_file_mom)[
                'sweep_' + str(int(sweep))].to_dataset()
            data_rho = dttree.open_datatree(nc_file_rho)[
                'sweep_' + str(int(sweep))].to_dataset()
            data_temp = dttree.open_datatree(nc_file_temp)[
                'sweep_' + str(int(sweep))].to_dataset()
            data_temp2 = data_temp.interp(
                coords=data.drop(['longitude', 'latitude',
                                  'altitude', 'elevation']).coords,
                method='nearest')
            data.RHOHV.values = data_rho.RHOHV_NC2P.values
            data = data.assign({'temp_beamtop': data_temp2.temp_beamtop})
            data = data.transpose('time', 'azimuth', 'range')
            data2 = zdr_at_zero_elev(data)
            data2 = data.assign({'ZDR': data2.ZDR0})
            lr_off, lr_nm = cal_zdr_lightrain(data, band='C',
                                              plot=[False, False],
                                              axes=None,
                                              colorbar=[False, False])
            sd_off, sd_nm = cal_zdr_smalldrops(data, band='C',
                                               plot=[False, False],
                                               axes=None,
                                               colorbar=[False, False])
            lr_off_ec, lr_nm_ec = cal_zdr_lightrain(data2, band='C',
                                                    plot=[False, False],
                                                    axes=None,
                                                    colorbar=[False, False])
            sd_off_ec, sd_nm_ec = cal_zdr_smalldrops(data2, band='C',
                                                     plot=[False, False],
                                                     axes=None,
                                                     colorbar=[False, False])
            path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
            if not overwrite and os.path.exists(path_out_nc):
                print('exists: ' + path_out_nc + ' -> continue')
            else:
                remo_var = list(data.data_vars.keys())
                remo_var.remove('ZDR')
                data = data.drop_vars(remo_var)
                data['zdr_off_lr'] = lr_off
                data['zdr_off_lr'].attrs["long_name"] = 'ZDR offset ' + \
                                                        'from light rain'
                data['zdr_off_lr'].attrs["short_name"] = 'ZDR off LR'
                data['zdr_off_lr'].attrs["units"] = 'dB'
                data['zdr_off_lr'].attrs["comment"] = 'to subtract ' + \
                                                      'from ZDR'
                data['zdr_off_lr_n'] = lr_nm
                data['zdr_off_lr_n'].attrs["long_name"] = 'number ' + \
                                                          'values for LR'
                data['zdr_off_lr_n'].attrs["short_name"] = 'nm for LR'
                data['zdr_off_lr_n'].attrs["units"] = '1'
                data['zdr_off_lr_ec'] = lr_off_ec
                data['zdr_off_lr_ec'].attrs["long_name"] = \
                    'ZDR offset from light rain and elevation ' + \
                    'corrected ZDR'
                data['zdr_off_lr_ec'].attrs["short_name"] = 'ZDR off LR ec'
                data['zdr_off_lr_ec'].attrs["units"] = 'dB'
                data['zdr_off_lr_ec'].attrs["comment"] = 'to subtract ' + \
                                                         'from ZDR'
                data['zdr_off_lr_n_ec'] = lr_nm
                data['zdr_off_lr_n_ec'].attrs["long_name"] = \
                    'number values for LR ec'
                data['zdr_off_lr_n_ec'].attrs["short_name"] = 'nm for ' + \
                                                              'LR ec'
                data['zdr_off_lr_n_ec'].attrs["units"] = '1'
                data['zdr_off_sd'] = sd_off
                data['zdr_off_sd'].attrs["long_name"] = 'ZDR offset ' + \
                                                        'from small drops'
                data['zdr_off_sd'].attrs["short_name"] = 'ZDR off SD'
                data['zdr_off_sd'].attrs["units"] = 'dB'
                data['zdr_off_sd'].attrs["comment"] = 'to subtract ' + \
                                                      'from ZDR'
                data['zdr_off_sd_n'] = sd_nm
                data['zdr_off_sd_n'].attrs["long_name"] = 'number' + \
                                                          ' values for SD'
                data['zdr_off_sd_n'].attrs["short_name"] = 'nm for SD'
                data['zdr_off_sd_n'].attrs["units"] = '1'
                data['zdr_off_sd_ec'] = sd_off_ec
                data['zdr_off_sd_ec'].attrs["long_name"] = \
                    'ZDR offset from small drops and elevation ' + \
                    'corrected ZDR'
                data['zdr_off_sd_ec'].attrs["short_name"] = 'ZDR off SD ec'
                data['zdr_off_sd_ec'].attrs["units"] = 'dB'
                data['zdr_off_sd_ec'].attrs["comment"] = 'to subtract ' + \
                                                         'from ZDR'
                data['zdr_off_sd_n_ec'] = sd_nm
                data['zdr_off_sd_n_ec'].attrs["long_name"] = \
                    'number values for SD ec'
                data['zdr_off_sd_n_ec'].attrs["short_name"] = 'nm for ' + \
                                                              'SD ec'
                data['zdr_off_sd_n_ec'].attrs["units"] = '1'
                dtree = dttree.DataTree(name="root")
                dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                parent=dtree)
                print('saving: ... ' + path_out_nc.split('/')[-1] + ' ...')
                dtree.load().to_netcdf(path_out_nc)

            data.close()
            data2.close()
            data_rho.close()
            data_temp.close()
            data_temp2.close()
    else:
        print('exists: ' + path_out_nc + ' -> continue')

    return


# J. Steinheuer
def calibrate_zdr_with_plot(date, location,
                            elevation_degs=np.array([5.5, 5.5, 4.5, 3.5, 2.5,
                                                     1.5, 0.5, 8., 12., 17.,
                                                     25., 5.5]),
                            modes=np.array(['pcp', 'vol', 'vol', 'vol', 'vol',
                                            'vol', 'vol', 'vol', 'vol', 'vol',
                                            'vol', '90grad']),
                            overwrite=False, pdf_or_png='png',
                            ):
    """
    calibrate zdr and plot it.
    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    elevation_degs : elevations in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    modes : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool; if *allmoms*-output and plot exists, it can be
                      overwritten
    pdf_or_png : 'png' or 'pdf'
    """
    index = 0
    n_rows = 4
    if (modes == '90grad').all():
        n_rows = 1

    n_cols = elevation_degs.size
    plt.figure(figsize=(3 * n_cols, 3 * n_rows))
    save = True
    # ----------------------------------------------------------------------- #
    # loop over all elevations:
    # ----------------------------------------------------------------------- #
    for elevation_deg, mode in zip(elevation_degs, modes):
        print(mode + ' ' +
              (str(elevation_deg) if (mode == 'vol') else '') +
              ' ' + location + ' ' + date)
        folder_plot = header.folder_plot
        sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                                   float(elevation_deg))[0][0])
        year = date[0:4]
        mon = date[4:6]
        day = date[6:8]
        folder_in = "/".join([header.dir_data_obs + '*', year,
                              year + '-' + mon,
                              year + '-' + mon + '-' + day,
                              location, mode + '*', sweep])
        nc_file_mom = glob.glob(folder_in + '/*allmoms*')
        if len(nc_file_mom) > 1:
            print('mom: too many files')
            return
        elif len(nc_file_mom) == 0:
            print('mom: no files')
            return
        else:
            nc_file_mom = nc_file_mom[0]

        path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
        if overwrite or not os.path.exists(path_out_nc) or plot:
            if elevation_deg == elevation_degs[0] and mode == modes[0]:
                file_in = nc_file_mom.split('/')[-1]
                file_out = file_in.replace('.hd5', '_' + str(n_rows) +
                                           'x' + str(n_cols) + '.' +
                                           pdf_or_png).replace('_' + sweep,
                                                               '_all')
                path_out_plot = folder_plot + 'ZDR_calibration/BB/' + \
                                              'ZDR_calibration_' + \
                                str(n_rows) + 'x' + str(n_cols) + '_' + \
                                file_out
                if os.path.isfile(path_out_plot) and not overwrite:
                    print(path_out_plot + ' exists;\n' + ' ... set: ' +
                          '> overwrite = True < for recalculation')
                    plt.close()
                    save = False
                    break

            nc_file_temp = glob.glob(folder_in + '/*temp*')
            if len(nc_file_temp) > 1:
                print('temp: too many files')
                if mode != '90grad':
                    return

            elif len(nc_file_temp) == 0:
                print('temp: no files')
                if mode != '90grad':
                    return
            else:
                nc_file_temp = nc_file_temp[0]

            nc_file_rho = glob.glob(folder_in + '/*rhohv_nc*')
            if len(nc_file_rho) > 1:
                print('rho: too many files')
                if mode != '90grad':
                    return

            elif len(nc_file_rho) == 0:
                print('rho: no files')
                if mode != '90grad':
                    return

            else:
                nc_file_rho = nc_file_rho[0]

            # --------------------------------------------------------------- #
            # vol/pcp or bird bad                                             #
            # --------------------------------------------------------------- #
            if mode == '90grad':
                data = dttree.open_datatree(nc_file_mom)[
                    'sweep_' + str(int(sweep))].to_dataset()
                index = index + 1
                ax = plt.subplot(n_rows, n_cols, index + 0 * n_cols)
                bb_off, bb_nm, bb_sd = cal_zdr_birdbath(data, plot=plot, ax=ax)
                ax.set_title('$90°$ ' + location.upper() + ' ' + date)
                path_out_nc = nc_file_mom.replace(
                    '_allmoms_', '_zdr_off_')
                if not overwrite and os.path.exists(path_out_nc):
                    print('exists: ' + path_out_nc + ' -> continue')
                else:
                    remo_var = list(data.data_vars.keys())
                    remo_var.remove('ZDR')
                    data = data.drop_vars(remo_var)
                    data['zdr_off_bb'] = bb_off
                    data['zdr_off_bb'].attrs["long_name"] = 'ZDR offset ' + \
                                                            'from bird bath'
                    data['zdr_off_bb'].attrs["short_name"] = 'ZDR off BB'
                    data['zdr_off_bb'].attrs["units"] = 'dB'
                    data['zdr_off_bb'].attrs["comment"] = 'to subtract ' + \
                                                          'from ZDR'
                    data['zdr_off_sd_bb'] = bb_sd
                    data['zdr_off_sd_bb'].attrs["long_name"] = \
                        'ZDR offset standard deviation from bird bath'
                    data['zdr_off_sd_bb'].attrs["short_name"] = 'ZDR off sd BB'
                    data['zdr_off_sd_bb'].attrs["units"] = 'dB'
                    data['zdr_off_bb_n'] = bb_nm
                    data['zdr_off_bb_n'].attrs["long_name"] = 'number ' + \
                                                              'values for BB'
                    data['zdr_off_bb_n'].attrs["short_name"] = 'nm for BB'
                    data['zdr_off_bb_n'].attrs["units"] = '1'
                    dtree = dttree.DataTree(name="root")
                    dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                    parent=dtree)
                    print('saving: ... ' + path_out_nc.split('/')[-1] + ' ...')
                    dtree.load().to_netcdf(path_out_nc)

                data.close()
            else:
                data = dttree.open_datatree(nc_file_mom)[
                    'sweep_' + str(int(sweep))].to_dataset()
                data_rho = dttree.open_datatree(nc_file_rho)[
                    'sweep_' + str(int(sweep))].to_dataset()
                data_temp = dttree.open_datatree(nc_file_temp)[
                    'sweep_' + str(int(sweep))].to_dataset()
                data_temp2 = data_temp.interp(
                    coords=data.drop(['longitude', 'latitude',
                                      'altitude', 'elevation']).coords,
                    method='nearest')
                data.RHOHV.values = data_rho.RHOHV_NC2P.values
                data = data.assign({'temp_beamtop': data_temp2.temp_beamtop})
                data = data.transpose('time', 'azimuth', 'range')
                data2 = zdr_at_zero_elev(data)
                data2 = data.assign({'ZDR': data2.ZDR0})
                colorbar = [True, True]
                index = index + 1
                axes = [plt.subplot(n_rows, n_cols, index + 0 * n_cols),
                        plt.subplot(n_rows, n_cols, index + 0 * n_cols)]
                lr_off, lr_nm = cal_zdr_lightrain(data, band='C',
                                                  plot=[plot, False],
                                                  axes=axes,
                                                  colorbar=colorbar)
                if mode == 'pcp0':
                    axes[0].set_title('PCP ' + location.upper() + ' ' +
                                      date + '    ')
                else:
                    axes[0].set_title(str(elevation_deg) + '° ' +
                                      location.upper() + ' ' + date +
                                      '    ')

                axes = [plt.subplot(n_rows, n_cols, index + 1 * n_cols),
                        plt.subplot(n_rows, n_cols, index + 1 * n_cols)]
                sd_off, sd_nm = cal_zdr_smalldrops(data, band='C',
                                                   plot=[plot, False],
                                                   axes=axes,
                                                   colorbar=colorbar)
                axes = [plt.subplot(n_rows, n_cols, index + 2 * n_cols),
                        plt.subplot(n_rows, n_cols, index + 2 * n_cols)]
                lr_off_ec, lr_nm_ec = cal_zdr_lightrain(data2, band='C',
                                                        plot=[plot, False],
                                                        axes=axes,
                                                        colorbar=colorbar)
                if mode == 'pcp':
                    axes[0].set_title(r'$Z_{DR}(PCP) \rightarrow Z_{DR}(0°)$')
                else:
                    axes[0].set_title(r'$Z_{DR}($' + str(elevation_deg) +
                                      r'$°) \rightarrow Z_{DR}(0°)$')

                axes = [plt.subplot(n_rows, n_cols, index + 3 * n_cols),
                        plt.subplot(n_rows, n_cols, index + 3 * n_cols)]
                sd_off_ec, sd_nm_ec = cal_zdr_smalldrops(data2, band='C',
                                                         plot=[plot, False],
                                                         axes=axes,
                                                         colorbar=colorbar)
                path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
                if not overwrite and os.path.exists(path_out_nc):
                    print('exists: ' + path_out_nc + ' -> continue')
                else:
                    remo_var = list(data.data_vars.keys())
                    remo_var.remove('ZDR')
                    data = data.drop_vars(remo_var)
                    data['zdr_off_lr'] = lr_off
                    data['zdr_off_lr'].attrs["long_name"] = 'ZDR offset ' + \
                                                            'from light rain'
                    data['zdr_off_lr'].attrs["short_name"] = 'ZDR off LR'
                    data['zdr_off_lr'].attrs["units"] = 'dB'
                    data['zdr_off_lr'].attrs["comment"] = 'to subtract ' + \
                                                          'from ZDR'
                    data['zdr_off_lr_n'] = lr_nm
                    data['zdr_off_lr_n'].attrs["long_name"] = 'number ' + \
                                                              'values for LR'
                    data['zdr_off_lr_n'].attrs["short_name"] = 'nm for LR'
                    data['zdr_off_lr_n'].attrs["units"] = '1'
                    data['zdr_off_lr_ec'] = lr_off_ec
                    data['zdr_off_lr_ec'].attrs["long_name"] = \
                        'ZDR offset from light rain and elevation ' + \
                        'corrected ZDR'
                    data['zdr_off_lr_ec'].attrs["short_name"] = 'ZDR off LR ec'
                    data['zdr_off_lr_ec'].attrs["units"] = 'dB'
                    data['zdr_off_lr_ec'].attrs["comment"] = 'to subtract ' + \
                                                             'from ZDR'
                    data['zdr_off_lr_n_ec'] = lr_nm
                    data['zdr_off_lr_n_ec'].attrs["long_name"] = \
                        'number values for LR ec'
                    data['zdr_off_lr_n_ec'].attrs["short_name"] = 'nm for ' + \
                                                                  'LR ec'
                    data['zdr_off_lr_n_ec'].attrs["units"] = '1'
                    data['zdr_off_sd'] = sd_off
                    data['zdr_off_sd'].attrs["long_name"] = 'ZDR offset ' + \
                                                            'from small drops'
                    data['zdr_off_sd'].attrs["short_name"] = 'ZDR off SD'
                    data['zdr_off_sd'].attrs["units"] = 'dB'
                    data['zdr_off_sd'].attrs["comment"] = 'to subtract ' + \
                                                          'from ZDR'
                    data['zdr_off_sd_n'] = sd_nm
                    data['zdr_off_sd_n'].attrs["long_name"] = 'number' + \
                                                              ' values for SD'
                    data['zdr_off_sd_n'].attrs["short_name"] = 'nm for SD'
                    data['zdr_off_sd_n'].attrs["units"] = '1'
                    data['zdr_off_sd_ec'] = sd_off_ec
                    data['zdr_off_sd_ec'].attrs["long_name"] = \
                        'ZDR offset from small drops and elevation ' + \
                        'corrected ZDR'
                    data['zdr_off_sd_ec'].attrs["short_name"] = 'ZDR off SD ec'
                    data['zdr_off_sd_ec'].attrs["units"] = 'dB'
                    data['zdr_off_sd_ec'].attrs["comment"] = 'to subtract ' + \
                                                             'from ZDR'
                    data['zdr_off_sd_n_ec'] = sd_nm
                    data['zdr_off_sd_n_ec'].attrs["long_name"] = \
                        'number values for SD ec'
                    data['zdr_off_sd_n_ec'].attrs["short_name"] = 'nm for ' + \
                                                                  'SD ec'
                    data['zdr_off_sd_n_ec'].attrs["units"] = '1'
                    dtree = dttree.DataTree(name="root")
                    dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                    parent=dtree)
                    print('saving: ... ' + path_out_nc.split('/')[-1] + ' ...')
                    dtree.load().to_netcdf(path_out_nc)

                data.close()
                data2.close()
                data_rho.close()
                data_temp.close()
                data_temp2.close()
        else:
            print('exists: ' + path_out_nc + ' -> continue')
            save = False

    # ----------------------------------------------------------------------- #
    # SAVE                                                                    #
    # ----------------------------------------------------------------------- #
    if save:
        plt.tight_layout()
        if not os.path.exists(folder_plot + 'ZDR_calibration/' +
                              location.upper() + '/'):
            os.makedirs(folder_plot + 'ZDR_calibration/' +
                        location.upper() + '/')

        plt.savefig(path_out_plot, format=pdf_or_png, transparent=True)
        plt.close()

    return


# --------------------------------------------------------------------------- #
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20210604",  # case01
    "20210620", "20210621",  # case02
    "20210628", "20210629",  # case03
    "20220519", "20220520",  # case04
    "20220623", "20220624", "20220625",  # case05
    "20220626", "20220627", "20220628",  # case06+07
    "20220630", "20220701",  # case08
    "20210714",  # case09
    "20221222",  # case10
    "20170725",  # caseX  # OLD CASE
]
LOCATIONS = [
    'asb', 'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
overwrite = False
# ----------------------------------- #
plot = True
pdf_or_png = 'png'
include_sweep = np.array([
    False,
    False, False, False, False, False, False,
    False, False, False, False,
    True,
])
elevation_degs = np.array([
    5.5,
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8., 12., 17., 25.,
    5.5,
])
modes = np.array([
    'pcp',
    'vol', 'vol', 'vol', 'vol', 'vol', 'vol',
    'vol', 'vol', 'vol', 'vol',
    '90grad'
])
elevation_degs = elevation_degs[include_sweep]
modes = modes[include_sweep]
# sorting volume
elevation_degs_2_sort = elevation_degs.copy()
elevation_degs_2_sort[modes == 'pcp'] = 0
elevation_degs_2_sort[modes == '90grad'] = 90
sort_i = elevation_degs_2_sort.argsort()
elevation_degs = elevation_degs[sort_i]
modes = modes[sort_i]
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        if plot:
            calibrate_zdr_with_plot(date, location,
                                    elevation_degs=elevation_degs,
                                    modes=modes, overwrite=overwrite,
                                    pdf_or_png=pdf_or_png
                                    )
        else:
            for elevation_deg, mode in zip(elevation_degs, modes):
                calibrate_zdr(date, location, elevation_deg=elevation_deg,
                              mode=mode, overwrite=overwrite)

