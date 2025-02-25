#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# PLOT_SYN_RADAR.py                                                           #
#                                                                             #
# Functions to plot QVPs, PPIs, pseudoRHIs,... from given synthetic           #
# (EMVORADO) QVP netcdf files.                                                #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
from PLOT_RADAR import d0_bringi, ice_retrieval_carlin
import datatree as dttree
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
from pathlib import Path
from statsmodels.stats.weightstats import DescrStatsW
import warnings

warnings.simplefilter('ignore')
import wradlib as wrl
import xarray as xr

xr.set_options(keep_attrs=True)


# --------------------------------------------------------------------------- #

# def d0_bringi(
#         zdr
# ):
#     d0 = np.full(zdr.shape, np.nan)
#     mask_thresh = zdr < 1.25
#     d0[mask_thresh] = \
#         0.0203 * zdr[mask_thresh] ** 4 - \
#         0.149 * zdr[mask_thresh] ** 3 + \
#         0.221 * zdr[mask_thresh] ** 2 + \
#         0.557 * zdr[mask_thresh] + 0.801
#     d0[~mask_thresh] = \
#         0.0355 * zdr[~mask_thresh] ** 3 - \
#         0.302 * zdr[~mask_thresh] ** 2 + \
#         1.06 * zdr[~mask_thresh] + 0.684
#     return d0


def plot_qvp_of_polarimetric_variable(
        mom,
        cmap,
        norm,
        levels,
        mom_cs=None,
        levels_cs=None,
        mom_cf=None,
        levels_cf=None,
        ax=None,
        cbar_title=None,
        title=None,
        scale_font=1.2,
        scale_numbers=1,
        y_coord='height',
        extend='both',
        xlabel='UTC [mm-dd hh]',
        ylabel='height [km]',
        top_height=10000,  # in m  # TODO check for TS obs
        mom_height_unit='m',
        save=False,
        save_name='test',
        save_path=''
):
    if ax is None:
        ax = plt.gca()

    fg = mom.plot.contourf(y=y_coord, cmap=cmap, norm=norm,
                           extend=extend, cbar_kwargs={'ticks': levels})
    if mom_cs is not None:
        cs = mom_cs.plot.contour(y=y_coord, colors='black',
                                 levels=levels_cs,
                                 extend=extend, linewidths=.5)
        ax.clabel(cs, **dict(fmt=r'%2.0f', inline=True,
                             fontsize=8 * scale_numbers))
    if mom_cf is not None:
        cf = mom_cf.plot.contour(y=y_coord, colors='dimgrey',
                                 levels=levels_cf,
                                 linestyles='dotted')
        ax.clabel(cf, **dict(fmt='%2.0f °C', inline=True,
                             fontsize=7 * scale_numbers))
    plt.xlabel(xlabel, fontsize=10 * scale_font)
    plt.xticks(fontsize=8 * scale_numbers)
    plt.ylabel(ylabel, fontsize=10 * scale_font)
    plt.ylim(0, top_height)
    if mom_height_unit == 'm':
        plt.yticks(ticks=plt.yticks()[0], labels=plt.yticks()[0] / 1000,
                   fontsize=10 * scale_numbers)
    elif mom_height_unit == 'km':
        plt.yticks(ticks=plt.yticks()[0], labels=plt.yticks()[0],
                   fontsize=10 * scale_numbers)

    if cbar_title is not None:
        fg.colorbar.set_label(cbar_title, fontsize=10 * scale_font)

    fg.colorbar.ax.tick_params(labelsize=10 * scale_numbers)
    if title is not None:
        plt.title(title, fontsize=7 * scale_font)

    plt.tight_layout()
    if save:
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)


def plot_syn_pseudoRHI(nc_file,
                       ax=None,
                       azimuth=0,
                       moment='KDP_NC',
                       cmap=None,  # 'jet',
                       levels=None,
                       norm=None,
                       title=None,
                       range_max=None,
                       # func='contourf',
                       func='pcolormesh',
                       ):
    # sweep = nc_file.split('/')[-2]
    vol = nc_file
    vol = vol.assign_coords(sweep_fixed_angle=vol.elevation)
    vol = vol.swap_dims({'elevation': 'sweep_fixed_angle'})
    vol = vol.assign_coords(prt_mode=(
        'sweep_fixed_angle',
        np.repeat('not set', vol.elevation.size)))
    vol = vol.assign_coords(sweep_mode=(
        'sweep_fixed_angle',
        np.repeat('azimuth_surveillance', vol.elevation.size)))
    vol["sweep_mode"] = vol["sweep_mode"].min()
    vol = vol.assign_coords(follow_mode=(
        'sweep_fixed_angle',
        np.repeat('not set', vol.elevation.size)))
    vol = vol.assign_coords(sweep_number=(
        'sweep_fixed_angle',
        np.arange(0, vol.elevation.size)))
    vol = vol.assign_coords(longitude=vol.station_longitude)
    vol = vol.assign_coords(latitude=vol.station_latitude)
    vol = vol.assign_coords(altitude=vol.station_height)
    vol.transpose('sweep_fixed_angle', 'azimuth', 'range')
    if type(moment) == list:
        for moment_i in moment:
            if moment_i in vol.variables:
                moment = moment_i
                break

    log = False
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
    elif 'rho' in moment.lower() or 'rhv' in moment.lower():
        if levels is None:
            levels = header.levels_rhohv
        if norm is None:
            norm = header.norm_rhohv
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['DBZH', 'zh']:
    elif 'zh' in moment.lower() or 'zr' in moment.lower():
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
    elif 'w' == moment.lower():
        if levels is None:
            levels = np.arange(-16, 16.1, 2)
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar2  # _smooth
    elif 'dm_' in moment.lower():
        if levels is None:
            levels = [i for i in np.arange(0, 8.5, .5)]
        if moment == 'Dm_cloud':
            levels = [i for i in np.arange(0, 0.17, 0.01)]
        if moment == 'Dm_hail':
            levels = [i for i in np.arange(0, 33, 2)]
        if moment == 'Dm_ice':
            levels = [i for i in np.arange(0, 1.7, .1)]
        if moment == 'Dm_graupel':
            levels = [i for i in np.arange(0, 17, 1)]
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar2  # _smooth
    elif 'qnr' in moment.lower() or 'qnc' in moment.lower() or \
            'qng' in moment.lower() or 'qnh' in moment.lower() or \
            'qni' in moment.lower() or 'qns' in moment.lower():
        if levels is None:
            # levels = [10.0**i for i in np.arange(-4, 8, 1)]
            levels = [i for i in np.arange(-2, 15, 1)]
            vol[moment] = np.log(vol[moment])
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar2  # _smooth
        log = True
    elif 'qr' in moment.lower() or 'qc' in moment.lower() or \
            'qg' in moment.lower() or 'qh' in moment.lower() or \
            'qi' in moment.lower() or 'qs' in moment.lower():
        if levels is None:
            # levels = [10.0**i for i in np.arange(-11, 3, 1)]
            levels = [i for i in np.arange(-20, -3, 1)]
            vol[moment] = np.log(vol[moment])
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar2  # _smooth
        log = True
    else:
        if cmap is None:
            cmap = 'jet'

    if ax is None:
        plt.figure(figsize=(5, 4))
        ax = plt.gca()

    rec_rhi = wrl.util.cross_section_ppi(vol, azimuth, method="nearest",
                                         bw=1.0)
    rec_rhi[moment].plot(x="gr", y="z", ax=ax, cmap=cmap, levels=levels,
                         norm=norm,
                         extend='both',
                         ylim=[0, 15000],
                         xlim=[0, range_max * 1000],
                         )
    # img = wrl.georef.create_xarray_dataarray(data=ppi[moment],
    #                                          r=ppi.range.values / 1000,
    #                                          phi=phi,
    #                                          theta=theta,
    #                                          # theta=ppi.elevation.values.flatten(),#[:],#.item(),
    #                                          site=[ppi.station_longitude.values,
    #                                                ppi.station_latitude.values,
    #                                                ppi.station_height.values]
    #                                          )
    # img = img.wrl.georef.georeference()

    # if range_max:
    #     img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
    #                      extend='both', func=func,
    #                      ylim=[-range_max, range_max],
    #                      xlim=[-range_max, range_max])
    # else:
    #     img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
    #                      extend='both', func=func, )

    if title:
        plt.title(title)

    if log:
        plt.text(1.1, -.05, '1e^',
                 transform=ax.transAxes)
    #
    # ax.set_xlabel("easting [km]")
    # ax.set_ylabel("northing [km]")


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def plot_syn_PPI(nc_file,
                 ax=None,
                 time_i=0,
                 moment='KDP_NC',
                 cmap=None,  # 'jet',
                 levels=None,
                 norm=None,
                 title=None,
                 range_max=None,
                 # func='contourf',
                 func='pcolormesh',
                 extend='both'
                 ):
    # sweep = nc_file.split('/')[-2]
    vol = nc_file
    ppi = vol.isel(time=time_i)
    ppi = ppi.transpose('azimuth', 'range', )

    if type(moment) == list:
        for moment_i in moment:
            if moment_i in ppi.variables:
                moment = moment_i
                break

    format = None
    log = False
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
    elif 'rho' in moment.lower() or 'rhv' in moment.lower():
        if levels is None:
            levels = header.levels_rhohv
        if norm is None:
            norm = header.norm_rhohv
        if cmap is None:
            cmap = header.cmap_radar
    # elif moment in ['DBZH', 'zh']:
    elif 'zh' in moment.lower() or 'zr' in moment.lower():
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
    elif 'w' == moment.lower():
        if levels is None:
            levels = np.arange(-16, 16.1, 2)
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar  # _smooth
    elif 'dm_' in moment.lower():
        if levels is None:
            levels = [i for i in np.arange(0, 8.5, .5)]
        if moment == 'Dm_cloud':
            levels = [i for i in np.arange(0, 0.17, 0.01)]
        if moment == 'Dm_hail':
            levels = [i for i in np.arange(0, 33, 2)]
        if moment == 'Dm_ice':
            levels = [i for i in np.arange(0, 1.7, .1)]
        if moment == 'Dm_graupel':
            levels = [i for i in np.arange(0, 17, 1)]
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar  # _smooth
    elif 'qnr' in moment.lower() or 'qnc' in moment.lower() or \
            'qng' in moment.lower() or 'qnh' in moment.lower() or \
            'qni' in moment.lower() or 'qns' in moment.lower():
        if levels is None:
            # levels = [10.0**i for i in np.arange(-4, 8, 1)]
            levels = [i for i in np.arange(-2, 15, 1)]
            ppi[moment] = np.log(ppi[moment])
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar  # _smooth
        log = True
    elif 'qr' in moment.lower() or 'qc' in moment.lower() or \
            'qg' in moment.lower() or 'qh' in moment.lower() or \
            'qi' in moment.lower() or 'qs' in moment.lower():
        if levels is None:
            # levels = [10.0**i for i in np.arange(-11, 3, 1)]
            levels = [i for i in np.arange(-20, -3, 1)]
            ppi[moment] = np.log(ppi[moment])
        if norm is None:
            norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
        if cmap is None:
            cmap = header.cmap_radar  # _smooth
        log = True
    # elif 'w' == moment.lower():
    #     if levels is None:
    #         levels = np.arange(-16, 16.1, 2)
    #     if norm is None:
    #         norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
    #     if cmap is None:
    #         cmap = header.cmap_radar_smooth
    # elif 'qnr' in moment.lower() or 'qnc' in moment.lower() or \
    #         'qng' in moment.lower() or 'qnh' in moment.lower() or \
    #         'qni' in moment.lower() or 'qns' in moment.lower():
    #     if levels is None:
    #         # levels = [10.0**i for i in np.arange(-4, 8, 1)]
    #         levels = [i for i in np.arange(-4, 16, 1)]
    #         ppi[moment] = np.log(ppi[moment])
    #     if norm is None:
    #         norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
    #     if cmap is None:
    #         cmap = header.cmap_radar_smooth
    #     log=True
    # elif 'qr' in moment.lower() or 'qc' in moment.lower() or \
    #         'qg' in moment.lower() or 'qh' in moment.lower() or \
    #         'qi' in moment.lower() or 'qs' in moment.lower():
    #     if levels is None:
    #         # levels = [10.0**i for i in np.arange(-11, 3, 1)]
    #         levels = [i for i in np.arange(-23, -3, 1)]
    #         ppi[moment] = np.log(ppi[moment])
    #     if norm is None:
    #         norm = mpl.colors.BoundaryNorm(levels, len(levels) - 1)
    #     if cmap is None:
    #         cmap = header.cmap_radar_smooth
    #     log=True
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
                                             site=[
                                                 ppi.station_longitude.values,
                                                 ppi.station_latitude.values,
                                                 ppi.station_height.values],
                                             format=format
                                             )
    img = img.wrl.georef.georeference()
    if ax is None:
        plt.figure(figsize=(5, 4))
        ax = plt.gca()

    if range_max:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         extend=extend, func=func,
                         ylim=[-range_max, range_max],
                         xlim=[-range_max, range_max])
    else:
        img.wrl.vis.plot(ax=ax, cmap=cmap, levels=levels, norm=norm,
                         extend=extend, func=func)

    if title:
        plt.title(title)

    ax.set_xlabel("easting [km]")
    ax.set_ylabel("northing [km]")
    if log:
        plt.text(1.1, -.05, '1e^',
                 transform=ax.transAxes)


# not working yet
def plot_syn_PPI_temp_ring(nc_file,
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
                                             site=[
                                                 ppi.station_longitude.values,
                                                 ppi.station_latitude.values,
                                                 ppi.station_height.values]
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


def mom_plot_dict(mom_name='zh'):
    if mom_name in ['zh', 'ZH', 'ZH_AC', 'zrsim']:
        return dict(mom_name=mom_name, mom_min=0, mom_max=40, bins_mom=30, )
    elif mom_name in ['zdr', 'zdrsim', 'ZDR', 'ZDR_AC_OC']:
        return dict(mom_name=mom_name, mom_min=-0.5, mom_max=3, bins_mom=40, )
    elif mom_name in ['kdpsim', 'KDP',  'KDP_NC',
                    'kdpsim_ml_corrected', 'KDP_ML_corrected']:
        return dict(mom_name=mom_name, mom_min=0, mom_max=0.3, bins_mom=40, )
    elif mom_name in ['rho', 'rhvsim', 'RHOHV', 'RHOHV_NC2P', ]:
        # return dict(mom_name=mom_name, mom_min=0.951, mom_max=1.051,
        return dict(mom_name=mom_name, mom_min=0.925, mom_max=1.025,
                    bins_mom=25, )
                    # bins_mom=40, )
    else:
        print('unknown moment: ' + mom_name)
        return dict(mom_name=mom_name, mom_min=-20, mom_max=20,
                    bins_mom=40, )


def plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=['PRO'],
        dates=['20170725'],
        hhmm_start='00:00',
        hhmm_end='23:55',
        elevation_deg=12,
        paths_in=None,
        title=None,
        da_icon_emvorado_run='ASS_2211/MAIN_2401.3/EMVO_00500000.2',
        spin_up_mm='60',
        moment='zrsim',
        mom_min=None,
        mom_max=None,
        bins_mom=None,
        vert_temp=True,  # CFTD
        temp_min=-20,
        temp_max=16,
        bins_temp=18,
        height_min=0,  # in km
        height_max=10,  # in km
        bins_height=20,
        vmax=15,
        filter_entr=False,
        filter_moms=True,
        ax=None,
        save=False,
        save_path=header.folder_plot + 'CFADs/',
        save_name='test_CFAD'
):
    # locations = ['PRO']
    # dates = ['20170725']
    # hhmm_start = '00:00'
    # hhmm_end = '23:55'
    # elevation_deg = 12
    # paths_in = ['/automount/data02/agradar/operation_hydrometeors/data/QVP/20170725/ASS_2211/MAIN_2401.3/EMVO_00500000.2/60min_spinup/QVP_12_Syn_PRO_201707250000_201707252355.nc']
    # paths_in = ['/automount/data02/agradar/operation_hydrometeors/data/obs_qvp/OpHymet2-caseX-20170725/2017/2017-07/2017-07-25/pro/vol5minng01/07/ras07-qvp5minng01_sweeph5onem_polmoms_nc_07-201707250003-201707252358-pro-10392.hd5']
    # title = None
    # da_icon_emvorado_run = 'ASS_2211/MAIN_2401.3/EMVO_00500000.2'
    # da_icon_emvorado_run = None
    # spin_up_mm = '60'
    # moment = 'zrsim'
    # moment = 'ZH_AC'
    # mom_min = 0
    # mom_max = 40
    # bins_mom = 40
    # vert_temp = True # CFTD
    # temp_min = -20
    # temp_max = 16
    # bins_temp = 18
    # height_min = 0  # in km
    # height_max = 10  # in km
    # bins_height = 20
    # vmax = None
    # filter_entr = False
    # ax = None
    # save = False
    # save_path = header.folder_plot + 'CFADs/'
    # save_name = 'test_CFAD'

    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if not isinstance(dates, list):
        dates = [dates]

    if not isinstance(locations, list):
        locations = [locations]

    if paths_in is None:
        if da_icon_emvorado_run is None:
            paths_in = glob.glob('/'.join([header.dir_obs_qvp + '*', '*',
                                           '*', '*', '*', '*', sweep, '*', ]))
        else:
            paths_in = glob.glob('/'.join([header.dir_data_qvp + '*',
                                           da_icon_emvorado_run + '/' +
                                           str(spin_up_mm) +
                                           'min_spinup/QVP_' +
                                           str(elevation_deg) + '_Syn_*', ]))
            if paths_in == []:
                print('nothing found in ' +
                      '/'.join([header.dir_data_qvp + '*',
                                da_icon_emvorado_run + '/' +
                                str(spin_up_mm) +
                                'min_spinup/QVP_' +
                                str(elevation_deg) + '_Syn_*', ]))
                return

    if not isinstance(paths_in, list):
        paths_in = [paths_in]

    if mom_max is None:
        mom_max = mom_plot_dict(moment)['mom_max']

    if mom_min is None:
        mom_min = mom_plot_dict(moment)['mom_min']

    if bins_mom is None:
        bins_mom = mom_plot_dict(moment)['bins_mom']

    mom_all = np.array([], )
    y_all = np.array([], )
    weights_all = np.array([], )
    n_cases = 0
    if da_icon_emvorado_run:  # Synthetic
        if not title:
            title = '-'.join([da_icon_emvorado_run[4:8],
                              da_icon_emvorado_run[14:20],
                              da_icon_emvorado_run[26:36],
                              spin_up_mm + 'min'])
    else:  # Observation
        if not title:
            title = 'C-band Observations'

    for path_in in paths_in:
        # check for proper files
        file_in = path_in.split('/')[-1]
        if da_icon_emvorado_run:  # Synthetic
            date = file_in.split('_')[4][:8]
            if (da_icon_emvorado_run + '/' + str(spin_up_mm) +
                'min_spinup/QVP_' + str(elevation_deg) + '_Syn_' in path_in) \
                    and (file_in.split('_')[1] == str(elevation_deg)) \
                    and (file_in.split('_')[3] in locations) \
                    and (date in dates)\
                    and (date == path_in.split('/')[-6]):
                year = date[0:4]
                mon = date[4:6]
                day = date[6:8]
                date_start = '-'.join([year, mon, day, hhmm_start])
                date_end = '-'.join([year, mon, day, hhmm_end])
            else:
                print('Not used: ' + path_in)
                continue

        else:  # Observation
            date = file_in.split('-')[2][:8]
            if ('/obs_qvp/' in path_in) \
                    and (path_in.split('/')[-2] == sweep) \
                    and (path_in.split('/')[-4].upper() in locations) \
                    and (date in dates) \
                    and ('_polmoms_nc' in file_in):
                year = date[0:4]
                mon = date[4:6]
                day = date[6:8]
                date_start = '-'.join([year, mon, day, hhmm_start])
                date_end = '-'.join([year, mon, day, hhmm_end])
            else:
                # print('Not used: ' + path_in)
                continue

        print('Included: ' + path_in)
        n_cases = n_cases + 1
        # OPEN
        file_nc = xr.open_dataset(path_in).transpose('time', ...)
        file_nc = file_nc.sel(time=slice(date_start, date_end))
        if filter_entr:
            file_nc = xr.where(file_nc['min_entropy'] > 0.8,
                               file_nc, np.nan)

        if filter_moms:
            if da_icon_emvorado_run:  # Synthetic
                zh = file_nc['zrsim']
                zdr = file_nc['zdrsim']
                kdp = file_nc['kdpsim']
                rho = file_nc['rhvsim']
            else:  # Observations
                zh = file_nc['ZH_AC']
                zdr = file_nc['ZDR_AC_OC']
                kdp = file_nc['KDP_NC']
                rho = file_nc['RHOHV_NC2P']

            file_nc = xr.where(zh > 0, file_nc, np.nan)
            file_nc = xr.where(zh < 80, file_nc, np.nan)
            file_nc = xr.where(zdr > -1, file_nc, np.nan)
            file_nc = xr.where(zdr < 8, file_nc, np.nan)
            file_nc = xr.where(rho > .7, file_nc, np.nan)
            file_nc = xr.where(rho < 1.1, file_nc, np.nan)
            file_nc = xr.where(kdp > -1, file_nc, np.nan)
            file_nc = xr.where(kdp < 8, file_nc, np.nan)

        if moment in list(file_nc.keys()):
            mom = file_nc[moment].values
        else:
            print(moment + ' not in ' + file_in)
            continue

        if vert_temp:
            y = file_nc.temp.values
            y_min, y_max, bins_y = temp_min, temp_max, bins_temp
            if 'units' in file_nc.temp.attrs:  # else: TS in °C
                if file_nc.temp.units == 'K':
                    y = y - 273.15
        else:
            y = np.array([file_nc.height.values] * file_nc.time.size) / 1000
            y_min, y_max, bins_y = height_min, height_max, bins_height
            if 'units' in file_nc.height.attrs:  # else: TS in m
                if file_nc.height.units == 'km':
                    y = y * 1000

        mom_min_outer = 2 * mom_min - mom_max  # TODO: necessary?
        mom_max_outer = 2 * mom_max - mom_min  # TODO: necessary?
        # mom_min_outer = mom_min   # TODO: necessary?
        # mom_max_outer = mom_max   # TODO: necessary?
        mask = (y >= y_min) & (y <= y_max) & \
               (mom >= -999) & (mom <= 999)
        mom = mom[mask]
        y = y[mask]
        i_sort = np.argsort(y)
        mom = mom[i_sort]
        y = y[i_sort]
        if ax is None:
            plt.figure(figsize=(6, 5))

        a = plt.hist(y, bins=bins_y, range=(y_min, y_max))
        weights = np.repeat(100 / a[0], np.int16(a[0]))
        mom_all = np.append(mom_all, mom)
        y_all = np.append(y_all, y)
        weights_all = np.append(weights_all, weights)
        # file_nc.close()

    # out of loop
    weights_all = weights_all / (n_cases)
    # PLOT
    if vmax:
        extend = 'max'
    else:
        extend = 'neither'

    h2d, mom2d, y2d, fg = plt.hist2d(mom_all, y_all,
                                     bins=[bins_mom * 3, bins_y],
                                     range=[[mom_min_outer, mom_max_outer],
                                            [y_min, y_max]], vmax=vmax,
                                     weights=weights_all, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]', extend=extend)
    if vert_temp:
        plt.gca().invert_yaxis()
        plt.ylabel('temperature (°C)')
    else:
        plt.ylabel(r'height (km)')

    plt.xlabel(moment + ' [' + file_nc[moment].units + ']')
    plt.title(title)
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
    plt.xlim([mom_min, mom_max])
    if vmax and np.max(h2d) > vmax:
        if np.max(h2d) > 100:
            print('above 100% :' + str(np.max(h2d)))

        plt.text(1.04, 1.03, min(np.round(np.max(h2d), 1), 100.0),
                 transform=ax.transAxes)
    file_nc.close()
    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)


def plot_CFAD_or_CFTD_from_QVP(
        dates=['20170725'],
        hhmm_start='00:00',
        hhmm_end='23:55',
        paths_in=None,
        title=None,
        folder_syn=header.dir_data_qvp,
        locations='PRO',
        elevation_deg=12,
        da_run='ASS_2211',
        icon_emvorado_run='MAIN_2401.3/EMVO_00500000.2',
        spin_up_mm='60',
        moment='zrsim',
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
        filter_entr_ML=False,
        filter_ML=None,
        filter2=False,
        ax=None,
        save=False,
        save_path=header.folder_plot + 'CFADs/',
        save_name='test_CFAD'
):
    if not isinstance(dates, list):
        dates = [dates]

    if not isinstance(locations, list):
        locations = [locations]

    if not isinstance(da_run, list):
        da_run = [da_run]

    if len(da_run) == 1:
        da_run = [da_run[0] for i in range(len(dates))]

    if not isinstance(icon_emvorado_run, list):
        icon_emvorado_run = [icon_emvorado_run]

    if len(icon_emvorado_run) == 1:
        icon_emvorado_run = [icon_emvorado_run[0] for i in range(len(dates))]

    if paths_in is not None and not isinstance(paths_in, list):
        paths_in = [paths_in]

    mom_all = np.array([], )
    y_all = np.array([], )
    weights_all = np.array([], )
    n_cases = 0
    for j, location in enumerate(locations):
        for i, date in enumerate(dates):
            year = date[0:4]
            mon = date[4:6]
            day = date[6:8]
            date_start = '-'.join([year, mon, day, hhmm_start])
            date_end = '-'.join([year, mon, day, hhmm_end])
            if paths_in:
                path_in = paths_in[j * len(dates) + i]
                if not title:
                    title = path_in.split('/')[-1]
            else:
                path_in = '/'.join([folder_syn + date, da_run[i],
                                    icon_emvorado_run[i],
                                    str(spin_up_mm) + 'min_spinup', 'QVP_' +
                                    str(elevation_deg) + '_Syn_' + location +
                                    '_' + date + '0000_' + date + '2355.nc'])
                if not title:
                    title = '-'.join([da_run[i][4:],
                                      icon_emvorado_run[i].split('/')[0][5:],
                                      icon_emvorado_run[i].split('/')[1][5:],
                                      spin_up_mm + 'min'])

            if not os.path.exists(path_in):
                print('Missing: ' + path_in)
                continue

            print('including: ' + path_in)
            n_cases = n_cases + 1
            # OPEN
            syn_nc = xr.open_dataset(path_in)
            syn_nc = syn_nc.sel(time=slice(date_start, date_end))
            if filter_entr_ML:
                min_entropy = syn_nc['min_entropy'].transpose('time', ...)
                syn_nc = xr.where(min_entropy > 0.8, syn_nc, np.nan)
                if filter_ML == None:
                    filter_ML = True

            if filter_ML:
                if 'mlh_top' in list(syn_nc.keys()):
                    ml_top = syn_nc['mlh_top'].transpose('time', ...)
                # elif 'height_ml' in list(syn_nc.keys()):
                elif 'height_ml' in list(syn_nc.coords):  # TS: coordinate
                    ml_top = syn_nc['height_ml'].transpose('time', ...)
                else:
                    print('no ML height found')
                    continue

                syn_nc = xr.where(ml_top > 0, syn_nc, np.nan)

            if 'zh' in list(syn_nc.keys()):
                zh = syn_nc['zh'].transpose('time', ...)
            elif 'zrsim' in list(syn_nc.keys()):
                zh = syn_nc['zrsim'].transpose('time', ...)
            elif 'ZH_AC' in list(syn_nc.keys()):
                zh = syn_nc['ZH_AC'].transpose('time', ...)
            else:
                print('no zh found')
                continue

            if 'zdr' in list(syn_nc.keys()):
                zdr = syn_nc['zdr'].transpose('time', ...)
            elif 'zdrsim' in list(syn_nc.keys()):
                zdr = syn_nc['zdrsim'].transpose('time', ...)
            elif 'ZDR_AC_OC' in list(syn_nc.keys()):
                zdr = syn_nc['ZDR_AC_OC'].transpose('time', ...)
            else:
                print('no zdr found')
                continue

            if 'kdp' in list(syn_nc.keys()):
                kdp = syn_nc['kdp'].transpose('time', ...)
            elif 'kdpsim' in list(syn_nc.keys()):
                kdp = syn_nc['kdpsim'].transpose('time', ...)
            elif 'KDP_NC' in list(syn_nc.keys()):
                kdp = syn_nc['KDP_NC'].transpose('time', ...)
            else:
                print('no kdp found')
                continue

            if 'rho' in list(syn_nc.keys()):
                rho = syn_nc['rho'].transpose('time', ...)
            elif 'rhvsim' in list(syn_nc.keys()):
                rho = syn_nc['rhvsim'].transpose('time', ...)
            elif 'RHOHV_NC2P' in list(syn_nc.keys()):
                rho = syn_nc['RHOHV_NC2P'].transpose('time', ...)
            else:
                print('no rho found')
                continue

            syn_nc = xr.where(zh > 0, syn_nc, np.nan)
            syn_nc = xr.where(zh < 80, syn_nc, np.nan)
            syn_nc = xr.where(zdr > -1, syn_nc, np.nan)
            syn_nc = xr.where(zdr < 8, syn_nc, np.nan)
            syn_nc = xr.where(rho > .7, syn_nc, np.nan)
            syn_nc = xr.where(rho < 1.1, syn_nc, np.nan)
            syn_nc = xr.where(kdp > -1, syn_nc, np.nan)
            syn_nc = xr.where(kdp < 8, syn_nc, np.nan)

            if moment == 'D0_r_Bringi_obs':
                mom = d0_bringi(syn_nc['zdr'].transpose('time', ...))[
                    'dm'].values
            elif moment == 'iwc_obs':
                mom = ice_retrieval_carlin(
                    syn_nc['zh'].transpose('time', ...),
                    syn_nc['zdr'].transpose('time', ...),
                    syn_nc['kdp'].transpose('time', ...))['iwc'].values
            elif moment == 'nt_obs':
                mom = ice_retrieval_carlin(
                    syn_nc['zh'].transpose('time', ...),
                    syn_nc['zdr'].transpose('time', ...),
                    syn_nc['kdp'].transpose('time', ...))['nt'].values
            elif moment == 'dm_obs':
                mom = ice_retrieval_carlin(
                    syn_nc['zh'].transpose('time', ...),
                    syn_nc['zdr'].transpose('time', ...),
                    syn_nc['kdp'].transpose('time', ...))['dm'].values
            elif moment == 'qtotice':
                mom = syn_nc['qg'].transpose('time', ...).values + \
                      syn_nc['qs'].transpose('time', ...).values + \
                      syn_nc['qi'].transpose('time', ...).values + \
                      syn_nc['qh'].transpose('time', ...).values
                mom = mom * 1000
            elif moment == 'D0_r_Bringi_syn':
                mom = d0_bringi(syn_nc['zdrsim'].transpose('time', ...))[
                    'dm'].values
                if filter2:
                    qvp_d0q_syn = syn_nc['D0_g'].transpose('time', ...).values
                    qvp_d0r_syn = syn_nc['D0_r'].transpose('time', ...).values
                    mom = np.where(qvp_d0q_syn < 1, mom, np.nan)
                    mom = np.where(qvp_d0r_syn > 0, mom, np.nan)
            else:
                mom = syn_nc[moment].transpose('time', ...).values

            if moment == 'vol_qntotice' or moment == 'nt_obs':
                mom = np.log10(mom) - 3

            if vert_temp:
                y = syn_nc.temp.transpose('time', ...).values
                y_min, y_max, bins_y = temp_min, temp_max, bins_temp
                if 'units' in syn_nc.temp.attrs:  # else: TS in °C
                    if syn_nc.temp.units == 'K':
                        y = y - 273.15
            else:
                y = np.array([syn_nc.height.values] * syn_nc.time.size) / 1000
                y_min, y_max, bins_y = height_min, height_max, bins_height
                if 'units' in syn_nc.height.attrs:  # else: TS in m
                    if syn_nc.height.units == 'km':
                        y = y * 1000

            mom_min_outer = 2 * mom_min - mom_max  # TODO
            mom_max_outer = 2 * mom_max - mom_min  # TODO
            # bins_y_3=bins_y*3 # TODO
            # mask = (y >= y_min) & (y <= y_max) & \
            #        (mom >= mom_min) & (mom <= mom_max)
            # mask = (y >= y_min) & (y <= y_max) & \
            #        (mom >= mom_min_outer) & (mom <= mom_max_outer) # TODO
            # TODO: realistic values might be set nan
            mask = (y >= y_min) & (y <= y_max) & \
                   (mom >= -999) & (mom <= 999)
            mom = mom[mask]
            y = y[mask]
            i_sort = np.argsort(y)
            mom = mom[i_sort]
            y = y[i_sort]
            if ax is None:
                plt.figure(figsize=(6, 5))

            a = plt.hist(y, bins=bins_y, range=(y_min, y_max))  # TODO
            # a = plt.hist(y, bins=bins_y-1, range=(y_min,y_max)) # TODO
            # a = plt.hist(y, bins=bins_y_3-1, )
            weights = np.repeat(100 / a[0], np.int16(a[0]))
            mom_all = np.append(mom_all, mom)
            y_all = np.append(y_all, y)
            weights_all = np.append(weights_all, weights)

    weights_all = weights_all / (n_cases)
    # PLOT
    if vmax:
        extend = 'max'
    else:
        extend = 'neither'

    # h2d, mom2d, y2d, fg = plt.hist2d(mom_all, y_all, bins=[bins_mom, bins_y],
    #                                  range=[[mom_min, mom_max],
    h2d, mom2d, y2d, fg = plt.hist2d(mom_all, y_all,
                                     bins=[bins_mom * 3, bins_y],  # TODO
                                     range=[[mom_min_outer, mom_max_outer],
                                            # TODO
                                            [y_min, y_max]], vmax=vmax,
                                     weights=weights_all, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]', extend=extend)
    if vert_temp:
        plt.gca().invert_yaxis()
        plt.ylabel('temperature (°C)')
    else:
        plt.ylabel(r'height (km)')

    if 'D0_r_Bringi' in moment:
        plt.xlabel('D0_r (mm)')
    elif 'vol_qntotice' in moment:
        plt.xlabel('$N_t\,(log_{10}(L^{-1}))$')
    elif 'qtotice' == moment:
        plt.xlabel('IWC$\,(g/m^{3})$')
    elif 'iwc_obs' in moment:
        plt.xlabel('IWC$\,(g/m^{3})$')
    elif 'dm_obs' in moment:
        plt.xlabel('$D_m$ (mm)')
    elif 'nt_obs' in moment:
        plt.xlabel('$N_t\,(log_{10}(L^{-1}))$')
    else:
        if 'standard_name' in syn_nc[moment].attrs:  # else: TS without any
            plt.xlabel(
                syn_nc[moment].standard_name + ' (' + syn_nc[
                    moment].units + ')')
        else:
            plt.xlabel(moment)

    plt.title(title)
    y_mid = y2d[1:] / 2 + y2d[:-1] / 2
    mom_mid = mom2d[1:] / 2 + mom2d[:-1] / 2
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid, weights=h2d[:, t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    plt.plot(quant_prof[0,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.2}$')
    plt.plot(quant_prof[1,], y_mid, color='red', ls='solid',
             linewidth=2, label='$Q_{0.5}$')
    plt.plot(quant_prof[2,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.8}$')
    plt.plot(mean_prof, y_mid, color='orange', ls='solid',
             linewidth=2, label='$\mu$')
    plt.legend()
    plt.xlim([mom_min, mom_max])  # TODO
    if vmax and np.max(h2d) > vmax:
        if np.max(h2d) > 100:
            print('above 100% :' + str(np.max(h2d)))

        plt.text(1.04, 1.03, min(np.round(np.max(h2d), 1), 100.0),
                 transform=ax.transAxes)
    syn_nc.close()
    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)


def calc_vradh(data):
    data['vradh'] = data.u * \
                    np.sin(np.deg2rad(data.azimuth)) * \
                    np.cos(np.deg2rad(data.elevation)) + \
                    data.v * \
                    np.cos(np.deg2rad(data.azimuth)) * \
                    np.cos(np.deg2rad(data.elevation)) + \
                    data.w * \
                    np.sin(np.deg2rad(data.elevation))
    return data
