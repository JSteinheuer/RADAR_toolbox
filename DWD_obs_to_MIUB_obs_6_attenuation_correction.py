#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 27.05.24                                                 #
# DWD_obs_to_MIUB_obs_6_attenuation_correction.py                             #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 6: attenuation correction.                                             #
#         see Ryzhkov and Zrnic (2019): 6.4 (pp. 162, 167)                    #
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


# J. Steinheuer
def attenuation_correction(swp_cf, alpha=0.08, beta=0.02, until_temp_bottom=4):
    """
    ...
    """
    swp_mask = swp_cf.where((swp_cf.DBZH > 0) &
                            (swp_cf.DBZH < 40) &
                            (swp_cf.RHOHV > 0.98) &
                            np.isnan(swp_cf.CMAP))
    zdroffset = np.nanmedian(swp_mask.ZDR)
    nm = np.sum(~np.isnan(swp_mask.ZDR.values))
    if plot:
        if ax is None:
            plt.figure(figsize=(4, 3))
            ax = plt.subplot(1, 1, 1)
        # ax.hist(swp_cf.ZDR.values.flatten(),
        #          bins=np.arange(-1, 3, .025), density=True, alpha=.5)
        ax.hist(swp_mask.ZDR.values.flatten(),
                 bins=np.arange(-1, 3, .025))
        ax.axvline(zdroffset, color='black')
        ax.set_title('Birdbath $90°$')
        ax.set_xlabel(r'$Z_{DR}$', fontsize=15)
        ax.set_ylabel(r'$counts$', fontsize=15)
        ax.grid(which='both', color='black', linestyle=':', alpha=0.5)
        ax.legend(title=r'$\Delta Z_{DR}$: ' + str(np.round(zdroffset, 3)) +
                        'dB\n' + r'$N$: ' + str(nm))
        plt.tight_layout()
        # plt.show()

    return zdroffset, nm


# --------------------------------------------------------------------------- #
DATES = [
    "20210714",  # case09
    "20210604",  # case01
    "20210620", "20210621",  # case02
    "20210628", "20210629",  # case03
    "20220519", "20220520",  # case04
    "20220623", "20220624", "20220625",  # case05
    "20220626", "20220627", "20220628",  # case06+07
    "20220630", "20220701",  # case08
    "20221222",  # case10
]
LOCATIONS = [
    'asb',
    'boo', 'drs', 'eis',
    'ess',
    'fbg', 'fld', 'hnr', 'isn',
    'mem', 'neu', 'nhb',
    'oft',
    'pro', 'ros', 'tur', 'umd',
]
# --------------------------------------------------------------------------- #
overwrite = False

# case (adjust!):
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# sorting volume
elevation_degs_2_sort = elevation_degs.copy()
elevation_degs_2_sort[modes == 'pcp'] = 0
elevation_degs_2_sort[modes == '90grad'] = 90
sort_i = elevation_degs_2_sort.argsort()
elevation_degs = elevation_degs[sort_i]
modes = modes[sort_i]
# date = DATES[0]
# location = LOCATIONS[0]
for date in DATES:
    for location in LOCATIONS:
        # plot parameters
        if plot:
            index = 0
            n_rows = 4
            n_cols = elevation_degs.size
            fig = plt.figure(figsize=(3 * n_cols, 3 * n_rows))
            save = True

        # ------------------------------------------------------------------- #
        # loop over all elevations:
        # ------------------------------------------------------------------- #
        for elevation_deg, mode in zip(elevation_degs, modes):
            print(elevation_deg)
            # ----------------------------------------------------------------#
            # folder and file search
            folder_plot = header.folder_plot
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

            path_out_nc = nc_file_mom.replace('_allmoms_', '_zdr_off_')
            if overwrite or not os.path.exists(path_out_nc) or plot:
                if elevation_deg == elevation_degs[0] and mode == modes[0]:
                    file_in = nc_file_mom.split('/')[-1]
                    if plot:
                        file_out = file_in.replace('.hd5', '_' + str(n_rows) +
                                                   'x' + str(n_cols) + '.' +
                                                   pdf_or_png).replace('_' +
                                                                       sweep,
                                                                       '_all')
                        path_out = folder_plot + 'ZDR_calibration/' + \
                                   location.upper() + '/ZDR_calibration_' + \
                                   str(n_rows) + 'x' + str(n_cols) + '_' + \
                                   file_out
                        if os.path.isfile(path_out) and not overwrite:
                            print(path_out + ' exists;\n' + ' ... set: ' +
                                  '> overwrite = True < for recalculation')
                            save = False
                            break

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
                # ----------------------------------------------------------- #
                # time
                time_start = date + nc_file_mom.split('-')[-4][-4:]
                time_end = date + nc_file_mom.split('-')[-3][-4:]
                dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
                dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
                dti = pd.date_range(dti_start, dti_end, freq="5min",
                                    inclusive='both')
                # ----------------------------------------------------------- #
                # ----------------------------------------------------------- #
                if mode == '90grad':
                    data = dttree.open_datatree(nc_file_mom)[
                        'sweep_' + str(int(sweep))].to_dataset()
                    # ------------------------------------------------------- #
                    # plot calibration
                    if plot:
                        index = index + 1
                        ax = plt.subplot(n_rows, n_cols, index + 0*n_cols)
                    else:
                        ax = None

                    bb_off, bb_nm = cal_zdr_birdbath(data, plot=plot, ax=ax)
                    # saving
                    path_out_nc = nc_file_mom.replace(
                        '_allmoms_', '_zdr_off_')
                    if not overwrite and os.path.exists(path_out_nc):
                        print('exists: ' + path_out_nc + ' -> continue')
                    else:
                        remo_var = list(data.data_vars.keys())
                        remo_var.remove('ZDR')
                        data = data.drop_vars(remo_var)
                        data['zdr_off_bb'] = bb_off
                        data['zdr_off_bb'].attrs["long_name"] = \
                            'ZDR offset from bird bath'
                        data['zdr_off_bb'].attrs["short_name"] = 'ZDR off BB'
                        data['zdr_off_bb'].attrs["units"] = 'dB'
                        data['zdr_off_bb'].attrs["comment"] = \
                            'to subtract from ZDR'

                        data['zdr_off_bb_n'] = bb_nm
                        data['zdr_off_bb_n'].attrs["long_name"] = \
                            'number values for BB'
                        data['zdr_off_bb_n'].attrs["short_name"] = 'nm for BB'
                        data['zdr_off_bb_n'].attrs["units"] = '1'

                        dtree = dttree.DataTree(name="root")
                        dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                        parent=dtree)
                        print('saving: ... ' + path_out_nc.split('/')[-1] +
                              ' ...')
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
                    data = data.assign({'temp_beamtop':
                                            data_temp2.temp_beamtop})
                    remo_var = list(data.data_vars.keys())
                    data = data.transpose('time', 'azimuth', 'range')
                    # elevation contribution included:
                    data2 = zdr_at_zero_elev(data)
                    data2 = data.assign({'ZDR': data2.ZDR0})
                    # ------------------------------------------------------- #
                    # plot calibration
                    colorbar = [True, True]
                    axes = None
                    if plot:
                        index = index + 1
                        axes = [plt.subplot(n_rows, n_cols, index + 0*n_cols),
                                plt.subplot(n_rows, n_cols, index + 0*n_cols)]

                    lr_off, lr_nm = cal_zdr_lightrain(data, band='C',
                                                      plot=[plot, False],
                                                      axes=axes,
                                                      colorbar=colorbar)
                    if mode == 'pcp' and plot:
                        axes[0].set_title('PCP ' + location.upper() + ' ' +
                                          date + '    ')
                    elif plot:
                        axes[0].set_title(str(elevation_deg) + '° ' +
                                          location.upper() + ' ' + date +
                                          '    ')

                    if plot:
                        axes = [plt.subplot(n_rows, n_cols, index + 1*n_cols),
                                plt.subplot(n_rows, n_cols, index + 1*n_cols)]

                    sd_off, sd_nm = cal_zdr_smalldrops(data, band='C',
                                                       plot=[plot, False],
                                                       axes=axes,
                                                       colorbar=colorbar)
                    if plot:
                        axes = [plt.subplot(n_rows, n_cols, index + 2*n_cols),
                                plt.subplot(n_rows, n_cols, index + 2*n_cols)]

                    lr_off_ec, lr_nm_ec = cal_zdr_lightrain(data2, band='C',
                                                            plot=[plot, False],
                                                            axes=axes,
                                                            colorbar=colorbar)
                    if mode == 'pcp' and plot:
                        axes[0].set_title(
                            r'$Z_{DR}(PCP) \rightarrow Z_{DR}(0°)$')
                    elif plot:
                        axes[0].set_title(r'$Z_{DR}($' + str(elevation_deg) +
                                          r'$°) \rightarrow Z_{DR}(0°)$')

                    if plot:
                        axes = [plt.subplot(n_rows, n_cols, index + 3*n_cols),
                                plt.subplot(n_rows, n_cols, index + 3*n_cols)]

                    sd_off_ec, sd_nm_ec = cal_zdr_smalldrops(data2, band='C',
                                                             plot=[plot,
                                                                   False],
                                                             axes=axes,
                                                             colorbar=colorbar)
                    # saving
                    path_out_nc = nc_file_mom.replace(
                        '_allmoms_', '_zdr_off_')
                    if not overwrite and os.path.exists(path_out_nc):
                        print('exists: ' + path_out_nc + ' -> continue')
                    else:
                        remo_var = list(data.data_vars.keys())
                        remo_var.remove('ZDR')
                        data = data.drop_vars(remo_var)

                        data['zdr_off_lr'] = lr_off
                        data['zdr_off_lr'].attrs["long_name"] = \
                            'ZDR offset from light rain'
                        data['zdr_off_lr'].attrs["short_name"] = 'ZDR off LR'
                        data['zdr_off_lr'].attrs["units"] = 'dB'
                        data['zdr_off_lr'].attrs["comment"] = \
                            'to subtract from ZDR'

                        data['zdr_off_lr_n'] = lr_nm
                        data['zdr_off_lr_n'].attrs["long_name"] = \
                            'number values for LR'
                        data['zdr_off_lr_n'].attrs["short_name"] = 'nm for LR'
                        data['zdr_off_lr_n'].attrs["units"] = '1'

                        data['zdr_off_lr_ec'] = lr_off_ec
                        data['zdr_off_lr_ec'].attrs["long_name"] = \
                            'ZDR offset from light rain and elevation ' \
                            'corrected ZDR'
                        data['zdr_off_lr_ec'].attrs["short_name"] = \
                            'ZDR off LR ec'
                        data['zdr_off_lr_ec'].attrs["units"] = 'dB'
                        data['zdr_off_lr_ec'].attrs["comment"] = \
                            'to subtract from ZDR'

                        data['zdr_off_lr_n_ec'] = lr_nm
                        data['zdr_off_lr_n_ec'].attrs["long_name"] = \
                            'number values for LR ec'
                        data['zdr_off_lr_n_ec'].attrs["short_name"] = \
                            'nm for LR ec'
                        data['zdr_off_lr_n_ec'].attrs["units"] = '1'

                        data['zdr_off_sd'] = sd_off
                        data['zdr_off_sd'].attrs["long_name"] = \
                            'ZDR offset from small drops'
                        data['zdr_off_sd'].attrs["short_name"] = 'ZDR off SD'
                        data['zdr_off_sd'].attrs["units"] = 'dB'
                        data['zdr_off_sd'].attrs["comment"] = \
                            'to subtract from ZDR'

                        data['zdr_off_sd_n'] = sd_nm
                        data['zdr_off_sd_n'].attrs["long_name"] = \
                            'number values for SD'
                        data['zdr_off_sd_n'].attrs["short_name"] = 'nm for SD'
                        data['zdr_off_sd_n'].attrs["units"] = '1'

                        data['zdr_off_sd_ec'] = sd_off_ec
                        data['zdr_off_sd_ec'].attrs["long_name"] = \
                            'ZDR offset from small drops and elevation ' \
                            'corrected ZDR'
                        data['zdr_off_sd_ec'].attrs["short_name"] = \
                            'ZDR off SD ec'
                        data['zdr_off_sd_ec'].attrs["units"] = 'dB'
                        data['zdr_off_sd_ec'].attrs["comment"] = \
                            'to subtract from ZDR'

                        data['zdr_off_sd_n_ec'] = sd_nm
                        data['zdr_off_sd_n_ec'].attrs["long_name"] = \
                            'number values for SD ec'
                        data['zdr_off_sd_n_ec'].attrs["short_name"] = \
                            'nm for SD ec'
                        data['zdr_off_sd_n_ec'].attrs["units"] = '1'

                        dtree = dttree.DataTree(name="root")
                        dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                        parent=dtree)
                        print('saving: ... ' + path_out_nc.split('/')[-1] +
                              ' ...')
                        dtree.load().to_netcdf(path_out_nc)
                        data.close()
            else:
                print('exists: ' + path_out_nc + ' -> continue')
                save = False

        # ------------------------------------------------------------------- #
        # SAVE                                                                #
        # ------------------------------------------------------------------- #
        if save and plot:
            plt.tight_layout()
            if not os.path.exists(folder_plot + 'ZDR_calibration/' +
                                  location.upper() + '/'):
                os.makedirs(folder_plot + 'ZDR_calibration/' +
                            location.upper() + '/')

            plt.savefig(path_out, format=pdf_or_png, transparent=True)
            plt.close()
