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
"""
@author: jgiles
Script for noise-correcting RHOHV.
# """

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
import os
import xarray as xr
from wradlib.dp import kdp_from_phidp

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils


# Function from radarmet by Kai
def xr_rolling(da, window, window2=None, method="mean",
               min_periods=2, **kwargs):
    """Apply rolling function `method` to 2D datasets

    Parameter
    ---------
    da : xarray.DataArray
        array with data to apply rolling function
    window : int
        size of window in range dimension

    Keyword Arguments
    -----------------
    window2 : int
        size of window in azimuth dimension
    method : str
        function name to apply
    min_periods : int
        minimum number of valid bins
    **kwargs : dict
        kwargs to feed to rolling function

    Return
    ------
    da_new : xarray.DataArray
        DataArray with applied rolling function
    """
    prng = window // 2
    srng = slice(prng, -prng)
    da_new = da.pad(range=prng, mode="reflect", reflect_type="odd")

    dim = dict(range=window)
    isel = dict(range=srng)

    if window2 is not None:
        paz = window2 // 2
        saz = slice(paz, -paz)
        da_new = da_new.pad(azimuth=paz, mode="wrap")
        dim.update(dict(azimuth=window2))
        isel.update(dict(azimuth=saz))

    rolling = da_new.rolling(dim=dim, center=True, min_periods=min_periods)

    da_new = getattr(rolling, method)(**kwargs)
    da_new = da_new.isel(**isel)
    return da_new


# Velibor Pejcic
def phase_offset(phioff, rng=3000):  # 3000.0):
    """Calculate Phase offset.

    Parameter
    ---------
    phioff : xarray.DataArray
        differential phase array

    Keyword Arguments
    -----------------
    rng : float
        range in m to calculate system phase offset

    Return
    ------
    start_range : xarray.DataArray
        DataArray with start range values
    off : xarray.DataArray
        DataArray with phase offset values
    """
    range_step = np.diff(phioff.range)[0]
    nprec = int(rng / range_step)
    if nprec % 2:
        nprec += 1

    # create binary array
    phib = xr.where(np.isnan(phioff), 0, 1)

    # take nprec range bins and calculate sum
    phib_sum = phib.rolling(range=nprec, center=True).sum(skipna=True)

    # get start range of first N consecutive precip bins
    start_range = phib_sum.idxmax(dim="range") - \
                  nprec // 2 * np.diff(phib_sum.range)[0]
    # get range of first non-nan value per ray
    # start_range = (~np.isnan(phioff)).idxmax(dim='range', skipna=True)
    # add range
    stop_range = start_range + rng
    # get phase values in specified range
    off = phioff.where(
        (phioff.range >= start_range) & (phioff.range <= stop_range),
        drop=True
    )
    # calculate nan median over range
    off = off.median(dim="range", skipna=True)
    return xr.Dataset(
        dict(PHIDP_OFFSET=off, start_range=start_range, stop_range=stop_range)
    )


# Velibor Pejcic
def proc_phidp_kdp(swp_cf, uh_tresh=0, rho_tresh=0.9, win_r=25,
                   win_azi=None, wkdp_light=9, wkdp_heavy=25):
    """
    Processing Phidp and KDP

    Input:
    ------
    swp_cf ::: Quality controlled sweep

    uh_tresh ::: ZH Threshold
    rho_trsh ::: RHOHV Threshold

    win_r ::: Window size for 2d Medianfilter in range
    win_azi ::: Window size for 2d Medianfilter in azimuth

    wkdp_light ::: Window for KDP derivation (light)
    wkdp_heavy  ::: Window for KDP derivation (heavy)

    Output:
    -------

    swp_cf ::: Sweep with filterd/smoothed and system offest corrected PHIDP
    and combined KDP (see Park et al. 2009)

    """

    # Thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) & (
            swp_cf.RHOHV > rho_tresh) & np.isnan(
        swp_cf.CMAP))  # TODO: rho_nc !!!!

    # Median filtering 2d
    phimed = swp_mask.UPHIDP.copy()
    # phimed = phimed.pipe(rm.filter_data, medwin)

    # Smoothing
    # gaussian convolution 1d - smoothing
    # gkern = rm.gauss_kernel(kwidth, sigma)
    # phiclean = phimed.pipe(rm.smooth_data, gkern)

    # Median 2D
    window = win_r
    window2 = win_azi
    phi_median = phimed.copy().pipe(xr_rolling, window, window2=window2,
                                    method="median", skipna=True,
                                    min_periods=3)

    phioff = swp_mask.UPHIDP.copy().where(
        (swp_cf.RHOHV >= 0.9))  # & (swp.DBZH>=0))
    off = phase_offset(phioff, 1000.0)
    phi_offset = off.PHIDP_OFFSET.load().median(skipna=True)

    phi_corr = phi_median - phi_offset
    # phidp/kdp
    kdp_light = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_light)
    kdp_heavy = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_heavy)

    # phisp1 = kdp1.pipe(rm.phidp_from_kdp)
    # phisp2, kdp2 = (phi_median -
    #                 phi_offset).wrl.dp.phidp_kdp_vulpiani(winlen=23)

    kdp_comb = kdp_heavy.where(swp_cf.DBZH < 40, kdp_light)

    swp_cf = swp_cf.assign(kdp_corr=kdp_comb)
    swp_cf = swp_cf.assign(phi_corr=phi_corr)

    return swp_cf


# --------------------------------------------------------------------------- #

DATES = ["20210604",  # case01
         "20210620", "20210621",  # case02
         "20210628", "20210629",  # case03
         "20220519", "20220520",  # case04
         "20220623", "20220624", "20220625",  # case05
         "20220626", "20220627", "20220628",  # case06+07
         "20220630", "20220701",  # case08
         "20210714",  # case09
         "20221222",  # case10
         ]
LOCATIONS = ['asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld', 'hnr', 'isn',
             'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd',
             ]
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']
# overwrite = True
overwrite = False
winlen = 13
# import time
# time.sleep(60*60*24)

# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# LOCATIONS = ['pro']
# ELEVATIONS = np.array([12])
# # MODE = ['pcp']

date = '20210604'
# date = '20210714'
location = 'pro'
elevation_deg = 5.5
mode = 'vol'

# for date in DATES:
#     for location in LOCATIONS:
#         for mode in MODE:
#             for elevation_deg in ELEVATIONS:
for i_t_a, i_t_b in zip([0, 144], [143, 288]):

    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        print('pcp only 00')
        # continue

    path_in = "/".join([header.dir_data_obs + '*' + date,
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_allmoms_*'])
    files = sorted(glob.glob(path_in))
    if not files:
        path_in = "/".join([header.dir_data_obs +
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep,
                            'ras*_allmoms_*'])
        files = sorted(glob.glob(path_in))
    if not files:
        print('no input data *_allmoms_*')
        # continue
    else:
        path_in = files[0]
        path_out = path_in.replace('_allmoms_', '_kdp_nc_'+ str(i_t_a))

    data = dttree.open_datatree(path_in)[
        # 'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    rem_var = list(data.data_vars.keys())
    rem_var.remove('CMAP')
    rem_var.remove('UPHIDP')
    # data = data.drop_vars(rem_var)
    # dummy_tra = np.empty(shape=[data.time.size, data.range.size,
    #                                 data.azimuth.size, ])
    data = data.transpose('time', 'azimuth', 'range')
    # dummy_tra[:] = np.nan
    data = data.isel(time=slice(i_t_a, i_t_b))
    kdp_new = kdp_from_phidp(data['UPHIDP'].values,
                             winlen=winlen,  # cband
                             # winlen=31, # xband ?!
                             min_periods=3)
    # data['KDP_NC'] = (['time', 'range', 'azimuth'], kdp_new,
    data['KDP_NC'] = (['time', 'azimuth', 'range'], kdp_new,
                      dict(standard_name='specific differential phase',
                           comments='KDP noise corrected with winlen=' +
                                    str(winlen),
                           units='deg/km'))

    # swp_cf = dttree.open_datatree(path_in)
    # swp_cf = dttree.open_datatree(path_in)[
    #     'sweep_' + str(int(sweep))].to_dataset().chunk('auto')

    # data_new=proc_phidp_kdp(swp_cf, uh_tresh=0, rho_tresh=0.9, win_r=25,
    #                         win_azi=None, wkdp_light=9, wkdp_heavy=25)

    data = proc_phidp_kdp(data, uh_tresh=0, rho_tresh=0.9, win_r=25,
                          win_azi=None, wkdp_light=9, wkdp_heavy=25)

    mom_use = [x for x in list(data.keys())]
    for mom in mom_use:
        # data[mom].encoding["coordinates"] = "time range azimuth"
        data[mom].encoding["coordinates"] = "time azimuth range"

    data = data.drop_vars(rem_var)
    dtree = dttree.DataTree(name="root")
    dttree.DataTree(data, name=f"sweep_{int(sweep)}", parent=dtree)
    print('saving: ... ' + path_out + ' ...')
    dtree.load().to_netcdf(path_out)
    data.close()
    print('saved:  ' + path_out + ' !')

    # mom_use = [x for x in list(data_new.keys())]
    # for mom in mom_use:
    #     data_new[mom].encoding["coordinates"] = "time range azimuth"

    # dtree2 = dttree.DataTree(name="root")
    # dttree.DataTree(data_new, name=f"sweep_{int(sweep)}", parent=dtree2)
    # print('saving: ... ' + path_out.replace('_kdp_nc_', '_kdp_nc_new_') + ' ...')
    # dtree2.load().to_netcdf(path_out.replace('_kdp_nc_', '_kdp_nc_new_'))
    # print('saved:  ' + path_out.replace('_kdp_nc_', '_kdp_nc_new_') + ' !')

    # neu=data_new.isel(time=1)
    # dtree3 = dttree.DataTree(name="root")
    # dttree.DataTree(neu, name=f"sweep_{int(sweep)}", parent=dtree3)
    # dtree3.load().to_netcdf(path_out.replace('_kdp_nc_', '_kdp_nc_neu3_'))
