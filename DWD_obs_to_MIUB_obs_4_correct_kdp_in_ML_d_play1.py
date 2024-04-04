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
from wradlib.dp import kdp_from_phidp
from xhistogram.xarray import histogram
# from scipy.ndimage.filters import uniform_filter1d
# from scipy.ndimage import uniform_filter1d
from scipy.ndimage import uniform_filter, gaussian_filter

# from scipy.signal import savgol_filter
# import scipy as sp

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils


# --------------------------------------------------------------------------- #

# Function from radarmet by Kai
def xr_rolling(da, window, window2=None, method="median",
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

    # TODO: new: throw out if too noisy
    rolling = da_new.rolling(dim=dim, center=True, min_periods=min_periods)
    da_md = getattr(rolling, 'median')(**kwargs)
    da_new = da_new.where((abs(da_md - da_new) < 5))
    # TODO: new: throw out if too noisy

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
def phase_offset(phioff, rng=3000):
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
def proc_phidp_kdp(swp_cf, uh_tresh=0, rho_tresh=0.8, snr_tresh=10,
                   win_r=25, win_azi=None, wkdp_light=9, wkdp_heavy=25,
                   rng=3000):
    """
    Processing Phidp and KDP

    Input:
    ------
    swp_cf ::: Quality controlled sweep

    uh_tresh ::: ZH Threshold
    rho_tresh ::: RHOHV Threshold
    snr_tresh ::: SNR Threshold

    win_r ::: Window size for 2d Medianfilter in range
    win_azi ::: Window size for 2d Medianfilter in azimuth

    wkdp_light ::: Window for KDP derivation (light)
    wkdp_heavy  ::: Window for KDP derivation (heavy)

    rng : float
        range in m to calculate system phase offset

    Output:
    -------

    swp_cf ::: Sweep with filterd/smoothed and system offest corrected PHIDP
    and combined KDP (see Park et al. 2009)

    """

    # 0: thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) &
                            (swp_cf.RHOHV > rho_tresh) &
                            (swp_cf.SNRH > snr_tresh) &
                            np.isnan(swp_cf.CMAP))

    # 1: centering PHI first guess (around modus; reduced by 100):
    phi_hist = np.histogram(swp_mask.UPHIDP.copy(),
                            bins=np.linspace(-180, 180, 361))
    phi_modus = np.argmax(phi_hist[0]) - 180
    phi_c = ((swp_mask.UPHIDP.copy() - xr.DataArray(
        phi_modus) - 100) % 360 + 180) % 360 - 180

    # 2: smoothing PHI along range (reduced by 100):
    phi_s = phi_c.copy().pipe(
        xr_rolling, window=win_r, window2=win_azi, method="median",
        skipna=True, min_periods=max(3, int((win_r - 1) / 4))  # TODO: 3 o 1/6?
    )

    # 3: offset Part 1 of 7: calculate offset 2d (reduced by 100):
    phi_off_2d = phase_offset(phi_s.copy().where(swp_cf.range > 1000),
                              rng).PHIDP_OFFSET.load()

    # 3: offset Part 2 of 7: filter outl. (>10 in neighbourhood; red. by 100):
    phi_off_2d_f = xr.where(np.isnan(phi_off_2d), -100, phi_off_2d)
    phi_off_2d_f = xr.DataArray(uniform_filter(
        phi_off_2d_f, size=5, mode='mirror'), dims=['time', 'azimuth'])
    phi_off_2d_f = phi_off_2d.where((abs(phi_off_2d - phi_off_2d_f)) < 10)

    # 3: offset Part 3 of 7:  interpolate na's (red. by 100):
    phi_off_2d_f_i = phi_off_2d_f.interpolate_na(dim='azimuth', limit=5)
    phi_off_2d_f_i = xr.where(np.isnan(phi_off_2d_f),
                              phi_off_2d_f_i, phi_off_2d_f)

    # 3: offset Part 4 of 7:  median offsets per time (red. by 100):
    phi_off_1d = phi_off_2d_f_i.median("azimuth", skipna=True)
    phi_off_0d = phi_off_1d.median("time", skipna=True)
    phi_off_1d = xr.where(np.isnan(phi_off_1d), phi_off_0d, phi_off_1d)

    # 3: offset Part 5 of 7:  filling gaps with 0 (red. by 100)
    phi_off_2d_f_i_f = phi_off_2d_f_i - phi_off_1d
    phi_off_2d_f_i_f = xr.where(np.isnan(phi_off_2d_f_i_f), 0,
                                phi_off_2d_f_i_f)
    phi_off_2d_f_i_f = phi_off_2d_f_i_f + phi_off_1d

    # 3: offset Part 6 of 7:  smoothing
    phi_off_2d_f_i_f_s = xr.DataArray(
        gaussian_filter(phi_off_2d_f_i_f, 3, mode='wrap'),
        dims=['time', 'azimuth'])

    # 3: offset Part 7 of 7:  noise corrected phi (NOT red. by 100 anymore!)
    phi_nc = phi_s - phi_off_2d_f_i_f_s

    # # 3b: check if phi needs to be reverted:
    phi_diffs = phi_nc.copy().diff('range', 1)
    phi_diffs = xr.where(abs(phi_diffs) > 1, np.nan, phi_diffs)
    phi_trues = phi_diffs.notnull().sum(["range", "azimuth"])
    phi_factor = phi_diffs.median(["range", "azimuth"], skipna=True)
    phi_factor = xr.where(phi_factor < 0, -1, 1)
    phi_factor = xr.where(phi_trues < 500, 1, phi_factor)  # TODO: 500???
    if sum(phi_factor.values) < phi_factor.size:
        print('REVERING PHI!')
        print(phi_trues.values)
        print(phi_factor.values)
        swp_cf["UPHIDP"] = swp_cf.UPHIDP * phi_factor
        return proc_phidp_kdp(swp_cf, uh_tresh, rho_tresh, snr_tresh,
                              win_r, win_azi, wkdp_light, wkdp_heavy, rng)

    # 4: KDP
    kdp_light = phi_nc.wrl.dp.kdp_from_phidp(winlen=wkdp_light)
    kdp_heavy = phi_nc.wrl.dp.kdp_from_phidp(winlen=wkdp_heavy)
    kdp_comb = kdp_heavy.where(swp_cf.DBZH < 40, kdp_light)

    # phi/kdp attributes
    phi_nc.attrs["long_name"] = 'Differential phase shift'
    phi_nc.attrs["short_name"] = 'PHI_DP'
    phi_nc.attrs["units"] = 'degrees'
    phi_nc.attrs["comments"] = 'PHI_DP smoothing with win_r=' + \
                               str(win_r) + ' and win_azi=' + str(win_azi)
    kdp_comb.attrs["comments"] = 'KDP noise corrected with winlen=' + \
                                 str(wkdp_heavy) + ' (DBZH<40) and ' + \
                                 'winlen=' + str(wkdp_light) + ' (DBZH>=40)'

    # assign!
    swp_cf = swp_cf.assign(KDP_NC=kdp_comb)
    swp_cf = swp_cf.assign(PHI_NC=phi_nc)

    # assign (maybe leave out)?
    swp_cf = swp_cf.assign(phi_c=phi_c + 100)
    swp_cf = swp_cf.assign(phi_s=phi_s + 100)

    swp_cf = swp_cf.assign(PHI_off_2d=phi_off_2d + 100)
    swp_cf = swp_cf.assign(phi_off_2d_f=phi_off_2d_f + 100)
    swp_cf = swp_cf.assign(phi_off_2d_f_i=phi_off_2d_f_i + 100)
    swp_cf = swp_cf.assign(phi_off_2d_f_i_f=phi_off_2d_f_i_f + 100)
    swp_cf = swp_cf.assign(phi_off_2d_f_i_f_s=phi_off_2d_f_i_f_s + 100)

    swp_cf = swp_cf.assign(phi_diffs=phi_diffs)

    return swp_cf


# --------------------------------------------------------------------------- #
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
# --------------------------------------------------------------------------- #
DATES = ["20210604",  # case01
         "20210714",  # case09
         # "20210620", "20210621",  # case02
         # "20210628", "20210629",  # case03
         # "20220519", "20220520",  # case04
         # "20220623", "20220624", "20220625",  # case05
         # "20220626", "20220627", "20220628",  # case06+07
         # "20220630", "20220701",  # case08
         # "20221222",  # case10
         ]
# LOCATIONS = ['asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld', 'hnr', 'isn',
#              'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd',
#              ]
LOCATIONS = ['ess', 'pro', 'umd', 'tur', 'asb', 'boo', 'drs', 'eis', 'fbg',
             # 'fld', 'hnr', 'isn',
             # 'mem', 'neu', 'nhb', 'oft', 'ros',
             ]
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']

merge = True
remove_parts = True
# overwrite = False
overwrite = True

snr_tresh = 15  # TODO: VP threshold
win_r = 25  # TODO: VP or PARK?!
# win_r = 31  # TODO: new
win_azi = None  # TODO: None, instead Gaussian filter for azimuth and time
# win_azi = 3  # TODO: None, instead Gaussian filter for azimuth and time
wkdp_light = 9  # TODO: PARK: light filtering
wkdp_heavy = 25  # TODO: PARK: heavy filtering
rng = 3000  # TODO: VP
rng = 1000  # TODO:  # next ckeck (together with range>1000) # before was good!!!!!!!!

# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# # DATES = ['20210714']
DATES = ['20210604', '20210714', ]
# # LOCATIONS = ['pro']
# # LOCATIONS = ['ess']
# # LOCATIONS = ['umd']
LOCATIONS = ['ess', 'tur', ]
# # LOCATIONS = ['tur']
# ELEVATIONS = np.array([12])
# MODE = ['vol']

# # TODO: 18.3.24 10:35
# DATES = ['20210604', ]
# ELEVATIONS = np.array([4.5, 12, 0.5, ])
# LOCATIONS = ['ess', ]
# MODE = ['vol', ]

# # # TODO: 21.3.24 15:40
# DATES = ['20210714', '20210604',]
DATES = ['20210604', ]
LOCATIONS = ['ess', ]
MODE = ['vol', ]
ELEVATIONS = np.array([12, 0.5])

# date = '20210604'
# # date = '20210714'
# location = 'pro'
# # location = 'ess'
# elevation_deg = 5.5
# mode = 'vol'

# # # TODO: 27.3.24 17:00
DATES = ['20210714', '20210604', ]
LOCATIONS = ['ess', 'tur', 'pro', ]
MODE = ['vol', ]
ELEVATIONS = np.array([12, 0.5])

# # # # TODO: 28.3.24 9:00
DATES = ['20210714', '20210604', ]
LOCATIONS = ['pro', 'ess', 'tur',]
MODE = ['vol', ]
ELEVATIONS = np.array([0.5, 12])

for date in DATES:
    for location in LOCATIONS:
        # if date == '20210604' and not location == 'ess':  #TODO: remove
        #     print(date + '_' + location)
        #     continue

        for mode in MODE:
            for elevation_deg in ELEVATIONS:
                print('start: ' + date + ' ' + location + ' ' +
                      str(elevation_deg))
                # parts = 288
                parts = 4
                merge_files = []
                for p in range(parts):
                # for p in [169,170,171]:
                # for p in [212,221]:
                    i_t_a = int(288 / parts * p)
                    i_t_b = int(288 / parts * (p + 1))
                    year = date[0:4]
                    mon = date[4:6]
                    day = date[6:8]
                    sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                                               float(elevation_deg))[0][0])
                    if mode == 'pcp' and sweep != '00':
                        print('pcp only 00')
                        break

                    path_in = "/".join([header.dir_data_obs + '*',
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
                        break
                    else:
                        path_in = files[0]
                        path_out = path_in.replace(
                            '_allmoms_', '_' +
                                         str(pd.Timestamp.today())[11:16] +
                                         '_kdp_nc_d_' + str(p) + '_')

                    if (os.path.isfile(path_out) and not overwrite) or \
                            (os.path.isfile(path_out.replace(
                                'kdp_nc_d_' + str(p), 'kdp_nc'))
                             and not overwrite):
                        print(path_out + ' exists;\n' + ' ... set: > ' +
                              'overwrite = True < for recalculation')
                        continue

                    path_rho_nc = path_in.replace('_allmoms_', '_rhohv_nc_')
                    data = dttree.open_datatree(path_in)[
                        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
                    data_rho = dttree.open_datatree(path_rho_nc)[
                        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
                    data.RHOHV.values = data_rho.RHOHV_NC2P.values
                    data = data.assign({'SNRH': data_rho.SNRH})
                    remo_var = list(data.data_vars.keys())
                    remo_var.remove('CMAP')
                    remo_var.remove('UPHIDP')
                    data = data.transpose('time', 'azimuth', 'range')
                    data = data.isel(time=slice(i_t_a, i_t_b))
                    if data.time.size == 0:
                        parts = p
                        break

                    merge_files.append(path_out)
                    data = proc_phidp_kdp(data, uh_tresh=0, rho_tresh=0.8,
                                          # TODO 0.9 is fine as well
                                          snr_tresh=snr_tresh,
                                          win_r=win_r, win_azi=win_azi,
                                          wkdp_light=wkdp_light,
                                          wkdp_heavy=wkdp_heavy,
                                          rng=rng)
                    mom_use = [x for x in list(data.keys())]
                    for mom in mom_use:
                        data[mom].encoding["coordinates"] = \
                            "time azimuth range"

                    data = data.drop_vars(remo_var)
                    dtree = dttree.DataTree(name="root")
                    dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                                    parent=dtree)
                    print('saving: ... ' + path_out + ' ...')
                    dtree.load().to_netcdf(path_out)
                    data.close()
                    print('saved:  ' + path_out + ' !')

                if merge and merge_files != []:
                    path_out_new = merge_files[0].replace(
                        'kdp_nc_d_0_', 'kdp_nc_'
                    ).replace('ras', str(pd.Timestamp.today())[
                                     5:10] + '-' +  # TODO: testing!
                              str(pd.Timestamp.today())[11:16] + '_ras'
                              # TODO: testing!
                              )
                    if os.path.isfile(path_out_new) and not overwrite:
                        print(path_out_new + ' exists;\n' + ' ... set: ' +
                              '> overwrite = True < for recalculation')
                    else:
                        data_merged = xr.merge([
                            dttree.open_datatree(merge_files[p])[
                                'sweep_' + str(int(sweep))].to_dataset(
                            ).chunk(-1) for p in range(parts)])
                        mom_use = [x for x in list(data_merged.keys())]
                        for mom in mom_use:
                            data_merged[mom].encoding["coordinates"] = \
                                "time azimuth range"

                        data_merged.attrs['processing_date'] = str(
                            pd.Timestamp.today())[:16]
                        print('saving: ... ' + path_out_new + ' ...')
                        dtree = dttree.DataTree(name="root")
                        dttree.DataTree(data_merged,
                                        name=f"sweep_{int(sweep)}",
                                        parent=dtree)
                        dtree.load().to_netcdf(path_out_new)
                        data_merged.close()
                        print('saved:  ' + path_out_new + ' !')
                        if remove_parts:
                            for file_part in merge_files:
                                os.remove(file_part)
