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

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils


# --------------------------------------------------------------------------- #

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
def phase_offset(phioff, rng=1000):  # 3000.0):
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
def proc_phidp_kdp(swp_cf, uh_tresh=0, rho_tresh=0.9, snr_tresh=10,
                   win_r=25, win_azi=None, wkdp_light=9, wkdp_heavy=25):
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

    Output:
    -------

    swp_cf ::: Sweep with filterd/smoothed and system offest corrected PHIDP
    and combined KDP (see Park et al. 2009)

    """

    # Thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) &
                            (swp_cf.RHOHV > rho_tresh) &
                            (swp_cf.DBSNRH > snr_tresh) &
                            np.isnan(swp_cf.CMAP))

    # Median filtering 2d
    phimed = swp_mask.UPHIDP.copy()

    # Median 2D
    window = win_r
    window2 = win_azi
    phi_median = phimed.copy().pipe(xr_rolling, window, window2=window2,
                                    method="median", skipna=True,
                                    min_periods=3)

    phioff = swp_mask.UPHIDP.copy().where(
        (swp_cf.RHOHV >= 0.9))  # & (swp.DBZH>=0))
    off = phase_offset(phioff, 3000.0)  # TODO: nochmal Wellen fragen
    phi_offset = off.PHIDP_OFFSET.load().median("azimuth",skipna=True) # neu
    phi_corr = phi_median - phi_offset

    # phidp/kdp
    kdp_light = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_light)
    kdp_heavy = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_heavy)
    kdp_comb = kdp_heavy.where(swp_cf.DBZH < 40, kdp_light)

    # attributes
    kdp_comb.attrs["comments"] = 'KDP noise corrected with winlen=' + \
                                 str(wkdp_heavy) + ' (DBZH<40) and ' + \
                                 'winlen=' + str(wkdp_light) + ' (DBZH>=40)'
    phi_corr.attrs["long_name"] = \
        'Differential phase shift'
    phi_corr.attrs["short_name"] = \
        'PHI_DP'
    phi_corr.attrs["units"] = \
        'degrees'
    phi_corr.attrs["comments"] = 'PHI_DP smoothing with win_r=' + \
                                 str(win_r) + ' and win_azi=' + str(win_azi)
    swp_cf = swp_cf.assign(KDP_NC=kdp_comb)
    swp_cf = swp_cf.assign(PHI_NC=phi_corr)

    return swp_cf


# --------------------------------------------------------------------------- #
ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                           8.0, 12.0, 17.0, 25.0])
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
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']

merge = True
remove_parts = False
overwrite = True

win_r = 25
win_azi = 1
wkdp_light = 9
wkdp_heavy = 25

snr_tresh = 15
win_rs = [25, 9, 15]
wkdp_light = 9
wkdp_heavy = 25

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

# i_t_a = 0
# i_t_b = 143
# part = 'a'
for win_r in win_rs:
    merge_files = []
    for i_t_a, i_t_b, part in \
            zip([0, 144], [144, 288], ['a', 'b']):

        year = date[0:4]
        mon = date[4:6]
        day = date[6:8]
        sweep = '0' + str(np.where(ELEVATIONS_ALL ==
                                   float(elevation_deg))[0][0])
        if mode == 'pcp' and sweep != '00':
            print('pcp only 00')
            continue

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
            continue
        else:
            path_in = files[0]
            path_out = path_in.replace('_allmoms_', '_kdp_nc_' +
                                       part + '_' # )
                                       + str(win_r))  # TODO: remove

        merge_files.append(path_out)
        if (os.path.isfile(path_out) and not overwrite) or \
                (os.path.isfile(path_out.replace(
                    'kdp_nc_' + str(part), 'kdp_nc'))
                 and not overwrite):
            print(path_out + ' exists;\n' + ' ... set: > ' +
                  'overwrite = True < for recalculation')
            continue

        path_rho_nc = path_in.replace('_allmoms_', '_rhohv_nc_')
        data = dttree.open_datatree(path_in)[
            'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
        # 'sweep_' + str(int(sweep))].to_dataset().chunk(-1)

        data_rho = dttree.open_datatree(path_rho_nc)[
            'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
        # 'sweep_' + str(int(sweep))].to_dataset().chunk(-1)

        data.RHOHV.values = data_rho.RHOHV_NC2P.values
        rem_var = list(data.data_vars.keys())
        rem_var.remove('CMAP')
        rem_var.remove('UPHIDP')
        data = data.transpose('time', 'azimuth', 'range')
        data = data.isel(time=slice(i_t_a, i_t_b))
        data = proc_phidp_kdp(data, uh_tresh=0, rho_tresh=0.9,
                              snr_tresh=snr_tresh,
                              win_r=win_r, win_azi=win_azi,
                              wkdp_light=wkdp_light,
                              wkdp_heavy=wkdp_heavy)
        mom_use = [x for x in list(data.keys())]
        for mom in mom_use:
            # data[mom].encoding["coordinates"] = "time range azimuth"
            data[mom].encoding["coordinates"] = \
                "time azimuth range"

        data = data.drop_vars(rem_var)
        dtree = dttree.DataTree(name="root")
        dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                        parent=dtree)
        print('saving: ... ' + path_out + ' ...')
        dtree.load().to_netcdf(path_out)
        data.close()
        print('saved:  ' + path_out + ' !')

    if merge:
        path_out_new = merge_files[0].replace(
            'kdp_nc_a_', 'kdp_nc_')
        if os.path.isfile(path_out_new) and not overwrite:
            print(path_out_new + ' exists;\n' + ' ... set: ' +
                  '> overwrite = True < for recalculation')
        else:
            data_merged = xr.merge([
                dttree.open_datatree(merge_files[0])[
                    'sweep_' + str(int(sweep))].to_dataset(
                ).chunk(-1),
                dttree.open_datatree(merge_files[1])[
                    'sweep_' + str(int(sweep))].to_dataset(
                ).chunk(-1),
            ])
            mom_use = [x for x in list(data_merged.keys())]
            for mom in mom_use:
                # data[mom].encoding["coordinates"] = "time range azimuth"
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
