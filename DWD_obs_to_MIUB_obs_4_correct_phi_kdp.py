#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 23.01.24                                                 #
# DWD_obs_to_MIUB_obs_4_correct_phi_kdp.py                                    #
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

warnings.filterwarnings("ignore")


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

    # new: throw out values that are too noisy
    rolling = da_new.rolling(dim=dim, center=True, min_periods=min_periods)
    da_md = getattr(rolling, 'median')(**kwargs)
    da_new = da_new.where((abs(da_md - da_new) < 5))

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

    range_step = abs(np.diff(phioff.range)[0])

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
        drop=False
    )

    # calculate nan median over range
    # off = off.median(dim="range", skipna=True)
    # better: 1.Quartile
    off = off.chunk(dict(range=-1))
    off = off.quantile(q=0.25, dim="range", skipna=True)

    return xr.Dataset(
        dict(PHIDP_OFFSET=off,
             start_range=start_range,
             stop_range=stop_range)
    )


# Velibor Pejcic
def proc_phidp_kdp(swp_cf, uh_tresh=0, rho_tresh=0.8, snr_tresh=15,
                   win_r=25, win_azi=None, wkdp_light=9, wkdp_heavy=25,
                   rng=3000, flip_default=1):
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

    # 1: thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) &
                            (swp_cf.RHOHV > rho_tresh) &
                            (swp_cf.SNRH > snr_tresh) &
                            np.isnan(swp_cf.CMAP))

    # 2: check if phi needs to be flipped:
    phi_c0 = swp_mask.UPHIDP.where(swp_cf.RHOHV > 0.95)
    phi_diffs = ((phi_c0.where((swp_cf.range > 3000) &
                               (swp_cf.RHOHV > 0.95) &
                               (swp_cf.DBZH > 20)).diff('range', 1) +
                  180) % 360) - 180
    phi_diffs = xr.where(abs(phi_diffs) > 10, np.nan, phi_diffs)
    phi_flip = phi_diffs.sum(["range", "azimuth"], skipna=True)
    # phi_flip_0 = phi_diffs.sum(["time", "range", "azimuth"], skipna=True)
    # if abs(phi_flip_0) < 1000:
    #     phi_flip_0 = flip_default

    # 3: flipper
    phi_flip_i = phi_flip.copy()
    phi_flip = xr.where(abs(phi_flip) < 1000, flip_default, phi_flip)
    phi_flip = xr.where(phi_flip < 0, -1, phi_flip)
    phi_flip = xr.where(phi_flip > 0, 1, phi_flip)
    # phi_flip_0 = xr.where(phi_flip_0 < 0, -1, 1)
    phi_flip = xr.where(np.isnan(phi_flip), flip_default, phi_flip)
    phi_c0 = phi_c0 * phi_flip
    phi_c1 = swp_mask.UPHIDP * phi_flip

    # 4: centering PHI first guess (around modus; reduced by 100):
    phi_hist = np.histogram(phi_c0, bins=np.linspace(-180, 180, 361))
    phi_modus = np.argmax(phi_hist[0]) - 180
    phi_c0 = ((phi_c0 - xr.DataArray(
        phi_modus) - 100) % 360 + 180) % 360 - 180
    phi_c1 = ((phi_c1 - xr.DataArray(
        phi_modus) - 100) % 360 + 180) % 360 - 180

    # 5: smoothing PHI along range (reduced by 100):
    phi_s = phi_c1.pipe(
        xr_rolling, window=win_r, window2=win_azi, method="median",
        skipna=True, min_periods=max(3, int((win_r - 1) / 4))
    )

    # 6: offset Part 1 of 6: calculate offset 2d (reduced by 100):
    phi_off_2d = phase_offset(phi_c0.where(swp_cf.range > 1000),
                              rng).PHIDP_OFFSET.load()

    # 6: offset Part 2 of 6: first filter (>5 in neighbourhood; red. by 100):
    phi_off_2d_f = xr.where(np.isnan(phi_off_2d), -100, phi_off_2d)
    phi_off_2d_f = xr.DataArray(uniform_filter(
        phi_off_2d_f, size=[0, 3], mode='mirror'), dims=['time', 'azimuth'])
    phi_off_2d_f = phi_off_2d.where((abs(phi_off_2d - phi_off_2d_f)) < 5)

    # 6: offset Part 3 of 6: median offsets per time (red. by 100):
    phi_off_1d = phi_off_2d_f.median("azimuth", skipna=True)
    phi_off_0d = phi_off_1d.median("time", skipna=True)
    phi_off_1d = xr.where(abs(phi_off_1d - phi_off_0d) > 10, np.nan,
                          phi_off_1d)
    phi_off_0d = phi_off_1d.median("time", skipna=True)
    phi_off_1d = xr.where(abs(phi_off_1d - phi_off_0d) > 10,
                          phi_off_0d, phi_off_1d)
    phi_off_1d = xr.where(np.isnan(phi_off_1d), phi_off_0d, phi_off_1d)

    # 6: offset Part 4 of 6: snd. filter (>10 to 1d median; red. by 100):
    phi_off_2d_f_f = phi_off_2d_f.where((abs(phi_off_2d_f - phi_off_1d)) < 10)
    phi_off_2d_f_f = xr.where(np.isnan(phi_off_2d_f_f), phi_off_1d,
                              phi_off_2d_f_f)

    # 6: offset Part 5 of 6:  smoothing
    phi_off_2d_f_f_s = xr.DataArray(
        gaussian_filter(phi_off_2d_f_f, sigma=[0, 5], mode='wrap'),
        dims=['time', 'azimuth'])

    # 6: offset Part 6 of 6:  noise corrected phi (NOT red. by 100 anymore!)
    phi_nc = phi_s - phi_off_2d_f_f_s

    # 7: KDP
    kdp_light = phi_nc.wrl.dp.kdp_from_phidp(winlen=wkdp_light)
    kdp_heavy = phi_nc.wrl.dp.kdp_from_phidp(winlen=wkdp_heavy)
    kdp_comb = kdp_heavy.where(swp_cf.DBZH < 40, kdp_light)

    # 8: assign: variables and attributes
    phi_nc.attrs["long_name"] = 'Differential phase shift'
    phi_nc.attrs["short_name"] = 'PHI_DP'
    phi_nc.attrs["units"] = 'degrees'
    phi_nc.attrs["comments"] = 'PHI_DP smoothing with win_r=' + \
                               str(win_r) + ' and win_azi=' + str(win_azi)
    swp_cf = swp_cf.assign(PHI_NC=phi_nc)

    kdp_comb.attrs["comments"] = 'KDP noise corrected with winlen=' + \
                                 str(wkdp_heavy) + ' (DBZH<40) and ' + \
                                 'winlen=' + str(wkdp_light) + ' (DBZH>=40)'
    swp_cf = swp_cf.assign(KDP_NC=kdp_comb)

    phi_c1 = phi_c1 + 100
    phi_c1.attrs[
        "long_name"] = 'raw differential phase shift centered at median'
    phi_c1.attrs["short_name"] = 'PHI_C'
    phi_c1.attrs["units"] = 'degrees'
    swp_cf = swp_cf.assign(PHI_C=phi_c1)  # + 100)

    phi_off_2d = phi_off_2d + 100
    phi_off_2d.attrs["long_name"] = 'instrument phase offset raw'
    phi_off_2d.attrs["short_name"] = 'PHI_DP offset raw'
    phi_off_2d.attrs["units"] = 'degrees'
    swp_cf = swp_cf.assign(PHI_OFFSET_raw=phi_off_2d)  # + 100)

    phi_off_2d_f_f_s = phi_off_2d_f_f_s + 100
    phi_off_2d_f_f_s.attrs[
        "long_name"] = 'instrument phase offset centered at 0'
    phi_off_2d_f_f_s.attrs["short_name"] = 'PHI_DP offset centered'
    phi_off_2d_f_f_s.attrs["units"] = 'degrees'
    swp_cf = swp_cf.assign(PHI_OFFSET_centered=phi_off_2d_f_f_s)  # + 100)

    phi_off_2d_f_f_s = phi_off_2d_f_f_s + phi_modus
    phi_off_2d_f_f_s.attrs["long_name"] = 'instrument phase offset'
    phi_off_2d_f_f_s.attrs["short_name"] = 'PHI_DP offset'
    phi_off_2d_f_f_s.attrs["units"] = 'degrees'
    swp_cf = swp_cf.assign(PHI_OFFSET=phi_off_2d_f_f_s)  # + 100 + phi_modus)

    phi_flip.attrs["long_name"] = 'PHI_DP flipped'
    phi_flip.attrs["short_name"] = 'PHI_DP flipped'
    phi_flip.attrs["units"] = '0,1'
    swp_cf = swp_cf.assign(PHI_FLIPPED=phi_flip)

    phi_flip_i.attrs["long_name"] = 'PHI_DP flipping index'
    phi_flip_i.attrs["short_name"] = 'PHI_DP flipping index'
    phi_flip_i.attrs["units"] = '°'
    swp_cf = swp_cf.assign(PHI_FLIPPED_INDEX=phi_flip_i)

    return swp_cf


# J. Steinheuer
def correct_phi_kdp(date, location, elevation_deg=5.5, mode='vol',
                    overwrite=False, dir_data_obs=header.dir_data_obs,
                    parts=6, merge=True, remove_parts=True,
                    uh_tresh=0, rho_tresh=0.8, snr_tresh=15,
                    win_r=25, win_azi=None, rng=3000,
                    wkdp_light=9, wkdp_heavy=25):
    """
    Smooth phi and kdp and calculate phi offset.

    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    elevation_deg : elevation in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    mode : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
    dir_data_obs : directory to search for input cases
                  (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
    dir_data_era5 : directory to search for era5 files.
    [...]
    """
    parts_current = parts
    time_a = time.time()
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp':
        if sweep != '00':
            return

        print('\nstart: ' + date + ' ' + location + ' pcp')
    else:
        print('\nstart: ' + date + ' ' + location + ' ' +
              str(elevation_deg))

    merge_files = []
    for p in range(parts):
        time_f = time.time()
        time_l = time.time()
        i_t_a = int(288 / parts * p)
        i_t_b = int(288 / parts * (p + 1))
        path_in = "/".join([dir_data_obs + '*',
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep,
                            'ras*_allmoms_*'])
        files = sorted(glob.glob(path_in))
        if not files:
            path_in = "/".join([dir_data_obs +
                                year, year + '-' + mon,
                                year + '-' + mon + '-' + day,
                                location, mode + '*', sweep,
                                'ras*_allmoms_*'])
            files = sorted(glob.glob(path_in))
        if not files:
            print('no input data *_allmoms_*')
            return
        else:
            path_in = files[0]
            path_out = path_in.replace(
                '_allmoms_', '_kdp_nc_' + str(p) + '_')

        if (os.path.isfile(path_out) and not overwrite) or \
                (os.path.isfile(path_out.replace(
                    '_kdp_nc_' + str(p), '_kdp_nc'))
                 and not overwrite):
            print(path_out + ' exists;\n' + ' ... set: > ' +
                  'overwrite = True < for recalculation')
            merge_files.append(path_out)
            continue

        path_rho_nc = path_in.replace('_allmoms_', '_rhohv_nc_')
        data = dttree.open_datatree(path_in)[
            'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
        data_rho = dttree.open_datatree(path_rho_nc)[
            'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
        data.RHOHV.values = data_rho.RHOHV_NC2P.values
        data = data.assign({'SNRH': data_rho.SNRH})
        remo_var = list(data.data_vars.keys())
        # remo_var.remove('CMAP')
        # remo_var.remove('UPHIDP')
        data = data.transpose('time', 'azimuth', 'range')
        data = data.isel(time=slice(i_t_a, i_t_b))
        if data.time.size == 0:
            parts_current = p
            print('no more time steps')
            break

        merge_files.append(path_out)
        if location == 'umd':
            flip_default = -1
        else:
            flip_default = 1

        time_l = time.time() - time_l
        time_p = time.time()
        data = proc_phidp_kdp(data,
                              uh_tresh=uh_tresh,
                              rho_tresh=rho_tresh,
                              snr_tresh=snr_tresh,
                              win_r=win_r,
                              win_azi=win_azi,
                              wkdp_light=wkdp_light,
                              wkdp_heavy=wkdp_heavy,
                              rng=rng,
                              flip_default=flip_default)
        time_p = time.time() - time_p
        time_s = time.time()
        mom_use = [x for x in list(data.keys())]
        for mom in mom_use:
            data[mom].encoding["coordinates"] = \
                "time azimuth range"
        data = data.drop_vars(remo_var)
        dtree = dttree.DataTree(name="root")
        dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                        parent=dtree)
        print('saving: ... ' + path_out.split('/')[-1] + ' ...')
        dtree.load().to_netcdf(path_out)
        data.close()
        print('saved (' + str(p + 1) + '/' + str(parts_current) +
              '): ' + path_out + ' !')
        time_s = time.time() - time_s
        time_f = time.time() - time_f
        print('full time: ' + str(round(time_f, 1)) +
              's: 1/3: loading time: ' + str(round(time_l, 1)) +
              's; 2/3: processing time: ' + str(round(time_p, 1)) +
              's; 3/3: saving time: ' + str(round(time_s, 1)) +
              's')

    if merge and merge_files != []:
        path_out_new = merge_files[0].replace(
            'kdp_nc_0_', 'kdp_nc_')
        if os.path.isfile(path_out_new) and not overwrite:
            print(path_out_new + ' exists;\n' + ' ... set: ' +
                  '> overwrite = True < for recalculation')
        else:
            data_merged = xr.merge([
                dttree.open_datatree(merge_files[p])[
                    'sweep_' + str(int(sweep))].to_dataset(
                ).chunk(-1) for p in range(parts_current)])
            mom_use = [x for x in list(data_merged.keys())]
            for mom in mom_use:
                data_merged[mom].encoding["coordinates"] = \
                    "time azimuth range"

            data_merged.attrs['processing_date'] = str(
                pd.Timestamp.today())[:16]
            print('saving: ... ' + path_out_new.split('/')[-1] +
                  ' ...')
            dtree = dttree.DataTree(name="root")
            dttree.DataTree(data_merged,
                            name=f"sweep_{int(sweep)}",
                            parent=dtree)
            dtree.load().to_netcdf(path_out_new)
            data_merged.close()
            print('combined:  ' + path_out_new + ' !')
            time_a = time.time() - time_a
            print(' -> full case: ' +
                  str(round(time_a, 1)) + 's\n')
            if remove_parts:
                for file_part in merge_files:
                    os.remove(file_part)


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
]
LOCATIONS = [
    'asb', 'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,
    4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,
])
MODE = [
    'pcp',
    'vol',
]
overwrite = False
# --------------------------------------------------------------------------- #
uh_tresh = 0
rho_tresh = 0.8
snr_tresh = 15
win_r = 25  # VP or PARK?!
win_azi = None
rng = 3000
wkdp_light = 9  # PARK: light filtering
wkdp_heavy = 25  # PARK: heavy filtering
parts = 6
merge = True
remove_parts = True
print('Departing into: ' + str(parts))
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                correct_phi_kdp(date=date, location=location,
                                elevation_deg=elevation_deg, mode=mode,
                                overwrite=overwrite,
                                dir_data_obs=header.dir_data_obs,
                                parts=parts, merge=merge,
                                remove_parts=remove_parts, uh_tresh=uh_tresh,
                                rho_tresh=rho_tresh, snr_tresh=snr_tresh,
                                win_r=win_r, win_azi=win_azi, rng=rng,
                                wkdp_light=wkdp_light, wkdp_heavy=wkdp_heavy)

# --------------------------------------------------------------------------- #
# OLD CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20170719",
    "20170725",
]
LOCATIONS = [
    'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,
    4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0,
])
MODE = [
    'pcp',
    'vol',
]
overwrite = False
# --------------------------------------------------------------------------- #
uh_tresh = 0
rho_tresh = 0.8
snr_tresh = 15
win_r = 25  # VP or PARK?!
win_azi = None
rng = 3000
wkdp_light = 9  # PARK: light filtering
wkdp_heavy = 25  # PARK: heavy filtering
parts = 6
merge = True
remove_parts = True
print('Departing into: ' + str(parts))
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                correct_phi_kdp(date=date, location=location,
                                elevation_deg=elevation_deg, mode=mode,
                                overwrite=overwrite,
                                dir_data_obs=header.dir_data_obs,
                                parts=parts, merge=merge,
                                remove_parts=remove_parts, uh_tresh=uh_tresh,
                                rho_tresh=rho_tresh, snr_tresh=snr_tresh,
                                win_r=win_r, win_azi=win_azi, rng=rng,
                                wkdp_light=wkdp_light, wkdp_heavy=wkdp_heavy)

# --------------------------------------------------------------------------- #
# CONTINUE?
# import DWD_obs_to_MIUB_obs_5_calibrate_zdr
