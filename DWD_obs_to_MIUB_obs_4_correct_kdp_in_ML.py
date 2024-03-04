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
from scipy.ndimage.filters import uniform_filter1d

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

    # Thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) &
                            (swp_cf.RHOHV > rho_tresh) &
                            (swp_cf.SNRH > snr_tresh) &
                            np.isnan(swp_cf.CMAP))

    # Median filtering 2d:
    # handel close 'overshooting' values, i.e. -179° and 179° so that in each
    # sweep most phi lays around 0 (phi_centered
    phi_hist = histogram(swp_mask.UPHIDP.copy(),
                         bins=np.linspace(-180, 180, 360),
                         dim=['azimuth', 'range', 'time'])
    # get modus of phi (minus 180 gives then degree instead of index)
    phi_modus = np.argmax(phi_hist.values) - 180
    # phi in [-180-modus, 180-modus]  ->  [0,360] -> [180, 540] ->
    # [180,360,0,180] -> [-180,new modus=0, 180]
    phi_c = ((swp_mask.UPHIDP.copy() - xr.DataArray(phi_modus))#new, dims='time'))
             % 360 + 180) % 360 - 180

    # Median 2D
    # TODO: many min_periods requested now!
    phi_smoothed = phi_c.copy().pipe(xr_rolling, window=win_r, window2=win_azi,
                                     method="median", skipna=True,
                                     min_periods=int((win_r - 1) / 2))

    # offset per azimuth and time
    phi_off_2d = phase_offset(phi_c.copy().where((swp_cf.RHOHV >= 0.9)), rng)

    # 1d per time first guess
    phi_off_fg = phi_off_2d.PHIDP_OFFSET.load().median("azimuth", skipna=True)

    # 2d (a and t) centered at 0
    phi_off_2d_c = phi_off_2d.PHIDP_OFFSET - phi_off_fg

    # throw out outliers
    phi_off_2d_c = phi_off_2d_c.where(
        (phi_off_2d_c > -20) & (phi_off_2d_c < 20))

    # mean or median?!
    phi_offset_mn = phi_off_2d_c.mean("azimuth", skipna=True)  # this=!
    phi_offset_md = phi_off_2d_c.median("azimuth", skipna=True)
    # phi_corr_mn = phi_smoothed - phi_offset_mn  # this=!
    # phi_corr_md = phi_smoothed - phi_offset_md
    # phi_corr = phi_corr_mn

    # new smoth:
    phi_offset_md_s = xr.DataArray(uniform_filter1d(phi_offset_md,
                                                    size=13, mode='mirror'),
                                   dims='time')
    phi_corr_md_s = phi_smoothed - phi_offset_md_s
    phi_corr = phi_corr_md_s

    # phidp/kdp
    kdp_light = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_light)
    kdp_heavy = phi_corr.wrl.dp.kdp_from_phidp(winlen=wkdp_heavy)
    kdp_comb = kdp_heavy.where(swp_cf.DBZH < 40, kdp_light)

    # attributes
    kdp_comb.attrs["comments"] = 'KDP noise corrected with winlen=' + \
                                 str(wkdp_heavy) + ' (DBZH<40) and ' + \
                                 'winlen=' + str(wkdp_light) + ' (DBZH>=40)'
    phi_corr.attrs["long_name"] = 'Differential phase shift'
    phi_corr.attrs["short_name"] = 'PHI_DP'
    phi_corr.attrs["units"] = 'degrees'
    phi_corr.attrs["comments"] = 'PHI_DP smoothing with win_r=' + \
                                 str(win_r) + ' and win_azi=' + str(win_azi)
    swp_cf = swp_cf.assign(KDP_NC=kdp_comb)
    swp_cf = swp_cf.assign(PHI_NC=phi_corr)
    swp_cf = swp_cf.assign(PHI_off_md=phi_offset_md)
    swp_cf = swp_cf.assign(PHI_off_mn=phi_offset_mn)
    swp_cf = swp_cf.assign(PHI_off_mn_s=phi_offset_mn_s)  # new smoth
    swp_cf = swp_cf.assign(PHI_off_md_s=phi_offset_md_s)  # new smoth
    swp_cf = swp_cf.assign(start_range=phi_off_2d.start_range)
    swp_cf = swp_cf.assign(stop_range=phi_off_2d.stop_range)
    swp_cf = swp_cf.assign(PHI_off_2d=phi_off_2d.PHIDP_OFFSET)
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
LOCATIONS = ['ess', 'pro', 'umd', 'tur', 'asb', 'boo', 'drs', 'eis', 'fbg',
             'fld', 'hnr', 'isn',
             'mem', 'neu', 'nhb', 'oft', 'ros',
             ]
ELEVATIONS = ELEVATIONS_ALL.copy()
MODE = ['pcp', 'vol']

merge = True
remove_parts = True
# remove_parts = False
# overwrite = False
overwrite = True

snr_tresh = 15
win_r = 25
win_azi = None
wkdp_light = 9
wkdp_heavy = 25
rng = 3000

# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# # LOCATIONS = ['pro']
# LOCATIONS = ['ess']
# ELEVATIONS = np.array([5.5])
# MODE = ['vol']

# date = '20210604'
# # date = '20210714'
# location = 'pro'
# # location = 'ess'
# elevation_deg = 5.5
# mode = 'vol'

for location in LOCATIONS:
    for date in DATES:
        for mode in MODE:
            for elevation_deg in ELEVATIONS:
                parts = 4
                merge_files = []
                for p in range(parts):
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
                        path_out = path_in.replace('_allmoms_', '_kdp_nc_' +
                                                   str(p) + '_')

                    if (os.path.isfile(path_out) and not overwrite) or \
                            (os.path.isfile(path_out.replace(
                                'kdp_nc_' + str(p), 'kdp_nc'))
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
                    data = proc_phidp_kdp(data, uh_tresh=0, rho_tresh=0.9,
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
                        'kdp_nc_0_', 'kdp_nc_')
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
