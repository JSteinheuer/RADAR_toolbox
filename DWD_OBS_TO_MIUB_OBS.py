#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 04.07.24                                                 #
# DWD_OBS_TO_MIUB_OBS.py                                                      #
#                                                                             #
# Functions for processing raw DWD C-Band data towards:                       #
#   STEP 1: oad all moments                                                   #
#   STEP 2: correct rho_hv                                                    #
#   STEP 3: add ERA5 temperature                                              #
#   STEP 4: correct phi, smooth kdp                                           #
#   STEP 5: calibrate zdr                                                     #
#   STEP 6: correct for attenuation                                           #
#   STEP 7: combine everything                                                #
# --------------------------------------------------------------------------- #

import cdsapi
import datatree as dttree
import glob
import matplotlib.pyplot as plt
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

import HEADER_RADAR_toolbox as header


# --------------------------------------------------------------------------- #
# STEP 1: Load all moments and all times (of day) in one file.                #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/build_radar_database/concat_dwd_data_to_d* #
#                                                                             #
# @author: jgiles                                                             #
# This script takes all dwd radar files from a folder (for one elevation) and #
# merges them into a single file combining all moments along all timesteps.   #
# Then saves the resulting dataset into a new file with the same naming       #
# style but with "allmoms" instead of the moment name. Additionally, it saves #
# either a true.txt or false.txt file alongside, if the data fulfills certain #
# condition, as an attempt to check if there is actually something            #
# interesting in that period of data.                                         #
# --------------------------------------------------------------------------- #


# J. Steinheuer adapted J. Gilles
def align(ds):
    """
    Reduce and align dataset coordinates

    Parameter
    ---------
    ds : xarray.DataArray or xarray.Dataset
    """
    # reduce time in the azimuth
    ds["time"] = ds["time"].load().min()
    # remove elevation in time
    ds["elevation"] = ds["elevation"].load().median()
    ds["azimuth"] = ds["azimuth"].load().round(1)
    # in case there are duplicate rays, remove them
    ds = xd.util.remove_duplicate_rays(ds)
    # in case there are duplicate times
    ds["time"] = np.unique(ds["time"])
    return ds.set_coords(["sweep_mode", "sweep_number", "prt_mode",
                          "follow_mode", "sweep_fixed_angle"])


# J. Steinheuer adapted J. Gilles
def align_pcp(ds):
    """
    Reduce and align dataset coordinates

    Parameter
    ---------
    ds : xarray.DataArray or xarray.Dataset
    """
    # reduce time in the azimuth
    ds["time"] = ds["time"].load().min()
    # remove elevation in time
    # ds["elevation"] = ds["elevation"].load().median('time')
    ds["azimuth"] = ds["azimuth"].load().round(1)
    # in case there are duplicate rays, remove them
    ds = xd.util.remove_duplicate_rays(ds)
    # in case there are duplicate times
    ds["time"] = np.unique(ds["time"])
    return ds.set_coords(["sweep_mode", "sweep_number", "prt_mode",
                          "follow_mode", "sweep_fixed_angle"])


# J. Steinheuer adapted J. Gilles
def reduce_duplicate_timesteps(ds):
    """
    Reduce duplicate time values by combining the data variables available at
    each timestep. In case there are duplicate data vars (same name in the
    duplicate time value) then remove the duplicates, otherwise keep all
    variables by renaming them with _n, with n = 2,... starting from the
    second one with the same name.

    Originally made for handling turkish radar datasets.

    Parameter
    ---------
    ds : xarray.Dataset
    """
    # Check if there are duplicate time values
    if (ds.time.diff("time").astype(int) == 0).any():
        # Create a new dummy dataset to concatenate the reduced data
        ds_reduced = ds.isel(time=[0]).copy()
        for tt in sorted(set(ds.time.values)):
            # Make a list for the selected variables
            reduced_vars = []
            # Select data for the current timestep
            ds_time = ds.sel(time=tt)
            # If time is not a dimension anymore
            if "time" not in ds_time.dims:
                # then just add the current ds
                ds_reduced = xr.concat([ds_reduced, ds_time], dim="time")
            else:
                # for each variable,
                for vv in ds_time.data_vars:
                    # start by removing NA
                    ds_time_nona = ds_time[vv].dropna("time", how="all")
                    if len(ds_time_nona["time"]) == 1:
                        # if the array got reduced to only 1 timestep,
                        # then attach it to reduced_vars
                        reduced_vars.append(ds_time_nona.copy())
                    elif len(ds_time_nona["time"]) == 0:
                        # if the array got reduced to 0 timestep, then skip
                        continue
                    else:
                        # if there are still more than 1 timestep, divide
                        # them into separate variables
                        for itt in range(len(ds_time.time)):
                            count = 0
                            if itt == 0:
                                # set the first one as the "pristine" variable
                                reduced_vars.append(
                                    ds_time_nona.isel(time=[itt]).copy())
                                count += 1
                            else:
                                for ivv in range(count):
                                    ivv += ivv  # we need to sum 1 because
                                    # we are going backwards
                                    if ds_time_nona.isel(time=[itt]).equals(
                                            reduced_vars[-ivv]):
                                        # if it is equal to any of the already
                                        # selected arrays, ignore it
                                        break
                                    else:
                                        # if it is not equal to previous
                                        # arrays, append it as a variant of
                                        # that variable
                                        reduced_vars.append(ds_time_nona.isel(
                                            time=[itt]).rename(
                                            vv + "_" + str(count)).copy())
                                        count += 1

                # Merge all variables into a single dataset
                reduced_vars_ds = xr.merge(reduced_vars)
                # Add the current timestep data to the new dataset
                ds_reduced = xr.concat([ds_reduced, reduced_vars_ds],
                                       dim="time")

        # delete the dummy data from the first timestep
        # Here I could also use .drop_isel but there is a bug and it does
        # not work https://github.com/pydata/xarray/issues/6605
        ds_reduced = ds_reduced.drop_duplicates("time", keep="last")
        return ds_reduced
    else:  # in case there are no duplicate timesteps, just return ds
        return ds


# J. Steinheuer adapted J. Gilles
def fix_time_in_coords(ds):
    """
    Fix time coord issues and reduce time dimension if present in dataset
    coordinates where it is not needed

    Parameter
    ---------
    ds : xarray.DataArray or xarray.Dataset
    """
    # It may happen that some time value is missing or that time values
    # are repeated, attempt to fix that using info in rtime
    if ds["time"].isnull().any() or (
            ds["time"].diff("time").compute().astype(int) <= 0).any():
        ds.coords["time"] = ds.rtime.min(dim="azimuth", skipna=True).compute()

    for coord in ["latitude", "longitude", "altitude", "elevation"]:
        # if some coord has dimension time, reduce using median
        if "time" in ds[coord].dims:
            ds.coords[coord] = ds.coords[coord].median("time")

    # in case there are still duplicate timesteps,
    # attempt to reduce the time dim
    ds = reduce_duplicate_timesteps(ds)
    return ds


# J. Steinheuer adapted J. Gilles
def load_dwd_raw(filepath, moments, pcp=False):
    """
    Load DWD raw data.

    Parameter
    ---------
    filepath : str, list
               Location of the file or path with wildcards to find files
               using glob or list of paths
    moments : str, list
              a list of requested pol. moments. For all set >*< or >any<
    """
    # collect files
    if type(filepath) is list:
        files = sorted(filepath)
    else:
        files = sorted(glob.glob(filepath))

    # extract list of moments
    moments_available = set(fp.split("_")[-2] for fp in files)
    if moments != '*':
        moments.append('ANY')
        moments_remove = [x for x in moments_available if x.upper()
                          not in moments]
        for mom_rem in moments_remove:
            moments_available.discard(mom_rem)

    # discard "allmoms" from the set if it exists
    moments_available.discard("allmoms")
    try:
        # for every moment, open all files in folder (all timesteps)
        # per moment into a dataset
        vardict = {}  # a dict for putting a dataset per moment
        for mom in moments_available:
            # open the odim files (single moment and elevation,
            # several timesteps)
            llmom = sorted([ff for ff in files if "_" + mom + "_" in ff])
            if pcp:
                vardict[mom] = xr.open_mfdataset(llmom, engine="odim",
                                                 combine="nested",
                                                 concat_dim="time",
                                                 preprocess=align_pcp)
            else:
                vardict[mom] = xr.open_mfdataset(llmom, engine="odim",
                                                 combine="nested",
                                                 concat_dim="time",
                                                 preprocess=align)

            vardict[mom] = fix_time_in_coords(vardict[mom])
    except OSError:
        pathparts = [xx if len(xx) == 10 and "20" in xx else None for xx in
                     llmom[0].split("/")]
        pathparts.sort(key=lambda e: (e is None, e))
        date = pathparts[0]
        print(date + " " + mom + ": Error opening files. Some file is "
                                 "corrupt or truncated.")
        sys.exit("Script terminated early. " + date + " " + mom +
                 ": Error opening files. Some file is corrupt or truncated.")

    # merge all moments
    return xr.merge(vardict.values())


# J. Steinheuer
def load_all_moms(date, location, elevation_deg=5.5, mode='vol',
                  moments=['CMAP', 'DBSNRH', 'DBZH',
                           'RHOHV', 'UPHIDP', 'ZDR', 'SNRHC'],
                  overwrite=False,
                  dir_data_obs=header.dir_data_obs,
                  dir_data_obs_realpep=header.dir_data_obs_realpep):
    """
    Load DWD data of one day and put the moments in one new *all_moms*-file.

    Parameter
    ---------
    date : 'yyyymmdd' datestring.
    location : 'rrr' 3-letter string for radar location.
    elevation_deg : elevation in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    mode : set 'vol' for volume and 'pcp' for precipitation.
    moments : a list of requested pol. moments. For all set >*< or >any<
    overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
    dir_data_obs : directory to search for input cases
                  (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
    dir_data_obs_realpep : a second directory to search for input.
                         (>dir_data_obs</yyyy/yyyy-mm/yyyy-mm-dd/).
    """
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        return

    if mode == '90grad' and sweep != '00':
        return

    path_in = "/".join([dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep, 'ras*'])
    files = sorted(glob.glob(path_in))
    files_temp = []
    for file in files:
        if 'allmoms' not in file:
            if 'rhohv_nc' not in file:
                if 'ERA5' not in file:
                    if 'kdp_nc' not in file:
                        if 'zdr_off' not in file:
                            if 'zh_zdr_ac' not in file:
                                if 'polmoms' not in file:
                                    if 'vradh' not in file:
                                        files_temp.append(file)

    files = files_temp
    if not files:
        print_a = 'No input: ' + path_in
        path_in = "/".join([dir_data_obs_realpep + '' +
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep, 'ras*'])
        files = sorted(glob.glob(path_in))
        if not files:
            print(print_a)
            print(' ... nor: ' + path_in + ' -> continue')
            return

        path_out = '/'.join((files[0].split('/'))[:-1])
        if date == '20210714':
            path_out = path_out.replace(dir_data_obs_realpep, dir_data_obs +
                                        'OpHymet2-case09-20210714/')
        else:
            path_out = path_out.replace(dir_data_obs_realpep, dir_data_obs +
                                        'OpHymet2-caseX-' + date + '/')

    else:
        path_out = '/'.join((files[0].split('/'))[:-1])

    files_temp = []
    for file in files:
        if 'allmoms' not in file:
            if 'rhohv_nc' not in file:
                if 'ERA5' not in file:
                    if 'kdp_nc' not in file:
                        if 'zdr_off' not in file:
                            if 'zh_zdr_ac' not in file:
                                if 'polmoms' not in file:
                                    if 'vradh' not in file:
                                        files_temp.append(file)

    files = files_temp
    name = files[0].split("/")[-1].split("_")
    t_start = files[0].split("/")[-1].split("-")[2][:12]
    t_end = files[-1].split("/")[-1].split("-")[2][:12]
    name[-2] = "allmoms"
    name_m1 = name[-1]
    name_m1 = name_m1.replace('-hd5', '.hd5')
    name_m1 = name_m1.replace('-h5', '.hd5')
    name_m1 = name_m1.split("-")
    name_m1[1] = t_start + '-' + t_end
    name[-1] = '-'.join(name_m1)
    name_out = ("_".join(name))
    file_out = '/'.join([path_out, name_out])
    Path(path_out).mkdir(parents=True, exist_ok=True)
    if not overwrite and os.path.exists(file_out):
        print('exists: ' + file_out + ' -> continue')
        return

    ds = load_dwd_raw(files, moments, pcp=(mode == 'pcp'))
    if moments != '*' or moments != 'any':
        mom_remove = [x for x in list(ds.keys()) if x not in moments]
        ds = ds.drop_vars(mom_remove)

    coords_needed = ['altitude', 'latitude', 'longitude',
                     'azimuth', 'elevation', 'range', 'time']
    coors_remove = [x for x in list(ds.coords) if x not in coords_needed]
    ds = ds.drop_vars(coors_remove)
    mom_use = [x for x in list(ds.keys())]
    for mom in mom_use:
        ds[mom].encoding["coordinates"] = "elevation azimuth range time"

    dtree = dttree.DataTree(name="root")
    dttree.DataTree(ds, name=f"sweep_{int(sweep)}", parent=dtree)
    print('saving: ... ' + file_out + ' ...')
    dtree.load().to_netcdf(file_out)
    ds.close()
    print('saved:  ' + file_out + ' !')
    return


# --------------------------------------------------------------------------- #
# STEP 2: correct for rho_hv.                                                 #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/correct_rhohv.py                           #
# --------------------------------------------------------------------------- #


# J. Steinheuer adapted J. Gilles
def noise_correction2(dbz, rho, noise_level):
    """
    Calculate SNR, apply to RHOHV
    Formula from Ryzhkov book page 187/203
    """
    # noise calculations
    snrh = dbz - 20 * np.log10(dbz.range * 0.001) - noise_level - \
           0.033 * dbz.range / 1000
    snrh = snrh.where(snrh >= 0).fillna(0)
    attrs = xd.model.sweep_vars_mapping['SNRH']
    attrs.pop('gamic', None)
    snrh = snrh.assign_attrs(attrs)
    snrh.name = "SNRH"
    rho_nc = rho * (1. + 1. / 10. ** (snrh * 0.1))
    rho_nc = rho_nc.assign_attrs(rho.attrs)
    rho_nc.name = "RHOHV_NC"
    return snrh, rho_nc


# J. Steinheuer adapted J. Gilles
def calculate_noise_level(dbz, rho, noise=(-40, -20, 1),
                          rho_bins=(0.9, 1.1, 0.005), snr_bins=(5., 30., .1)):
    """
    This functions calculates the noise levels and noise corrections for
    RHOHV, for a range of noise values. It returns a list of signal-to-noise
    and corrected rhohv arrays, as well as histograms involved in the
    calculations, a list of standard deviations for every result and the noise
    value with the minimum std. The final noise correction should be chosen
     based on the rn value (minumum std).

    The default noise range is based on BoXPol data, it may be good to extend
    it a bit for C-Band.

    The final noise level (rn) should be used with noise_correction2 one more
    time to get the final result. It may happen that the correction is too
    strong and we get some RHOHV values over 1. We should check this for some
    days of data and if that is the case, then select a noise level that is
    slightly less (about 2% less).
    """
    noise = np.arange(*noise)
    rho_bins = np.arange(*rho_bins)
    snr_bins = np.arange(*snr_bins)
    corr = [noise_correction2(dbz, rho, n) for n in noise]
    hist = [histogram(rho0, snr0, bins=[rho_bins, snr_bins],
                      block_size=rho.time.size) for snr0, rho0 in corr]
    std = [np.std(r.idxmax('RHOHV_NC_bin')).values for r in hist]
    rn = noise[np.argmin(std)]
    return corr, hist, std, rn


# J. Steinheuer
def correct_rho_hv(date, location, elevation_deg=5.5, mode='vol',
                   overwrite=False, dir_data_obs=header.dir_data_obs):
    """
    Load *all_moms*-file of one day and correct for rho_hv.

    Parameter
    ---------
    date : 'yyyymmdd' datestring.
    location : 'rrr' 3-letter string for radar location.
    elevation_deg : elevation in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    mode : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool;, if *allmoms*-output exists, it can be overwritten.
    dir_data_obs : directory to search for input cases
                  (>dir_data_obs</*/yyyy/yyyy-mm/yyyy-mm-dd).
    """
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        return

    path_in = "/".join([dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_allmoms_*'])
    files = sorted(glob.glob(path_in))
    if not files:
        print('No input: ' + path_in + ' -> continue')
        return
    else:
        path_in = files[0]
        path_out = path_in.replace('_allmoms_', '_rhohv_nc_')

    if not overwrite and os.path.exists(path_out):
        print('exists: ' + path_out + ' -> continue')
        return

    dbzh_names = ["DBZH"]  # names to look for the DBZH variable
    rhohv_names = ["RHOHV"]  # same but for RHOHV
    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    # get RHOHV name
    for X_RHO in rhohv_names:
        if X_RHO in data.data_vars:
            break

    # get DBZH name
    for X_DBZH in dbzh_names:
        if X_DBZH in data.data_vars:
            break

    # check that the variables actually exist, otherwise continue
    if X_DBZH not in data.data_vars:
        print("DBZH not found in data")
        data.close()
        return
    if X_RHO not in data.data_vars:
        print("RHOHV not found in data")
        sys.exit("RHOHV not found in data.")
        data.close()
        return

    rho_nc = calculate_noise_level(
        data[X_DBZH], data[X_RHO], noise=(-45, -15, 1))
    # print('lin.fit:' + path_out + ' ...')
    fits = []
    for nn, rhon in enumerate(rho_nc[0]):
        merged = xr.merge(rhon)
        rhonc_snrh = xr.DataArray(
            merged.RHOHV_NC.values.flatten(),
            coords={"SNRH": merged.SNRH.values.flatten()})
        try:
            fits.append(float(rhonc_snrh.where(
                (0 < rhonc_snrh.SNRH) &
                (rhonc_snrh.SNRH < 20) &
                (rhonc_snrh > 0.7)
            ).polyfit("SNRH", deg=1, skipna=True
                      ).polyfit_coefficients[0].values))
        except:
            # if it does not work, just attach nan
            fits.append(np.nan)

    # checking which fit has the slope closest to zero
    try:
        bci = np.nanargmin(np.abs(np.array(fits)))
    except ValueError:
        # if all slopes are nan, no correction good enough, abort
        print("Could not calculate noise correction (possibly "
              "due to not enough data points). Aborting...")
        return

    # get the best noise correction level according to bci
    ncl = np.arange(-45, -15, 1)[bci]
    # merge into a single array
    rho_nc_out = xr.merge(rho_nc[0][bci])
    # add noise correction level as attribute
    rho_nc_out.attrs["noise correction level"] = ncl
    rho_nc_out['SNRH'].encoding["coordinates"] = \
        "elevation azimuth range time"
    rho_nc_out["RHOHV_NC"].encoding = data[X_RHO].encoding
    dtree = dttree.DataTree(name="root")
    # Just in case, calculate again for a NCL slightly lower (2%),
    # in case the automatically-selected one is too strong
    rho_nc2 = noise_correction2(
        data[X_DBZH], data[X_RHO], ncl * 1.02)
    rho_nc_out2 = xr.merge(rho_nc2)
    rho_nc_out2.attrs["noise correction level"] = ncl * 1.02
    rho_nc_out2['SNRH'].encoding["coordinates"] = \
        "elevation azimuth range time"
    rho_nc_out2["RHOHV_NC"].encoding = data[X_RHO].encoding
    rho_nc_out["RHOHV_NC2P"] = rho_nc_out2["RHOHV_NC"]
    rho_nc_out["RHOHV_NC2P"].attrs["comments"] = \
        "for a NCL slightly lower (2%), i.e. NCL_2P=NCL*1.02"
    dttree.DataTree(rho_nc_out, name=f"sweep_{int(sweep)}", parent=dtree)
    print('saving: ... ' + path_out + ' ...')
    dtree.load().to_netcdf(path_out)
    data.close()
    rho_nc_out.close()
    rho_nc_out2.close()
    print('saved:  ' + path_out + ' !')
    return


# --------------------------------------------------------------------------- #
# STEP 3: write hourly temperature values to each radar and its grid. ERA5 3D #
#         temperature fields from ERA5_download_data_DE.py                    #
# --------------------------------------------------------------------------- #


def download_ERA5_temp(date, overwrite=False, dir_out=header.dir_data_era5):
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]

    c = cdsapi.Client()
    file_out1 = dir_out + str(year) + str(mon) + str(day) + "-3D-T-q-ml.grib"
    if not overwrite and os.path.exists(file_out1):
        print('exists: ' + file_out1 + ' -> continue')
    elif not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out1.replace('grib', 'nc' + ' -> continue'))
    else:
        c.retrieve("reanalysis-era5-complete", {
            "class": "ea",
            "date": "%s-%s-%s" % (year, mon, day),
            "expver": "1",
            "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/"
                        "21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/"
                        "38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/"
                        "55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/"
                        "72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/"
                        "89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/"
                        "104/105/106/107/108/109/110/111/112/113/114/115/"
                        "116/117/118/119/120/121/122/123/124/125/126/127/"
                        "128/129/130/131/132/133/134/135/136/137",
            "levtype": "ml",
            "param": "130/133",  # Temperature (t) and specific humidity (q)
            "stream": "oper",
            "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/"
                    "06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/"
                    "12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/"
                    "18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type": "an",
            "grid": [0.25, 0.25],
            'area': [56, 2, 46.25, 17],
        }, file_out1)
    file_out2 = dir_out + str(year) + str(mon) + str(day) + \
                "-2D-z-lnsp-ml.grib"
    if not overwrite and os.path.exists(file_out2):
        print('exists: ' + file_out2 + ' -> continue')
    elif not overwrite and os.path.exists(file_out2.replace('grib', 'nc')):
        print('exists: ' + file_out2.replace('grib', 'nc' + ' -> continue'))
    else:
        c.retrieve("reanalysis-era5-complete", {
            "class": "ea",
            "date": "%s-%s-%s" % (year, mon, day),
            "expver": "1",
            "levelist": "1",
            "levtype": "ml",
            "param": "129/152",  # Geopotential (z) and ln of s. pres. (lnsp)
            "stream": "oper",
            "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/"
                    "06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/"
                    "12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/"
                    "18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type": "an",
            "grid": [0.25, 0.25],
            'area': [56, 2, 46.25, 17],
        }, file_out2)

    file_out3 = file_out1.replace('T-q-ml', 'z')
    if not overwrite and os.path.exists(file_out3):
        print('exists: ' + file_out3 + ' -> continue')
    elif not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out3.replace('grib', 'nc' + ' -> continue'))
    else:
        os.system('python /user/s6justei/PyCharm/PyCharmProjects/' +
                  'RADAR_toolbox/ERA5_compute_geopotential_on_ml.py ' +
                  file_out1 +
                  ' ' + file_out2 + ' -o ' + file_out3)

    if not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out1.replace('grib', 'nc') + ' -> continue')
    else:
        os.system('cdo -f nc copy ' + file_out1 + ' ' +
                  file_out1.replace('grib', 'nc'))
        # os.system('rm ' + file_out1)

    # if not overwrite and os.path.exists(file_out2.replace('grib', 'nc')):
    #     print(
    #         'exists: ' + file_out2.replace('grib', 'nc') + ' -> continue')
    # else:
    #     os.system('cdo -f nc copy ' + file_out2 + ' ' +
    #               file_out2.replace('grib', 'nc'))
    #     # os.system('rm ' + file_out2)

    if not overwrite and os.path.exists(file_out3.replace('grib', 'nc')):
        print(
            'exists: ' + file_out3.replace('grib', 'nc') + ' -> continue')
    else:
        os.system('cdo -f nc copy ' + file_out3 + ' ' +
                  file_out3.replace('grib', 'nc'))
        # os.system('rm ' + file_out3)

    # os.system('cdo -f nc copy ' + file_out + ' ' +
    #           file_out.replace('grib', 'nc'))
    # os.system('rm ' + file_out)

    # file_out4 = dir_out + str(year) + str(mon) + str(day) + "-2D-SP-T2m.nc"
    # if not overwrite and os.path.exists(file_out4):
    #     print('exists: ' + file_out4 + ' -> continue')
    # else:
    #     c.retrieve('reanalysis-era5-single-levels', {
    #         'product_type': 'reanalysis',
    #         'variable': ['2m_temperature', 'surface_pressure',
    #                      'geopotential'],  # TODO: new
    #         # 'variable': ['2m_temperature', 'surface_pressure', ],  # old
    #         'year': year,
    #         'month': mon,
    #         'day': day,
    #         'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
    #                  '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
    #                  '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
    #                  '18:00', '19:00', '20:00', '21:00', '22:00', '23:00', ],
    #         'area': [56, 2, 46.25, 17],
    #         'format': 'netcdf',
    #     }, file_out4)


# V. Pejcic (K. Mühlbauer) / from J.Steinheuer SYN_RADAR_1_CREATE_VOLUME_SCAN
def get_lon_lat_alt(r, az, el, sitecoords):
    """
    Get lon/lat/alt from range, azimuth, elevation, and site-coordinates.

    A wradlib.georef wrapper.

     Args:
        r: array-like radar range.
        az: array-like radar azimuth.
        el: radar elevation.
        sitecoords: list of lon,lat, alt of radar site.

    Returns:
        lat: latitude array of dims(r.shape,az.shape)
        lon: longitude array of dims(r.shape,az.shape)
        alt: altidude array of dims(r.shape,az.shape)
    """

    proj_wgs84 = wrl.georef.epsg_to_osr(4326)
    cent_coords = wrl.georef.spherical_to_centroids(r, az, el, sitecoords,
                                                    crs=proj_wgs84)
    cent_coords = np.squeeze(cent_coords)
    lon = cent_coords[..., 0]
    lat = cent_coords[..., 1]
    alt = cent_coords[..., 2]
    return lon, lat, alt


# V. Pejcic (K. Mühlbauer) / from J.Steinheuer SYN_RADAR_1_CREATE_VOLUME_SCAN
def ipol_fc_to_radgrid(mod_lon, mod_lat, mod_z, rad_lon, rad_lat, rad_alt,
                       method='Nearest'):
    """
    Interpolate forecast grid to the radar grid.

    Search for the models lon/lat/alt indices, that are closest to the
    given radar lon/lat/alt. The return ist a >mask< array, that shrinks the
    model domain and an instance of wradlib.inpol.Nearest class. Later, for
    a given variable shaped in model grid, i.e. >n_heights x n_model_cells<,
    one needs to apply 1) the masking and 2) call the instance func_ipol to
    get the nearest variables on the radar grid, i.e. shaped in
    [flattened n_lo x n_la x n_al].

     Args:
        mod_lon: array of model longitudes [n_heights x n_model_cells].
        mod_lat: array of model latitudes [n_heights x n_model_cells].
        mod_z: array of model level altitudes [n_heights x n_model_cells].
        rad_lon: radar longitude array [flattened n_lo x n_la x n_al].
        rad_lat: radar latitude array [flattened n_lo x n_la x n_al].
        rad_alt: altitude array of [flattened n_lo x n_la x n_al].

    Returns:
        func_ipol: wradlib.inpol.Nearest class that is initialized and
        mask: necessary for model_fields to reduce their dim
            ([n_heights x n_model_cells] towards [n_mask]) in order to
            call func_ipol.
    """

    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)
    mod_x, mod_y = wrl.georef.reproject(mod_lon,
                                        mod_lat,
                                        projection_target=proj_stereo,
                                        projection_source=proj_wgs)
    rad_x, rad_y = wrl.georef.reproject(rad_lon,
                                        rad_lat,
                                        projection_target=proj_stereo,
                                        projection_source=proj_wgs)

    # # only those model data that are in radar domain
    # mask = (mod_x >= rad_x.min()) & (mod_x <= rad_x.max()) & (
    #         mod_y >= rad_y.min()) & (mod_y <= rad_y.max()) & (
    #         mod_z >= rad_alt.min()) & (mod_z <= rad_alt.max())

    # only those model data that are in radar domain + bordering volume
    outer_x = max(0.3 * (rad_x.max() - rad_x.min()), 1)
    outer_y = max(0.3 * (rad_y.max() - rad_y.min()), 1)
    lower_z = 50
    upper_z = 2000
    mask = (mod_x >= rad_x.min() - outer_x) & (
            mod_x <= rad_x.max() + outer_x) & (
                   mod_y >= rad_y.min() - outer_y) & (
                   mod_y <= rad_y.max() + outer_y) & (
                   mod_z >= rad_alt.min() - lower_z) & (
                   mod_z <= rad_alt.max() + upper_z)
    mod_x = mod_x[mask]
    mod_y = mod_y[mask]
    mod_alt = mod_z[mask]

    # source coordinates model and target coordinates radar
    src = np.vstack((mod_x.ravel(),
                     mod_y.ravel(),
                     mod_alt.ravel() / 1e3)).T
    trg = np.vstack((rad_x.ravel(),
                     rad_y.ravel(),
                     rad_alt.ravel() / 1e3)).T
    if method == 'Linear':
        func_ipol = wrl.ipol.Linear(src, trg)
    else:
        func_ipol = wrl.ipol.Nearest(src, trg)

    return func_ipol, mask


# J. Steinheuer
def era5_temp(date, location, elevation_deg=5.5, mode='vol',
              overwrite=False, dir_data_obs=header.dir_data_obs,
              dir_data_era5=header.dir_data_era5):
    """
    Load one *all_moms*-file (of one day and sweep) and create on the
    radar grid all temperature values taken from ERA5.

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
    """
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    if mode == 'pcp' and sweep != '00':
        return

    path_in = "/".join([dir_data_obs + '*',
                        year, year + '-' + mon,
                        year + '-' + mon + '-' + day,
                        location, mode + '*', sweep,
                        'ras*_allmoms_*'])
    files = sorted(glob.glob(path_in))
    if not files:
        print('No input: ' + path_in + ' -> continue')
        return
    else:
        path_in = files[0]
        path_out = path_in.replace('_allmoms_', '_ERA5_temp_')
        if len(files) > 1:
            print('More than 1 input -> take files[0]: ' + path_in)

    if not overwrite and os.path.exists(path_out):
        print('exists: ' + path_out + ' -> continue')
        return

    path_in_any = '/'.join(path_in.split('/')[:-1]) + '/*any*'
    files_any = sorted(glob.glob(path_in_any))
    if not files_any:
        path_in_any = "/".join([header.dir_data_obs_realpep + '' +
                                year, year + '-' + mon,
                                year + '-' + mon + '-' + day,
                                location, mode + '*', sweep, 'ras*'])
        files_any = sorted(glob.glob(path_in_any))
        if not files_any:
            print('nothing found -> bw=1')
            bw = 1
        else:
            path_in_any = files_any[0]
            try:
                bw = dttree.open_datatree(path_in_any)['how'].attrs['beamwidth']
            except:
                print('nothing in *any* found -> bw=1')
                bw = 1

    else:
        path_in_any = files_any[0]
        try:
            bw = dttree.open_datatree(path_in_any)['how'].attrs['beamwidth']
        except:
            print('nothing in *any* found -> bw=1')
            bw = 1

    bw = round(bw, 3)
    print(date + ' ' + location + ' bw=' + str(bw))
    path_in_t = dir_data_era5 + date + '-3D-T-q-ml.nc'
    path_in_z = dir_data_era5 + date + '-3D-z.nc'
    if not os.path.exists(path_in_t):
        print('missing ERA5: ' + path_in_t + ' -> continue')
        return

    if not os.path.exists(path_in_z):
        print('missing ERA5: ' + path_in_z + ' -> continue')
        return

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data = data.drop_vars(list(data.data_vars.keys()))
    data_t_era5 = dttree.open_datatree(path_in_t).to_dataset()
    data_z_era5 = dttree.open_datatree(path_in_z).to_dataset()
    r = data.range
    az = data.azimuth
    el = data.elevation
    sitecoords = [data.longitude.values,
                  data.latitude.values,
                  data.altitude.values]
    lon, lat, alt = get_lon_lat_alt(r, az, el, sitecoords)
    # variables
    data['lat'] = (['range', 'azimuth'],
                   lat.transpose(),
                   dict(standard_name='latitude',
                        units='degrees_north'))
    data['lon'] = (['range', 'azimuth'],
                   lon.transpose(),
                   dict(standard_name='longitude',
                        units='degrees_east'))
    if mode == 'pcp':
        data['alt'] = (['range', 'azimuth'],
                       alt.transpose(),
                       dict(standard_name='altitude',
                            comments='height above mean sea level',
                            units='m'))
    else:
        data['alt'] = (['range'],
                       alt[0, :],
                       dict(standard_name='altitude',
                            comments='height above mean sea level',
                            units='m'))

    data['time'] = data_t_era5['time']
    dummy_tra = np.empty(shape=[data.time.size, data.range.size,
                                data.azimuth.size, ])
    dummy_tra[:] = np.nan
    data['temp'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                    dict(standard_name='air temperature', units='K'))
    data['temp_beamtop'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                            dict(standard_name='air temperature at beamtop',
                                 units='K',
                                 comment='beambroadening induced temperature of '
                                         'bin volume top (approximated for pw=0)'))
    data['temp_beambottom'] = (['time', 'range', 'azimuth'], dummy_tra.copy(),
                               dict(
                                   standard_name='air temperature at beambottom',
                                   units='K', comment='beambroadening induced '
                                                      'temperature of bin volume '
                                                      'bottom (approximated for pw=0)'))
    shape_ra = (data.range.size, data.azimuth.size)
    for t_i in range(data['time'].size):
        # grid for searching later the closest ERA5 cells
        rad_lon = data['lon'].data.flatten()
        rad_lat = data['lat'].data.flatten()
        if mode == 'pcp':
            rad_alt = data['alt'].data.flatten()
        else:
            rad_alt = np.repeat(data['alt'].data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        era5_lon, era5_lat = np.meshgrid(data_z_era5.lon, data_z_era5.lat)
        era5_alt = data_z_era5.z.isel(time=t_i) / 9.81
        era5_lon = era5_lon.flatten()
        era5_lat = era5_lat.flatten()
        era5_z = era5_alt.data.reshape(era5_alt.shape[0],
                                       era5_alt.shape[1] * era5_alt.shape[2])
        # temp center
        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
            # lon with shape (137, 2013):
            np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
            # alt with shape (137, 2013):
            era5_z,
            rad_lon, rad_lat, rad_alt,  # (259200,)
            method='Linear',
        )  # mask.shape (137, 2013)
        data['temp'][t_i, :, :] = func_ipol(
            data_t_era5['t'].data[t_i, :, :, :].reshape(
                era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
            [mask]).reshape(shape_ra)

        # temp top
        beamradius = wrl.util.half_power_radius(data.range, bw)
        if mode == 'pcp':
            rad_alt = (data['alt'] + beamradius).data.flatten()
        else:
            rad_alt = np.repeat((data['alt'] + beamradius).data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
            # lon with shape (137, 2013):
            np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
            # alt with shape (137, 2013):
            era5_z,
            rad_lon, rad_lat, rad_alt,  # (259200,)
            method='Linear',
        )  # mask.shape (137, 2013)
        data['temp_beamtop'][t_i, :, :] = func_ipol(
            data_t_era5['t'].data[t_i, :, :, :].reshape(
                era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
            [mask]).reshape(shape_ra)

        # temp bottom
        if mode == 'pcp':
            rad_alt = (data['alt'] - beamradius).data.flatten()
        else:
            rad_alt = np.repeat((data['alt'] - beamradius).data[:, np.newaxis],
                                data['azimuth'].shape, axis=1).flatten()

        func_ipol, mask = ipol_fc_to_radgrid(
            # lon with shape (137, 2013):
            np.repeat(era5_lon[np.newaxis, :], era5_z.shape[0], axis=0),
            # lon with shape (137, 2013):
            np.repeat(era5_lat[np.newaxis, :], era5_z.shape[0], axis=0),
            # alt with shape (137, 2013):
            era5_z,
            rad_lon, rad_lat, rad_alt,  # (259200,)
            method='Linear',
        )  # mask.shape (137, 2013)
        data['temp_beambottom'][t_i, :, :] = func_ipol(
            data_t_era5['t'].data[t_i, :, :, :].reshape(
                era5_alt.shape[0], era5_alt.shape[1] * era5_alt.shape[2])
            [mask]).reshape(shape_ra)

    mom_use = [x for x in list(data.keys())]
    for mom in mom_use:
        data[mom].encoding["coordinates"] = "time range azimuth"

    dtree = dttree.DataTree(name="root")
    dttree.DataTree(data, name=f"sweep_{int(sweep)}", parent=dtree)
    print('saving: ... ' + path_out + ' ...')
    dtree.load().to_netcdf(path_out)
    data.close()
    print('saved:  ' + path_out + ' !')


# --------------------------------------------------------------------------- #
# STEP 4: correct KDP in ML.                                                  #
# --------------------------------------------------------------------------- #


# Function from radarmet by K. Mühlbauer
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
# observations towards MIUB 'standard'.                                       #
# STEP 5: calibrate ZDR.                                                      #
#         Adapted from Velibor Pejcic                                         #
# --------------------------------------------------------------------------- #


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
    swp_cf = swp_cf.chunk(chunks=-1)
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
    swp_cf = swp_cf.chunk(chunks=-1)
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
    swp_cf = swp_cf.chunk(chunks=-1)
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
    swp_cf = swp_cf.chunk(chunks=-1)
    swp_mask = swp_cf.where((swp_cf.DBZH > 0) &
                            (swp_cf.DBZH < 40) &
                            (swp_cf.RHOHV > 0.98) &
                            np.isnan(swp_cf.CMAP))
    zdroffset = np.nanmedian(swp_mask.ZDR)
    zdroffset_sd = np.nanstd(swp_mask.ZDR)
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
    if overwrite or not os.path.exists(path_out_nc):
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
            bb_off, bb_nm, bb_sd = cal_zdr_birdbath(data, plot=False,
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
        if elevation_deg == elevation_degs[0] and mode == modes[0]:
            file_in = nc_file_mom.split('/')[-1]
            file_out = file_in.replace('.hd5', '_' + str(n_rows) +
                                       'x' + str(n_cols) + '.' +
                                       pdf_or_png).replace('_' + sweep,
                                                           '_all')
            path_out_plot = folder_plot + 'ZDR_calibration/' + \
                            location.upper() + '/ZDR_calibration_' + \
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
            bb_off, bb_nm, bb_sd = cal_zdr_birdbath(data, plot=True, ax=ax)
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
                                              plot=[True, False],
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
                                               plot=[True, False],
                                               axes=axes,
                                               colorbar=colorbar)
            axes = [plt.subplot(n_rows, n_cols, index + 2 * n_cols),
                    plt.subplot(n_rows, n_cols, index + 2 * n_cols)]
            lr_off_ec, lr_nm_ec = cal_zdr_lightrain(data2, band='C',
                                                    plot=[True, False],
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
                                                     plot=[True, False],
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
# STEP 6: attenuation correction.                                             #
#         see Ryzhkov and Zrnic (2019): 6.4 (pp. 162, 167)                    #
# --------------------------------------------------------------------------- #


# J. Steinheuer
def correct_zh_zdr(swp_cf, uh_tresh=0,
                   alpha=0.08, beta=0.02,
                   # until_temp_beamtop=273.15+4,# TODO if diff in ML is need
                   ML_until_temp_beambottom=273.15 + 0,
                   # alpha_ML_multipl=1,  # TODO if alpha in ML is to change
                   # beta_ML_multipl=1,  # TODO if beta in ML is to change
                   ):
    """
    ...
    """
    # 1: thresholding
    swp_mask = swp_cf.where((swp_cf.DBZH > uh_tresh) &
                            np.isnan(swp_cf.CMAP)
                            )
    swp_mask['PHI_NC'] = xr.where(swp_mask.PHI_NC < 0, 0, swp_mask.PHI_NC)
    swp_last_l = swp_mask.where(
        (swp_mask.temp_beambottom >= ML_until_temp_beambottom) &
        (swp_mask.temp_beambottom < ML_until_temp_beambottom + 1))
    phi_const = swp_last_l.PHI_NC.median(dim="range", skipna=True)
    phi_const = xr.where(np.isnan(phi_const), swp_mask.where(
        swp_mask.temp_beambottom >= ML_until_temp_beambottom).PHI_NC.max(
        dim="range", skipna=True), phi_const)
    phi_const = xr.where(np.isnan(phi_const), swp_mask.where(
        swp_mask.temp_beambottom < ML_until_temp_beambottom).PHI_NC.min(
        dim="range", skipna=True), phi_const)

    phi_4ac = xr.where(swp_mask.temp_beambottom >= ML_until_temp_beambottom,
                       swp_mask.PHI_NC, phi_const)

    zh_ac = swp_mask.DBZH + phi_4ac * alpha
    zdr_ac = swp_mask.ZDR + phi_4ac * beta
    zh_ac.attrs["long_name"] = 'reflectivity factor attenuation corrected'
    zh_ac.attrs["short_name"] = 'ZH ac'
    zh_ac.attrs["units"] = 'dBZ'
    swp_cf = swp_cf.assign(ZH_AC=zh_ac)
    zdr_ac.attrs["long_name"] = 'Log differential reflectivity ' + \
                                'attenuation corrected'
    zdr_ac.attrs["short_name"] = 'ZDR ac'
    zdr_ac.attrs["units"] = 'dB'
    swp_cf = swp_cf.assign(ZDR_AC=zdr_ac)
    swp_cf = swp_cf.assign(PHI_4AC=phi_4ac)
    return swp_cf


# J. Steinheuer
def attenuation_correction(date, location, elevation_deg=5.5, mode='vol',
                           overwrite=False, dir_data_obs=header.dir_data_obs):
    """
    .

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
    """
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

    path_out = path_in.replace('_allmoms_', '_zh_zdr_ac_')
    if os.path.isfile(path_out) and not overwrite:
        print(path_out + ' exists;\n' + ' ... set: > ' +
              'overwrite = True < for recalculation')
        return

    path_kdp = path_in.replace('_allmoms_', '_kdp_nc_')
    if not os.path.exists(path_kdp):
        print('not exists: ' + path_kdp + ' -> continue')
        return

    path_rho = path_in.replace('_allmoms_', '_rhohv_nc_')
    if not os.path.exists(path_rho):
        print('not exists: ' + path_rho + ' -> continue')
        return

    path_temp = path_in.replace('_allmoms_', '_ERA5_temp_')
    if not os.path.exists(path_temp):
        print('not exists: ' + path_temp + ' -> continue')
        return

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_rho = dttree.open_datatree(path_rho)[
        'sweep_' + str(int(sweep))].to_dataset()
    data_kdp = dttree.open_datatree(path_kdp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_temp = dttree.open_datatree(path_temp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk('auto')
    data_temp2 = data_temp.interp(
        coords=data.drop(['longitude', 'latitude',
                          'altitude', 'elevation']).coords,
        method='nearest')

    data.RHOHV.values = data_rho.RHOHV_NC2P.values
    data = data.assign({'SNRH': data_rho.SNRH})
    data = data.assign({'PHI_NC': data_kdp.PHI_NC})
    data = data.assign({'temp_beambottom': data_temp2.temp_beambottom})
    data = data.assign({'temp_beamtop': data_temp.temp_beamtop})
    remo_var = list(data.data_vars.keys())
    data = data.transpose('time', 'azimuth', 'range')
    data = correct_zh_zdr(data)
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
    data_kdp.close()
    data_temp.close()


# --------------------------------------------------------------------------- #
# STEP 7: combine all noise corrected polarimetric moments                    #
# --------------------------------------------------------------------------- #


# J. Steinheuer
def combine_pol_mom_nc(date, location, elevation_deg=5.5, mode='vol',
                       overwrite=False, dir_data_obs=header.dir_data_obs,
                       method_zdr_priorities=['BB', 'SD_V', 'LR_V'
                                                            'SD_I', 'LR_I', ],
                       other_zdr_off_day='',
                       n_zdr_lowest=1000, std_zdr_highest=2):
    """
    .

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
    """
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

    if other_zdr_off_day:
        year_of = other_zdr_off_day[0:4]
        mon_of = other_zdr_off_day[4:6]
        day_of = other_zdr_off_day[6:8]
    else:
        year_of = date[0:4]
        mon_of = date[4:6]
        day_of = date[6:8]

    path_in_bb = "/".join([dir_data_obs + '*',
                           year_of, year_of + '-' + mon_of,
                           year_of + '-' + mon_of + '-' + day_of,
                           location, '90grad*', '00',
                           'ras*_zdr_off_*'])
    files = sorted(glob.glob(path_in_bb))
    if not files:
        path_in_bb = "/".join([dir_data_obs +
                               year_of, year_of + '-' + mon_of,
                               year_of + '-' + mon_of + '-' + day_of,
                               location, '90grad*', '00',
                               'ras*_zdr_off_*'])
        files = sorted(glob.glob(path_in_bb))
    if not files:
        print('no birdbath input data *90grad*_zdr_off_*')
        if 'BB' in method_zdr_priorities:
            method_zdr_priorities.remove('BB')

    else:
        path_in_bb = files[0]

    paths_zdr_off = []
    path_zdr_off = "/".join([dir_data_obs + '*',
                             year_of, year_of + '-' + mon_of,
                             year_of + '-' + mon_of + '-' + day_of,
                             location, 'pcp*', '00',
                             'ras*_zdr_off_*'])
    files = sorted(glob.glob(path_zdr_off))
    if not files:
        path_zdr_off = "/".join([dir_data_obs +
                                 year_of, year_of + '-' + mon_of,
                                 year_of + '-' + mon_of + '-' + day_of,
                                 location, 'pcp*', '00',
                                 'ras*_zdr_off_*'])
        files = sorted(glob.glob(path_zdr_off))
    if not files:
        print('no pcp input data *_zdr_off_*')
        # return
    else:
        paths_zdr_off.append(files[0])
        if sweep == '00' and mode == 'pcp':
            path_zdr_off = files[0]

    for swe in ['00', '01', '02', '03', '04',
                '05', '06', '07', '08', '09']:
        path_zdr_off_i = "/".join([dir_data_obs + '*',
                                   year_of, year_of + '-' + mon_of,
                                   year_of + '-' + mon_of + '-' + day_of,
                                   location, 'vol*', swe,
                                   'ras*_zdr_off_*'])
        files = sorted(glob.glob(path_zdr_off_i))
        if not files:
            path_zdr_off_i = "/".join([dir_data_obs +
                                       year_of, year_of + '-' + mon_of,
                                       year_of + '-' + mon_of + '-' + day_of,
                                       location, 'vol*', swe,
                                       'ras*_zdr_off_*'])
            files = sorted(glob.glob(path_zdr_off_i))
        if not files:
            print('no ' + swe + ' input data *_zdr_off_*')
            # return
        else:
            paths_zdr_off.append(files[0])
            if swe == sweep and mode == 'vol':
                path_zdr_off = files[0]

    if not paths_zdr_off:
        if 'SD_V' in method_zdr_priorities:
            method_zdr_priorities.remove('SD_V')

        if 'LR_V' in method_zdr_priorities:
            method_zdr_priorities.remove('LR_V')

    path_out = path_in.replace('_allmoms_', '_polmoms_nc_')
    if os.path.isfile(path_out) and not overwrite:
        print(path_out + ' exists;\n' + ' ... set: > ' +
              'overwrite = True < for recalculation')
        return

    path_kdp = path_in.replace('_allmoms_', '_kdp_nc_')
    if not os.path.exists(path_kdp):
        print('not exists: ' + path_kdp + ' -> continue')
        return

    path_rho = path_in.replace('_allmoms_', '_rhohv_nc_')
    if not os.path.exists(path_rho):
        print('not exists: ' + path_rho + ' -> continue')
        return

    path_temp = path_in.replace('_allmoms_', '_ERA5_temp_')
    if not os.path.exists(path_temp):
        print('not exists: ' + path_temp + ' -> continue')
        return

    if not os.path.exists(path_zdr_off):
        print('not exists: ' + path_zdr_off + ' -> continue')
        if 'SD_I' in method_zdr_priorities:
            method_zdr_priorities.remove('SD_I')

        if 'LR_I' in method_zdr_priorities:
            method_zdr_priorities.remove('LR_I')

    path_zhzdr_ac = path_in.replace('_allmoms_', '_zh_zdr_ac_')
    if not os.path.exists(path_zhzdr_ac):
        print('not exists: ' + path_zhzdr_ac + ' -> continue')
        return

    data = dttree.open_datatree(path_in)[
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    data_rho = dttree.open_datatree(path_rho)[
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    data_kdp = dttree.open_datatree(path_kdp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    data_ac = dttree.open_datatree(path_zhzdr_ac)[
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    data_temp = dttree.open_datatree(path_temp)[
        'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
    data_temp2 = data_temp.interp(
        coords=data.drop(['longitude', 'latitude',
                          'altitude', 'elevation']).coords,
        method='nearest')

    # combine:
    remo_var = list(data_kdp.data_vars.keys())
    remo_var.remove('PHI_NC')
    remo_var.remove('KDP_NC')
    data_combined = data_kdp.drop_vars(remo_var)
    if 'VRADH' in data.variables.keys():
        data_combined = data_combined.assign({'VRADH': data.VRADH})

    data_combined = data_combined.assign({'RHOHV_NC2P': data_rho.RHOHV_NC2P})
    data_combined = data_combined.assign({'ZH_AC': data_ac.ZH_AC.reindex_like(
        other=data_rho.RHOHV_NC2P, method='nearest')})
    data_combined.ZH_AC.attrs[
        "long_name"] = 'reflectivity factor attenuation corrected'
    data_combined.ZH_AC.attrs["short_name"] = 'ZH ac'
    data_combined.ZH_AC.attrs["units"] = 'dBZ'
    data_combined = data_combined.assign({'ZDR_AC_OC':
                                              data_ac.ZDR_AC.reindex_like(
                                                  other=data_rho.RHOHV_NC2P,
                                                  method='nearest')})
    xr.set_options(keep_attrs=True)
    data_combined = data_combined.where(np.isnan(data.CMAP))
    data_combined = data_combined.where(data.DBZH >= 0)
    data_combined = data_combined.assign(
        {'PHI_offset_median': np.nanmedian(data_kdp.PHI_OFFSET)})
    data_combined.PHI_offset_median.attrs["long_name"] = \
        'median instrument phase offset'
    data_combined.PHI_offset_median.attrs["short_name"] = \
        'median PHI_DP offset'
    data_combined.PHI_offset_median.attrs["units"] = 'degrees'
    for method_zdr in method_zdr_priorities + ['00']:
        if method_zdr == 'BB':
            data_bb = dttree.open_datatree(path_in_bb)[
                'sweep_0'].to_dataset().chunk(-1)
            n = data_bb.zdr_off_bb_n.data
            sd = data_bb.zdr_off_sd_bb.data
            of = data_bb.zdr_off_bb.data
            print('BB:' + other_zdr_off_day + ' n=' + str(n) +
                  ' offset=' + str(of) + ' std=' + str(sd))
            if n >= n_zdr_lowest and sd <= std_zdr_highest:
                data_combined.ZDR_AC_OC.data = \
                    data_combined.ZDR_AC_OC.data - of
                data_combined = data_combined.assign(
                    {'ZDR_offset': data_bb.zdr_off_bb})
                data_bb.close()
                break
            else:
                data_bb.close()

        elif method_zdr == 'SD_I':
            data_zdr_off = dttree.open_datatree(path_zdr_off)[
                'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
            n = data_zdr_off.zdr_off_sd_n.data
            of = data_zdr_off.zdr_off_sd.data
            print('SD_I:' + other_zdr_off_day + ' n=' + str(n) +
                  ' offset=' + str(of))
            if n >= n_zdr_lowest:
                data_combined.ZDR_AC_OC.data = \
                    data_combined.ZDR_AC_OC.data - of
                data_combined = data_combined.assign(
                    {'ZDR_offset': data_zdr_off.zdr_off_sd})
                data_zdr_off.close()
                break
            else:
                data_zdr_off.close()

        elif method_zdr == 'LR_I':
            data_zdr_off = dttree.open_datatree(path_zdr_off)[
                'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
            n = data_zdr_off.zdr_off_lr_n.data
            of = data_zdr_off.zdr_off_lr.data
            print('LR_I:' + other_zdr_off_day + ' n=' + str(n) +
                  ' offset=' + str(of))
            if n >= n_zdr_lowest:
                data_combined.ZDR_AC_OC.data = \
                    data_combined.ZDR_AC_OC.data - of
                data_combined = data_combined.assign(
                    {'ZDR_offset': data_zdr_off.zdr_off_lr})
                data_zdr_off.close()
                break
            else:
                data_zdr_off.close()

        elif method_zdr == 'SD_V':
            n_vol = []
            of_vol = []
            for path_zdr_off_i, swe in zip(paths_zdr_off,
                                           [0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]):
                data_zdr_off = dttree.open_datatree(path_zdr_off_i)[
                    'sweep_' + str(int(swe))].to_dataset().chunk(-1)
                n_vol.append(data_zdr_off.zdr_off_sd_n.data.item())
                of_vol.append(data_zdr_off.zdr_off_sd.data.item())
                if path_zdr_off_i == path_zdr_off:
                    print('SD_I:' + other_zdr_off_day + ' n=' +
                          str(n_vol[-1]) + ' offset=' + str(of_vol[-1]))

                data_zdr_off.close()

            n = sum(n_vol)
            of = sum(np.array(n_vol) * np.array(of_vol)) / n
            print('SD_all: n=' + str(np.round(np.array(n_vol) / 1000)) +
                  'k \n   offset=' + str(np.round(np.array(of_vol), 2)))
            print('SD_V:' + other_zdr_off_day + ' n=' + str(n) +
                  ' offset=' + str(of))
            if n >= n_zdr_lowest:
                data_zdr_off = dttree.open_datatree(path_zdr_off_i)[
                    'sweep_' + str(int(swe))].to_dataset().chunk(-1)
                data_combined.ZDR_AC_OC.data = \
                    data_combined.ZDR_AC_OC.data - of
                data_combined = data_combined.assign(
                    {'ZDR_offset': data_zdr_off.zdr_off_sd})
                data_combined.ZDR_offset.attrs['comment'] = \
                    data_combined.ZDR_offset.comment + \
                    ' (combined offset from all sweeps)'
                data_combined.ZDR_offset.data = of
                data_zdr_off.close()
                break

        elif method_zdr == 'LR_V':
            n_vol = []
            of_vol = []
            for path_zdr_off_i, swe in zip(paths_zdr_off,
                                           [0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]):
                data_zdr_off = dttree.open_datatree(path_zdr_off_i)[
                    'sweep_' + str(int(swe))].to_dataset().chunk(-1)
                n_vol.append(data_zdr_off.zdr_off_lr_n.data.item())
                of_vol.append(data_zdr_off.zdr_off_lr.data.item())
                if path_zdr_off_i == path_zdr_off:
                    print('LR_I:' + other_zdr_off_day + ' n=' +
                          str(n_vol[-1]) + ' offset=' + str(of_vol[-1]))

                data_zdr_off.close()

            n = sum(n_vol)
            of = sum(np.array(n_vol) * np.array(of_vol)) / n
            print('LR_all: n=' + str(np.round(np.array(n_vol) / 1000)) +
                  'k \n   offset=' + str(np.round(np.array(of_vol), 2)))
            print('LR_V:' + other_zdr_off_day + ' n=' + str(n) +
                  ' offset=' + str(of))
            if n >= n_zdr_lowest:
                data_zdr_off = dttree.open_datatree(path_zdr_off_i)[
                    'sweep_' + str(int(swe))].to_dataset().chunk(-1)
                data_combined.ZDR_AC_OC.data = \
                    data_combined.ZDR_AC_OC.data - of
                data_combined = data_combined.assign(
                    {'ZDR_offset': data_zdr_off.zdr_off_lr})
                data_combined.ZDR_offset.attrs['comment'] = \
                    data_combined.ZDR_offset.comment + \
                    ' (combined offset from all sweeps)'
                data_combined.ZDR_offset.data = of
                data_zdr_off.close()
                break

        elif method_zdr == '00':
            print('no offset could be estimated')
            data_combined = data_combined.assign({'ZDR_offset': 0})
            data_combined['ZDR_offset'].attrs["short_name"] = 'ZDR off'
            data_combined['ZDR_offset'].attrs["units"] = 'dB'
            data_combined['ZDR_offset'].attrs["comment"] = \
                'no proper offset found from methods!'
            break

        else:
            print(method_zdr + ' is no method that is implemented')

    if other_zdr_off_day:
        data_combined['ZDR_offset'].attrs["comment"] = \
            data_combined.ZDR_offset.comment + \
            ' (from day = ' + other_zdr_off_day + ')'

    data_combined = data_combined.assign({'temp': data_temp2.temp})
    data_combined = data_combined.assign({'temp_beambottom':
                                              data_temp2.temp_beambottom})
    data_combined.temp_beambottom.attrs["long_name"] = \
        'air temperature at beambottom'
    data_combined = data_combined.assign({'temp_beamtop':
                                              data_temp2.temp_beamtop})
    data_combined.temp_beambottom.attrs["long_name"] = \
        'air temperature at beamtop'
    data_combined = data_combined.assign({'lat': data_temp.lat})
    data_combined = data_combined.assign({'lon': data_temp.lon})
    data_combined = data_combined.assign({'alt': data_temp.alt})
    data_combined.elevation.attrs["standard_name"] = 'ray_elevation_angle'
    data_combined.elevation.attrs["long_name"] = 'elevation_angle_from_surface'
    data_combined.elevation.attrs["units"] = 'degrees'
    data_combined = data_combined.drop_vars('quantile')
    data_combined = data_combined.transpose('time', 'azimuth', 'range')
    data_combined = data_combined.chunk(chunks=-1)
    data_combined.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
    mom_use = [x for x in list(data_combined.keys())]
    for mom in mom_use:
        data_combined[mom].encoding["coordinates"] = \
            "time azimuth range lat lon"

    data_combined['alt'].encoding["coordinates"] = "azimuth range lat lon"
    data_combined['lat'].encoding["coordinates"] = "azimuth range"
    data_combined['lon'].encoding["coordinates"] = "azimuth range"
    data_combined['time'].encoding["units"] = "seconds since " + \
                                              year + "-" + mon + "-" + day
    data_combined['time'].attrs["comment"] = "UTC"
    dtree = dttree.DataTree(name="root")
    dttree.DataTree(data_combined, name=f"sweep_{int(sweep)}",
                    parent=dtree)
    dtree.load().to_netcdf(path_out)
    print('Done: ... ' + path_out.split('/')[-1] + ' ...')
    data.close()
    data_ac.close()
    data_combined.close()
    data_kdp.close()
    data_rho.close()
    data_temp.close()
    data_temp2.close()
    return

# --------------------------------------------------------------------------- #
