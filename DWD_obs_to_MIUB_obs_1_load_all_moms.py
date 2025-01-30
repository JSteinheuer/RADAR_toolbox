#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 10.01.24                                                 #
# DWD_obs_to_MIUB_obs_1_load_all_moms.py                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 1: Load all moments and all times (of day) in one file.                #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/build_radar_database/concat_dwd_data_to_d* #
# --------------------------------------------------------------------------- #
"""
@author: jgiles
This script takes all dwd radar files from a folder (for one elevation) and
merges them into a single file combining all moments along all timesteps.
Then saves the resulting dataset into a new file with the same naming
style but with "allmoms" instead of the moment name. Additionally, it saves
either a true.txt or false.txt file alongside, if the data fulfills certain
condition, as an attempt to check if there is actually something interesting
in that period of data.
"""

import datatree as dttree
import numpy as np
import sys
import glob
import HEADER_RADAR_toolbox as header
from pathlib import Path
import os
import xarray as xr
import xradar as xd
import time as time_p
import datetime as dt

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).

import HEADER_RADAR_toolbox as header


# --------------------------------------------------------------------------- #
# FUNCTIONS: J. Steinheuer modified functions from J. Gilles


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
                                    files_temp.append(file)
                                    # if 'vradh' not in file:
                                    #     files_temp.append(file) # TODO: why?

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
                                    files_temp.append(file)
                                    # if 'vradh' not in file:
                                    #     files_temp.append(file)  # TODO: why

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
    if type(overwrite) == str and os.path.isfile(file_out):
        out_of_date = dt.datetime.strptime(overwrite, '%Y-%m-%d')
        file_date = dt.datetime.strptime(
            time_p.strftime("%Y-%m-%d", time_p.localtime(
                os.path.getctime(file_out))), '%Y-%m-%d')
        if out_of_date > file_date:
            print('exists: ' + file_out + '\n' +
                  ' ... but out-of-date as ' +
                  out_of_date.strftime("%Y-%m-%d") + ' > ' +
                  file_date.strftime("%Y-%m-%d"))
            overwrite = True
        else:
            overwrite = False
    else:
        overwrite = False

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
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    # "20210604",  # case01
    # "20210620", "20210621",  # case02
    # "20210628", "20210629",  # case03
    # "20220519", "20220520",  # case04
    # "20220623", "20220624", "20220625",  # case05
    # "20220626", "20220627", "20220628",  # case06+07
    # "20220630", "20220701",  # case08
    "20210713",  # case09
    "20210714",  # case09
    # "20221222",  # case10
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
    '90grad',
]
moments = ['CMAP', 'DBSNRH', 'DBZH', 'RHOHV', 'UPHIDP', 'ZDR', 'SNRHC',
           'VRADH', ]
overwrite = False
overwrite = True
overwrite = '2025-01-28'
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                load_all_moms(date, location, elevation_deg,
                              mode, moments, overwrite)

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
    '90grad',
]
moments = ['CMAP', 'DBSNRH', 'DBZH', 'RHOHV', 'UPHIDP', 'ZDR', 'SNRHC',
           'VRADH', ]
overwrite = False
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
# for date in DATES:
#     for location in LOCATIONS:
#         for elevation_deg in ELEVATIONS:
#             for mode in MODE:
#                 load_all_moms(date, location, elevation_deg,
#                               mode, moments, overwrite)

# --------------------------------------------------------------------------- #
# CONTINUE?
# import DWD_obs_to_MIUB_obs_2_correct_rho_hv
