#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 18.01.24                                                 #
# DWD_obs_to_MIUB_obs_2_correct_rho_hv.py                                     #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 2: correct for rho_hv.                                                 #
#         Adapted from Julian Giles:                                          #
#         radar_processing_scripts/correct_rhohv.py                           #
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

sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).

# from radar_processing_scripts import utils
from xhistogram.xarray import histogram
import xradar as xd


# --------------------------------------------------------------------------- #
# FUNCTIONS: J. Steinheuer modified functions from J. Gilles


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
def correct_rho_hv(date, location , elevation_deg=5.5, mode='vol',
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
    ELEVATIONS_ALL = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
                               8.0, 12.0, 17.0, 25.0])
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    sweep = '0' + str(np.where(ELEVATIONS_ALL ==
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
# SET PARAMS:

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
LOCATIONS = ['asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld',  'hnr', 'isn',
             'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd', ]
ELEVATIONS = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0])
MODE = ['pcp', 'vol']
overwrite = False

# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:

# # DATES = ['20210604']
# DATES = ['20210714']
# LOCATIONS = ['pro']
# ELEVATIONS = np.array([5.5])
# MODE = ['vol']
# overwrite = True

for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                correct_rho_hv(date, location, elevation_deg,
                               mode, overwrite)


# --------------------------------------------------------------------------- #
# OLD CASES                                                                   #
# --------------------------------------------------------------------------- #
# go to ags!
# header.dir_data_vol = '/automount/ags/operation_hydrometeors/data/Syn_vol/'
# header.dir_data_qvp = '/automount/ags/operation_hydrometeors/data/QVP/'
# header.dir_data_mod = '/automount/ags/operation_hydrometeors/data/mod/'
# header.dir_data_era5 = '/automount/ags/operation_hydrometeors/data/ERA5/'
# header.dir_projects = '/automount/user/s6justei/PyCharm/PyCharmProjects/'
# header.dir_data_obs = '/automount/ags/operation_hydrometeors/data/obs/'
# header.dir_data_obs_realpep = '/automount/realpep/upload/RealPEP-SPP/DWD-CBand/'
# header.folder_plot = '/automount/ags/operation_hydrometeors/plots/'
# header.folder_qvp_plot = '/automount/ags/operation_hydrometeors/plots/QVPs/'
# header.folder_ppi_plot = '/automount/ags/operation_hydrometeors/plots/PPIs/'

# --------------------------------------------------------------------------- #
# SET PARAMS:

DATES = ["20170719",
         ]
LOCATIONS = ['pro', 'umd', 'nhb', 'fld',
             # 'asb', 'boo', 'drs', 'eis', 'ess', 'fbg', 'fld',  'hnr', 'isn',
             # 'mem', 'neu', 'nhb', 'oft', 'pro', 'ros', 'tur', 'umd',
             ]
# ELEVATIONS = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0, 25.0])
ELEVATIONS = np.array([5.5, 12.0, ])
MODE = ['pcp', 'vol']  # TODO: '90grad' Birth Bath ?!
moments = ['CMAP', 'DBSNRH', 'DBZH', 'RHOHV', 'UPHIDP', 'ZDR', 'SNRHC']
overwrite = False
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:

for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                correct_rho_hv(date, location, elevation_deg,
                               mode, overwrite)

# import DWD_obs_to_MIUB_obs_3_ERA5_temp_to_RAD
