# !/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 22.05.24                                                 #
# DWD_obs_to_MIUB_obs_6_calibrate_zdr.py                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 5: calibrate ZDR.                                                      #
#         Adapted from Velibor Pejcic                                         #
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
def calibrate_zdr_with_QVP(date, location, elevation_deg=5.5, mode='pcp',
                           overwrite=False):
    """
    calibrate zdr.
    Parameter
    ---------
    date : 'yyyymmdd' date string.
    location : 'rrr' 3-letter string for radar location.
    elevation_deg : elevations in degrees, set to 5.5 for precipitation scan
                    (as this is the sweep 0 for the volume).
    modes : set 'vol' for volume and 'pcp' for precipitation.
    overwrite : Bool;, if *zdr_off*-output exists, it can be
                       overwritten
    """
    print(mode + ' QVP ' +
          (str(elevation_deg) if (mode == 'vol') else '') +
          ' ' + location + ' ' + date)
    sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                               float(elevation_deg))[0][0])
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    folder_in = "/".join([header.dir_obs_qvp + '*', year,
                          year + '-' + mon,
                          year + '-' + mon + '-' + day,
                          location, mode + '*', sweep])
    nc_file_mom = glob.glob(folder_in + '/*_polmoms_nc_*')
    if len(nc_file_mom) > 1:
        print('mom: too many files')
        return
    elif len(nc_file_mom) == 0:
        print('mom: no files')
        return
    else:
        nc_file_mom = nc_file_mom[0]

    path_out_nc = nc_file_mom.replace('_polmoms_', '_zdroffQVPpolmoms_')
    if overwrite or not os.path.exists(path_out_nc):
        # --------------------------------------------------------------- #
        # vol/pcp or bird bad                                             #
        # --------------------------------------------------------------- #
        if mode == '90grad':
            print('todo')
        else:
            data = dttree.open_datatree(nc_file_mom).to_dataset()
            swp_cf = data.copy()

        swp_cf = swp_cf.chunk(chunks=-1)

        # determine melting layer
        swp_ML = swp_cf.where((swp_cf.height > swp_cf.mlh_bottom)&
                              (swp_cf.height < swp_cf.mlh_top))

        # calc rho min in ML
        rho_min = swp_ML.RHOHV_NC2P.min(dim="height", skipna=True)

        # determine ML top
        top_height_ML = swp_cf.mlh_top

        # check aggregation gradient:
        swp_agg_g = swp_cf.where((swp_cf.ZH_AC > 0) &
                                 (swp_cf.RHOHV_NC2P > 0.7) &
                                 (swp_cf.PHI_NC < 12) &  # Veli
                                 (swp_cf.KDP_NC < 0.1) &
                                 (rho_min < 0.95) &
                                 (swp_cf.height > top_height_ML) &
                                 (swp_cf.height < top_height_ML + 3))
        # negativ for increasing zh with height:
        agg_diff = swp_agg_g.ZH_AC.diff(dim='height')
        # negativ and already per km:
        vert_resol = (swp_agg_g.height[1] - swp_agg_g.height[2])
        agg_diff = agg_diff.median(dim=['height'], skipna=True) / vert_resol

        # agg_diff.compute()
        # print(agg_diff[:].data)
        # print(agg_diff[:].data.compute())
        # print(vert_resol.data)

        # check aggregation size:
        swp_agg_s = swp_cf.where((swp_cf.ZH_AC > 0) &
                                 (swp_cf.RHOHV_NC2P > 0.7) &
                                 (rho_min < 0.95) &
                                 (swp_cf.PHI_NC < 12) &  # Veli
                                 (swp_cf.KDP_NC < 0.1) &
                                 (swp_cf.height > top_height_ML) &
                                 (swp_cf.height < top_height_ML + 1))
        agg_med = swp_agg_s.ZH_AC.median(dim=['height'], skipna=True)

        # used for DAS:
        swp_mask = swp_cf.where((swp_cf.ZH_AC > 0) &  # HM at all
                                # (swp_cf.RHOHV_NC2P > 0.98) &
                                (swp_cf.RHOHV_NC2P > 0.7) &
                                (rho_min < 0.95) &
                                (swp_cf.PHI_NC < 12) &  # low attenuation
                                (swp_cf.KDP_NC < 0.1) &
                                (agg_diff > 3) &  # aggregation gradient down
                                (agg_med > 15) &  # aggregation
                                (swp_cf.height > top_height_ML) &
                                (swp_cf.height < top_height_ML + 1))

        # zdr_1 = np.nanmedian(swp_mask.ZDR_AC_OC)
        zdr_qvp = swp_mask.ZDR_AC_OC.median(dim=['height'])
        mask = np.isnan(zdr_qvp.data.compute())
        idx = np.where(~mask, np.arange(mask.shape[0]), 0)  # 0 or zdr_1?
        np.maximum.accumulate(idx, axis=0, out=idx)

        # print(idx)

        out = zdr_qvp.data[idx]

        # print(out.compute())

        zdr_qvp.values = out

        # print(zdr_qvp.values)
        # zh_ppi = swp_mask.ZH_AC.median(dim=['height'])

        zdr_zh_qvp_n = swp_mask.ZDR_AC_OC.count(dim=['height'])

        # zdroffset = zdr_1 - 0.35 + 0.01 * zh_1  # S band?!
        # zdroffset_ppi = zdr_ppi - 0.35 + 0.01 * zh_ppi  # S band?!
        # zdroffset = zdr_1 - 0.15  # S band?!

        zdroffset_qvp = zdr_qvp - 0.15  # S band?!
        nm = np.sum(~np.isnan(swp_mask.ZDR_AC_OC.values))

        # zdr_off_das_qvp
        data['zdr_off_das_qvp'] = zdroffset_qvp
        data['zdr_off_das_qvp'].attrs["short_name"] = 'ZDR off DAS QVP'
        data['zdr_off_das_qvp'].attrs["long_name"] = \
            'ZDR offset from dry aggregated snow after QVP'
        data['zdr_off_das_qvp'].attrs["units"] = 'dB'
        data['zdr_off_das_qvp'].attrs["comment"] = 'to subtract from ZDR'

        # zdr_off_das_n_qvp
        data['zdr_off_das_n_qvp'] = zdr_zh_qvp_n
        data['zdr_off_das_n_qvp'].attrs["short_name"] = 'nm for DAS QVP'
        data['zdr_off_das_n_qvp'].attrs["long_name"] = \
            'number values for dry aggregated snow after QVP'
        data['zdr_off_das_n_qvp'].attrs["units"] = '1'

        # zdr_2
        data['ZDR_AC_OC_OCDAS'] = data.ZDR_AC_OC - zdroffset_qvp
        data['ZDR_AC_OC_OCDAS'].attrs["short_name"] = 'ZDR_AC_OC_OCDAS'
        data['ZDR_AC_OC_OCDAS'].attrs["long_name"] = \
            'ZDR twice offset corrected dfrom DAS after QVP'
        data['ZDR_AC_OC_OCDAS'].attrs["units"] = 'db'
        data['ZDR_AC_OC_OCDAS'].attrs["coordinates"] = 'time height'

        dtree = dttree.DataTree(name="root")
        dttree.DataTree(data, name=f"sweep_{int(sweep)}",
                        parent=dtree)
        print('saving: ... ' + path_out_nc.split('/')[-1] + ' ...')
        dtree.load().to_netcdf(path_out_nc)
        data.close()
    else:
        print('exists: ' + path_out_nc + ' -> continue')

    return


# --------------------------------------------------------------------------- #
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20210714",  # case09
]
LOCATIONS = [
    'ess',
]
overwrite = False

include_sweep = np.array([
    False,
    False, False, False, False, False, False,
    False, True, False, False,
    False,
])
elevation_degs = np.array([
    5.5,
    5.5, 4.5, 3.5, 2.5, 1.5, 0.5,
    8., 12., 17., 25.,
    5.5,
])
modes = np.array([
    'pcp',
    'vol', 'vol', 'vol', 'vol', 'vol', 'vol',
    'vol', 'vol', 'vol', 'vol',
    '90grad'
])
elevation_degs = elevation_degs[include_sweep]
modes = modes[include_sweep]
# sorting volume
elevation_degs_2_sort = elevation_degs.copy()
elevation_degs_2_sort[modes == 'pcp'] = 0
elevation_degs_2_sort[modes == '90grad'] = 90
sort_i = elevation_degs_2_sort.argsort()
elevation_degs = elevation_degs[sort_i]
modes = modes[sort_i]
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg, mode in zip(elevation_degs, modes):
            calibrate_zdr_with_QVP(date, location, elevation_deg=elevation_deg,
                                   mode=mode, overwrite=overwrite)

