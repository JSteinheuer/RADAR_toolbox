#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 04.06.24                                                 #
# DWD_obs_to_MIUB_obs_7_combine_pol_mom_nc.py                                 #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# STEP 7: combine all noise corrected polarimetric moments                    #
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
        return
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
            return
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
        data_ac.ZDR_AC.reindex_like(other=data_rho.RHOHV_NC2P,
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
# NEW CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    "20210604",  # case01  # TODO: redo because of time and temp!
    "20210620", "20210621",  # case02  # TODO: redo because of time and temp!
    "20210628", "20210629",  # case03  # TODO: redo because of time and temp!
    "20220519",  # "20220520",  # case04  # TODO: redo because of time and temp!
    "20220623", "20220624", "20220625",  # case05 # TODO: redo because of tnt!
    "20220626", "20220627", "20220628",  # case06+07  # TODO: redo because of t
    "20220630", "20220701",  # case08  # TODO: redo because of temp!
    "20210714",  # case09 # TODO: redo because of temp!
    "20221222",  # case10 # TODO: redo because of temp!
]
LOCATIONS = [
    'asb', 'boo', 'drs', 'eis', 'ess', 'fbg',
    'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    'oft', 'pro', 'ros', 'tur', 'umd',
]
ELEVATIONS = np.array([
    5.5,
    4.5, 3.5, 2.5, 1.5, 0.5, 8.0, 12.0, 17.0,
    25.0,
])
MODE = [
    'pcp',
    'vol',
]
overwrite = True
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                other_zdr_off_day = ''
                n_zdr_lowest = 1000
                std_zdr_highest = 2
                method_zdr_priorities = ['BB', 'SD_V', 'LR_V' 'SD_I', 'LR_I']
                if (date == '20210604') & (location in ['asb']):
                    other_zdr_off_day = '20210621'
                if (date == '20210604') & (location in ['boo']):
                    method_zdr_priorities = ['LR_V']
                if (date == '20210620') & \
                        (location in ['asb', 'boo', 'hnr', 'ros']):
                    other_zdr_off_day = '20210621'
                if (date == '20210621') & \
                        (location in ['isn', 'tur']):
                    other_zdr_off_day = '20210620'
                if (date == '20210628') & \
                        (location in ['boo']):
                    method_zdr_priorities = ['SD_V']
                if (date == '20210628') & \
                        (location in ['drs', 'neu']):
                    other_zdr_off_day = '20210629'
                if (date == '20210628') & \
                        (location in ['pro']):
                    std_zdr_highest = 4.1
                if (date == '20210629') & \
                        (location in ['boo']):
                    method_zdr_priorities = ['SD_V']
                    other_zdr_off_day = '20210628'
                if (date == '20210629') & \
                        (location in ['pro']):
                    other_zdr_off_day = '20210628'
                    std_zdr_highest = 4.1
                if (date == '20220519') & \
                        (location in ['pro', 'mem']):
                    other_zdr_off_day = '20220520'
                if (date == '20220623') & \
                        (location in ['asb', 'boo', 'drs', 'eis', 'fld',
                                      'hnr', 'neu', 'oft', 'pro', 'umd']):
                    other_zdr_off_day = '20220624'
                if (date in ['20220623', '20220625']) & (location == 'ros'):
                    other_zdr_off_day = '20220624'
                    std_zdr_highest = 2.5
                if (date in ['20220624']) & (location == 'ros'):
                    std_zdr_highest = 2.5
                if (date == '20220625') & \
                        (location in ['asb', 'boo', 'ess', 'fbg', 'fld', 'hnr',
                                      'mem', 'nhb', 'oft', 'tur', 'umd']):
                    other_zdr_off_day = '20220624'
                if (date == '20220626') & \
                        (location in ['boo', 'drs', 'eis', 'isn',
                                      'mem', 'neu', 'umd']):
                    other_zdr_off_day = '20220627'
                if (date == '20220626') & \
                        (location in ['pro']):
                    other_zdr_off_day = '20220627'
                    std_zdr_highest = 2.1
                if (date == '20220627') & \
                        (location in ['pro']):
                    std_zdr_highest = 2.1
                if (date == '20220628') & \
                        (location in ['asb', 'boo', 'ess', 'fld',
                                      'hnr', 'isn', 'nhb', 'oft']):
                    other_zdr_off_day = '20220627'
                if (date == '20220630') & \
                        (location in ['boo', 'drs', 'eis', 'neu',
                                      'pro', 'ros', 'umd']):
                    other_zdr_off_day = '20220701'

                combine_pol_mom_nc(date=date, location=location,
                                   elevation_deg=elevation_deg, mode=mode,
                                   overwrite=overwrite,
                                   n_zdr_lowest=n_zdr_lowest,
                                   std_zdr_highest=std_zdr_highest,
                                   other_zdr_off_day=other_zdr_off_day,
                                   method_zdr_priorities=method_zdr_priorities,
                                   dir_data_obs=header.dir_data_obs)

# --------------------------------------------------------------------------- #
# OLD CASES                                                                   #
# --------------------------------------------------------------------------- #
# SET PARAMS:
DATES = [
    # "20170719",
    "20170725",
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
overwrite = True
# --------------------------------------------------------------------------- #
# START: Loop over cases, dates, and radars:
for date in DATES:
    for location in LOCATIONS:
        for elevation_deg in ELEVATIONS:
            for mode in MODE:
                other_zdr_off_day = ''
                n_zdr_lowest = 1000
                std_zdr_highest = 2
                method_zdr_priorities = ['BB', 'SD_V', 'LR_V' 'SD_I', 'LR_I']
                if date == '20170719':
                    method_zdr_priorities = ['SD_V', 'LR_V' 'SD_I', 'LR_I']

                combine_pol_mom_nc(date=date, location=location,
                                   elevation_deg=elevation_deg, mode=mode,
                                   overwrite=overwrite,
                                   n_zdr_lowest=n_zdr_lowest,
                                   std_zdr_highest=std_zdr_highest,
                                   other_zdr_off_day=other_zdr_off_day,
                                   method_zdr_priorities=method_zdr_priorities,
                                   dir_data_obs=header.dir_data_obs)

# --------------------------------------------------------------------------- #
