#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# SET_SYN_RADAR.py                                                            #
#                                                                             #
# Functions to calculate synthetic volume scans from EMVORADO and ICON.       #
# --------------------------------------------------------------------------- #

import multiprocessing as mp
import numpy as np
import glob
import pandas as pd
from pathlib import Path
import os
from osgeo import osr
import wradlib as wrl
import xarray as xr
import pyinterp
import time as time_p
import datetime as dt

import HEADER_RADAR_toolbox as header


def rad_dict(xband_res=None):
    """Return dictionary with radar_locations and ids.

    Args:
        xband_res: If not None, it needs to be an integer below 1000 (m)
        of xband resolution.
    Returns:
        radar_dict: dictionary with 16 Cband radars (NNN: nnnnnn) and
        potentially the two Xbands (BOX:66nnn, JUX:67nnn).
    """

    radar_dict = {
        'BOO': '010132',
        'DRS': '010488',
        'EIS': '010780',
        'ESS': '010410',
        'FBG': '010908',
        'FLD': '010440',
        'HNR': '010339',
        'ISN': '010873',
        'MEM': '010950',
        'NEU': '010557',
        'NHB': '010605',
        'OFT': '010629',
        'PRO': '010392',
        'ROS': '010169',
        'TUR': '010832',
        'UMD': '010356',
    }
    if xband_res is not None:
        if int(xband_res) < 1000:
            radar_dict['BOX'] = '066' + str(int(xband_res)).zfill(3)
            radar_dict['JUX'] = '067' + str(int(xband_res)).zfill(3)

    return radar_dict


def get_path_syn_volume(date, time, spin_up_mm,
                        radar_id, dir_data, da_run='',
                        icon_run='', icon_emvorado_run='',
                        icon_folder='ICONdata',
                        naming2026=False):
    """Get paths to ICON/EMVORADO outputs for synthetic volume scans.

    Get the right path to the ICON forecast file and the synthetic C band
    volume file for specific date, spin-up time, and radar. All data needs
    to lay correct in subfolders of >dir_data<.

    Args:
        date: Date string in yyyymmdd.
        time: Time string in HHMMSS.
        spin_up_mm: Lowest spin-up time since last DA as string in MM.
        radar_id: Radar identification number.
        dir_data: Path to folder containing data in within.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        icon_folder: ICON folder string. Default 'ICONdata'
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)

    Returns:
        A dict with corresponding folder names for the forecast of ICON and
        secondly EMVORADO. Further the two specific netcdf-files of the ICON -
        and EMVORADO forecast. No valid data results in empty entries in dict.
        For example:

        {dir_of_fc: '/20170725/ASS_2211/MAIN_2203.0/ICONdata/' + \
            '20170725210000/main0200/',
         file_fc: 'fc_R19B07.20170725230000_55min.nc',
         dir_of_vol: '/20170725/ASS_2211/MAIN_2211.0/EMVO_00000000.2/' + \
            20170725210000/main0200/det/' + \
            'cdfin_volscan_201707252300-201707260100/',
         'file_vol': 'cdfin_allsim_id-010392_201707252355_201707252355' + \
            '_volscan.nc'}
    """

    t_spin = pd.Timedelta(spin_up_mm + 'min')
    dti_fc = pd.Timestamp(date + time)

    # time of data assimilation (full hour with at least spin-up time since DA)
    dti_da = (pd.Timestamp(date + time) - t_spin - pd.Timedelta('1799s')
              ).round('1h')

    # last full 3hourly time from DA needed for time folder
    dti_3hda = (dti_da - pd.Timedelta('5399s')).round('3h')

    # forecast time (difference from DA to time and at least spin-up)
    t_fc = dti_fc - dti_da

    # forecast scan
    # dir_of_fc = dir_data + dti_fc.strftime('%Y%m%d') + '/' + \
    #             da_run + '/' + icon_run + '/ICONdata/' + \
    #             dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
    #             'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_fc = dir_data + '*/' + \
                da_run + '/' + icon_run + '/' + icon_folder + '/' + \
                dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_fc = glob.glob(dir_of_fc + '*')
    if len(dir_of_fc) > 0:
        dir_of_fc = dir_of_fc[0] + '/'
    # else:  # previous day folder (if DA at 21:00)?
    #     dir_of_fc = dir_data + \
    #                 dti_3hda.strftime('%Y%m%d') + '/' + \
    #                 da_run + '/' + icon_run + '/ICONdata/' + '/' + \
    #                 dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
    #                 'main0' + str(int(dti_da.strftime('%H')) % 3) + '00' + \
    #                 icon_run + '/'
    #     dir_of_fc = glob.glob(dir_of_fc + '*')
    #     if len(dir_of_fc) > 0:
    #         dir_of_fc = dir_of_fc[0] + '/'

    if not dir_of_fc:
        dir_of_fc = ''
        file_fc = ''
    else:
        fc_minutes = str(int(t_fc.seconds / 60 % 60))
        fc_hours = str(int(t_fc.seconds / 3600))
        fc_time_str = '0000'[:-len(fc_hours)] + fc_hours + \
                      '00'[:-len(fc_minutes)] + fc_minutes + '00'
        file_fc = 'fc_R19B07.' + dti_da.strftime('%Y%m%d%H%M%S') + '_' + \
                  fc_time_str + '.RADOLAN.nc'  # .RADOLAN or .RADOLAN.nc?!
        if not os.path.isfile(dir_of_fc + file_fc):
            file_fc = 'fc_R19B07.' + dti_da.strftime('%Y%m%d%H%M%S') + \
                      '_' + str(int(t_fc.seconds / 60)) + 'min.nc'  # or nc
            if not os.path.isfile(dir_of_fc + file_fc):
                file_fc = 'fc_R19B07.' + dti_da.strftime(
                    '%Y%m%d%H%M%S') + '_' + fc_time_str + '.RADOLAN'
                if not os.path.isfile(dir_of_fc + file_fc):
                    file_fc = 'fc_R19B07.' + dti_da.strftime(
                        '%Y%m%d%H%M%S') + '_' + fc_time_str + '.SW'
                    if not os.path.isfile(dir_of_fc + file_fc):
                        dir_of_fc = ''

    # volume scan
    if naming2026:
        vol_name_ending='_PPI0080.nc'
    else:
        vol_name_ending='_volscan.nc'

    # dir_of_vol = dir_data + dti_fc.strftime('%Y%m%d') + '/' + \
    #              da_run + '/' + icon_emvorado_run + '/' + \
    #              dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
    #              'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_vol = dir_data + '*/' + \
                 da_run + '/' + icon_emvorado_run + '/' + \
                 dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                 'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_vol = glob.glob(dir_of_vol + '*')
    if len(dir_of_vol) > 0:
        dir_of_vol_tmp = dir_of_vol[0] + '/det' + '/cdfin_volscan_' + \
                         dti_da.strftime('%Y%m%d%H%M') + '-' + \
                         (dti_da +
                          pd.Timedelta('2h')).strftime('%Y%m%d%H%M') + '/'
        if not os.path.isdir(dir_of_vol_tmp):
            dir_of_vol_tmp = dir_of_vol[0] + '/det' + '/cdfin_volscan_' + \
                             dti_da.strftime('%Y%m%d%H%M') + '-' + \
                             (dti_da +
                              pd.Timedelta('3h')).strftime('%Y%m%d%H%M') + '/'

        dir_of_vol = dir_of_vol_tmp

    # else:  # previous day folder (if DA at 21:00)?
    #     dir_of_vol = dir_data + dti_3hda.strftime('%Y%m%d') + '/' + \
    #                  da_run + '/' + icon_emvorado_run + '/' + \
    #                  dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
    #                  'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    #     dir_of_vol = glob.glob(dir_of_vol + '*')
    #     if len(dir_of_vol) > 0:
    #         dir_of_vol = dir_of_vol[0] + '/det' + '/cdfin_volscan_' + \
    #                      dti_da.strftime('%Y%m%d%H%M') + '-' + \
    #                      (dti_da + pd.Timedelta('2h')).strftime(
    #                          '%Y%m%d%H%M') + '/'

    if not dir_of_vol:
        dir_of_vol = ''
        file_vol = ''
    else:
        file_vol = 'cdfin_allsim_id-' + radar_id + '_' + \
                   dti_fc.strftime('%Y%m%d%H%M') + '_' + \
                   dti_fc.strftime('%Y%m%d%H%M') + vol_name_ending
        if not os.path.isfile(dir_of_vol + file_vol):
            file_vol = ''

    return dict(dir_of_fc=dir_of_fc, file_fc=file_fc,
                dir_of_vol=dir_of_vol, file_vol=file_vol)


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
        method: 'Nearest' or 'Linear' for method to interpolate mod fields to
                rad grid.

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


def create_vol_nc_old(time_start='2017072500', time_end='2017072506',
                      dir_data_in=header.dir_data_mod,
                      dir_data_out=header.dir_data_vol,
                      radar_loc='PRO', radar_id='010392', spin_up_mm=30,
                      da_run='', icon_run='', icon_emvorado_run='',
                      overwrite=False, include_icon=True, include_emv=True,
                      method='Nearest',naming2026=False):
    """
    Create a synthetic volume scan from EMVORADO and ICON data.

     Args:
        time_start: start time in >yyyymmddhh<.
        time_end: end time in >yyyymmddhh<.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs
        dir_data_out: directory for the output.
        radar_loc: string  >RRR< naming the radar
        radar_id: string (a number) >010NNN< of that radar.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        overwrite: If True, then process only if output is not existing.
        include_icon: If True, ICON variables are included.
        include_emv: If True, synthetic pol. var from EVMORADO are included.
        method: 'Nearest' or 'Linear' for method to interpolate mod fields to
                rad grid.
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)
    Returns:
    """

    spin_up_mm = str(spin_up_mm)
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='left')
    dir_out = dir_data_out + dti[0].strftime('%Y%m%d') + '/' + \
              da_run + '/' + icon_emvorado_run + '/' + \
              str(spin_up_mm) + 'min_spinup/'
    if not include_icon:
        file_out = 'EMV_Vol_'
        if not include_emv:
            print('Nothing to do. Please include ICON or EMVORADO!')
            return
    elif not include_emv:
        file_out = 'ICON_Vol_'
        dir_out = dir_out.replace(icon_emvorado_run, icon_run + '/ICONdata')
    else:
        file_out = 'Syn_Vol_'

    file_out = file_out + radar_loc + '_' + dti[0].strftime('%Y%m%d%H%M') + \
               '_' + dti[-1].strftime('%Y%m%d%H%M') + '.nc'
    if type(overwrite) == str and os.path.isfile(dir_out + file_out):
        out_of_date = dt.datetime.strptime(overwrite, '%Y-%m-%d')
        file_date = dt.datetime.strptime(
            time_p.strftime("%Y-%m-%d", time_p.localtime(
                os.path.getctime(dir_out + file_out))), '%Y-%m-%d')
        if out_of_date > file_date:
            print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
            print(file_out + ' exists;\n' +
                  ' ... but out-of-date as ' +
                  out_of_date.strftime("%Y-%m-%d") + ' > ' +
                  file_date.strftime("%Y-%m-%d"))
            # print('___________________________')
            overwrite = True
        else:
            overwrite = False
    else:
        overwrite = False

    if os.path.isfile(dir_out + file_out) and not overwrite:
        print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
        print(file_out + ' exists;\n' +
              ' ... set: > overwrite = True < for recalculation')
        # print('___________________________')
        return

    ncells = 0
    for t_i in range(len(dti)):
        date = dti[t_i].strftime('%Y%m%d')
        time = dti[t_i].strftime('%H%M%S')
        print(radar_loc, ' - ', date, ' - ', time)
        dir_of_fc, file_fc, dir_of_vol, file_vol = \
            get_path_syn_volume(date, time, spin_up_mm, radar_id, dir_data_in,
                                da_run, icon_run, icon_emvorado_run,
                                naming2026).values()

        if ('vol_scan' not in locals() and include_icon
            and not dir_of_fc + file_fc == '') or \
                ('vol_scan' not in locals() and include_emv
                 and not dir_of_vol + file_vol == ''):
            if not os.path.isfile(dir_of_vol + file_vol):
                print('No EMVORADO input for time')
                continue

            single_scan = xr.open_dataset(dir_of_vol + file_vol, lock=False)
            single_scan = single_scan.transpose('n_range', 'n_azimuth',
                                                'records')
            print('Initialization of Variables')
            station_latitude = np.array(single_scan['station_latitude'][:])
            station_longitude = np.array(single_scan['station_longitude'][:])
            station_height = np.array(single_scan['station_height'][:])
            range_resolution = np.array(single_scan['range_resolution'][:])
            ppi_elevation = np.array(single_scan['ppi_elevation'][:])
            n_range_bins = np.array(single_scan['n_range_bins'][:])
            range_start = np.array(single_scan['range_start'][:])
            azimuth_start = np.array(single_scan['azimuth_start'][:])
            ray_azimuth = np.array(single_scan['ray_azimuth'][:])

            # range
            r = np.arange(range_start[0],
                          (range_resolution[0] * n_range_bins[0]) +
                          range_resolution[0],
                          range_resolution[0])

            # azimuth
            az = ray_azimuth[:, 0] - azimuth_start[0]

            # sitecoords
            sitecoords = [station_longitude[0],
                          station_latitude[0],
                          station_height[0]]

            # new dataset
            vol_scan = xr.Dataset()

            # initialize dims
            vol_scan = vol_scan.expand_dims(dim=dict(azimuth=az))
            vol_scan.azimuth.attrs = dict(
                standard_name='azimuth',
                comments='increasing clockwise with 0°/360° pointing ' +
                         'northwards', units='degrees')
            vol_scan = vol_scan.expand_dims(dim=dict(range=r))
            vol_scan.range.attrs = dict(standard_name='range',
                                        units='m')
            vol_scan = vol_scan.expand_dims(dim=dict(elevation=ppi_elevation))
            vol_scan.elevation.attrs = dict(standard_name='elevation',
                                            units='degrees')
            vol_scan = vol_scan.expand_dims(dim=dict(time=dti), axis=0)

            # initialize variables
            dummy_trae = np.empty(
                shape=[vol_scan.time.size, vol_scan.range.size,
                       vol_scan.azimuth.size,
                       vol_scan.elevation.size])
            dummy_trae[:] = np.nan

            # grid variables
            vol_scan['lat'] = (['range', 'azimuth', 'elevation'],
                               dummy_trae[0, :].copy(),
                               dict(standard_name='latitude',
                                    units='degrees_north'))
            vol_scan['lon'] = (['range', 'azimuth', 'elevation'],
                               dummy_trae[0, :].copy(),
                               dict(standard_name='longitude',
                                    units='degrees_east'))
            vol_scan['alt'] = (['range', 'elevation'],
                               dummy_trae[0, :, 0, :].copy(),
                               dict(standard_name='altitude',
                                    comments='height above mean sea level',
                                    units='m'))
            if include_icon:
                # ICON variables
                vol_scan['temp'] = (['time', 'range', 'azimuth', 'elevation'],
                                    dummy_trae.copy(),
                                    dict(standard_name='air temperature',
                                         units='K'))
                vol_scan['pres'] = (['time', 'range', 'azimuth', 'elevation'],
                                    dummy_trae.copy(),
                                    dict(standard_name='pressure',
                                         units='Pa'))

                # ICON specific contents
                vol_scan['qc'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(), dict(
                    standard_name='specific cloud water content',
                    units='kg kg-1'))
                vol_scan['qg'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(), dict(
                    standard_name='specific graupel content',
                    units='kg kg-1'))
                vol_scan['qh'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(),
                                  dict(standard_name='specific hail content',
                                       units='kg kg-1'))
                vol_scan['qi'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(), dict(
                    standard_name='specific cloud ice content',
                    units='kg kg-1'))
                vol_scan['qr'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(),
                                  dict(standard_name='rain mixing ratio',
                                       units='kg kg-1'))
                vol_scan['qs'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(),
                                  dict(standard_name='snow mixing ratio',
                                       units='kg kg-1'))
                vol_scan['qv'] = (['time', 'range', 'azimuth', 'elevation'],
                                  dummy_trae.copy(),
                                  dict(standard_name='specific humidity',
                                       units='kg kg-1'))

                # ICON number concentrations
                vol_scan['qnc'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration cloud droplets',
                    units='kg-1'))
                vol_scan['qng'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration graupel',
                    units='kg-1'))
                vol_scan['qnh'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration hail', units='kg-1'))
                vol_scan['qni'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration cloud ice',
                    units='kg-1'))
                vol_scan['qnr'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration rain droplet',
                    units='kg-1'))
                vol_scan['qns'] = (['time', 'range', 'azimuth', 'elevation'],
                                   dummy_trae.copy(), dict(
                    standard_name='number concentration snow',
                    units='kg-1'))

                # ICON wind variables
                vol_scan['u'] = (['time', 'range', 'azimuth', 'elevation'],
                                 dummy_trae.copy(),
                                 dict(standard_name='eastward wind',
                                      comments='zonal wind',
                                      units='m s-1'))
                vol_scan['v'] = (['time', 'range', 'azimuth', 'elevation'],
                                 dummy_trae.copy(),
                                 dict(standard_name='northward wind',
                                      comments='meridional wind',
                                      units='m s-1'))
                vol_scan['w'] = (['time', 'range', 'azimuth', 'elevation'],
                                 dummy_trae.copy(),
                                 dict(standard_name='vertical velocity',
                                      comments='upward air movement',
                                      units='m s-1'))
                vol_scan.attrs['icon_run'] = icon_run

            if include_emv:
                # EMVORADO pol. variables
                vol_scan['zrsim'] = (['time', 'range', 'azimuth', 'elevation'],
                                     dummy_trae.copy(), dict(
                    standard_name='horizontal reflectivity',
                    comments='simulated radar reflectivity', units='dBZ'))
                vol_scan['zdrsim'] = (
                    ['time', 'range', 'azimuth', 'elevation'],
                    dummy_trae.copy(), dict(
                        standard_name='differential reflectivity',
                        comments='simulated differential reflectivity',
                        units='dB'))
                vol_scan['rhvsim'] = (
                    ['time', 'range', 'azimuth', 'elevation'],
                    dummy_trae.copy(), dict(
                        standard_name='co-polar correlation coefficient',
                        comments='simulated rhohv',
                        units='1'))
                vol_scan['kdpsim'] = (
                    ['time', 'range', 'azimuth', 'elevation'],
                    dummy_trae.copy(),
                    dict(standard_name='specific differential phase',
                         comments='simulated KDP',
                         units='deg/km'))
                vol_scan.attrs['icon_emvorado_run'] = icon_emvorado_run

            # spin-up time and location of radar
            vol_scan['spin_up_time'] = ([], spin_up_mm, dict(
                standard_name='spin-up time',
                comments='lower threshold for the time that must have' +
                         ' elapsed since the last DA of the pol. variable',
                units='min'))
            vol_scan['station_latitude'] = (
                [], single_scan.station_latitude.data[0], dict(
                    standard_name='station latitude',
                    comments='latitude of the radar location',
                    units='degrees_north'))
            vol_scan['station_longitude'] = (
                [], single_scan.station_longitude.data[0], dict(
                    standard_name='station longitude',
                    comments='longitude of the radar location',
                    units='degrees_east'))
            vol_scan['station_height'] = (
                [], single_scan.station_height.data[0], dict(
                    standard_name='station height',
                    comments='height above means sea level of the radar' +
                             ' location', units='m'))

            # write grid variables
            for sweep_i in range(len(ppi_elevation)):
                el = ppi_elevation[sweep_i]
                lon, lat, alt = get_lon_lat_alt(r, az, el, sitecoords)
                vol_scan['lat'][:, :, sweep_i] = lat.transpose()
                vol_scan['lon'][:, :, sweep_i] = lon.transpose()
                vol_scan['alt'][:, sweep_i] = alt[0,
                                              :]  # TODO: check!!!! alt per azimuth constant?

            # grid for searching later the closest ICON cells
            rad_lon = vol_scan['lon'].data.flatten()
            rad_lat = vol_scan['lat'].data.flatten()
            rad_alt_3d = np.repeat(vol_scan['alt'].data[:, np.newaxis, :],
                                   vol_scan['azimuth'].shape, axis=1)
            rad_alt = np.repeat(vol_scan['alt'].data[:, np.newaxis, :],
                                vol_scan['azimuth'].shape, axis=1).flatten()

            # global attributes
            vol_scan.attrs['title'] = 'Synthetic C-band radar variables'
            vol_scan.attrs['institution'] = 'University of Bonn'
            vol_scan.attrs['history'] = 'DWD: ICON + EMVORADO'
            vol_scan.attrs['author'] = \
                'Julian Steinheuer, J.Steinheuer@uni-bonn.de'
            vol_scan.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
            vol_scan.attrs['creation_date'] = str(
                pd.Timestamp(single_scan.Creation_date))[:16]
            vol_scan.attrs['station_ID_national'] = radar_id
            vol_scan.attrs['station_name'] = radar_loc
            vol_scan.attrs['data_assimilation_run'] = da_run

            # finishing initialization
            single_scan.close()

        if include_emv:
            if not os.path.isfile(dir_of_vol + file_vol):
                print('No EMVORADO input file for time')
                continue

            single_scan = xr.open_dataset(dir_of_vol + file_vol, lock=False)
            single_scan = single_scan.transpose('n_range', 'n_azimuth',
                                                'records')
            # copy pol. variables
            zrsim = single_scan['zrsim'].data
            zdrsim = single_scan['zdrsim'].data
            rhvsim = single_scan['rhvsim'].data
            kdpsim = single_scan['kdpsim'].data

            # only applied filter is:
            mask_rho = rhvsim < 0

            zrsim[mask_rho] = np.nan
            zdrsim[mask_rho] = np.nan
            rhvsim[mask_rho] = np.nan
            kdpsim[mask_rho] = np.nan
            vol_scan['zrsim'][t_i, :, :, :] = zrsim
            vol_scan['zdrsim'][t_i, :, :, :] = zdrsim
            vol_scan['rhvsim'][t_i, :, :, :] = rhvsim
            vol_scan['kdpsim'][t_i, :, :, :] = kdpsim

            # single_scan done!
            single_scan.close()

        # read direct ICON fc:
        if include_icon:
            if not os.path.isfile(dir_of_fc + file_fc):
                print(' No input file for t=' + time)
                continue

            single_fc = xr.open_dataset(dir_of_fc + file_fc, lock=False)

            # check for ICON grid
            if not single_fc.ncells.size == ncells:
                grid_fc = xr.open_dataset(
                    dir_data_in + '/grid/hgrd_R19B07.RADOLAN.nc',
                    lock=False)
                if grid_fc.cell.size == single_fc.ncells.size:
                    fc_lon = np.rad2deg(grid_fc.clon.data)
                    fc_lat = np.rad2deg(grid_fc.clat.data)
                    ncells = grid_fc.cell.size
                    grid_fc.close()
                    levels_fc = xr.open_dataset(
                        dir_data_in + '/grid/vgrd_R19B07.RADOLAN.nc',
                        lock=False)
                    fc_alt = levels_fc.z_mc.data[0, :, :]
                    levels_fc.close()
                else:
                    grid_fc.close()
                    grid_fc = xr.open_dataset(
                        dir_data_in + '/grid/hgrd_R19B07.ICON-D2.nc',
                        lock=False)
                    if grid_fc.cell.size == single_fc.ncells.size:
                        fc_lon = np.rad2deg(grid_fc.clon.data)
                        fc_lat = np.rad2deg(grid_fc.clat.data)
                        ncells = grid_fc.cell.size
                        grid_fc.close()
                        levels_fc = xr.open_dataset(
                            dir_data_in + '/grid/vgrd_R19B07.ICON-D2.nc',
                            lock=False)
                        fc_alt = levels_fc.z_mc.data[:, :]
                        levels_fc.close()
                    else:
                        grid_fc.close()
                        print('No forcast data, as differing ncells ' +
                              'were found in ' + dir_data_in + '/grid/')
                        continue

                # rad_alt=np.repeat(sitecoords[2], rad_lon.size)
                # rad_alt=1rad_alt  # TODO: check altitudal error
                func_ipol, mask = ipol_fc_to_radgrid(
                    # lon with shape (65,258775):
                    np.repeat(fc_lon[np.newaxis, :], fc_alt.shape[0], axis=0),
                    # lat with shape (65,258775)
                    np.repeat(fc_lat[np.newaxis, :], fc_alt.shape[0], axis=0),
                    fc_alt,
                    rad_lon, rad_lat, rad_alt, method
                )
                shape_rae = vol_scan['lon'].shape

            # copy fc variables
            vol_scan['temp'][t_i, :, :, :] = \
                func_ipol(single_fc['temp'].data[0, :, :][mask]).reshape(
                    shape_rae)
            try:
                vol_scan['pres'][t_i, :, :, :] = \
                    func_ipol(single_fc['pres'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qv'][t_i, :, :, :] = \
                    func_ipol(single_fc['qv'].data[0, :, :][mask]).reshape(
                        shape_rae)

                vol_scan['qc'][t_i, :, :, :] = \
                    func_ipol(single_fc['qc'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qr'][t_i, :, :, :] = \
                    func_ipol(single_fc['qr'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qi'][t_i, :, :, :] = \
                    func_ipol(single_fc['qi'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qs'][t_i, :, :, :] = \
                    func_ipol(single_fc['qs'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qg'][t_i, :, :, :] = \
                    func_ipol(single_fc['qg'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qh'][t_i, :, :, :] = \
                    func_ipol(single_fc['qh'].data[0, :, :][mask]).reshape(
                        shape_rae)

                vol_scan['qnc'][t_i, :, :, :] = \
                    func_ipol(single_fc['qnc'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qnr'][t_i, :, :, :] = \
                    func_ipol(single_fc['qnr'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qni'][t_i, :, :, :] = \
                    func_ipol(single_fc['qni'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qns'][t_i, :, :, :] = \
                    func_ipol(single_fc['qns'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qng'][t_i, :, :, :] = \
                    func_ipol(single_fc['qng'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['qnh'][t_i, :, :, :] = \
                    func_ipol(single_fc['qnh'].data[0, :, :][mask]).reshape(
                        shape_rae)

                vol_scan['u'][t_i, :, :, :] = \
                    func_ipol(single_fc['u'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['v'][t_i, :, :, :] = \
                    func_ipol(single_fc['v'].data[0, :, :][mask]).reshape(
                        shape_rae)
                vol_scan['w'][t_i, :, :, :] = \
                    ((func_ipol(single_fc['w'].data[0, 1:, :][mask]) +
                      func_ipol(single_fc['w'].data[0, :-1, :][mask])) / 2
                     ).reshape(
                        shape_rae)
            except:
                print('some variable not in ICON fc')

            # single_fc done!
            single_fc.close()

    if 'vol_scan' in locals():
        Path(dir_out).mkdir(parents=True, exist_ok=True)
        print('start saving')
        vol_scan.to_netcdf(dir_out + file_out, unlimited_dims='time')
        print('saved')
        vol_scan.close()
        print('    ! Case done now !      ')
        print('___________________________')
    else:
        if include_icon and not include_emv:
            print('   ! No ICON for case !    ')
            print('___________________________')
        elif not include_icon and include_emv:
            print('   ! No EMVO for case !    ')
            print('___________________________')
        else:
            print('   ! No data for case !    ')
            print('___________________________')


def load_emvorado_to_radar_volume(path_or_data, rename=False):
    """
    Load and reorganize EMVORADO output into a xarray.Dataset in the same
    flavor as DWD data. Optimized for EMVORADO output with one file per
    timestep containing all elevations and variables.
    WARNING: The resulting volume has its elevations ordered from lower to
             higher and not according to the scan strategy.

    Parameters
    ----------
    path_or_data : str or nested sequence of paths or xarray.Dataset
        – Either a string glob in the form "path/to/my/files/*.nc" or an
        explicit list of files to open. Paths can be given as strings or
        as pathlib Paths. Feeds into xarray.open_mfdataset. Alternatively,
        already-loaded data in the form of an xarray.Dataset can be passed.
    rename : bool. If True, then rename the variables to DWD-like naming.

    Returns
    -------
    data_vol : xarray.Dataset

    """
    if type(path_or_data) is xr.Dataset:
        data_emvorado_xr = path_or_data
    else:
        data_emvorado_xr = xr.open_mfdataset(path_or_data,
                                             concat_dim="time",
                                             combine="nested",
                                             lock=False)

    try:
        data = data_emvorado_xr.rename_dims({"n_range": "range",
                                             "n_azimuth": "azimuth"})
    except ValueError:
        data = data_emvorado_xr

    # we make the coordinate arrays
    if "time" in data.dims and "time" not in data.coords:
        range_coord = \
            np.array([np.arange(rs, rr * rb + rs, rr) for rr, rs, rb in
                      zip(data.range_resolution[0],
                          data.range_start[0],
                          data.n_range_bins[0])])[0]
        azimuth_coord = np.array([np.arange(azs, azr * azb + azs, azr)
                                  for azr, azs, azb in
                                  zip(data.azimuthal_resolution[0],
                                      data.azimuth_start[0],
                                      data.azimuth.shape * np.ones_like(
                                          data.records))])[0]

        # create time coordinate
        time_coord = xr.DataArray([dt.datetime(int(yy), int(mm), int(dd),
                                               int(hh), int(mn), int(ss))
                                   for yy, mm, dd, hh, mn, ss in
                                   zip(data.year.isel(records=0),
                                       data.month.isel(records=0),
                                       data.day.isel(records=0),
                                       data.hour.isel(records=0),
                                       data.minute.isel(records=0),
                                       data.second.isel(records=0),
                                       )
                                   ], dims=["time"])

    elif "time" not in data.coords:
        range_coord = \
            np.array([np.arange(rs, rr * rb + rs, rr) for rr, rs, rb in
                      zip(data.range_resolution,
                          data.range_start,
                          data.n_range_bins)])[0]
        azimuth_coord = np.array([np.arange(azs, azr * azb + azs, azr)
                                  for azr, azs, azb in
                                  zip(data.azimuthal_resolution,
                                      data.azimuth_start,
                                      data.azimuth.shape * np.ones_like(
                                          data.records))])[0]

        # create time coordinate
        time_coord = xr.DataArray([dt.datetime(int(yy), int(mm), int(dd),
                                               int(hh), int(mn), int(ss))
                                   for yy, mm, dd, hh, mn, ss in
                                   zip(data.year,
                                       data.month,
                                       data.day,
                                       data.hour,
                                       data.minute,
                                       data.second
                                       )], dims=["time"])[0]

    # add coordinates for range, azimuth, time, latitude, longitude,
    # altitude, elevation, sweep_mode
    try:
        data.coords["range"] = ("range", range_coord)
    except NameError:
        pass
    try:
        data.coords["azimuth"] = ("azimuth", azimuth_coord)
    except NameError:
        pass
    try:
        data.coords["time"] = time_coord
    except NameError:
        pass
    data.coords["latitude"] = float(
        data["station_latitude"].values.flatten()[0])
    data.coords["longitude"] = float(
        data["station_longitude"].values.flatten()[0])
    try:
        data.coords["altitude"] = float(
            [ss for ss in data.attrs["Data_description"].split(" ")
             if "radar_alt_msl_mod" in ss][0].split("=")[1]
        )
    except KeyError:
        data.coords["altitude"] = float(
            data["station_height"].values.flatten()[0]
        )
    if "elevation" not in data.coords:
        if "time" in data['ray_elevation'].dims:
            data.coords["elevation"] = data["ray_elevation"].mean(
                ("azimuth", "time"))
        else:
            data.coords["ray_elevation"] = data["ray_elevation"].mean(
                "azimuth")

    data.coords["sweep_mode"] = 'azimuth_surveillance'

    # move some variables to attributes
    vars_to_attrs = ["station_name", "country_ID", "station_ID_national",
                     "station_longitude", "station_height",
                     "station_latitude", "range_resolution",
                     "azimuthal_resolution", "range_start", "azimuth_start",
                     "extended_nyquist", "high_nyquist", "dualPRF_ratio",
                     "range_gate_length", "n_ranges_averaged",
                     "n_pulses_averaged", "DATE", "TIME",
                     "year", "month", "day", "hour", "minute", "second",
                     "ppi_azimuth", "ppi_elevation", "n_range_bins"
                     ]
    for vta in vars_to_attrs:
        try:
            data[vta] = data[vta].isel(records=0, time=0)
            tmp = data[vta]
            data = data.drop_vars(vta, errors="ignore")
            data.attrs[vta] = tmp
        except KeyError:
            pass

    # add attribute "fixed_angle"
    try:
        # if one timestep
        data.attrs["fixed_angle"] = float(data.attrs["ppi_elevation"])
    except TypeError:
        # if multiple timesteps
        data.attrs["fixed_angle"] = float(
            data.attrs["ppi_elevation"].values.flatten()[0])
    except KeyError:
        data.attrs["fixed_angle"] = float(
            data["elevation"].values.flatten()[0])

    # for each remaining variable add "long_name" and "units" attribute
    for vv in data.data_vars.keys():
        try:
            data[vv].attrs["long_name"] = data[vv].attrs["Description"]
        except:
            pass
            # print("no long_name attribute in " + vv)

        try:
            data[vv].attrs["units"] = data[vv].attrs["Unit"]
        except:
            pass
            # print("no long_name attribute in " + vv)

        return data


def create_vol_nc(time_start='2021071412', time_end='2021071418',
                  dir_data_in=header.dir_data_mod,
                  dir_data_out=header.dir_data_vol,
                  radar_loc='ESS', radar_id='010410', spin_up_mm=120,
                  da_run='ASS_2411', icon_run='MAIN_2411',
                  icon_emvorado_run='MAIN_2411.1/EMVO_00510000.2',
                  overwrite=False, include_icon=True, include_emv=False,
                  icon_folder='ICONdata', naming2026=False):
    """
    Create a synthetic volume scan from EMVORADO and ICON data.

     Args:
        time_start: start time in >yyyymmddhh<.
        time_end: end time in >yyyymmddhh<.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs
        dir_data_out: directory for the output.
        radar_loc: string  >RRR< naming the radar
        radar_id: string (a number) >010NNN< of that radar.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        overwrite: If True, then process only if output is not existing, if
                   false process not if output existing. If overwrite is a
                   string of type 'yyyy-mm-dd', then a potentioal file is
                   overwritten if its creation date is older than yyyy-mm-dd.
        include_icon: If True, ICON variables are included.
        include_emv: If True, synthetic pol. var from EVMORADO are included.
        icon_folder: ICON folder string. Default 'ICONdata'
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)

    Returns:
    """

    # start
    current_time = time_p.time()

    # time
    spin_up_mm = str(spin_up_mm)
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='left')

    # output
    dir_out = dir_data_out + dti[0].strftime('%Y%m%d') + '/' + \
              da_run + '/' + icon_emvorado_run + '/' + \
              str(spin_up_mm) + 'min_spinup/'
    if not include_icon:
        file_out = 'EMV_Vol_'
        if not include_emv:
            print('Nothing to do. Please include ICON or EMVORADO!')
            return
    elif not include_emv:
        file_out = 'ICON_Vol_'
        dir_out = dir_out.replace(icon_emvorado_run,
                                  icon_run + '/ICONdata')
    else:
        file_out = 'Syn_Vol_'

    file_out = file_out + radar_loc + '_' + dti[0].strftime('%Y%m%d%H%M') + \
               '_' + dti[-1].strftime('%Y%m%d%H%M') + '.nc'
    if type(overwrite) == str and os.path.isfile(dir_out + file_out):
        out_of_date = dt.datetime.strptime(overwrite, '%Y-%m-%d')
        file_date = dt.datetime.strptime(
            time_p.strftime("%Y-%m-%d", time_p.localtime(
                os.path.getctime(dir_out + file_out))), '%Y-%m-%d')
        if out_of_date > file_date:
            print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
            print(file_out + ' exists;\n' +
                  ' ... but out-of-date as ' +
                  out_of_date.strftime("%Y-%m-%d") + ' > ' +
                  file_date.strftime("%Y-%m-%d"))
            print('___________________________')
            overwrite = True
        else:
            overwrite = False
    else:
        overwrite = False

    if os.path.isfile(dir_out + file_out) and not overwrite:
        print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
        print(file_out + ' exists;\n' +
              ' ... set: > overwrite = True < for recalculation')
        print('___________________________')
        return dir_out + file_out

    # loop over time
    ncells = 0  # dummy to start grid-loading
    for t_i in range(len(dti)):

        # time
        date = dti[t_i].strftime('%Y%m%d')
        time = dti[t_i].strftime('%H%M%S')

        # file names
        dir_of_fc, file_fc, dir_of_vol, file_vol = \
            get_path_syn_volume(date, time, spin_up_mm, radar_id, dir_data_in,
                                da_run, icon_run, icon_emvorado_run,
                                icon_folder, naming2026).values()

        # initialization of volume necessary?
        if ('vol_scan' not in locals() and include_icon
            and not dir_of_fc + file_fc == '') or \
                ('vol_scan' not in locals() and include_emv
                 and not dir_of_vol + file_vol == ''):
            if not os.path.isfile(dir_of_vol + file_vol):
                print('No EMVORADO input for time')
                continue

            print('Initialization of Variables')

            # load RADAR Vol (EMVORADO) even if its data should not be included
            radar_volume = load_emvorado_to_radar_volume(dir_of_vol + file_vol)
            radar_volume = radar_volume.transpose(
                'time', 'records', 'azimuth', 'range', ...)

            # RADAR coords
            sitecoords = [float(radar_volume.station_longitude),
                          float(radar_volume.station_latitude),
                          float(radar_volume.station_height)]

            # # write grid variables
            xyz, proj_aeqd = wrl.georef.spherical_to_centroids(
                radar_volume["range"] +
                radar_volume["range"].diff("range").mean().values / 2,
                radar_volume["azimuth"],
                radar_volume["elevation"],
                sitecoords,
                crs=None,
            )

            # xyz needs to be included
            if "x" not in radar_volume.coords:
                radar_volume = wrl.georef.georeference(radar_volume)

            # RADAR target grid in shape=raveled(elevation x azimuth x range)
            trg = np.vstack((radar_volume.x.values.ravel(),
                             radar_volume.y.values.ravel(),
                             radar_volume.z.values.ravel())).T

            # new dataset
            vol_scan = xr.Dataset()

            # initialize dims
            vol_scan = vol_scan.expand_dims(
                dim=dict(elevation=radar_volume.elevation.values))
            vol_scan.elevation.attrs = dict(standard_name='elevation',
                                            units='degrees')
            vol_scan = vol_scan.expand_dims(
                dim=dict(range=radar_volume.range.values))
            vol_scan.range.attrs = dict(standard_name='range',
                                        units='m')
            vol_scan = vol_scan.expand_dims(
                dim=dict(azimuth=radar_volume.azimuth.values))
            vol_scan.azimuth.attrs = dict(
                standard_name='azimuth',
                comments='increasing clockwise with 0°/360° pointing ' +
                         'northwards', units='degrees')
            vol_scan = vol_scan.expand_dims(dim=dict(time=dti), axis=0)
            vol_scan = vol_scan.transpose(
                'time', 'elevation', 'azimuth', 'range', )

            # initialize variables
            dummy_tear = np.empty(
                shape=[vol_scan.time.size,
                       vol_scan.elevation.size,
                       vol_scan.azimuth.size,
                       vol_scan.range.size])
            dummy_tear[:] = np.nan

            # grid variables
            vol_scan['lat'] = (['elevation', 'azimuth', 'range'],
                               dummy_tear[0, :].copy(),
                               dict(standard_name='latitude',
                                    units='degrees_north'))
            vol_scan['lon'] = (['elevation', 'azimuth', 'range'],
                               dummy_tear[0, :].copy(),
                               dict(standard_name='longitude',
                                    units='degrees_east'))
            vol_scan['alt'] = (['elevation', 'range', ],
                               dummy_tear[0, :, 0, :].copy(),
                               dict(standard_name='altitude',
                                    comments='height above mean sea level',
                                    units='m'))
            if include_icon:
                # ICON variables
                vol_scan['temp'] = (['time', 'elevation', 'azimuth', 'range'],
                                    dummy_tear.copy(),
                                    dict(standard_name='air temperature',
                                         units='K'))
                vol_scan['pres'] = (['time', 'elevation', 'azimuth', 'range'],
                                    dummy_tear.copy(),
                                    dict(standard_name='pressure',
                                         units='Pa'))

                # ICON specific contents
                vol_scan['qc'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(), dict(
                    standard_name='specific cloud water content',
                    units='kg kg-1'))
                vol_scan['qg'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(), dict(
                    standard_name='specific graupel content',
                    units='kg kg-1'))
                vol_scan['qh'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(),
                                  dict(standard_name='specific hail content',
                                       units='kg kg-1'))
                vol_scan['qi'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(), dict(
                    standard_name='specific cloud ice content',
                    units='kg kg-1'))
                vol_scan['qr'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(),
                                  dict(standard_name='rain mixing ratio',
                                       units='kg kg-1'))
                vol_scan['qs'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(),
                                  dict(standard_name='snow mixing ratio',
                                       units='kg kg-1'))
                vol_scan['qv'] = (['time', 'elevation', 'azimuth', 'range'],
                                  dummy_tear.copy(),
                                  dict(standard_name='specific humidity',
                                       units='kg kg-1'))

                # ICON number concentrations
                vol_scan['qnc'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration cloud droplets',
                    units='kg-1'))
                vol_scan['qng'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration graupel',
                    units='kg-1'))
                vol_scan['qnh'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration hail', units='kg-1'))
                vol_scan['qni'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration cloud ice',
                    units='kg-1'))
                vol_scan['qnr'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration rain droplet',
                    units='kg-1'))
                vol_scan['qns'] = (['time', 'elevation', 'azimuth', 'range'],
                                   dummy_tear.copy(), dict(
                    standard_name='number concentration snow',
                    units='kg-1'))

                # ICON wind variables
                vol_scan['u'] = (['time', 'elevation', 'azimuth', 'range'],
                                 dummy_tear.copy(),
                                 dict(standard_name='eastward wind',
                                      comments='zonal wind',
                                      units='m s-1'))
                vol_scan['v'] = (['time', 'elevation', 'azimuth', 'range'],
                                 dummy_tear.copy(),
                                 dict(standard_name='northward wind',
                                      comments='meridional wind',
                                      units='m s-1'))
                vol_scan['w'] = (['time', 'elevation', 'azimuth', 'range'],
                                 dummy_tear.copy(),
                                 dict(standard_name='vertical velocity',
                                      comments='upward air movement',
                                      units='m s-1'))
                vol_scan.attrs['icon_run'] = icon_run

            if include_emv:
                # EMVORADO pol. variables
                vol_scan['zrsim'] = (['time', 'elevation', 'azimuth', 'range'],
                                     dummy_tear.copy(), dict(
                    standard_name='horizontal reflectivity',
                    comments='simulated radar reflectivity', units='dBZ'))
                vol_scan['zdrsim'] = (
                    ['time', 'elevation', 'azimuth', 'range'],
                    dummy_tear.copy(), dict(
                        standard_name='differential reflectivity',
                        comments='simulated differential reflectivity',
                        units='dB'))
                vol_scan['rhvsim'] = (
                    ['time', 'elevation', 'azimuth', 'range'],
                    dummy_tear.copy(), dict(
                        standard_name='co-polar correlation coefficient',
                        comments='simulated rhohv',
                        units='1'))
                vol_scan['kdpsim'] = (
                    ['time', 'elevation', 'azimuth', 'range'],
                    dummy_tear.copy(),
                    dict(standard_name='specific differential phase',
                         comments='simulated KDP',
                         units='deg/km'))
                vol_scan.attrs['icon_emvorado_run'] = icon_emvorado_run

            # to include always: spin-up time and location of radar
            vol_scan['spin_up_time'] = ([], spin_up_mm, dict(
                standard_name='spin-up time',
                comments='lower threshold for the time that must have' +
                         ' elapsed since the last DA of the pol. variable',
                units='min'))
            vol_scan['station_latitude'] = (
                [], radar_volume.station_latitude.data, dict(
                    standard_name='station latitude',
                    comments='latitude of the radar location',
                    units='degrees_north'))
            vol_scan['station_longitude'] = (
                [], radar_volume.station_longitude.data, dict(
                    standard_name='station longitude',
                    comments='longitude of the radar location',
                    units='degrees_east'))
            vol_scan['station_height'] = (
                [], radar_volume.station_height.data, dict(
                    standard_name='station height',
                    comments='height above means sea level of the radar' +
                             ' location', units='m'))

            # lat lon alt per scan (could be improved, but works)
            ra = radar_volume["range"].values.copy()
            az = radar_volume["azimuth"].values.copy()
            for sweep_i in range(radar_volume.elevation.size):
                el = radar_volume["elevation"].values[sweep_i].copy()
                lon, lat, alt = get_lon_lat_alt(ra,
                                                az,
                                                el, sitecoords)
                vol_scan['lat'][sweep_i, :, :] = lat
                vol_scan['lon'][sweep_i, :, :] = lon
                vol_scan['alt'][sweep_i, :] = alt[0, :]

            # global attributes
            vol_scan.attrs['title'] = 'Synthetic C-band radar variables'
            vol_scan.attrs['institution'] = 'University of Bonn'
            vol_scan.attrs['history'] = 'DWD: ICON + EMVORADO'
            vol_scan.attrs['author'] = \
                'Julian Steinheuer, J.Steinheuer@uni-bonn.de'
            vol_scan.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
            vol_scan.attrs['creation_date'] = str(
                pd.Timestamp(radar_volume.Creation_date))[:16]
            vol_scan.attrs['station_ID_national'] = radar_id
            vol_scan.attrs['station_name'] = radar_loc
            vol_scan.attrs['data_assimilation_run'] = da_run

            # finishing initialization
            radar_volume.close()

        print(radar_loc, ' - ', date, ' - ', time)

        # per every timestep
        if include_emv:
            if not os.path.isfile(dir_of_vol + file_vol):
                print('No EMVORADO input file for time')
                continue

            radar_volume = load_emvorado_to_radar_volume(dir_of_vol + file_vol)
            radar_volume = radar_volume.transpose(
                'time', 'records', 'azimuth', 'range', ...)

            # copy pol. variables
            zrsim = radar_volume['zrsim'].data
            zdrsim = radar_volume['zdrsim'].data
            rhvsim = radar_volume['rhvsim'].data
            kdpsim = radar_volume['kdpsim'].data

            # only applied filter is:
            mask_rho = rhvsim < 0

            zrsim[mask_rho] = np.nan
            zdrsim[mask_rho] = np.nan
            rhvsim[mask_rho] = np.nan
            kdpsim[mask_rho] = np.nan
            vol_scan['zrsim'][t_i, :, :, :] = zrsim[0, :]
            vol_scan['zdrsim'][t_i, :, :, :] = zdrsim[0, :]
            vol_scan['rhvsim'][t_i, :, :, :] = rhvsim[0, :]
            vol_scan['kdpsim'][t_i, :, :, :] = kdpsim[0, :]

            # single_scan done!
            radar_volume.close()

        # read direct ICON fc:
        if include_icon:
            if not os.path.isfile(dir_of_fc + file_fc):
                print(' No input file for t=' + time)
                continue

            single_fc = xr.open_dataset(dir_of_fc + file_fc,
                                        chunks="auto", lock=False)
            single_fc = single_fc.transpose('time', 'height', 'ncells', ...)

            # get rid of height_2 (half level center as in w)
            single_fc = single_fc.interp(
                height_2=np.arange(1.5, single_fc.height_2.size))
            for vv in single_fc.data_vars:
                if "height_2" in single_fc[vv].dims:
                    vv_alone = single_fc[vv]
                    vv_alone = vv_alone.assign_coords(
                        height_2=np.arange(1., 66))
                    vv_alone = vv_alone.rename({"height_2": "height"})
                    single_fc[vv] = vv_alone

            single_fc = single_fc.drop_vars({'height_2'})

            # check for ICON grid
            if not single_fc.ncells.size == ncells:
                if single_fc.ncells.size == 258775:
                    ncells = 258775  # RADOLAN
                    grid_fc = xr.open_dataset(
                        dir_data_in + '/grid/hgrd_R19B07.RADOLAN.nc',
                        lock=False)
                    if grid_fc.clon.units == 'radian':
                        grid_fc.clon.data = np.rad2deg(grid_fc.clon.data)
                        grid_fc.clon.attrs['units'] = 'degree'
                    if grid_fc.clat.units == 'radian':
                        grid_fc.clat.data = np.rad2deg(grid_fc.clat.data)
                        grid_fc.clat.attrs['units'] = 'degree'

                    levels_fc = xr.open_dataset(
                        dir_data_in + '/grid/vgrd_R19B07.RADOLAN.nc',
                        lock=False)
                    levels_fc = levels_fc.isel(time=0)
                elif single_fc.ncells.size == 542040:
                    ncells = 542040  # ICON-D2
                    grid_fc = xr.open_dataset(
                        dir_data_in + '/grid/hgrd_R19B07.ICON-D2.nc',
                        lock=False)
                    if grid_fc.cell.size == single_fc.ncells.size:
                        if grid_fc.clon.units == 'radian':
                            grid_fc.clon.data = np.rad2deg(grid_fc.clon.data)
                            grid_fc.clon.attrs['units'] = 'degree'
                        if grid_fc.clat.units == 'radian':
                            grid_fc.clat.data = np.rad2deg(grid_fc.clat.data)

                        levels_fc = xr.open_dataset(
                            dir_data_in + '/grid/vgrd_R19B07.ICON-D2.nc',
                            lock=False)
                elif single_fc.ncells.size == 100122:
                    ncells = 100122  # ICON-SW LIFT
                    grid_fc = xr.open_dataset(
                        dir_data_in + '/grid/hgrd_R19B07.SW.nc', lock=False)
                    if grid_fc.cell.size == single_fc.ncells.size:
                        if grid_fc.clon.units == 'radian':
                            grid_fc.clon.data = np.rad2deg(grid_fc.clon.data)
                            grid_fc.clon.attrs['units'] = 'degree'
                        if grid_fc.clat.units == 'radian':
                            grid_fc.clat.data = np.rad2deg(grid_fc.clat.data)

                        levels_fc = xr.open_dataset(
                            dir_data_in + '/grid/vgrd_R19B07.SW.nc',
                            lock=False)
                else:
                    print('No forcast data, as differing ncells ' +
                          'were found in ' + dir_data_in + '/grid/')
                    print('continue')

                # if flipped:
                if levels_fc.height.size == levels_fc.height_2.size + 1:
                    levels_fc = levels_fc.rename(
                        {'height_2': 'height', 'height': 'height_2',
                         'height_2_bnds': 'height_bnds',
                         'height_bnds': 'height_2_bnds',
                         })

                levels_fc = levels_fc.transpose('height', 'ncells', ...)

                # shrink ICON variables
                vars_to_compute_t = []
                vars_to_compute_const = []
                vars_to_drop = []
                single_fc = single_fc.transpose('time', 'height', 'ncells',
                                                ...)
                for vv in single_fc.data_vars:
                    if single_fc[vv].dims == ('time', 'height', 'ncells'):
                        vars_to_compute_t.append(vv)
                    elif single_fc[vv].dims == ('height', 'ncells'):
                        vars_to_compute_const.append(vv)
                    else:
                        vars_to_drop.append(vv)
                        # print(vv + ' is not processed.')

                single_fc = single_fc.drop_vars(vars_to_drop)

                # reproject ICON into RADAR grid
                proj_wgs = osr.SpatialReference()
                proj_wgs.ImportFromEPSG(4326)
                mod_x, mod_y = wrl.georef.reproject(grid_fc["clon"].values,
                                                    grid_fc["clat"].values,
                                                    trg_crs=proj_aeqd,
                                                    src_crs=proj_wgs)

                # ICON source grid in shape=raveled(ncells_shrunken)
                src = np.vstack(
                    (np.repeat(mod_x[np.newaxis, :],
                               levels_fc.height.size, axis=0).ravel(),
                     np.repeat(mod_y[np.newaxis, :],
                               levels_fc.height.size, axis=0).ravel(),
                     levels_fc.z_mc.values.ravel())).T

                # calculate indices of ICON field that are nearest to RADAR
                mesh = pyinterp.RTree(ecef=True)
                stacked_single_fc = single_fc.isel(time=0).stack(
                    stacked=['height', 'ncells', ])
                data = stacked_single_fc[vars_to_compute_t[0]]
                indices_all_fc = np.arange(0, data.shape[0])
                mesh.packing(src, indices_all_fc)
                neighbors, indices_near_fc = mesh.value(trg, within=False, k=1)
                indices_near_fc = indices_near_fc.flatten()
                indices_near_fc = indices_near_fc.astype(int)

                shape_ear = vol_scan['lon'].shape

            # select fc variables
            stacked_single_fc = single_fc.isel(time=0).stack(
                stacked=['height', 'ncells', ])
            stacked_single_fc = stacked_single_fc.isel(stacked=indices_near_fc)

            # copy fc variables
            vol_scan['temp'][t_i, :, :, :] = \
                stacked_single_fc['temp'].values.reshape(shape_ear)
            try:
                vol_scan['pres'][t_i, :, :, :] = \
                    stacked_single_fc['pres'].values.reshape(shape_ear)
                vol_scan['qv'][t_i, :, :, :] = \
                    stacked_single_fc['qv'].values.reshape(shape_ear)
                vol_scan['qc'][t_i, :, :, :] = \
                    stacked_single_fc['qc'].values.reshape(shape_ear)
                vol_scan['qr'][t_i, :, :, :] = \
                    stacked_single_fc['qr'].values.reshape(shape_ear)
                vol_scan['qi'][t_i, :, :, :] = \
                    stacked_single_fc['qi'].values.reshape(shape_ear)
                vol_scan['qs'][t_i, :, :, :] = \
                    stacked_single_fc['qs'].values.reshape(shape_ear)
                vol_scan['qg'][t_i, :, :, :] = \
                    stacked_single_fc['qg'].values.reshape(shape_ear)
                vol_scan['qh'][t_i, :, :, :] = \
                    stacked_single_fc['qh'].values.reshape(shape_ear)
                vol_scan['qnc'][t_i, :, :, :] = \
                    stacked_single_fc['qnc'].values.reshape(shape_ear)
                vol_scan['qnr'][t_i, :, :, :] = \
                    stacked_single_fc['qnr'].values.reshape(shape_ear)
                vol_scan['qni'][t_i, :, :, :] = \
                    stacked_single_fc['qni'].values.reshape(shape_ear)
                vol_scan['qns'][t_i, :, :, :] = \
                    stacked_single_fc['qns'].values.reshape(shape_ear)
                vol_scan['qng'][t_i, :, :, :] = \
                    stacked_single_fc['qng'].values.reshape(shape_ear)
                vol_scan['qnh'][t_i, :, :, :] = \
                    stacked_single_fc['qnh'].values.reshape(shape_ear)
                vol_scan['u'][t_i, :, :, :] = \
                    stacked_single_fc['u'].values.reshape(shape_ear)
                vol_scan['v'][t_i, :, :, :] = \
                    stacked_single_fc['v'].values.reshape(shape_ear)
                vol_scan['w'][t_i, :, :, :] = \
                    stacked_single_fc['w'].values.reshape(shape_ear)
            except:
                print('some variable not in ICON fc')

            # single_fc done!
            single_fc.close()
            stacked_single_fc.close()

    if 'vol_scan' in locals():
        Path(dir_out).mkdir(parents=True, exist_ok=True)
        # transpose for panoply reasons (why do not know)
        vol_scan = vol_scan.transpose(
            'time', 'range', 'azimuth', 'elevation', )
        print('start saving')
        vol_scan.to_netcdf(dir_out + file_out, unlimited_dims='time')
        print('saved')
        vol_scan.close()
        print('    ! Case done now !      ')
        print(f"... which took "
              f"{(time_p.time() - current_time) / 60:.2f} min ..." 
              f" ({time_p.strftime('%d/%m %H:%M', time_p.gmtime())} UTC)"
              )
        print('___________________________')
        return dir_out + file_out
    else:
        if include_icon and not include_emv:
            print('   ! No ICON for case !    ')
            print('___________________________')
        elif not include_icon and include_emv:
            print('   ! No EMVO for case !    ')
            print('___________________________')
        else:
            print('   ! No data for case !    ')
            print('___________________________')


def create_8_vol_nc_of_day(day='20170725', da_run='ASS_2211',
                           icon_run='MAIN_2211.0',
                           icon_emvorado_run='MAIN_2211.0/EMVO_00000000.2',
                           spin_up_mm=30,
                           radar_locs=list(rad_dict().keys()),
                           dir_data_in=header.dir_data_mod,
                           dir_data_out=header.dir_data_vol,
                           overwrite_EMV=False,
                           overwrite_ICON=False,
                           naming2026=False,
                           ):
    """
    Create for day 8 synthetic volume scans from EMVORADO and ICON data.

     Args:
        day: day in >yyyymmdd<.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        radar_locs: list of strings  >RRR< naming the radars.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs.
        dir_data_out: directory for the output.
        overwrite_EMV=overwrite EMVORADO? True, False, <yyyy-mm-dd>,
        overwrite_ICON=overwrite ICON? True, False, <yyyy-mm-dd>,
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)

    Returns:
    """
    for radar_loc in radar_locs:
        time_start = pd.to_datetime(day, format="%Y%m%d")
        for i in range(4):
            time_end = time_start + pd.Timedelta('6h')
            print('________________________________________')
            print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
                  str(spin_up_mm) + '_spinup/')
            create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
                          time_end=time_end.strftime('%Y%m%d%H'),
                          spin_up_mm=spin_up_mm, da_run=da_run,
                          icon_run=icon_run,
                          icon_emvorado_run=icon_emvorado_run,
                          dir_data_in=dir_data_in,
                          dir_data_out=dir_data_out,
                          radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                          overwrite=overwrite_ICON,
                          include_icon=True, include_emv=False,
                          naming2026=naming2026)
            print('________________________________________')
            print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
                  str(spin_up_mm) + '_spinup/')
            create_vol_nc(time_start=time_start.strftime('%Y%m%d%H'),
                          time_end=time_end.strftime('%Y%m%d%H'),
                          spin_up_mm=spin_up_mm, da_run=da_run,
                          icon_run=icon_run,
                          icon_emvorado_run=icon_emvorado_run,
                          dir_data_in=dir_data_in,
                          dir_data_out=dir_data_out,
                          radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                          overwrite=overwrite_EMV,
                          include_icon=False, include_emv=True,
                          naming2026=naming2026)
            time_start = time_end


def create_8_vol_nc_of_day_cdo(day='20170725', da_run='ASS_2211',
                               icon_run='MAIN_2211.0',
                               icon_emvorado_run='MAIN_2211.0/EMVO_00000000.2',
                               spin_up_mm=30,
                               radar_locs=list(rad_dict().keys()),
                               dir_data_in=header.dir_data_mod,
                               dir_data_out=header.dir_data_vol,
                               overwrite_EMV=False,
                               overwrite_ICON=False,
                               naming2026=False,
                               ):
    """
    Create for day 8 synthetic volume scans from EMVORADO and ICON data.

     Args:
        day: day in >yyyymmdd<.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        radar_locs: list of strings  >RRR< naming the radars.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs.
        dir_data_out: directory for the output.
        overwrite_EMV=overwrite EMVORADO? True, False, <yyyy-mm-dd>,
        overwrite_ICON=overwrite ICON? True, False, <yyyy-mm-dd>,
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)

    Returns:
    """
    for radar_loc in radar_locs:
        time_start = pd.to_datetime(day, format="%Y%m%d")
        for i in range(4):
            time_end = time_start + pd.Timedelta('6h')
            print('________________________________________')
            print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
                  str(spin_up_mm) + '_spinup/')
            current_time = time_p.time()
            list_icon = []
            for ii in range(6):
                time_start_ii = time_start + pd.Timedelta(str(ii) + 'h')
                time_end_ii = time_start_ii + pd.Timedelta('1h')
                list_icon.append(create_vol_nc(
                    time_start=time_start_ii.strftime('%Y%m%d%H'),
                    time_end=time_end_ii.strftime('%Y%m%d%H'),
                    spin_up_mm=spin_up_mm, da_run=da_run,
                    icon_run=icon_run,
                    icon_emvorado_run=icon_emvorado_run,
                    dir_data_in=dir_data_in,
                    dir_data_out=dir_data_out,
                    radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                    overwrite=overwrite_ICON,
                    include_icon=True, include_emv=False,
                    naming2026=naming2026)
                )

            new = list_icon[0][:-15] + 'neu' + list_icon[0][-15:]
            print('cdo merge ' + ' '.join(list_icon) + ' ' + new)
            os.system('cdo merge ' + ' '.join(list_icon) + ' ' + new)
            print(f"... which took "
                  f"{(time_p.time() - current_time) / 60:.2f} min ..."
                  f" ({time_p.strftime('%d/%m %H:%M', time_p.gmtime())} UTC)"
                  )
            print('________________________________________')
            print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
                  str(spin_up_mm) + '_spinup/')
            list_emv = []
            for ii in range(6):
                time_start_ii = time_start + pd.Timedelta(str(ii) + 'h')
                time_end_ii = time_start_ii + pd.Timedelta(+'1h')
                list_emv.append(create_vol_nc(
                    time_start=time_start_ii.strftime('%Y%m%d%H'),
                    time_end=time_end_ii.strftime('%Y%m%d%H'),
                    spin_up_mm=spin_up_mm, da_run=da_run,
                    icon_run=icon_run,
                    icon_emvorado_run=icon_emvorado_run,
                    dir_data_in=dir_data_in,
                    dir_data_out=dir_data_out,
                    radar_loc=radar_loc, radar_id=rad_dict()[radar_loc],
                    overwrite=overwrite_EMV,
                    include_icon=False, include_emv=True,
                    naming2026=naming2026)
                )

            new = list_emv[0][:-15] + list_emv[0][-15:]
            os.system('cdo merge ' + ' '.join(list_emv) + ' ' + new)
            time_start = time_end


def create_8_vol_nc_of_day_paralell(day='20170725', da_run='ASS_2211',
                                    icon_run='MAIN_2211.0',
                                    icon_emvorado_run=
                                    'MAIN_2211.0/EMVO_00000000.2',
                                    spin_up_mm=30,
                                    radar_locs=list(rad_dict().keys()),
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol,
                                    naming2026=False,
                                    ):
    """
    Create for day 8 synthetic volume scans from EMVORADO and ICON data.
    Do it parallel for the 4 times 6h blocks.

     Args:
        day: day in >yyyymmdd<.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        radar_locs: list of strings  >RRR< naming the radars.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs.
        dir_data_out: directory for the output.
        method: 'Nearest' or 'Linear' for method to interpolate mod fields to
                rad grid.
        naming2026: True for new naming convention of 2026 (PPI instead
                       of volscan)

    Returns:
    """
    for radar_loc in radar_locs:
        time_start = pd.to_datetime(day, format="%Y%m%d")
        dti = pd.date_range(time_start, time_start + pd.Timedelta('24h'),
                            freq="6h", inclusive='both').strftime('%Y%m%d%H')
        times_start = list(dti[:4])
        times_end = list(dti[1:])
        pool = mp.Pool(min(1, mp.cpu_count()))
        # pool = mp.Pool(min(2, mp.cpu_count()))
        # pool = mp.Pool(min(4, mp.cpu_count()))
        print('________________________________________')
        print(day + '/' + da_run + '/' + icon_run + '/ICONdata/' +
              str(spin_up_mm) + '_spinup/')
        pool.starmap(create_vol_nc,
                     zip(times_start,
                         times_end,
                         np.repeat(dir_data_in, 4),
                         np.repeat(dir_data_out, 4),
                         np.repeat(radar_loc, 4),
                         np.repeat(rad_dict()[radar_loc], 4),
                         np.repeat(spin_up_mm, 4),
                         np.repeat(da_run, 4),
                         np.repeat(icon_run, 4),
                         np.repeat(icon_emvorado_run, 4),
                         np.repeat(False, 4),
                         np.repeat(True, 4),
                         np.repeat(False, 4),
                         np.repeat('ICONdata', 4),
                         np.repeat(naming2026, 4))
                     )
        print('________________________________________')
        print(day + '/' + da_run + '/' + icon_emvorado_run + '/' +
              str(spin_up_mm) + '_spinup/')
        pool.starmap(create_vol_nc,
                     zip(times_start,
                         times_end,
                         np.repeat(dir_data_in, 4),
                         np.repeat(dir_data_out, 4),
                         np.repeat(radar_loc, 4),
                         np.repeat(rad_dict()[radar_loc], 4),
                         np.repeat(spin_up_mm, 4),
                         np.repeat(da_run, 4),
                         np.repeat(icon_run, 4),
                         np.repeat(icon_emvorado_run, 4),
                         np.repeat(False, 4),
                         np.repeat(False, 4),
                         np.repeat(True, 4),
                         np.repeat('ICONdata', 4),
                         np.repeat(naming2026, 4))
                     )
