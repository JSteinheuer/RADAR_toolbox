#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.10.23                                                 #
# SYN_RADAR_1_CREATE_VOLUME_SCAN.py                                           #
#                                                                             #
# Functions to calculate synthetic volume scans from EMVORADO and ICON.       #
# --------------------------------------------------------------------------- #

import xarray as xr
import numpy as np
import wradlib as wrl
import pandas as pd
import os
from pathlib import Path
import wradlib.georef as georef
from osgeo import osr
import glob
import multiprocessing as mp
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

    radar_dict = {'OFT': '010629',
                  'PRO': '010392',
                  'NEU': '010557',
                  'FBG': '010908',
                  'UMD': '010356',
                  'ISN': '010873',
                  'ROS': '010169',
                  'MEM': '010950',
                  'EIS': '010780',
                  'ESS': '010410',
                  'HNR': '010339',
                  'FLD': '010440',
                  'NHB': '010605',
                  'DRS': '010488',
                  'TUR': '010832',
                  'BOO': '010132',
                  }
    if xband_res is not None:
        if int(xband_res) < 1000:
            radar_dict['BOX'] = '066' + str(int(xband_res)).zfill(3)
            radar_dict['JUX'] = '067' + str(int(xband_res)).zfill(3)

    return radar_dict


def get_path_syn_volume(date, time, spin_up_mm,
                        radar_id, dir_data, da_run='',
                        icon_run='', icon_emvorado_run=''):
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
    dir_of_fc = dir_data + dti_fc.strftime('%Y%m%d') + '/' + \
                da_run + '/' + icon_run + '/ICONdata/' + \
                dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_fc = glob.glob(dir_of_fc + '*')
    if len(dir_of_fc) > 0:
        dir_of_fc = dir_of_fc[0] + '/'
    else:  # previous day folder (if DA at 21:00)?
        dir_of_fc = dir_data + \
                    dti_3hda.strftime('%Y%m%d') + '/' + \
                    da_run + '/' + icon_run + '/ICONdata/' + '/' + \
                    dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                    'main0' + str(int(dti_da.strftime('%H')) % 3) + '00' + \
                    icon_run + '/'
        dir_of_fc = glob.glob(dir_of_fc + '*')
        if len(dir_of_fc) > 0:
            dir_of_fc = dir_of_fc[0] + '/'

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
                dir_of_fc = ''

    # volume scan
    dir_of_vol = dir_data + dti_fc.strftime('%Y%m%d') + '/' + \
                 da_run + '/' + icon_emvorado_run + '/' + \
                 dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                 'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
    dir_of_vol = glob.glob(dir_of_vol + '*')
    if len(dir_of_vol) > 0:
        dir_of_vol = dir_of_vol[0] + '/det' + '/cdfin_volscan_' + \
                     dti_da.strftime('%Y%m%d%H%M') + '-' + \
                     (dti_da + pd.Timedelta('2h')).strftime('%Y%m%d%H%M') + '/'

    else:  # previous day folder (if DA at 21:00)?
        dir_of_vol = dir_data + dti_3hda.strftime('%Y%m%d') + '/' + \
                     da_run + '/' + icon_emvorado_run + '/' + \
                     dti_3hda.strftime('%Y%m%d%H%M%S') + '/' + \
                     'main0' + str(int(dti_da.strftime('%H')) % 3) + '00'
        dir_of_vol = glob.glob(dir_of_vol + '*')
        if len(dir_of_vol) > 0:
            dir_of_vol = dir_of_vol[0] + '/det' + '/cdfin_volscan_' + \
                         dti_da.strftime('%Y%m%d%H%M') + '-' + \
                         (dti_da + pd.Timedelta('2h')).strftime(
                             '%Y%m%d%H%M') + '/'

    if not dir_of_vol:
        dir_of_vol = ''
        file_vol = ''
    else:
        file_vol = 'cdfin_allsim_id-' + radar_id + '_' + \
                   dti_fc.strftime('%Y%m%d%H%M') + '_' + \
                   dti_fc.strftime('%Y%m%d%H%M') + '_volscan.nc'
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

    proj_wgs84 = georef.epsg_to_osr(4326)
    cent_coords = georef.spherical_to_centroids(r, az, el, sitecoords,
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
    outer_x = 0.2 * (rad_x.max() - rad_x.min())
    outer_y = 0.2 * (rad_y.max() - rad_y.min())
    lower_z = 50
    upper_z = 1000
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


def create_vol_nc(time_start='2017072500', time_end='2017072506',
                  dir_data_in=header.dir_data_mod,
                  dir_data_out=header.dir_data_vol,
                  radar_loc='PRO', radar_id='010392', spin_up_mm=30,
                  da_run='', icon_run='', icon_emvorado_run='',
                  overwrite=False, include_icon=True, include_emv=True):
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
    if os.path.isfile(dir_out + file_out) and not overwrite:
        print(radar_loc, '   -   ', time_start, '-', time_end[-2:])
        print(file_out + ' exists;\n' +
              ' ... set: > overwrite = True < for recalculation')
        print('___________________________')
        return

    ncells = 0
    for t_i in range(len(dti)):
        date = dti[t_i].strftime('%Y%m%d')
        time = dti[t_i].strftime('%H%M%S')
        print(radar_loc, ' - ', date, ' - ', time)
        dir_of_fc, file_fc, dir_of_vol, file_vol = \
            get_path_syn_volume(date, time, spin_up_mm, radar_id, dir_data_in,
                                da_run, icon_run, icon_emvorado_run).values()

        if 'vol_scan' not in locals() and not dir_of_fc + dir_of_vol == '':
            if not os.path.isfile(dir_of_vol + file_vol):
                print('No EMVORADO input for time')
                continue

            single_scan = xr.open_dataset(dir_of_vol + file_vol)
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
                vol_scan['alt'][:, sweep_i] = alt[0, :]

            # grid for searching later the closest ICON cells
            rad_lon = vol_scan['lon'].data.flatten()
            rad_lat = vol_scan['lat'].data.flatten()
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

            single_scan = xr.open_dataset(dir_of_vol + file_vol)
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

            single_fc = xr.open_dataset(dir_of_fc + file_fc)

            # check for ICON grid
            if not single_fc.ncells.size == ncells:
                grid_fc = xr.open_dataset(
                    dir_data_in + '/grid/hgrd_R19B07.RADOLAN.nc')
                if grid_fc.cell.size == single_fc.ncells.size:
                    fc_lon = np.rad2deg(grid_fc.clon.data)
                    fc_lat = np.rad2deg(grid_fc.clat.data)
                    ncells = grid_fc.cell.size
                    grid_fc.close()
                    levels_fc = xr.open_dataset(
                        dir_data_in + '/grid/vgrd_R19B07.RADOLAN.nc')
                    fc_alt = levels_fc.z_mc.data[0, :, :]
                    levels_fc.close()
                else:
                    grid_fc.close()
                    grid_fc = xr.open_dataset(
                        dir_data_in + '/grid/hgrd_R19B07.ICON-D2.nc')
                    if grid_fc.cell.size == single_fc.ncells.size:
                        fc_lon = np.rad2deg(grid_fc.clon.data)
                        fc_lat = np.rad2deg(grid_fc.clat.data)
                        ncells = grid_fc.cell.size
                        grid_fc.close()
                        levels_fc = xr.open_dataset(
                            dir_data_in + '/grid/vgrd_R19B07.ICON-D2.nc')
                        fc_alt = levels_fc.z_mc.data[:, :]
                        levels_fc.close()
                    else:
                        grid_fc.close()
                        print('No forcast data, as differing ncells ' +
                              'were found in ' + dir_data_in + '/grid/')
                        continue

                func_ipol, mask = ipol_fc_to_radgrid(
                    # lon with shape (65,258775):
                    np.repeat(fc_lon[np.newaxis, :], fc_alt.shape[0], axis=0),
                    # lat with shape (65,258775)
                    np.repeat(fc_lat[np.newaxis, :], fc_alt.shape[0], axis=0),
                    fc_alt, rad_lon, rad_lat, rad_alt
                )
                shape_rae = vol_scan['lon'].shape

            # copy fc variables
            vol_scan['temp'][t_i, :, :, :] = \
                func_ipol(single_fc['temp'].data[0, :, :][mask]).reshape(
                    shape_rae)
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

            # single_fc done!
            single_fc.close()

    if 'vol_scan' in locals():
        Path(dir_out).mkdir(parents=True, exist_ok=True)
        vol_scan.to_netcdf(dir_out + file_out, unlimited_dims='time')
        vol_scan.close()
        print('    ! Case done now !      ')
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
                          include_icon=True, include_emv=False)
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
                          include_icon=False, include_emv=True)
            time_start = time_end


def create_8_vol_nc_of_day_paralell(day='20170725', da_run='ASS_2211',
                                    icon_run='MAIN_2211.0',
                                    icon_emvorado_run=
                                    'MAIN_2211.0/EMVO_00000000.2',
                                    spin_up_mm=30,
                                    radar_locs=list(rad_dict().keys()),
                                    dir_data_in=header.dir_data_mod,
                                    dir_data_out=header.dir_data_vol,
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

    Returns:
    """
    for radar_loc in radar_locs:
        time_start = pd.to_datetime(day, format="%Y%m%d")
        dti = pd.date_range(time_start, time_start + pd.Timedelta('24h'),
                            freq="6h", inclusive='both').strftime('%Y%m%d%H')
        times_start = list(dti[:4])
        times_end = list(dti[1:])
        pool = mp.Pool(min(2, mp.cpu_count()))
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
                         np.repeat(False, 4))
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
                         np.repeat(True, 4))
                     )
