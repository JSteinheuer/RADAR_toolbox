#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 17.06.24                                                 #
#                                                                             #
# separate script later to:                                                   #
# ?PROCESS_RADAR.py                                                           #
# ?process_RADAR_QVP_from_volume_scan.py                                      #
#                                                                             #
# Functions to calculate QVPs from given volume scans.                        #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# NOT WORKING YET                                                             #
# --------------------------------------------------------------------------- #

import datatree as dttree
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import os
import scipy
import wradlib as wrl
import xarray as xr
import datetime

xr.set_options(keep_attrs=True)
import HEADER_RADAR_toolbox as header


def entropy_of_vars(pol_vars, dim='azimuth', n_lowest=30):
    """
    Calculate the information Entropy, to estimate the homogeneity
    from a sector- or whole 360° PPi for each timestep

    Args:
      pol_vars: list of xarrays with polarimetric variables,
                e.g.: [10**(0.1*Z_H),10**(0.1*Z_DR), RHO_HV, K_DP]
      dim: dimension on which entropy is calculated (azimuth for QVPs)
      n_lowest: least nr. of non-nan values for returning non-nan entropy.

    Returns:
      Output: List of xarrays Entropy and the minimum of all.

     """

    def calc_entropy(var, dim, n_valid_values):
        minimum = np.nanmin(var.data)
        if minimum <= 0:
            var = xr.where(var.data <= 0, np.nan, var)
            print('variable has values <= 0 which are set to '
                  'nan')

        var_normed = var / var.sum(dim, skipna=True)
        values = -((var_normed * np.log10(var_normed)).sum(dim)
                   ) / np.log10(var_normed[dim].size)
        values = values.where(var_normed.count(dim=dim) >= n_valid_values)
        return values

    entropy = calc_entropy(pol_vars[0], dim=dim, n_valid_values=n_lowest)
    entropy.name = 'entropy'
    entropy = entropy.expand_dims(dim={'n_all': len(pol_vars) + 1}).copy()
    for i in range(1, len(pol_vars)):
        entropy[i, :] = calc_entropy(pol_vars[i], dim=dim,
                                     n_valid_values=n_lowest)

    entropy[-1, :] = entropy[:-1, :, :].min(dim='n_all', skipna=False)
    return entropy


def normalise(da, dim):
    damin = da.min(dim=dim, skipna=True, keep_attrs=True)
    damax = da.max(dim=dim, skipna=True, keep_attrs=True)
    da = (da - damin) / (damax - damin)
    return da


def sobel(da, dim):
    axis = da.get_axis_num(dim)
    func = lambda x, y: scipy.ndimage.sobel(x, y)
    return xr.apply_ufunc(func, da, axis, dask='parallelized')


def ml_normalise(ds, moments=dict(DBZH=(10., 60.),
                                  RHOHV=(0.65, 1.),
                                  PHIDP=(-0.5, 360)), dim='height'):
    assign = {}
    for m, span in moments.items():
        v = ds[m]
        v = v.where((v >= span[0]) & (v <= span[1]))
        v = v.pipe(normalise, dim=dim)
        assign.update({f'{m}_norm': v})
    return xr.Dataset(assign)


def ml_clean(ds, zh, dim='height'):
    # removing profiles with too few data (>93% of array is nan)
    good = ds[zh + "_norm"].count(dim=dim) / ds[zh + "_norm"].sizes[dim] * 100
    return ds.where(good > 7)


def ml_combine(ds, zh, rho):
    comb = (1 - ds[rho + "_norm"]) * ds[zh + "_norm"]
    return comb


def ml_gradient(ds, dim='height', ywin=None):
    assign = {}
    for m in ds.data_vars:
        # step 3 sobel filter
        dy = ds[m].pipe(sobel, dim)
        # step 3 smoothing
        # currently smoothing is only in y-direction (height)
        dy = dy.rolling({dim: ywin}, center=True).mean()
        assign.update({f'{m}_dy': dy})
    return ds.assign(assign)


def ml_noise_reduction(ds, thres=0.02):
    assign = {}
    for m in ds.data_vars:
        if '_dy' in m:
            dy = ds[m].where((ds[m] > thres) | (ds[m] < -thres))
            assign.update({f"{m}": dy})
    return ds.assign(assign)


def ml_interpolate2(ds, ppi, moment, **kwargs):
    # min_val = ds_top.argmin(dim='height', skipna=True)
    # select range and interpolate
    nan = np.isnan(ppi[moment])
    mom = ppi[moment].where((ppi.height < ds.mlh_bottom) |
                            (ppi.height > ds.mlh_top)).interpolate_na(
        dim='range', **kwargs)
    mom = xr.where(nan, np.nan, mom)
    mom = xr.where(np.isnan(mom), ppi[moment], mom)
    mom.attrs = ppi[moment].attrs
    mom.attrs['comments'] = mom.attrs['comments'] + '; Melting layer corrected'
    return mom


def ml_height_bottom_new(ds, moment='comb_dy', dim='height',
                         skipna=True):  # yes
    hgt = ds[moment].idxmax(dim=dim, skipna=skipna, fill_value=np.nan).load()

    ds = ds.assign(dict(
        mlh_bottom=hgt))
    return ds


def ml_height_top_new(ds, moment='comb_dy', dim='height', skipna=True):
    hgt = ds[moment].idxmin(dim=dim, skipna=skipna, fill_value=np.nan).load()
    ds = ds.assign(dict(mlh_top=hgt))
    return ds


def melting_layer_qvp_X_new(ds, moments=dict(ZH_AC=(10., 60.),
                                             RHOHV_NC2P=(0.65, 1.),
                                             PHI_NC=(-0.5, 360)),
                            dim='height',
                            thres=0.02, xwin=5, ywin=5, fmlh=0.3,
                            all_data=False):  # yes
    # step 0
    zh = [k for k in moments if ("zh" in k.lower()) or
          ("th" in k.lower()) or ("zh" in k.lower())][0]
    rho = [k for k in moments if ("rho" in k.lower()) or
           ("rho" in k.lower()) or ("cc" in k.lower())][0]
    phi = [k for k in moments if "phi" in k.lower()][0]

    # step 1
    ds0 = ml_normalise(ds, moments=moments, dim=dim)

    # step 1a: removing profiles with too few data (>93% of array is nan)
    good = ml_clean(ds0, zh, dim=dim)

    # step 2
    comb = ml_combine(good, zh, rho)
    ds0 = ds0.assign(dict(comb=comb))

    # step 3 (and 8)
    ds0 = ml_gradient(ds0, dim=dim, ywin=ywin)

    # step 4
    ds1 = ml_noise_reduction(ds0, thres=thres)

    # step 5
    ds2 = ml_height_bottom_new(ds1, dim=dim)
    ds2 = ml_height_top_new(ds2, dim=dim)

    # step 6
    med_mlh_bot = ds2.mlh_bottom.rolling(
        time=xwin, min_periods=xwin // 2, center=True).median(skipna=True)
    med_mlh_top = ds2.mlh_top.rolling(
        time=xwin, min_periods=xwin // 2, center=True).median(skipna=True)

    # step 7: step 5 again
    above = (1 + fmlh) * med_mlh_top
    below = (1 - fmlh) * med_mlh_bot
    ds3 = ds1.where((ds1.height >= below) & (ds1.height <= above), drop=False)
    ds3 = ml_height_bottom_new(ds3, dim=dim)
    ds3 = ml_height_top_new(ds3, dim=dim)
    # ds3.mlh_bottom and ds3.ml_bottom derived according Wolfensberger et al.

    # step 8: already done at step 3

    # step 9: ml top refinement: cut below ml_top
    ds4 = ds0.where(
        (ds0.height > ds3.mlh_top))  # This! # & (ds0.height < above))
    # cut above first local maxima (via differentiation & difference of signs)
    ds4 = ds4.chunk('auto')
    p4 = np.sign(ds4[zh + "_norm_dy"].differentiate(dim)).diff(dim).idxmin(
        dim, skipna=True)
    ds5 = ds4.where(ds4.height < p4)
    ds5 = ml_height_top_new(ds5, moment=phi + "_norm", dim=dim)
    # ds5.mlh_top, ds5.ml_top, derived according Wolfensberger et al.
    # but used PHIDP_norm instead of DBZH_norm

    # step 9b: ml bottom refinement
    ds6 = ds0.where((ds0.height < ds3.mlh_bottom) &
                    (ds0.height > below))  # This, i.e., both conditions!
    ds6 = ds6.chunk('auto')
    # cut below first local maxima (via differentiation & difference of signs)
    p4a = np.sign(ds6[phi + "_norm_dy"].differentiate(dim)).diff(dim).sortby(
        dim, ascending=False).idxmin(dim, skipna=True)
    ds7 = ds6.where(ds6.height > p4a)

    # step 10
    hgt = ds7[phi + "_norm"].idxmin(dim=dim, skipna=True, fill_value=np.nan)

    ds7 = ds7.assign(dict(mlh_bottom=hgt))
    # ds7.mlh_bottom ds7.ml_bottom, derived similar to Wolfensberger et al.
    # but for bottom uses PHIDP_norm

    # step 11: return
    ds = ds.assign(dict(mlh_top=ds5.mlh_top,
                        mlh_bottom=ds7.mlh_bottom, ), )
    ds['mlh_top'] = ds['mlh_top'].assign_attrs(dict(
        standard_name='melting layer height (top)',
        units=ds.height.units,
        comments='Height above radar. Derived from Wolfensberger et al. 2015. '
                 'Adapted to synthetic data.'))
    ds['mlh_bottom'] = ds['mlh_bottom'].assign_attrs(dict(
        standard_name='melting layer height (bottom)',
        units=ds.height.units,
        comments='Height above radar. Derived from Wolfensberger et al. 2015. '
                 'Adapted to synthetic data, but without much value '
                 '(use mlh_top).'))
    if all_data is True:
        ds = ds.assign({"comb": ds0.comb,
                        rho + "_norm": ds0[rho + "_norm"],
                        zh + "_norm": ds0[zh + "_norm"],
                        phi + "_norm": ds0[phi + "_norm"],
                        rho + "_norm_dy": ds0[rho + "_norm_dy"],
                        zh + "_norm_dy": ds0[zh + "_norm_dy"],
                        phi + "_norm_dy": ds0[phi + "_norm_dy"],
                        })

    return ds



date = '20170725'
elevation_deg = 12
location = 'pro'
# overwrite = False
overwrite = True
lowest_rhv = 0.7  # None#0.7
lowest_kdp = None
highest_kdp = None
n_lowest = 30
lowest_zh = 0  # None
highest_zh = None
lowest_snr = 10

# qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg, location=location,
#                     overwrite=overwrite, lowest_rhv=lowest_rhv,
#                     lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
#                     n_lowest=n_lowest,
#                     lowest_zh=lowest_zh, highest_zh=highest_zh)

DATES = [
    # "20210604",  # case01
    # "20210620", "20210621",  # case02
    # "20210628", "20210629",  # case03
    # "20220519", "20220520",  # case04
    # "20220623", "20220624", "20220625",  # case05
    # "20220626", "20220627", "20220628",  # case06+07
    # "20220630", "20220701",  # case08
    # "20210713",  # case09
    # "20210714",  # case09
    # "20221222",  # case10
    # "20170719",  # caseX -> old OP HM 1 case
    "20170725",  # caseX -> old OP HM 1 case
    # "20181223",  "20181224",  # case -> PRISTINE
]
LOCATIONS = [
    # 'ess',
    # 'asb', 'boo', 'drs', 'eis', 'fbg',
    # 'fld', 'hnr', 'isn', 'mem', 'neu', 'nhb',
    # 'oft',
    'pro',
    # 'ros', 'tur',
    # 'umd',
]
for location in LOCATIONS:
    for date in DATES:
        for elevation_deg in [12]:
            print(location + ' ' + date + '' + str(elevation_deg))
            # qvp_from_radar_PPIs(date=date, elevation_deg=elevation_deg,
            #                     location=location,
            #                     overwrite=overwrite, lowest_rhv=lowest_rhv,
            #                     lowest_snr=lowest_snr,
            #                     lowest_kdp=lowest_kdp, highest_kdp=highest_kdp,
            #                     lowest_zh=lowest_zh, highest_zh=highest_zh,
            #                     n_lowest=n_lowest)


dir_data_obs=header.dir_data_obs
dir_data_out=header.dir_obs_qvp

# def qvp_from_radar_PPIs(date='20170725', elevation_deg=12, location='pro',
#                         dir_data_obs=header.dir_data_obs,
#                         dir_data_out=header.dir_obs_qvp,
#                         overwrite=False,
#                         lowest_rhv=0.7, lowest_snr=10,
#                         lowest_kdp=None, highest_kdp=None,
#                         lowest_zh=0, highest_zh=None,
#                         n_lowest=30):
#     """
#     Create QVP for one day out of one elevation.
#
#      Args:
#         date: day in >yyyymmdd<.
#         elevation_deg: elevation angle to be choosen from Syn_vol for QVP
#                        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0].
#         location: string >rrr< naming the radar.
#         dir_data_obs: directory with folder >yyyymmdd< of the day outputs.
#         dir_data_out: directory for the output.
#         overwrite: If True, then process only if output is not existing.
#         lowest_rhv: If not None, it is the lowest acceptable value not
#                     filtered.
#         lowest_snr: If not None, it is the lowest acceptable value not
#                     filtered.
#         lowest_kdp: If not None, it is the lowest acceptable value not
#                     filtered.
#         highest_kdp: If not None, it is the highest acceptable value not
#                      filtered.
#         lowest_zh: If not None, it is the lowest acceptable value not
#                    filtered.
#         highest_zh: If not None, it is the highest acceptable value not
#                     filtered.
#         n_lowest: lowest nr of valid values for returning a non nan value per
#                   360 azimuths.
#
#     Returns:
#     """
# initialization
year = date[0:4]
mon = date[4:6]
day = date[6:8]
mode = 'vol'
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
path_in = "/".join([dir_data_obs + '*',
                    year, year + '-' + mon,
                    year + '-' + mon + '-' + day,
                    location, mode + '*', sweep,
                    'ras*_polmoms_nc_*'])
files = sorted(glob.glob(path_in))
if not files:
    print('no input data *_allmoms_*')
    # return
else:
    path_in = files[0]

# load
data = dttree.open_datatree(path_in)[
    'sweep_' + str(int(sweep))].to_dataset().chunk(-1)
data = data.transpose('time', 'azimuth', 'range')
time_start = date + '0000'
time_end = date + '2355'
dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
# TODO
dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
filter_comments = ''
path_qvp = path_in.replace(dir_data_obs, dir_data_out).replace('-vol',
                                                               '-qvp')
# TODO
path_qvp = path_in.replace(dir_data_obs, dir_data_out).replace(
    '-vol','-qvp'+datetime.datetime.now().strftime("%Y%m%d%H%M"))
if os.path.isfile(path_qvp) and not overwrite:
    print(path_qvp + ' already exists;\n' +
          'set >overwrite=True< for recalculation')
    # return

if lowest_rhv is not None:
    data['RHOHV_NC2P'] = \
        data['RHOHV_NC2P'].where(data['RHOHV_NC2P'] >= lowest_rhv)
    filter_comments = filter_comments + 'RHOHV_NC2P >= ' + \
                      str(lowest_rhv) + '. '
if lowest_kdp is not None:
    data['KDP_NC'] = \
        data['KDP_NC'].where(data['KDP_NC'] >= lowest_kdp)
    filter_comments = filter_comments + 'KDP_NC >= ' + \
                      str(lowest_kdp) + '. '
if highest_kdp is not None:
    data['KDP_NC'] = \
        data['KDP_NC'].where(data['KDP_NC'] <= highest_kdp)
    filter_comments = filter_comments + 'KDP_NC <= ' + \
                      str(highest_kdp) + ' . '
if lowest_zh is not None:
    data['ZH_AC'] = \
        data['ZH_AC'].where(data['ZH_AC'] >= lowest_zh)
    filter_comments = filter_comments + 'ZH_AC >= ' + \
                      str(lowest_zh) + ' . '
if highest_zh is not None:
    data['ZH_AC'] = \
        data['ZH_AC'].where(data['ZH_AC'] <= highest_zh)
    filter_comments = filter_comments + 'ZH_AC <= ' + \
                      str(highest_zh) + ' . '
if lowest_snr is not None:
    data['SNRH'] = \
        data['SNRH'].where(data['SNRH'] >= lowest_snr)
    filter_comments = filter_comments + 'SNR >= ' + \
                      str(lowest_snr) + ' . '

mask = ~ np.isnan(data['ZH_AC'])
for var in list(data.data_vars):
    if var not in ['KDP_NC', 'PHI_NC']:
        if data[var].dims == ('time', 'azimuth', 'range'):
            print(var)
            mask = mask * (~ np.isnan(data[var]))

for var in list(data.data_vars):
    if data[var].dims == ('time', 'azimuth', 'range'):
        data[var] = data[var].where(mask)

data = data.assign_coords({"height": ('range', data['range'].data *
                                      np.sin(elevation_deg * np.pi / 180.)
                                      / 1000.)})
data['height'] = data['height'].assign_attrs(dict(
    standard_name='height above radar',
    units='km'))
# ------------------------------------------------------------------- #
# QVP: KDP Melting layer correction                                   #
# ------------------------------------------------------------------- #
# Wolfensberger
qvp_first_nc = data.copy()
for var in ['ZH_AC', 'PHI_NC', 'RHOHV_NC2P']:
    count_values = qvp_first_nc[var].count(dim='azimuth')
    qvp_first_nc[var] = qvp_first_nc[var].median(
        dim='azimuth', skipna=True, keep_attrs=True)
    qvp_first_nc[var] = qvp_first_nc[var].where(
        count_values >= n_lowest)

filter_comments = filter_comments + 'QVP needs at least ' + \
                  str(n_lowest) + ' of 360 valid values per median ' \
                                  'calculation on azimuth. '
qvp_first_nc = qvp_first_nc.swap_dims({'range': 'height'})
qvp_first_nc = melting_layer_qvp_X_new(ds=qvp_first_nc, all_data=False)
data['mlh_top'] = qvp_first_nc.mlh_top
data['mlh_bottom'] = qvp_first_nc.mlh_bottom
data['PHI_NC_ml_corrected'] = \
    ml_interpolate2(qvp_first_nc, data, 'PHI_NC', method='linear')
data = data.transpose('time', 'azimuth', 'range')

# ----------------------------------------------------------------------- #
# NOT WORKING YET                                                         #
# ----------------------------------------------------------------------- #
kdp_new = wrl.dp.kdp_from_phidp(data['PHI_NC_ml_corrected'].data,  # TODO!?!?!?!?
                                winlen=13)#,  # cband  # TODO!?!?!?!?
                                # winlen=31, # xband ?!  # TODO!?!?!?!?
                               # in_periods=3)  # TODO!?!?!?!?
kdp_new[np.where(np.isnan(data['KDP_NC'].data))] = np.nan  # TODO!?!?!?!?
data['KDP_ml_corrected'] = (['time', 'azimuth', 'range'],  # TODO!?!?!?!?
                                 kdp_new, dict(
    standard_name=data['KDP_NC'].standard_name,
    units=data['KDP_NC'].units,
    comments=data['KDP_NC'].comments + '; Melting layer corrected'))
# Giagrande
cut_ml = qvp_first_nc.where(qvp_first_nc.height < qvp_first_nc.mlh_top)
cut_ml = cut_ml.where(qvp_first_nc.height > qvp_first_nc.mlh_bottom)
min_height_ML = cut_ml.rhvsim.idxmin(dim="height")
new_cut_below_min_ML = \
    qvp_first_nc.where(qvp_first_nc.height > min_height_ML)
new_cut_above_min_ML = \
    qvp_first_nc.where(qvp_first_nc.height < min_height_ML)
new_cut_below_min_ML_filter = new_cut_below_min_ML.rhvsim.where(
    (new_cut_below_min_ML.rhvsim >= 0.97) &
    (new_cut_below_min_ML.rhvsim <= 1))
new_cut_above_min_ML_filter = new_cut_above_min_ML.rhvsim.where(
    (new_cut_above_min_ML.rhvsim >= 0.97) &
    (new_cut_above_min_ML.rhvsim <= 1))
filter_comments = filter_comments + \
                  'ML refinement from Giagrande with rhvsim >= 0.97. '
panda_below_min = new_cut_below_min_ML_filter.to_pandas().transpose()
first_valid_height_ml = pd.DataFrame(panda_below_min).apply(
    pd.Series.first_valid_index)
first_valid_height_ml = first_valid_height_ml.to_xarray()
first_valid_height_ml = xr.where(
    first_valid_height_ml == None, data['mlh_top'].data,
    first_valid_height_ml)
# ML BOTTOM Giagrande refinement
panda_above_min = new_cut_above_min_ML_filter.to_pandas().transpose()
last_valid_height_ml = pd.DataFrame(panda_above_min).apply(
    pd.Series.last_valid_index)
last_valid_height_ml = last_valid_height_ml.to_xarray()
last_valid_height_ml = xr.where(last_valid_height_ml.data == None,
                                data['mlh_bottom'].data,
                                last_valid_height_ml)
data['mlh_top_gia'] = (["time"], first_valid_height_ml.data, dict(
    standard_name=data['mlh_top'].standard_name + ' refined',
    units=data['mlh_top'].units,
    comments=data['mlh_top'].comments + ' Refined with Giagrande.'))
data['mlh_bottom_gia'] = (["time"], last_valid_height_ml.data, dict(
    standard_name=data['mlh_bottom'].standard_name + ' refined',
    units=data['mlh_bottom'].units,
    comments=data['mlh_bottom'].comments + ' Refined with Giagrande.'))

# ------------------------------------------------------------------- #
# QVP: EMV                                                            #
# ------------------------------------------------------------------- #
pol_vars = [10 ** (0.1 * data.ZH_AC), 10 ** (0.1 * data.ZDR_AC_OC),
            data.RHOHV_NC2P, data.KDP_NC]  # TODO: ML correction
qvp_entropy = entropy_of_vars(pol_vars=pol_vars, dim='azimuth',
                              n_lowest=n_lowest)
for var in list(data.data_vars):
    if 'azimuth' in data[var].dims:
        count_values = data[var].count(dim='azimuth')
        data[var] = data[var].median(dim='azimuth',
                                     skipna=True, keep_attrs=True)
        data[var] = data[var].where(count_values >= n_lowest)

# TODO: Ice mycrophysical retrieval
# search for lambda
path_in_any = '/'.join(path_in.split('/')[:-1]) + '/*any*'
files_any = sorted(glob.glob(path_in_any))
if not files_any:
    path_in_any = "/".join([header.dir_data_obs_realpep + '' +
                            year, year + '-' + mon,
                            year + '-' + mon + '-' + day,
                            location, mode + '*', sweep, 'ras*'])
    files_any = sorted(glob.glob(path_in_any))

if not files_any:
    print('nothing found -> lamb=50')
    lamb = 50
else:
    path_in_any = files_any[0]
    try:
        lamb = dttree.open_datatree(  # we need it in mm, i.e. *10
            path_in_any)['how'].attrs['wavelength'] * 10
    except:
        print('nothing in *any* found -> lambda=5')
        lamb = 50

zh_lin = 10 ** (0.1 * data.ZH_AC)
zdr_lin = 10 ** (0.1 * data.ZDR_AC_OC)
# zdp_lin = zh_lin * (1 - zdr_lin ** (-1))
# data = data.where(data.KDP_NC > 0.01)  # TODO: really?

# zh=np.array([20,20,40,40])
# zh_lin = 10 ** (0.1 * zh)
# kdp=np.array([0.05,0.05,0.05,0.05])
# zdr=np.array([.3,.5,.3,5])
# zdr_lin = 10 ** (0.1 * zdr)
# iwc = xr.where(
#     zdr < 0.4,
#     0.033 * (kdp * lamb) ** 0.67 * zh_lin ** 0.33,
#     0.004 * kdp * lamb / (1 - zdr_lin ** (-1))
# )
# iwcb = xr.where(
#     zdr < 0.4,
#     0.033 * (kdp * lamb) ** 0.66 * zh_lin ** 0.28,
#     0.004 * kdp * lamb / (1 - zdr_lin ** (-1))
# )
# iwc_jst = xr.where(
#     zdr < 0.4,
#     0.033 * (kdp * lamb) ** 0.67 * zdr_lin ** 0.33,
#     0.004 * kdp * lamb / (1 - zdr_lin ** (-1))
# )
data['IWC'] = xr.where(
    data.ZDR_AC_OC < 0.4,
    0.033 * (data.KDP_NC * lamb) ** 0.67 * zh_lin ** 0.33,
    0.004 * data.KDP_NC * lamb / (1 - zdr_lin ** (-1))
)
data['IWC'].values = np.where(data.temp < 273.15 + 4,
                              data['IWC'].values, np.nan)
data.IWC.attrs["long_name"] = 'ice water content'
data.IWC.attrs["short_name"] = 'IWC'
data.IWC.attrs["standard_name"] = 'IWC'
data.IWC.attrs["units"] = 'g m^-3'
data.IWC.attrs["comments"] = 'hybrid retrieval cf. Carlin et al. (2021)'
# data.IWC.attrs["coordinates"] = 'time height'

# TODO: 6.69 or 3.
data['Nt_totice'] = 6.69 - 3 + 2 * np.log10(
    data.IWC.values) - 0.1 * data.ZH_AC

data['Nt_totice'].values = np.where(data.temp < 273.15 + 4,
                                    data['Nt_totice'].values, np.nan)
data['Nt_totice'].values = np.where(data['IWC'] < 0,
                                    np.nan, data['Nt_totice'].values)
data.Nt_totice.attrs["long_name"] = 'total number concentration of ice'
data.Nt_totice.attrs["short_name"] = 'Nt_totice'
data.Nt_totice.attrs["standard_name"] = 'Nt_totice'
data.Nt_totice.attrs["units"] = 'log_10 (L^-1)'  # L -> 3.36, m^3 -> 6.69
data.Nt_totice.attrs["comments"] = 'hybrid retrieval cf. Carlin et al. ' \
                                   '(2021)'
# data.Nt_totice.attrs["coordinates"] = 'time height'

data['Dm_totice'] = 0.67 * (zh_lin / (data.KDP_NC * lamb)) ** (1 / 3)

#     xr.where(
#     data.ZDR_AC_OC < 0.4,
#     -0.1 + 2 * (zdp_lin / (data.KDP_NC * lamb)) ** (1 / 2),
#     0.67 * (zh_lin / (data.KDP_NC*lamb)) ** (1/3)
# )
data['Dm_totice'].values = np.where(data.temp < 273.15 + 4,
                                    data['Dm_totice'].values, np.nan)
data.Dm_totice.attrs["long_name"] = 'mean volume diameter of total ice'
data.Dm_totice.attrs["short_name"] = 'Dm_totice'
data.Dm_totice.attrs["standard_name"] = 'Dm_totice'
data.Dm_totice.attrs["units"] = 'mm'
data.Dm_totice.attrs[
    "comments"] = 'hybrid retrieval cf. Carlin et al. (2021)'
# data.Dm_totice.attrs["coordinates"] = 'time height'

# TODO: Ice mycrophysical retrieval

# do the qvp median!
data = data.median('azimuth')
data['min_entropy'] = qvp_entropy[-1, :]
data['min_entropy'] = data['min_entropy'].assign_attrs(dict(
    units='1', standard_name='minimum of Entropy',
    comments='Minimum of Shaonon information Entropy, calculated for '
             'each of the four pol. variables (lin zh, lin zdr, rho, '
             'K_DP) on azimuth.  Note that K_DP<0 is filtered.'))

data = data.swap_dims({'range': 'height'})
for var in list(data.data_vars):
    data[var].encoding['coordinates'] = 'time height'

data.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
data.attrs['filter_comments'] = filter_comments
# ------------------------------------------------------------------- #
# QVP: SAVE                                                           #
# ------------------------------------------------------------------- #
Path('/'.join(path_qvp.split('/')[:-1])).mkdir(parents=True, exist_ok=True)
data.to_netcdf(path_qvp, unlimited_dims='time')
data.close()
# return


# --------------------------------------------------------------------------- #