#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# QVP_FROM_VOLUME_SCAN.py                                                     #
#                                                                             #
# Functions to calculate QVPs from given synthetic (EMVORADO) volume scans.   #
# --------------------------------------------------------------------------- #

import pandas as pd
import numpy as np
import xarray as xr
import os
from pathlib import Path
import scipy
from wradlib.dp import kdp_from_phidp
import HEADER_RADAR_toolbox as header


def phidp_from_kdp(da):
    """
    Derive PHIDP from KDP.

    Parameter
    ---------
        da: xarray.DataArray
            array with specific differential phase data

    Return
    ------
        phi: xarray.DataArray
             array with differential phase values
    """
    dr = da.range.diff('range').median('range').values / 1000.
    return xr.apply_ufunc(scipy.integrate.cumtrapz,
                          da,
                          input_core_dims=[["range"]],
                          output_core_dims=[["range"]],
                          dask='parallelized',
                          kwargs=dict(dx=dr, initial=0.0, axis=-1),
                          ) * 2


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

    def calc_entropy(var, dim=dim, n_valid_values=n_lowest):
        minimum = np.nanmin(var.data)
        if minimum <= 0:
            var = xr.where(var.data <= 0, np.nan, var)
            print(var.name + ' has values <= 0 which are set to '
                             'nan (minimum=' + str(minimum) + ')')

        var_normed = var / var.sum(dim, skipna=True)
        values = -((var_normed * np.log10(var_normed)).sum(dim)) \
                 / np.log10(var_normed[dim].size)
        values = values.where(var_normed.count(dim=dim) >= n_valid_values)
        return values

    entropy = calc_entropy(pol_vars[0], dim=dim)
    entropy.name = 'entropy'
    entropy = entropy.expand_dims(dim={'n_all': len(pol_vars) + 1}).copy()
    for i in range(1, len(pol_vars)):
        entropy[i, :] = calc_entropy(pol_vars[i], dim=dim)

    entropy[-1, :] = entropy[:-1, :, :].min(dim='n_all', skipna=False)
    return entropy


def icon_hydromets():
    """
    Preparing dict with ICON 2mom micro-physics (original PSD- and m-D-params)
    """
    hydromets = {}
    hydromets['cloud'] = seifert2general(
        {'nu': 1.0, 'mu': 1.0, 'xmax': 2.6e-10, 'xmin': 4.2e-15,
         'a': 1.24e-01, 'b': 0.333333})
    hydromets['rain'] = seifert2general(
        {'nu': 0.0, 'mu': 0.333333, 'xmax': 3.0e-6, 'xmin': 2.6e-10,
         'a': 1.24e-01, 'b': 0.333333})
    hydromets['ice'] = seifert2general(
        {'nu': 0.0, 'mu': 0.333333, 'xmax': 1.0e-5, 'xmin': 1.0e-12,
         'a': 0.835, 'b': 0.39})
    hydromets['snow'] = seifert2general(
        {'nu': 0.0, 'mu': 0.5, 'xmax': 2.0e-5, 'xmin': 1.0e-10,
         'a': 5.13, 'b': 0.5})
    hydromets['graupel'] = seifert2general(
        {'nu': 1.0, 'mu': 0.333333, 'xmax': 5.3e-4, 'xmin': 4.19e-9,
         'a': 1.42e-1, 'b': 0.314})
    hydromets['hail'] = seifert2general(
        {'nu': 1.0, 'mu': 0.333333, 'xmax': 5.0e-3, 'xmin': 2.6e-9,
         'a': 0.1366, 'b': 0.333333})
    return hydromets


def seifert2general(hymet_seif):
    """
    Convert original ICON (mass-based, Seifert notation) to Dmax-based,
    "general" notation micro-physical parameters.
    (in = COSMO/ICON) N(m) = N0 * m^nu_x * exp(-lam*m^mu_x)
                      D(m) = a * m^b
    (out = EMVORADO)  N(D) = N0 * D^mu_D * exp(-lam*D^nu_D)
                      m(D) = a * D^b   [ -> D = (m/a)^(1/b)
    """
    hymet_gen = {}
    hymet_gen['mu'] = (1.0 / hymet_seif['b']) * (hymet_seif['nu'] + 1.0) - 1.0
    hymet_gen['nu'] = (1.0 / hymet_seif['b']) * hymet_seif['mu']
    hymet_gen['xmax'] = hymet_seif['xmax']
    hymet_gen['xmin'] = hymet_seif['xmin']
    hymet_gen['a'] = (1.0 / hymet_seif['a']) ** (1.0 / hymet_seif['b'])
    hymet_gen['b'] = 1.0 / hymet_seif['b']
    return hymet_gen


def adjust_icon_fields(infields, hymets=icon_hydromets(), spec2dens=1):
    """
    Take hydrometeor fields with qx and qnx, convert them to mass/number
    densities (from specific mass/number contents to densities) if necessary
    (spec2dens=1), and adjust qn to meet ICON min/max bounds of mean
    hydrometeor masses.
    Required input:
        Hydrometeor qx and qnx (all 6), qv (water vapor content), temp
        (temperature), pres (pressure).
    Output:
        volume specific qx, qnx.
    """
    r_d = 287.05  # dry air gas constant
    r_v = 461.51  # water vapor gas constant
    q = {}
    qn = {}

    q['cloud'] = infields['qc'][:].data
    qn['cloud'] = infields['qnc'][:].data
    q['rain'] = infields['qr'][:].data
    qn['rain'] = infields['qnr'][:].data
    q['ice'] = infields['qi'][:].data
    qn['ice'] = infields['qni'][:].data
    q['snow'] = infields['qs'][:].data
    qn['snow'] = infields['qns'][:].data
    q['graupel'] = infields['qg'][:].data
    qn['graupel'] = infields['qng'][:].data
    q['hail'] = infields['qh'][:].data
    qn['hail'] = infields['qnh'][:].data

    qv = infields['qv'][:].data
    t = infields['temp'][:].data
    p = infields['pres'][:].data
    qx = q['cloud'] + q['rain'] + q['ice'] + \
         q['snow'] + q['graupel'] + q['hail']
    rvd_m_o = r_v / r_d - 1.0
    if spec2dens:
        rho = p / (r_d * t * (1.0 + rvd_m_o * qv - qx))
        for key in q:
            q[key] = q[key] * rho
            qn[key] = qn[key] * rho

    for key in q:
        x = q[key] / (qn[key] + 1e-38)
        x[x < hymets[key]['xmin']] = hymets[key]['xmin']
        x[x > hymets[key]['xmax']] = hymets[key]['xmax']
        qn[key] = q[key] / x
        q[key][q[key] < 1e-7] = 0.
        qn[key][q[key] < 1e-7] = 1.  # for numeric stability.
    return q, qn


def gfct(x):
    """
    Gamma function (as implemented and used in EMVORADO)
    """
    c1 = 76.18009173
    c2 = -86.50532033
    c3 = 24.01409822
    c4 = -1.231739516
    c5 = 0.120858003e-2
    c6 = -0.536382e-5
    stp = 2.50662827465

    tmp = x + 4.5
    p = stp * (1.0 + c1 / x + c2 / (x + 1.0) + c3 / (x + 2.0) +
               c4 / (x + 3.0) + c5 / (x + 4.0) + c6 / (x + 5.0))
    g_fct = p * np.exp((x - 0.5) * np.log(tmp) - tmp)
    return g_fct


def mgdparams(q, qn, hymets=icon_hydromets()):
    """
    Derive non-constant MGD-PSD parameters N0 and lam (mu & nu set constant
    from ICON microphysics)
    """
    mgd = {}
    for key in q:
        mgd[key] = {}
        tmp1 = (hymets[key]['mu'] + 1.0) / hymets[key]['nu']
        tmp2 = (hymets[key]['mu'] + hymets[key]['b'] + 1.0) / hymets[key]['nu']
        gfct_tmp1 = gfct(tmp1)
        gfct_tmp2 = gfct(tmp2)
        mgd[key]['lam'] = np.ones_like(q[key])
        mgd[key]['n0'] = np.zeros_like(q[key])
        mgd[key]['lam'][q[key] > 0.] = \
            ((hymets[key]['a'] * qn[key][q[key] > 0.] * gfct_tmp2) /
             (q[key][q[key] > 0.] * gfct_tmp1)) ** (hymets[key]['nu'] /
                                                    hymets[key]['b'])
        mgd[key]['n0'][q[key] > 0.] = \
            qn[key][q[key] > 0.] * hymets[key]['nu'] * \
            mgd[key]['lam'][q[key] > 0.] ** tmp1 / gfct_tmp1
    return mgd


def calc_moments(mgd, k=[0, 3, 4], hymets=icon_hydromets(), adjust=0):
    """
    Calculate (fields of) (Dmax-based) PSD moments from given constant (mu,nu)
    and qx/qnx derived fields of variable (N0,lam) MGD-parameters.
    Usage:
      k: list of moments to be computed
      adjust: Flag whether to zero out (all) k moments if one of them is 0.
            Beside for qx=0, mom[k]=0 might occur due to number accuracy
            limitations. So, to avoid different moments contributing
            differently to multi-PSD statistics, zero all remaining if one of
            them has number issues.
    """
    Mk = {}
    for key in mgd:
        Mk[key] = {}
        for kk in k:
            c1 = (hymets[key]['mu'] + kk + 1) / hymets[key]['nu']
            c2 = gfct(c1) / hymets[key]['nu']
            # initialize Mk[kk]
            Mk[key][kk] = np.zeros_like(mgd[key]['n0'])
            # calculate Mk[kk] only if q, ie n0, is non-zero
            Mk[key][kk][mgd[key]['n0'] > 0.] = \
                c2 * mgd[key]['n0'][mgd[key]['n0'] > 0.] / \
                mgd[key]['lam'][mgd[key]['n0'] > 0.] ** c1
        if adjust:
            for kk in k:
                for ik in k:
                    Mk[key][kk][Mk[key][ik] == 0.] = 0.
    return Mk


def calc_multimoments(Mk, inhm=['ice', 'snow', 'graupel', 'hail'],
                      outhm='totice'):
    Mk[outhm] = {}
    for hm in inhm:
        assert (hm in Mk.keys()), \
            'Class %s to accumulate moments not in input moments dict' % hm
    for kk in Mk[inhm[0]]:
        Mk[outhm][kk] = np.zeros_like(Mk[inhm[0]][kk])
        for hm in inhm:
            Mk[outhm][kk] += Mk[hm][kk]
    return Mk


def calc_Dmean(Mk, meantype='vol'):
    k_needed = {'vol': [4, 3], 'mass': [3, 0]}
    Dmean = {}
    for hm in Mk.keys():
        Dmean[hm] = np.zeros_like(Mk[hm][k_needed[meantype][0]])
        Dmean[hm][:] = np.nan
        Dmean[hm] = np.where(Mk[hm][k_needed[meantype][0]] > 0,
                             (Mk[hm][k_needed[meantype][0]] /
                              Mk[hm][k_needed[meantype][1]]) **
                             (1. / (k_needed[meantype][0] -
                                    k_needed[meantype][1])), 0)

    #
    return Dmean


def normalise(da, dim):  # yes
    damin = da.min(dim=dim, skipna=True, keep_attrs=True)
    damax = da.max(dim=dim, skipna=True, keep_attrs=True)
    da = (da - damin) / (damax - damin)
    return da


def sobel(da, dim):  # yes
    axis = da.get_axis_num(dim)
    func = lambda x, y: scipy.ndimage.sobel(x, y)
    return xr.apply_ufunc(func, da, axis, dask='parallelized')


def ml_normalise(ds, moments=dict(DBZH=(10., 60.), RHOHV=(0.65, 1.),
                                  PHIDP=(-0.5, 360)), dim='height'):  # yes
    assign = {}
    for m, span in moments.items():
        v = ds[m]
        v = v.where((v >= span[0]) & (v <= span[1]))
        v = v.pipe(normalise, dim=dim)
        assign.update({f'{m}_norm': v})
    return xr.Dataset(assign)


def ml_clean(ds, zh, dim='height'):  # yes
    # removing profiles with too few data (>93% of array is nan)
    good = ds[zh + "_norm"].count(dim=dim) / ds[zh + "_norm"].sizes[dim] * 100
    return ds.where(good > 7)


def ml_combine(ds, zh, rho):  # yes
    comb = (1 - ds[rho + "_norm"]) * ds[zh + "_norm"]
    return comb


def ml_gradient(ds, dim='height', ywin=None):  # yes
    assign = {}
    for m in ds.data_vars:
        # step 3 sobel filter
        dy = ds[m].pipe(sobel, dim)
        # step 3 smoothing
        # currently smoothes only in y-direction (height)
        dy = dy.rolling({dim: ywin}, center=True).mean()
        assign.update({f'{m}_dy': dy})
    return ds.assign(assign)


def ml_noise_reduction(ds, thres=0.02):  # yes
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


def melting_layer_qvp_X_new(ds, moments=dict(zrsim=(10., 60.),
                                             rhvsim=(0.65, 1.),
                                             phidpsim=(-0.5, 360)),
                            dim='height',
                            thres=0.02, xwin=5, ywin=5, fmlh=0.3,
                            all_data=False):  # yes
    # step 0
    zh = [k for k in moments if ("zh" in k.lower()) or
          ("th" in k.lower()) or ("zr" in k.lower())][0]
    rho = [k for k in moments if ("rho" in k.lower()) or
           ("rhv" in k.lower()) or ("cc" in k.lower())][0]
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
    p4 = np.sign(ds4[zh + "_norm_dy"].differentiate(dim)).diff(dim).idxmin(
        dim, skipna=True)
    ds5 = ds4.where(ds4.height < p4)
    ds5 = ml_height_top_new(ds5, moment=phi + "_norm", dim=dim)
    # ds5.mlh_top, ds5.ml_top, derived according Wolfensberger et al.
    # but used PHIDP_norm instead of DBZH_norm

    # step 9b: ml bottom refinement
    ds6 = ds0.where((ds0.height < ds3.mlh_bottom) &
                    (ds0.height > below))  # This, i.e., both conditions!
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


def qvp_from_syn_vol(day='20170725', da_run='ASS_2211',
                     icon_run='MAIN_2211.0',
                     icon_emvorado_run='MAIN_2211.0/EMVO_00000000.2',
                     spin_up_mm=30, elevation_deg=12, radar_loc='PRO',
                     dir_data_in=header.dir_data_vol,
                     dir_data_out=header.dir_data_qvp,
                     overwrite=False, merge=True, lowest_rhv=0.7,
                     lowest_kdp=None, highest_kdp=None, n_lowest=30):
    """
    Create QVP for one day out of 8 synthetic volume scans from EMVORADO and
    ICON data (each 4times 6hourly).

     Args:
        day: day in >yyyymmdd<.
        da_run: subfolder specifying the data_assimilation run.
        icon_run: subfolder specifying the ICON run.
        icon_emvorado_run: subfolders specifying the EMVORADO run.
        spin_up_mm: lower threshold for the time that must have elapsed
                    since the last DA.
        elevation_deg: elevation angle to be choosen from Syn_vol for QVP
                       [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0].
        radar_loc: string >RRR< naming the radar.
        dir_data_in: directory with folder >yyyymmdd< of the day outputs.
        dir_data_out: directory for the output.
        overwrite: If True, then process only if output is not existing.
        merge: If True, the 4 QVPs of the day will be merged.
        lowest_rhv: If not None, it is the lowest acceptable value not
                    filtered.
        lowest_kdp: If not None, it is the lowest acceptable value not
                    filtered.
        highest_kdp: If not None, it is the highest acceptable value not
                     filtered.
        n_lowest: lowest nr of valid values for returning a non nan value per
                  360 azimuths.

    Returns:
    """
    time_start = day + '0000'
    time_end = day + '2355'
    spin_up_mm = str(spin_up_mm)
    dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
    dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
    dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
    dir_qvp = dir_data_out + dti[0].strftime('%Y%m%d') + '/' + \
              da_run + '/' + icon_emvorado_run + '/' + \
              str(spin_up_mm) + 'min_spinup/'
    file_qvp_4 = 'QVP_' + str(elevation_deg) + '_Syn_' + radar_loc + '_' + \
                 dti[0].strftime('%Y%m%d%H%M') + '_' + \
                 dti[-1].strftime('%Y%m%d%H%M') + '.nc'
    files_qvp = []
    if os.path.isfile(dir_qvp + file_qvp_4) and not overwrite:
        print(dir_qvp + file_qvp_4 + ' already exists;\n' +
              'set >overwrite=True< for recalculation')
        return

    for hhmm_start in ['0000', '0600', '1200', '1800']:
        time_start = day + hhmm_start
        time_end = day + str(int(hhmm_start) + 555).zfill(4)
        filter_comments = 'PPI ICON data are only included, if ' \
                          'corresponding EMVORADO PPI data exists. '
        # ------------------------------------------------------------------- #
        # open nc files                                                       #
        # ------------------------------------------------------------------- #
        dti_start = pd.to_datetime(time_start, format="%Y%m%d%H%M")
        dti_end = pd.to_datetime(time_end, format="%Y%m%d%H%M")
        dti = pd.date_range(dti_start, dti_end, freq="5min", inclusive='both')
        dir_emv = dir_data_in + dti[0].strftime('%Y%m%d') + '/' + \
                  da_run + '/' + icon_emvorado_run + '/' + \
                  str(spin_up_mm) + 'min_spinup/'
        dir_icon = dir_data_in + dti[0].strftime('%Y%m%d') + '/' + \
                   da_run + '/' + icon_run + '/ICONdata/' + \
                   str(spin_up_mm) + 'min_spinup/'
        file_emv = 'EMV_Vol_' + radar_loc + '_' + \
                   dti[0].strftime('%Y%m%d%H%M') + '_' + \
                   dti[-1].strftime('%Y%m%d%H%M') + '.nc'
        file_icon = 'ICON_Vol_' + radar_loc + '_' + \
                    dti[0].strftime('%Y%m%d%H%M') + '_' + \
                    dti[-1].strftime('%Y%m%d%H%M') + '.nc'
        file_qvp = 'QVP_' + str(elevation_deg) + '_Syn_' + radar_loc + '_' + \
                   dti[0].strftime('%Y%m%d%H%M') + '_' + \
                   dti[-1].strftime('%Y%m%d%H%M') + '.nc'
        files_qvp.append(file_qvp)
        if os.path.isfile(dir_qvp + file_qvp) and not overwrite:
            print(dir_qvp + file_qvp + ' already exists;\n' +
                  'set >overwrite=True< for recalculation')
            continue

        if not os.path.isfile(dir_emv + file_emv):
            print(dir_emv + file_emv + ' not found')
            merge = False
            continue

        if not os.path.isfile(dir_icon + file_icon):
            print(dir_emv + file_emv + ' found but not ' +
                  dir_icon + file_icon)
            merge = False
            continue

        print('------------------------------------------------------')
        print('QVP for ' + dir_emv + file_emv)
        emv_nc = xr.open_dataset(dir_emv + file_emv)
        icon_nc = xr.open_dataset(dir_icon + file_icon)

        # ------------------------------------------------------------------- #
        # PPI: common EMV and ICON mask                                       #
        # ------------------------------------------------------------------- #
        emv_nc = emv_nc.sel(elevation=elevation_deg, drop=True)
        icon_nc = icon_nc.sel(elevation=elevation_deg, drop=True)
        emv_nc.attrs['elevation'] = str(elevation_deg) + ' degrees'
        emv_nc.attrs['icon_run'] = icon_nc.icon_run
        if 'phidpsim' not in list(emv_nc.data_vars):
            kdpsim = emv_nc.kdpsim.fillna(0)
            emv_nc['phidpsim'] = (['time', 'range', 'azimuth', ],
                phidp_from_kdp(kdpsim).transpose('time', 'range', 'azimuth'
                                                 ).data, dict(
                standard_name='Differential propagation phase shift',
                comments='PHIDP from KDP',
                units='deg/km'))

        emv_nc = emv_nc.transpose('time', 'range', 'azimuth')
        if lowest_rhv:
            emv_nc['rhvsim'] = \
                emv_nc['rhvsim'].where(emv_nc['rhvsim'] >= lowest_rhv)
            filter_comments = filter_comments + 'rhvsim >= ' + \
                              str(lowest_rhv) + '. '
        if lowest_kdp:
            emv_nc['kdpsim'] = \
                emv_nc['kdpsim'].where(emv_nc['kdpsim'] >= lowest_kdp)
            filter_comments = filter_comments + 'kdpsim >= ' + \
                              str(lowest_kdp) + ' . '
        if highest_kdp:
            emv_nc['kdpsim'] = \
                emv_nc['kdpsim'].where(emv_nc['kdpsim'] <= highest_kdp)
            filter_comments = filter_comments + 'kdpsim <= ' + \
                              str(highest_kdp) + ' . '

        mask = ~ np.isnan(emv_nc['zrsim'])
        for var in list(emv_nc.data_vars):
            if emv_nc[var].dims == ('time', 'range', 'azimuth'):
                mask = mask * (~ np.isnan(emv_nc[var]))

        # # should be all non nan:
        # for var in list(icon_nc.data_vars):
        #     if icon_nc[var].dims == ('time', 'range', 'azimuth'):
        #         mask = mask * (~ np.isnan(icon_nc[var]))
        #
        for var in list(emv_nc.data_vars):
            if emv_nc[var].dims == ('time', 'range', 'azimuth'):
                emv_nc[var] = emv_nc[var].where(mask)

        # ------------------------------------------------------------------- #
        # QVP: KDP Melting layer correction                                   #
        # ------------------------------------------------------------------- #
        # Wolfensberger
        qvp_first_nc = emv_nc.copy()
        for var in ['zrsim', 'phidpsim', 'rhvsim']:
            count_values = qvp_first_nc[var].count(dim='azimuth')
            qvp_first_nc[var] = qvp_first_nc[var].median(
                dim='azimuth', skipna=True, keep_attrs=True)
            qvp_first_nc[var] = qvp_first_nc[var].where(
                count_values >= n_lowest)

        filter_comments = filter_comments + 'QVP needs at least ' + \
            str(n_lowest) + ' of 360 valid values per median calculation on ' \
                            'azimuth. '
        qvp_first_nc['height'] = (['range'], qvp_first_nc['range'].data *
                                  np.sin(elevation_deg * np.pi / 180.) / 1000.,
                                  dict(standard_name='height above radar',
                                       units='km'))
        qvp_first_nc = qvp_first_nc.swap_dims({'range': 'height'})
        qvp_first_nc = melting_layer_qvp_X_new(ds=qvp_first_nc, all_data=False)
        emv_nc['mlh_top'] = qvp_first_nc.mlh_top
        emv_nc['mlh_bottom'] = qvp_first_nc.mlh_bottom
        emv_nc['height'] = (['range'], emv_nc['range'].data *
                            np.sin(elevation_deg * np.pi / 180.) / 1000.,
                            dict(standard_name='height above radar',
                                 units='km'))
        emv_nc['phidpsim_ml_corrected'] = \
            ml_interpolate2(qvp_first_nc, emv_nc, 'phidpsim', method='linear')
        emv_nc = emv_nc.transpose('time', 'azimuth', 'range')
        kdp_new = kdp_from_phidp(emv_nc['phidpsim_ml_corrected'].data,
                           winlen=13,  # cband
                           # winlen=31, # xband ?!
                           min_periods=3)
        kdp_new[np.where(np.isnan(emv_nc['kdpsim'].data))]=np.nan
        emv_nc['kdpsim_ml_corrected'] = (['time', 'azimuth', 'range'],
                                         kdp_new, dict(
            standard_name=emv_nc['kdpsim'].standard_name,
            units=emv_nc['kdpsim'].units,
            comments=emv_nc['kdpsim'].comments + '; Melting layer corrected'))
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
            first_valid_height_ml == None, emv_nc['mlh_top'].data,
            first_valid_height_ml)
        # ML BOTTOM Giagrande refinement
        panda_above_min = new_cut_above_min_ML_filter.to_pandas().transpose()
        last_valid_height_ml = pd.DataFrame(panda_above_min).apply(
            pd.Series.last_valid_index)
        last_valid_height_ml = last_valid_height_ml.to_xarray()
        last_valid_height_ml = xr.where(last_valid_height_ml.data == None,
                                        emv_nc['mlh_bottom'].data,
                                        last_valid_height_ml)
        emv_nc['mlh_top_gia'] = (["time"], first_valid_height_ml.data, dict(
            standard_name=emv_nc['mlh_top'].standard_name + ' refined',
            units=emv_nc['mlh_top'].units,
            comments=emv_nc['mlh_top'].comments + ' Refined with Giagrande.'))
        emv_nc['mlh_bottom_gia'] = (["time"], last_valid_height_ml.data, dict(
            standard_name=emv_nc['mlh_bottom'].standard_name + ' refined',
            units=emv_nc['mlh_bottom'].units,
            comments=emv_nc[
                         'mlh_bottom'].comments + ' Refined with Giagrande.'))

        # ------------------------------------------------------------------- #
        # QVP: EMV                                                            #
        # ------------------------------------------------------------------- #
        pol_vars = [10 ** (0.1 * emv_nc.zrsim), 10 ** (0.1 * emv_nc.zdrsim),
                    emv_nc.rhvsim, emv_nc.kdpsim_ml_corrected]
        qvp_entropy = entropy_of_vars(pol_vars, dim='azimuth',
                                      n_lowest=n_lowest)
        for var in list(emv_nc.data_vars):
            if 'azimuth' in emv_nc[var].dims:
                count_values = emv_nc[var].count(dim='azimuth')
                emv_nc[var] = emv_nc[var].median(dim='azimuth',
                                                 skipna=True, keep_attrs=True)
                emv_nc[var] = emv_nc[var].where(count_values >= n_lowest)

        emv_nc['azimuth'] = emv_nc['azimuth'].median(dim='azimuth',
                                                     skipna=True,
                                                     keep_attrs=True)
        emv_nc['azimuth'].data = np.nan

        # ------------------------------------------------------------------- #
        # PPI: qi,qni -> diameters                                            #
        # ------------------------------------------------------------------- #
        q_dens, qn_dens = adjust_icon_fields(icon_nc)
        multi_params = mgdparams(q_dens, qn_dens)
        moments = calc_moments(mgd=multi_params)
        multimoments = calc_multimoments(moments)
        mean_volume_diameter = calc_Dmean(multimoments)
        for hm in ['graupel', 'ice', 'rain', 'hail', 'cloud', 'snow']:
            icon_nc['vol_q' + hm[0]] = (
                ['time', 'range', 'azimuth', ], q_dens[hm], dict(
                    standard_name='volume ' + icon_nc[
                        'q' + hm[0]].standard_name,
                    units='m3 m-3'))
            icon_nc['vol_qn' + hm[0]] = (
                ['time', 'range', 'azimuth', ], qn_dens[hm], dict(
                    standard_name=icon_nc['qn' + hm[0]].standard_name +
                                  ' per volume', units='m-3'))
            icon_nc['D0_' + hm[0]] = (
                ['time', 'range', 'azimuth', ],
                mean_volume_diameter[hm] * 1000,
                dict(standard_name='mean volume diameter of ' +
                                   icon_nc['qn' + hm[0]].standard_name[21:],
                     units='mm'))

        vol_qtotice = xr.concat(
            [icon_nc.vol_qi, icon_nc.vol_qs, icon_nc.vol_qh,
             icon_nc.vol_qg], dim="ice_hydrometeors")
        vol_qtotice = vol_qtotice.sum(dim="ice_hydrometeors", skipna=False)
        icon_nc['vol_qtotice'] = (['time', 'range', 'azimuth', ],
                                  vol_qtotice.data, dict(
            standard_name='volume specific total ice water content',
            comments='vol_qg + vol_qh + vol_qi + vol_qs',
            units='m3 m-3'))
        vol_qntotice = xr.concat([icon_nc.vol_qni, icon_nc.vol_qns,
                                  icon_nc.vol_qnh, icon_nc.vol_qng],
                                 dim="ice_hydrometeors")
        vol_qntotice = vol_qntotice.sum(dim="ice_hydrometeors", skipna=False)
        icon_nc['vol_qntotice'] = (['time', 'range', 'azimuth', ],
                                   vol_qntotice.data, dict(
            standard_name='number concentration total ice water content',
            comments='vol_qng + vol_qnh + vol_qni + vol_qns',
            units='m-3'))
        icon_nc['D0_totice'] = (['time', 'range', 'azimuth', ],
                                mean_volume_diameter['totice'] * 1000, dict(
            standard_name='mean volume diameter of total ice',
            units='mm'))
        for var in list(icon_nc.data_vars):
            if icon_nc[var].dims == ('time', 'range', 'azimuth'):
                icon_nc[var] = icon_nc[var].where(mask)

        # ------------------------------------------------------------------- #
        # QVP: ICON                                                           #
        # ------------------------------------------------------------------- #
        for var, da in icon_nc.data_vars.items():
            if 'azimuth' in icon_nc[var].dims:
                count_values = icon_nc[var].count(dim='azimuth')
                icon_nc[var] = icon_nc[var].median(dim='azimuth',
                                                   skipna=True,
                                                   keep_attrs=True)
                icon_nc[var] = icon_nc[var].where(count_values >= n_lowest)

        icon_nc['azimuth'] = icon_nc['azimuth'].median(dim='azimuth',
                                                       skipna=True,
                                                       keep_attrs=True)
        icon_nc['azimuth'] = np.nan

        # ------------------------------------------------------------------- #
        # QVP: MERGE                                                          #
        # ------------------------------------------------------------------- #
        qvp_nc = xr.merge([emv_nc, icon_nc])
        qvp_nc.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
        qvp_nc['min_entropy'] = qvp_entropy[-1, :]
        qvp_nc['min_entropy'] = qvp_nc['min_entropy'].assign_attrs(dict(
            units='1', standard_name='minimum of Entropy',
            comments='Minimum of Shaonon information Entropy, calculated for '
                     'each of the four pol. variables (lin zh, lin zdr, rho, '
                     'K_DP) on azimuth.  Note that K_DP<0 is filtered.'))
        qvp_nc.attrs['filter_comments'] = filter_comments
        qvp_nc = qvp_nc.swap_dims({'range': 'height'})
        qvp_nc = qvp_nc.drop_vars('range')
        qvp_nc = qvp_nc.drop_vars('azimuth')

        # ------------------------------------------------------------------- #
        # QVP: SAVE                                                           #
        # ------------------------------------------------------------------- #
        emv_nc.close()
        icon_nc.close()
        Path(dir_qvp).mkdir(parents=True, exist_ok=True)
        qvp_nc.to_netcdf(dir_qvp + file_qvp, unlimited_dims='time')
        qvp_nc.close()

    # ----------------------------------------------------------------------- #
    if merge:
        qvp_nc = xr.merge([xr.open_dataset(dir_qvp + files_qvp[0]),
                           xr.open_dataset(dir_qvp + files_qvp[1]),
                           xr.open_dataset(dir_qvp + files_qvp[2]),
                           xr.open_dataset(dir_qvp + files_qvp[3]),
                           ])
        qvp_nc.attrs['processing_date'] = str(pd.Timestamp.today())[:16]
        qvp_nc.to_netcdf(dir_qvp + file_qvp_4, unlimited_dims='time')
        qvp_nc.close()
        for file_qvp in files_qvp:
            os.remove(dir_qvp + file_qvp)

    return
