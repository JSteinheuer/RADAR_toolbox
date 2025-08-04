#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# PROCESS_SYN_RADAR.py                                                        #
#                                                                             #
# Functions to calculate QVPs from these synthetic (EMVORADO) volume scans.   #
# --------------------------------------------------------------------------- #

import numpy as np
import pandas as pd
from pathlib import Path
import os
import scipy
import wradlib as wrl
import xarray as xr
import time as time_p
import datetime as dt

import HEADER_RADAR_toolbox as header



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

    return Dmean


# ------------------------------------------------------------------- #
# ------------------------------------------------------------------- #
icon_nc = xr.open_dataset('PATH_TO_ICON_VOLUME_OR_FIELD')
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

# ------------------------------------------------------------------- #