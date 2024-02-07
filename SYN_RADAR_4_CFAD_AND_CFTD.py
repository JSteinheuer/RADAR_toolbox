#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.02.24                                                 #
# SYN_RADAR_4_CFAD_AND_CFTD.py                                                #
#                                                                             #
# plotting routine for Continuous-frequency-altitude-diagram (CFAD) and       #
# Continuous-frequency-temperatur-diagram (CFTD).                             #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import warnings
from statsmodels.stats.weightstats import DescrStatsW
from pathlib import Path

warnings.simplefilter('ignore')


def plot_CFAD_or_CFTD_from_QVP(
        date='20170725',
        hhmm_start='00:00',
        hhmm_end='23:55',
        path_in=None,
        title=None,
        folder_syn=header.dir_data_qvp,
        location='PRO',
        elevation_deg=12,
        da_run='ASS_2211',
        icon_emvorado_run='MAIN_2401.3/EMVO_00500000.2',
        spin_up_mm='60',
        moment='zrsim',
        mom_min=0,
        mom_max=40,
        bins_mom=40,
        vert_temp=True,  # CFTD
        temp_min=-20,
        temp_max=16,
        bins_temp=18,
        height_min=0,  # in km
        height_max=10,  # in km
        bins_height=20,
        ax=None,
        save=False,
        save_path=header.folder_plot + 'CFADs/',
        save_name='test_CFAD'
):
    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]
    date_start = '-'.join([year, mon, day, hhmm_start])
    date_end = '-'.join([year, mon, day, hhmm_end])
    if path_in:
        if not title:
            title = path_in.split('/')[-1]
    else:
        path_in = '/'.join([folder_syn, date, da_run, icon_emvorado_run,
                            str(spin_up_mm) + 'min_spinup', 'QVP_' +
                            str(elevation_deg) + '_Syn_' + location + '_' +
                            date + '0000_' + date + '2355.nc'])
        if not title:
            title = '-'.join([da_run[4:],
                              icon_emvorado_run.split('/')[0][5:],
                              icon_emvorado_run.split('/')[1][5:],
                              spin_up_mm + 'min'])

    # OPEN
    syn_nc = xr.open_dataset(path_in)
    syn_nc = syn_nc.sel(time=slice(date_start, date_end))
    mom = syn_nc[moment].transpose('time', ...).values
    if vert_temp:
        y = syn_nc.temp.transpose('time', ...).values
        y_min, y_max, bins_y = temp_min, temp_max, bins_temp
        if 'units' in syn_nc.temp.attrs:  # else: TS in °C
            if syn_nc.temp.units == 'K':
                y = y - 273.15
    else:
        y = np.array([syn_nc.height.values] * syn_nc.time.size) / 1000
        y_min, y_max, bins_y = height_min, height_max, bins_height
        if 'units' in syn_nc.height.attrs:  # else: TS in m
            if syn_nc.height.units == 'km':
                y = y * 1000

    mask = (y >= y_min) & (y <= y_max) & \
           (mom >= mom_min) & (mom <= mom_max)  # TODO: realistic values might
    # TODO: ... be set nan
    mom = mom[mask]
    y = y[mask]
    i_sort = np.argsort(y)
    mom = mom[i_sort]
    y = y[i_sort]

    # PLOT
    if ax is None:
        plt.figure(figsize=(6, 5))

    a = plt.hist(y, bins=bins_y, )
    weights = np.repeat(100 / a[0], np.int16(a[0]))
    h2d, mom2d, y2d, fg = plt.hist2d(mom, y, bins=[bins_mom, bins_y],
                                     range=[[mom_min, mom_max],
                                            [y_min, y_max]],
                                     weights=weights, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]')
    if vert_temp:
        plt.gca().invert_yaxis()
        plt.ylabel('temperature (°C)')
    else:
        plt.ylabel(r'height (km)')

    if 'standard_name' in syn_nc[moment].attrs:  # else: TS without any
        plt.xlabel(
            syn_nc[moment].standard_name + ' (' + syn_nc[moment].units + ')')
    else:
        plt.xlabel(moment)

    plt.title(title)
    y_mid = y2d[1:] / 2 + y2d[:-1] / 2
    mom_mid = mom2d[1:] / 2 + mom2d[:-1] / 2
    quant_prof = np.zeros([3, len(y_mid)])
    mean_prof = np.zeros(len(y_mid))
    for t_i in range(len(y_mid)):
        wq = DescrStatsW(data=mom_mid, weights=h2d[:, t_i])
        quant_prof[:, t_i] = wq.quantile(probs=np.array([0.2, 0.5, 0.8]),
                                         return_pandas=False)
        mean_prof[t_i] = wq.mean

    plt.plot(quant_prof[0,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.2}$')
    plt.plot(quant_prof[1,], y_mid, color='red', ls='solid',
             linewidth=2, label='$Q_{0.5}$')
    plt.plot(quant_prof[2,], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.8}$')
    plt.plot(mean_prof, y_mid, color='orange', ls='solid',
             linewidth=2, label='$\mu$')
    plt.legend()
    syn_nc.close()
    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)

