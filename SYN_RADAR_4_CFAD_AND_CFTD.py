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
xr.set_options(keep_attrs=True)
import warnings
from statsmodels.stats.weightstats import DescrStatsW
from pathlib import Path

warnings.simplefilter('ignore')


def plot_CFAD_or_CFTD_from_QVP(
        dates=['20170725'],
        hhmm_start='00:00',
        hhmm_end='23:55',
        paths_in=None,
        title=None,
        folder_syn=header.dir_data_qvp,
        locations='PRO',
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
        vmax=None,
        filter_entr_ML=False,
        ax=None,
        save=False,
        save_path=header.folder_plot + 'CFADs/',
        save_name='test_CFAD'
):
    if not isinstance(dates, list):
        dates = [dates]

    if not isinstance(locations, list):
        locations = [locations]

    if paths_in is not None and not isinstance(paths_in, list):
        paths_in = [paths_in]

    mom_all = np.array([], )
    y_all = np.array([], )
    weights_all = np.array([], )
    for j, location in enumerate(locations):
        for i, date in enumerate(dates):
            year = date[0:4]
            mon = date[4:6]
            day = date[6:8]
            date_start = '-'.join([year, mon, day, hhmm_start])
            date_end = '-'.join([year, mon, day, hhmm_end])
            if paths_in:
                path_in = paths_in[j*len(dates)+i]
                if not title:
                    title = path_in.split('/')[-1]
            else:
                path_in = '/'.join([folder_syn, date, da_run,
                                    icon_emvorado_run,
                                    str(spin_up_mm) + 'min_spinup', 'QVP_' +
                                    str(elevation_deg) + '_Syn_' + location +
                                    '_' + date + '0000_' + date + '2355.nc'])
                if not title:
                    title = '-'.join([da_run[4:],
                                      icon_emvorado_run.split('/')[0][5:],
                                      icon_emvorado_run.split('/')[1][5:],
                                      spin_up_mm + 'min'])

            # OPEN
            syn_nc = xr.open_dataset(path_in)
            syn_nc = syn_nc.sel(time=slice(date_start, date_end))
            if filter_entr_ML:
                if 'mlh_top' in list(syn_nc.keys()):
                    ml_top = syn_nc['mlh_top'].transpose('time', ...)
                # elif 'height_ml' in list(syn_nc.keys()):
                elif 'height_ml' in list(syn_nc.coords):  # TS: coordinate
                    ml_top = syn_nc['height_ml'].transpose('time', ...)
                else:
                    print('no ML height found')
                    continue

                syn_nc = xr.where(ml_top > 0, syn_nc, np.nan)
                min_entropy = syn_nc['min_entropy'].transpose('time', ...)
                syn_nc = xr.where(min_entropy > 0.8, syn_nc, np.nan)

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

            mom_min_outer=2*mom_min-mom_max # TODO
            mom_max_outer=2*mom_max-mom_min # TODO
            # bins_y_3=bins_y*3 # TODO
            # mask = (y >= y_min) & (y <= y_max) & \
            #        (mom >= mom_min) & (mom <= mom_max)
            mask = (y >= y_min) & (y <= y_max) & \
                   (mom >= mom_min_outer) & (mom <= mom_max_outer) # TODO
            # TODO: realistic values might be set nan
            mom = mom[mask]
            y = y[mask]
            i_sort = np.argsort(y)
            mom = mom[i_sort]
            y = y[i_sort]
            if ax is None:
                plt.figure(figsize=(6, 5))

            a = plt.hist(y, bins=bins_y, range=(y_min,y_max)) # TODO
            # a = plt.hist(y, bins=bins_y-1, range=(y_min,y_max)) # TODO
            # a = plt.hist(y, bins=bins_y_3-1, )
            weights = np.repeat(100 / a[0], np.int16(a[0]))
            mom_all = np.append(mom_all, mom)
            y_all = np.append(y_all, y)
            weights_all = np.append(weights_all, weights)

    weights_all = weights_all/(len(dates) * len(locations))
    # PLOT
    if vmax:
        extend = 'max'
    else:
        extend = 'neither'

    # h2d, mom2d, y2d, fg = plt.hist2d(mom_all, y_all, bins=[bins_mom, bins_y],
    #                                  range=[[mom_min, mom_max],
    h2d, mom2d, y2d, fg = plt.hist2d(mom_all, y_all, bins=[bins_mom*3, bins_y], # TODO
                                     range=[[mom_min_outer, mom_max_outer], # TODO
                                            [y_min, y_max]], vmax=vmax,
                                     weights=weights_all, cmap='YlGnBu')
    plt.colorbar(label='frequency [%]', extend=extend)
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

    plt.plot(quant_prof[0, ], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.2}$')
    plt.plot(quant_prof[1, ], y_mid, color='red', ls='solid',
             linewidth=2, label='$Q_{0.5}$')
    plt.plot(quant_prof[2, ], y_mid, color='red', ls='dashed',
             linewidth=1, label='$Q_{0.8}$')
    plt.plot(mean_prof, y_mid, color='orange', ls='solid',
             linewidth=2, label='$\mu$')
    plt.legend()
    plt.xlim([mom_min, mom_max]) # TODO
    if vmax and np.max(h2d) > vmax:
        if np.max(h2d)>100:
            print('above 100% :' + str(np.max(h2d)))

        plt.text(1.04, 1.03, min(np.round(np.max(h2d), 1), 100.0),
                 transform=ax.transAxes)
    syn_nc.close()
    plt.tight_layout()
    if save:
        Path(save_path).mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)
