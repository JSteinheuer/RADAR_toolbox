#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 01.11.23                                                 #
# SYN_RADAR_3_PLOT_QVP.py                                                     #
#                                                                             #
# Functions to plot QVPs from given synthetic (EMVORADO) QVP netcdf files.    #
# --------------------------------------------------------------------------- #

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def d0_bringi(
        zdr
):
    d0 = np.full(zdr.shape, np.nan)
    mask_thresh = zdr < 1.25
    d0[mask_thresh] = \
        0.0203 * zdr[mask_thresh] ** 4 - \
        0.149 * zdr[mask_thresh] ** 3 + \
        0.221 * zdr[mask_thresh] ** 2 + \
        0.557 * zdr[mask_thresh] + 0.801
    d0[~mask_thresh] = \
        0.0355 * zdr[~mask_thresh] ** 3 - \
        0.302 * zdr[~mask_thresh] ** 2 + \
        1.06 * zdr[~mask_thresh] + 0.684
    return d0


def plot_qvp_of_polarimetric_variable(
        mom,
        cmap,
        norm,
        levels,
        mom_cs=None,
        levels_cs=None,
        mom_cf=None,
        levels_cf=None,
        ax=None,
        cbar_title=None,
        title=None,
        scale_font=1.2,
        scale_numbers=1,
        y_coord='height',
        extend='both',
        xlabel='UTC [mm-dd hh]',
        ylabel='height [km]',
        top_height=10000,  # in m  # TODO check for TS obs
        mom_height_unit='m',
        save=False,
        save_name='test',
        save_path=''
):
    if ax is None:
        ax = plt.gca()

    fg = mom.plot.contourf(y=y_coord, cmap=cmap, norm=norm,
                           extend=extend, cbar_kwargs={'ticks': levels})
    if mom_cs is not None:
        cs = mom_cs.plot.contour(y=y_coord, colors='black',
                                 levels=levels_cs,
                                 extend=extend, linewidths=.5)
        ax.clabel(cs, **dict(fmt=r'%2.0f', inline=True,
                             fontsize=8 * scale_numbers))
    if mom_cf is not None:
        cf = mom_cf.plot.contour(y=y_coord, colors='dimgrey',
                                 levels=levels_cf,
                                 linestyles='dotted')
        ax.clabel(cf, **dict(fmt='%2.0f °C', inline=True,
                             fontsize=7 * scale_numbers))
    plt.xlabel(xlabel, fontsize=10 * scale_font)
    plt.xticks(fontsize=8 * scale_numbers)
    plt.ylabel(ylabel, fontsize=10 * scale_font)
    plt.ylim(0, top_height)
    if mom_height_unit == 'm':
        plt.yticks(ticks=plt.yticks()[0], labels=plt.yticks()[0] / 1000,
                   fontsize=10 * scale_numbers)
    elif mom_height_unit == 'km':
        plt.yticks(ticks=plt.yticks()[0], labels=plt.yticks()[0],
                   fontsize=10 * scale_numbers)

    if cbar_title is not None:
        fg.colorbar.set_label(cbar_title, fontsize=10 * scale_font)

    fg.colorbar.ax.tick_params(labelsize=10 * scale_numbers)
    if title is not None:
        plt.title(title, fontsize=7 * scale_font)

    plt.tight_layout()
    if save:
        plt.savefig(save_path + save_name + '.pdf',
                    format='pdf', transparent=True)


# --------------------------------------------------------------------------- #
# Thresholds                                                                  #
# --------------------------------------------------------------------------- #
# TODO: necessary here?

d0_min = 0.4
d0_max = 2.5
zh = 0
rho = 0.7
kdp = 0.01
zdr = -1

# --------------------------------------------------------------------------- #
# Data                                                                        #
# --------------------------------------------------------------------------- #
# TODO: necessary here?

# 00 : 5.5
# 01 : 4.5
# 02 : 3.5
# 03 : 2.5
# 04 : 1.5
# 05 : 0.5
# 06 : 8.0
# 07 : 12.0   <- use that
# 08 : 17.0
# 09 : 25.0
scan = '07'
elevations = [5.5, 4.5, 3.5, 2.5, 1.5, 0.5, 8, 12, 17, 25]
elevation = elevations[int(scan)]
