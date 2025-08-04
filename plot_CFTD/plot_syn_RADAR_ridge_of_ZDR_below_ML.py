#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 13.02.25                                                 #
# plot_syn_RADAR_QVP_of_4_polmoms.py                                          #
#                                                                             #
# Run the QVP functions in PLOT_SYN_RADAR.py for generating specific QVP plot.#
# --------------------------------------------------------------------------- #

import os
import xarray as xr
from scipy.stats import wasserstein_distance

import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list, mom_plot_dict
from SET_SYN_RADAR import rad_dict
from statsmodels.stats.weightstats import DescrStatsW
import pandas as pd
import seaborn as sns
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #


# ------------------------------------ #
# CFTDs                                #
# ------------------------------------ #
hhmm_start_cftds = '00:00'
hhmm_end_cftds = '23:59'

# ------------------------------------ #
# full: ------------------------------ #
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
data_max = 31000
testing = False
# testing: --------------------------- #
# locations = ['ESS']  # TODO: remove
# dates = ['20210714']  # TODO: remove
# data_max = 2200  # TODO: remove
# testing = True
# ------------------------------------ #

# CFADs ? ---------------------------- #
vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ! ------------------------- #
vert_temp = True
temp_min = 2
temp_max = 16
bins_temp = 7
# ------------------------------------ #

# ------------------------------------ #
# QVPs & CFTDs                         #
# ------------------------------------ #
elevation_deg = 12
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
letters='abcdefghijklmnopqrstuvwxyz\u03B1\u03B2'
letters_i=0
# with entropy filter: --------------- #
filter_entr = True
filter_entr_at = 0.7
filter_moms = False
# without - not for paper!: ----- -----#
# filter_entr = False  # TODO  # only testing
# filter_entr_at = 0  # TODO  # only testing
# filter_moms = True  # TODO  # only testing
# ------------------------------------ #

# ------------------------------------ #
# MODELS                               #
# ------------------------------------ #
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
short_names = []
colors = []
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00010000.2')
# spin_up_mms.append('120')
# short_names.append('R0E1')
# # short_names.append(icon_emvorado_runs[-1][5:-15]+icon_emvorado_runs[-1][-10:] + '/' + spin_up_mms[-1])
# colors.append('red')
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
# spin_up_mms.append('120')
# short_names.append('R0E2')
# # short_names.append(icon_emvorado_runs[-1][5:-15]+icon_emvorado_runs[-1][-10:] + '/' + spin_up_mms[-1])
# colors.append('orange')
# ------------------------------------ #
# SYN data row 3                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R0E3')
# short_names.append(icon_emvorado_runs[-1][5:-15]+icon_emvorado_runs[-1][-10:] + '/' + spin_up_mms[-1])
colors.append('green')
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.1/EMVO_00510000.2')
spin_up_mms.append('120')
short_names.append('R1E3')
# short_names.append(icon_emvorado_runs[-1][5:-15]+icon_emvorado_runs[-1][-10:] + '/' + spin_up_mms[-1])
colors.append('magenta')
# ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
# da_runs.append('ASS_2411')
# icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
# spin_up_mms.append('120')
# short_names.append('R2E3')
# # short_names.append(icon_emvorado_runs[-1][5:-15]+icon_emvorado_runs[-1][-10:] + '/' + spin_up_mms[-1])
# colors.append('blue')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# CFTDs                                                                       #
# --------------------------------------------------------------------------- #
# ------------------------------------ #
# Ridge plot parameters                #
# ------------------------------------ #

if len(dates) == 1:
    dates_str = dates[0] + '_'
elif dates == ['20210714', '20210713'] or dates == ['20210713', '20210714']:
    dates_str = '20210713+14_'
else:
    dates_str = 'All_t_'

if len(locations) == 1:
    locations_str = locations[0] + '_'
else:
    locations_str = 'all_radars_'

folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
# --------------------------------------------------------------------------- #
# CFTDs OBS row 1                                                             #
# --------------------------------------------------------------------------- #
color='black'
# ------------------------------------ #
x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start_cftds,
    hhmm_end=hhmm_end_cftds,
    elevation_deg=elevation_deg,
    da_icon_emvorado_run=None,
    moment='ZDR_AC_OC',
    vert_temp=vert_temp,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    filter_entr=filter_entr,
    filter_entr_at=filter_entr_at,
    filter_moms=filter_moms,
    ax=None,
    save=False,
    color=color,
    plot_legend=True,
    plot_data=False,
    data_max=data_max,
)
plt.close()
# ------------------------------------ #

# Create the data
temp_rd=(np.round(y/2,0)*2)
df = pd.DataFrame(dict(mom=x, temp=temp_rd))
x_obs=x[:]
temp_obs=temp_rd[:]

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(19, rot=-.25, light=.7)
gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15,
               height=.4, palette=pal,ylim=[0,.6], xlim=[-1.99,6.1]) # zdr obs

# Draw the densities in a few steps
gg.map(sns.kdeplot, "mom",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
gg.map(sns.kdeplot, "mom", clip_on=False, color="lightgrey", lw=2, bw_adjust=.5)

# passing color=None to refline() uses the hue mapping
gg.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


gg.map(label, "mom")

# Set the subplots to overlap
gg.figure.subplots_adjust(hspace=-.25)

# Remove axes details that don't play well with overlap


# Remove axes details that don't play well with overlap
gg.set_titles("")
gg.set(yticks=[], ylabel="")
gg.despine(bottom=True, left=True)
gg.set(xlabel='$Z_{DR}$ [dB]')
gg.set(xlabel=None)
ax = plt.gca()
ax.text(.5, 11, '66', transform=ax.transAxes)

# ax.text(.5, 14, short_name, transform=ax.transAxes)
# # Create the data
# temp_rd=(np.round(y/2,0)*2)
# df = pd.DataFrame(dict(mom=x, temp=temp_rd))
#
# # Initialize the FacetGrid object
# pal = sns.cubehelix_palette(19, rot=-.25, light=.7)
# gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15,
#                    height=.35, palette=pal,ylim=[0,3], xlim=[-1,5.5]) # zdr obs
#
# # Draw the densities in a few steps
# gg.map(sns.kdeplot, "mom",
#       bw_adjust=.5, clip_on=False,
#       fill=True, alpha=1, linewidth=1.5)
# gg.map(sns.kdeplot, "mom", clip_on=False, color="lightgrey", lw=2, bw_adjust=.5)
#
# # passing color=None to refline() uses the hue mapping
# gg.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
#
# # Define and use a simple function to label the plot in axes coordinates
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color=color,
#             ha="left", va="center", transform=ax.transAxes)
#
# gg.map(label, "mom")
#
# # Set the subplots to overlap
# gg.figure.subplots_adjust(hspace=-.25)
#
# # Remove axes details that don't play well with overlap
# gg.set_titles("")
# gg.set(yticks=[], ylabel="")
# gg.despine(bottom=True, left=True)
# gg.set(xlabel='$Z_{DR}$ [dB]')
# ax = plt.gca()
# ax.text(.5, 7, letters[letters_i] + ') obs',  transform=ax.transAxes)
# letters_i=letters_i+1
# plt.savefig(
#     folder_plot +
#     '/Ridge_ZDR_' +
#     str(elevation_deg) + '°_' + dates_str + locations_str +
#     ['', 'entr_'][filter_entr] +
#     ['', str(filter_entr_at)+'_'][filter_entr] +
#     ['', 'mom_'][filter_moms] + 'obs_' +
#     ['', '_2e14'][testing] +
#     '.pdf', format='pdf', transparent=True, bbox_inches='tight')
plt.savefig(
    folder_plot +
    '/Ridge_ZDR_below_ML_' +
    str(elevation_deg) + '°_' + dates_str + locations_str +
    ['', 'entr_'][filter_entr] +
    ['', str(filter_entr_at)+'_'][filter_entr] +
    ['', 'mom_'][filter_moms] + 'obs_' +
    ['', '_2e14'][testing] +
    '.png', format='png', dpi=300, bbox_inches='tight')

plt.close()

# --------------------------------------------------------------------------- #
# CFTDs CBAND SYN row i                                                       #
# --------------------------------------------------------------------------- #
for da_run, icon_emvorado_run, spin_up_mm, color, short_name in zip(
        da_runs, icon_emvorado_runs, spin_up_mms, colors, short_names):
    da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
    model_name_file = '-'.join([icon_emvorado_run.split('/')[0][9:],
                                icon_emvorado_run.split('/')[1][5:]])
    mod_names = '_'.join([mod_names, model_name_file])
    model_name = '-'.join([da_run[4:],
                           icon_emvorado_run.split('/')[0][5:],
                           icon_emvorado_run.split('/')[1][5:],
                           spin_up_mm + 'min'])
    # ------------------------------------ #
    ax = None
    x ,y = plot_CFAD_or_CFTD_from_QVP_with_list(
        locations=locations,
        dates=dates,
        hhmm_start=hhmm_start_cftds,
        hhmm_end=hhmm_end_cftds,
        elevation_deg=elevation_deg,
        da_icon_emvorado_run=da_icon_emvorado_run,
        spin_up_mm=spin_up_mm,
        moment='zdrsim',
        vert_temp=vert_temp,
        temp_min=temp_min,
        temp_max=temp_max,
        bins_temp=bins_temp,
        height_min=height_min,  # in km
        height_max=height_max,  # in km
        bins_height=bins_height,
        filter_entr=filter_entr,
        filter_entr_at=filter_entr_at,
        filter_moms=filter_moms,
        ax=ax,
        save=False,
        color=color,
        plot_legend=True,
        plot_data=False,
        data_max=data_max,
    )
    plt.close()
    # ------------------------------------ #

    # Create the data
    temp_rd=(np.round(y/2,0)*2)
    df = pd.DataFrame(dict(mom=x, temp=temp_rd))

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(19, rot=-.25, light=.7)
    gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15,
                   height=.4, palette=pal,ylim=[0,.6], xlim=[-1.99,6.1]) # zdr obs

    # Draw the densities in a few steps
    gg.map(sns.kdeplot, "mom",
          bw_adjust=.5, clip_on=False,
          fill=True, alpha=1, linewidth=1.5)
    gg.map(sns.kdeplot, "mom", clip_on=False, color="lightgrey", lw=2, bw_adjust=.5)

    # passing color=None to refline() uses the hue mapping
    gg.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    for temp in np.sort(np.unique(temp_rd)):
        print(temp)
        print(wasserstein_distance(x[temp_rd==temp], x_obs[temp_obs==temp]))

    print('___________________')
    gg.map(label, "mom")

    # Set the subplots to overlap
    gg.figure.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    gg.set_titles("")
    gg.set(yticks=[], ylabel="")
    gg.despine(bottom=True, left=True)
    gg.set(xlabel='$Z_{DR}$ [dB]')
    gg.set(xlabel=None)
    ax = plt.gca()
    ax.text(.5, 11, '66', transform=ax.transAxes)
    letters_i = letters_i + 1
    # ax.text(.5, 7, letters[letters_i] + ') '+ short_name, transform=ax.transAxes)
    # letters_i = letters_i + 1

    # plt.savefig(
    #     folder_plot +
    #     '/Ridge_ZDR_' +
    #     str(elevation_deg) + '°_' + dates_str + locations_str +
    #     ['', 'entr_'][filter_entr] +
    #     ['', str(filter_entr_at)+'_'][filter_entr] +
    #     ['', 'mom_'][filter_moms] + model_name + '_' +
    #     ['', '_2e14'][testing] +
    #     '.pdf', format='pdf', transparent=True, bbox_inches='tight')
    plt.savefig(
        folder_plot +
        '/Ridge_ZDR_below_ML_' +
        str(elevation_deg) + '°_' + dates_str + locations_str +
        ['', 'entr_'][filter_entr] +
        ['', str(filter_entr_at)+'_'][filter_entr] +
        ['', 'mom_'][filter_moms] + model_name + '_' +
        ['', '_2e14'][testing] +
        '.png', format='png', dpi=300,  bbox_inches='tight')

    plt.close()

