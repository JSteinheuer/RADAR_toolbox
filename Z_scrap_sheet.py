#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

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
import HEADER_RADAR_toolbox as header
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import numpy as np
from PLOT_SYN_RADAR import plot_qvp_of_polarimetric_variable
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP_with_list
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
location = 'ESS'
date = '20210714'
locations = [location]
dates = [date]
hhmm_start = '00:00'
hhmm_end = '21:00'
elevation_deg = 12
top_height = 8
da_runs = []
icon_emvorado_runs = []
spin_up_mms = []
# ------------------------------------ #
filter_entr = True  # TODO  # only for CFTDS
filter_moms = False  # TODO  # only for CFTDS
# ------------------------------------ #
# not for paper!:
# ------------------------------------ #
# filter_entr = False  # TODO  # only for CFTDS
# filter_moms = True  # TODO  # only for CFTDS
# ------------------------------------ #
# CFADs ?
# vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
# or CFTDs ?
vert_temp = True
temp_min = -20
temp_max = 16
bins_temp = 18
# ------------------------------------ #
# SYN data row 1                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00010000.2')
spin_up_mms.append('120')
# ------------------------------------ #
# SYN data row 2                       #
# ------------------------------------ #
# # TODO: missing
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00410000.2')
spin_up_mms.append('120')
# ------------------------------------ #
# SYN data row 3                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.0/EMVO_00510000.2')
spin_up_mms.append('120')
# ------------------------------------ #
# SYN data row 4                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.1/EMVO_00510000.2')
spin_up_mms.append('120')
# # ------------------------------------ #
# SYN data row 5                       #
# ------------------------------------ #
da_runs.append('ASS_2411')
icon_emvorado_runs.append('MAIN_2411.3/EMVO_00510000.2')
spin_up_mms.append('120')

# ------------------------------------ #
year = date[0:4]
mon = date[4:6]
day = date[6:8]
date_start = '-'.join([year, mon, day, hhmm_start])
date_end = '-'.join([year, mon, day, hhmm_end])
mode = 'vol'
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# INITIALIZATION                                                              #
# --------------------------------------------------------------------------- #
locations = list(rad_dict().keys())
dates = ['20210714', '20210713']
hhmm_start = '00:00'
hhmm_end = '23:55'
elevation_deg = 12
top_height = 8
# ------------------------------------ #

# --------------------------------------------------------------------------- #

mode = 'vol'
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
folder_plot = header.folder_plot + 'CFADs'
mod_names = ''
n_rows = 1
n_cols = 2
plt.figure(figsize=(20, 10))  # TODO
n_i = 0
current_row = 0
current_col = 0

# --------------------------------------------------------------------------- #
# CBAND OBS row 1                                                             #
# --------------------------------------------------------------------------- #
# current_row = 1
# current_col = 1
# n_i = (current_row - 1) * n_cols + current_col
# print(n_i)
# ax = plt.subplot(n_rows, n_cols, n_i)
# mom ,temp = plot_CFAD_or_CFTD_from_QVP_with_list(
#     locations=locations,
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     elevation_deg=elevation_deg,
#     da_icon_emvorado_run=None,
#     moment='ZDR_AC_OC',
#     vert_temp=vert_temp,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     filter_entr=filter_entr,
#     filter_moms=filter_moms,
#     ax=ax,
#     save=False,
# )
# # ------------------------------------ #

da_run, icon_emvorado_run, spin_up_mm = (
    da_runs[-1], icon_emvorado_runs[-1], spin_up_mms[-1])
da_icon_emvorado_run = da_run + '/' + icon_emvorado_run
current_row = 1
current_col = 1
n_i = (current_row - 1) * n_cols + current_col
print(n_i)
ax = plt.subplot(n_rows, n_cols, n_i)
mom ,temp = plot_CFAD_or_CFTD_from_QVP_with_list(
    locations=locations,
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
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
    filter_moms=filter_moms,
    ax=ax,
    save=False,
)
# ------------------------------------ #

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Create the data
temp_rd=(np.round(temp/2,0)*2)
df = pd.DataFrame(dict(mom=mom, temp=temp_rd))

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(19, rot=-.25, light=.7)
# pal = sns.cubehelix_palette(19, rot=.5, light=.7)
# gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15, height=.5, palette=pal)
# gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15, height=.5, palette=pal,ylim=[0,0.05]) # zh obs
gg = sns.FacetGrid(df, row="temp", hue="temp", aspect=15, height=.5, palette=pal,ylim=[0,1]) # zdr obs

# Draw the densities in a few steps
gg.map(sns.kdeplot, "mom",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
gg.map(sns.kdeplot, "mom", clip_on=False, color="w", lw=2, bw_adjust=.5)

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
gg.set_titles("")
gg.set(yticks=[], ylabel="")
gg.despine(bottom=True, left=True)


###



# import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#
# # Create the data
# rs = np.random.RandomState(1979)
# x = rs.randn(500)
# g = np.tile(list("ABCDEFGHIJ"), 50)
# df = pd.DataFrame(dict(x=x, g=g))
# m = df.g.map(ord)
# df["x"] += m
#
# # Initialize the FacetGrid object
# pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
# g = sns.FacetGrid(df, row="g", hue="g", aspect=15, height=.5, palette=pal)
#
# # Draw the densities in a few steps
# g.map(sns.kdeplot, "x",
#       bw_adjust=.5, clip_on=False,
#       fill=True, alpha=1, linewidth=1.5)
# g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw_adjust=.5)
#
# # passing color=None to refline() uses the hue mapping
# g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
#
#
# # Define and use a simple function to label the plot in axes coordinates
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color=color,
#             ha="left", va="center", transform=ax.transAxes)
#
#
# g.map(label, "x")
#
# # Set the subplots to overlap
# g.figure.subplots_adjust(hspace=-.25)
#
# # Remove axes details that don't play well with overlap
# g.set_titles("")
# g.set(yticks=[], ylabel="")
# g.despine(bottom=True, left=True)



