#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 25.03.24                                                 #
# plot_syn_RADAR_CFTDs.py                                                     #
#                                                                             #
# Run the function in PLOT_SYN_RADAR.py for generating specific CFTD.         #
# --------------------------------------------------------------------------- #

# TODO: change plot_CFAD_or_CFTD_from_QVP towards:
# TODO:     - include entropy filter

import HEADER_RADAR_toolbox as header
import matplotlib.pyplot as plt
import numpy as np
import warnings
from pathlib import Path
import glob

warnings.simplefilter('ignore')
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# PARAMS
n_rows = 2
n_cols = 2
plt.figure(figsize=(n_cols * 4.5, n_rows * 3.6))
n_i = 0
current_row = 0
current_col = 0

mod_names = ''

# All
dates = ['20210714', '20210713']
hhmm_start = '00:00'
hhmm_end = '23:55'

# vert_temp = False
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20

vert_temp = True
temp_min = 0
temp_max = 16
bins_temp = 8

# adjust ###################### adjust ###################### adjust ##########
# filter_entr_ML = True  # as email 15.2.24
filter_entr_ML = False

# locations = list(rad_dict().keys())
locations = list(rad_dict().keys())
# locations = ['ESS', 'NHB', 'OFT', 'FBG', 'FLD']
# adjust ###################### adjust ###################### adjust ##########

mode = 'vol'
elevation_deg = 12
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
vmax = 15

# MOMENTS
moment_1o = 'D0_r_Bringi_obs'
moment_1s = 'D0_r_Bringi_syn'
mom_min_1 = 0
mom_max_1 = 2.5
bins_mom_1 = 40

moment_2s = 'D0_r'
mom_min_2 = 0
mom_max_2 = 2.5
bins_mom_2 = 40

moment_3s = 'D0_g'
mom_min_3 = 0
mom_max_3 = 2.5
bins_mom_3 = 40


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Obs row 1
years = [date[0:4] for date in dates]
mons = [date[4:6] for date in dates]
days = [date[6:8] for date in dates]
folder_ins = ["/".join([header.dir_obs_qvp + '*', year,
                          year + '-' + mon,
                          year + '-' + mon + '-' + day,
                          location.lower(), mode + '*', sweep])
              for year, mon, day in zip(years, mons, days)
              for location in locations]

nc_file_mom = [glob.glob(folder_in + '/ras*qvp*_polmoms_nc_*')
               for folder_in in  folder_ins]
paths_in = []
for qvp_dateloc in nc_file_mom:
    paths_in=paths_in + qvp_dateloc
    if qvp_dateloc.__len__() > 1:
        print('ERROR: qvp ambiguous:' + qvp_dateloc[0] + qvp_dateloc[1])

title = 'C-band observation'
current_row = current_row + 1
# --------------------------------------------------------------------------- #
# Obs 1
current_col = 1
current_row = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_1o,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_1,
    mom_min=mom_min_1,
    bins_mom=bins_mom_1,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# CBAND SYN 1 row 2
da_run = 'ASS_2411'
icon_emvorado_run = 'MAIN_2411.3/EMVO_00510000.2/'
spin_up_mm = '120'
folder_syn = header.dir_data_qvp
current_row = current_row + 1
model_name = '-'.join([da_run[4:],
                        icon_emvorado_run.split('/')[0][5:],
                        icon_emvorado_run.split('/')[1][5:],
                        spin_up_mm + 'min'])
mod_names = '-'.join([mod_names, model_name])
# --------------------------------------------------------------------------- #
# Syn 1
current_col = 1
current_row = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_1s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_1,
    mom_min=mom_min_1,
    bins_mom=bins_mom_1,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)
# Syn 2
current_col = 2
current_row = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    da_run=da_run,
    icon_emvorado_run=icon_emvorado_run,
    spin_up_mm=spin_up_mm,
    moment=moment_2s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_2,
    mom_min=mom_min_2,
    bins_mom=bins_mom_2,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)
# Syn 3
# current_col = 2
# current_row = 1
# plot_CFAD_or_CFTD_from_QVP(
#     dates=dates,
#     hhmm_start=hhmm_start,
#     hhmm_end=hhmm_end,
#     locations=locations,
#     elevation_deg=elevation_deg,
#     folder_syn=folder_syn,
#     da_run=da_run,
#     icon_emvorado_run=icon_emvorado_run,
#     spin_up_mm=spin_up_mm,
#     moment=moment_3s,
#     vert_temp=vert_temp,  # CFTD
#     mom_max=mom_max_3,
#     mom_min=mom_min_3,
#     bins_mom=bins_mom_3,
#     temp_min=temp_min,
#     temp_max=temp_max,
#     bins_temp=bins_temp,
#     height_min=height_min,  # in km
#     height_max=height_max,  # in km
#     bins_height=bins_height,
#     vmax=vmax,
#     filter_entr_ML=filter_entr_ML,
#     ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
# )
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

save_path = header.folder_plot + 'CFADs/'
Path(save_path).mkdir(parents=True, exist_ok=True)

locs_name = '_all'
if  locations.__len__()<16:
    locs_name= '_'.join(locations)

if vert_temp:
    save_name = 'PrePEP_CFTD_'
else:
    save_name = 'PrePEP_CFAD_'
plt.savefig(
    save_path + save_name + str(elevation_deg) + '_1314_07_21' +
    '_' + hhmm_start + '-' + hhmm_end + locs_name +
    mod_names +
    ['', '_filtered'][filter_entr_ML] +
    '_rain.pdf', format='pdf', transparent=True)
plt.close()
