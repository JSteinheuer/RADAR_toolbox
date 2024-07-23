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
import glob
import matplotlib.pyplot as plt
import numpy as np
import warnings
from pathlib import Path

warnings.simplefilter('ignore')
from PLOT_SYN_RADAR import plot_CFAD_or_CFTD_from_QVP
from SET_SYN_RADAR import rad_dict

# --------------------------------------------------------------------------- #
# PARAMS

n_rows = 4
n_cols = 4
plt.figure(figsize=(n_cols * 5, n_rows * 4))
n_i = 0
current_row = 0
current_col = 0
mod_names = ''
filter_entr_ML = True
filter_entr_ML = False
filter_ML = False

# All
# date = '20170725'
dates = [
    '20181223',
    '20181224',
]
hhmm_start = '00:00'
hhmm_end = '23:55'
vert_temp = True
temp_min = -28
temp_max = 8
bins_temp = 18
height_min = 0  # in km
height_max = 10  # in km
bins_height = 20
locations = list(rad_dict().keys())
elevation_deg = 12
vmax = 15

# MOMENTS
moment_1s = 'zrsim'
moment_1o = 'ZH_AC'
mom_min_1 = 0
mom_max_1 = 40
bins_mom_1 = 40

moment_2s = 'zdrsim'
moment_2o = 'ZDR_AC_OC'
mom_min_2 = -0.2
mom_max_2 = 2.2
bins_mom_2 = 40

moment_3s = 'kdpsim_ml_corrected'
moment_3o = 'KDP_NC'
mom_min_3 = -.1
mom_max_3 = 0.3
bins_mom_3 = 40

moment_4s = 'rhvsim'
moment_4o = 'RHOHV_NC2P'
mom_min_4 = 0.799
mom_max_4 = 1.051
bins_mom_4 = 30

# Obs row 1
# #TODO: hardcoded more or les!
# folder_obs = header.dir_obs_qvp + 'OpHymet2-caseX-20181223/2018/2018-12/*/*/vol*/07/*qvp*07*.hd5'
# files = glob.glob(folder_obs)
sweep = '0' + str(np.where(header.ELEVATIONS_ALL ==
                           float(elevation_deg))[0][0])
folder_obs = [header.dir_obs_qvp + 'OpHymet*/' + date[0:4] + '/' +
              date[0:4] + '-' + date[4:6] + '/' + date[0:4] + '-' +
              date[4:6] + '-' + date[6:8] + '/' + location.lower() + '/vol*/' +
              sweep + '/*qvp*07*.hd5'
              for location in locations for date in dates]
files_list = [glob.glob(folder) for folder in folder_obs]
for j in [i for i, x in enumerate(files_list) if not x]:
    files_list[j] = ['dummy']

paths_in = [x for xs in files_list for x in xs]

title = 'C-band observation'
current_row = current_row + 1

# Obs 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_1o,
    # moment=moment_1s,
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
    filter_ML=filter_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 2
current_col = 2
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_2o,
    # moment=moment_2s,
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
    filter_ML=filter_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 3
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_3o,
    # moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    filter_ML=filter_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# Obs 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    locations=locations,
    elevation_deg=elevation_deg,
    paths_in=paths_in,
    title=title,
    moment=moment_4o,
    # moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
    temp_min=temp_min,
    temp_max=temp_max,
    bins_temp=bins_temp,
    height_min=height_min,  # in km
    height_max=height_max,  # in km
    bins_height=bins_height,
    vmax=vmax,
    filter_entr_ML=filter_entr_ML,
    filter_ML=filter_ML,
    ax=plt.subplot(n_rows, n_cols, (current_row - 1) * n_cols + current_col),
)

# --------------------------------------------------------------------------- #
# CBAND SYN: 3

da_run1 = 'ASS_2407'
icon_emvorado_run1 = 'MAIN_2405.3/EMVO_00510000.2'

da_run2 = 'ASS_2407'
icon_emvorado_run2 = 'MAIN_2405.3/EMVO_00510200.2'

da_run3 = 'ASS_2407'
icon_emvorado_run3 = 'MAIN_2405.3/EMVO_00510300.2'

spin_up_mm = '120'
folder_syn = header.dir_data_qvp
current_row = current_row + 1

# --------------------------------------------------------------------------- #
# CBAND SYN 1 row 3

model_name = '-'.join([da_run1[4:],# adjust
                       icon_emvorado_run1.split('/')[0][5:],
                       icon_emvorado_run1.split('/')[1][5:],
                       spin_up_mm + 'min'])
path_in = ['/'.join([header.dir_data_qvp + date,
                     da_run1, icon_emvorado_run1,  # adjust
                     str(spin_up_mm) + 'min_spinup', 'QVP_' +
                     str(elevation_deg) + '_Syn_' + location + '_' +
                     date + '*' + date + '*.nc'])
           for location in locations for date in dates]
# -------------------------------------- #
mod_names = '-'.join([mod_names, model_name])
files_list = [glob.glob(folder) for folder in path_in]
for j in [i for i, x in enumerate(files_list) if not x]:
    files_list[j] = ['dummy']

paths_in = [x for xs in files_list for x in xs]
dates_in = [path[-28:-20] for path in paths_in]
locations_in = [path[-32:-29] for path in paths_in]
locations_now = locations
dates_now = ['20181223', '20181224', '20181224', ]
# --------------------------------------------------------------------------- #
# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    spin_up_mm=spin_up_mm,
    title=model_name,
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
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
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
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
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

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
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

current_row = current_row + 1

## --------------------------------------------------------------------------- #
# CBAND SYN 2 row 3

model_name = '-'.join([da_run2[4:],# adjust
                       icon_emvorado_run2.split('/')[0][5:],
                       icon_emvorado_run2.split('/')[1][5:],
                       spin_up_mm + 'min'])
path_in = ['/'.join([header.dir_data_qvp + date,
                     da_run2, icon_emvorado_run2,  # adjust
                     str(spin_up_mm) + 'min_spinup', 'QVP_' +
                     str(elevation_deg) + '_Syn_' + location + '_' +
                     date + '*' + date + '*.nc'])
           for location in locations for date in dates]
# -------------------------------------- #
mod_names = '-'.join([mod_names, model_name])
files_list = [glob.glob(folder) for folder in path_in]
for j in [i for i, x in enumerate(files_list) if not x]:
    files_list[j] = ['dummy']

paths_in = [x for xs in files_list for x in xs]
dates_in = [path[-28:-20] for path in paths_in]
locations_in = [path[-32:-29] for path in paths_in]
locations_now = locations
dates_now = ['20181223', '20181224', '20181224', ]
# --------------------------------------------------------------------------- #
# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    spin_up_mm=spin_up_mm,
    title=model_name,
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
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
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
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
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

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
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

current_row = current_row + 1

# --------------------------------------------------------------------------- #
# CBAND SYN 3 row 3

model_name = '-'.join([da_run3[4:],# adjust
                       icon_emvorado_run3.split('/')[0][5:],
                       icon_emvorado_run3.split('/')[1][5:],
                       spin_up_mm + 'min'])
path_in = ['/'.join([header.dir_data_qvp + date,
                     da_run3, icon_emvorado_run3,  # adjust
                     str(spin_up_mm) + 'min_spinup', 'QVP_' +
                     str(elevation_deg) + '_Syn_' + location + '_' +
                     date + '*' + date + '*.nc'])
           for location in locations for date in dates]
# -------------------------------------- #
mod_names = '-'.join([mod_names, model_name])
files_list = [glob.glob(folder) for folder in path_in]
for j in [i for i, x in enumerate(files_list) if not x]:
    files_list[j] = ['dummy']

paths_in = [x for xs in files_list for x in xs]
dates_in = [path[-28:-20] for path in paths_in]
locations_in = [path[-32:-29] for path in paths_in]
locations_now = locations
dates_now = ['20181223', '20181224', '20181224', ]
# --------------------------------------------------------------------------- #
# Syn 1
current_col = 1
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    spin_up_mm=spin_up_mm,
    title=model_name,
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
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
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
current_col = 3
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_3s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_3,
    mom_min=mom_min_3,
    bins_mom=bins_mom_3,
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

# Syn 4
current_col = 4
plot_CFAD_or_CFTD_from_QVP(
    dates=dates_now,
    hhmm_start=hhmm_start,
    hhmm_end=hhmm_end,
    paths_in=paths_in,
    locations=locations_now,
    elevation_deg=elevation_deg,
    folder_syn=folder_syn,
    title=model_name,
    spin_up_mm=spin_up_mm,
    moment=moment_4s,
    vert_temp=vert_temp,  # CFTD
    mom_max=mom_max_4,
    mom_min=mom_min_4,
    bins_mom=bins_mom_4,
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

current_row = current_row + 1

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# SAVE                                                                        #
# --------------------------------------------------------------------------- #

save_path = header.folder_plot + 'CFADs/'
Path(save_path).mkdir(parents=True, exist_ok=True)

if vert_temp:
    save_name = 'CFTD_'
else:
    save_name = 'CFAD_'

plt.savefig(
    save_path + save_name + str(elevation_deg) + '_all' +
    '_' + dates[0] + '-' + dates[1] + '_all' +
    mod_names +
    ['', '_filtered'][filter_entr_ML]
    + '.pdf', format='pdf', transparent=True)
plt.close()
