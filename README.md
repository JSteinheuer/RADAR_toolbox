# RADAR_toolbox

**Julian Steinheuer**

This repository contains code to process RADAR data.

**How to use:**

- Capital letter scripts _RADAR_SCRIPT.py_ shall contain the functions
- similar named scripts in lower case shall use that functions for specific cases with some name extensions for hints, as _radar_script_ess_210714_pcp.py_
- hard coded paths which are regular used shall be set in _HEADER_RADAR_toolbox.py_
- workflows, i.e. for processing radar observations step by steps, have different scripts with numbers indicating the step, as _radar_script_1_load_all_moments.py_, *radar_script_2_correct_rho_hv.py*

<br>
**Workflows:**

_DWD_OBS_TO_MIUB_OBS.py_
<br>
take raw DWD Cband radar data and process them in 'miub style' towards rho_hv corrected, phi_dp smoothed and offset corrected, k_dp smoothed, zdr calibrated, attenuation corrected with help of ERA5 temperatures, combined polarimetric moments netcdfs.

_ERA5....py_
<br>
download ERA5 temperatures for broad Germany (Radolan domain) from ECMWF climate data store

_PLOT_RADAR.py_
<br>
functions to plot the real DWD radar observations: PPIs, pseudoRHIs, CFTDs, CFADs, QVPs

_PLOT_SYN_RADAR.py_
<br>
functions to plot the synthetic ICON/EMVORADO radar data: PPIs, pseudoRHIs, CFTDs, CFADs, QVPs

_PROCESS_SYN_RADAR.py_
<br>
functions to process ready synthetic radar data: QVPs,

_SET_SYN_RADAR.py_
<br>
rearrange synthetic radar data from EMVORADO and ICON



<br>
**some Side-scripts**

_change_folder_structure.py_
<br>
mv and cp commands moving TS into JS structure

_rename_ass_icon_emv.py_
<br>
 new convention from JM: meaningful names for the old runs

_Z_scrap_sheet.py_
<br>
play around!


<br>
_July 11, 2024_

This program is a free software distributed under the terms of the GNU General 
Public License as published by the Free Software Foundation, version 3 
(GNU-GPLv3). You can redistribute and/or modify by citing the mentioned 
publication, but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
