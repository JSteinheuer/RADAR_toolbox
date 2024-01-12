#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 10.01.24                                                 #
# DWD_obs_to_MIUB_obs.py                                                      #
#                                                                             #
# Processing script to quality check, calibrate, and correct the DWD C-band   #
# observations towards MIUB 'standard'.                                       #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import sys
sys.path.insert(0, header.dir_projects +
                'RADAR_toolbox/radar_processing_scripts/')
# maybe awkward because >import utils< would work now, but the following
# additionally enables  pycharm to word-completing (i.e., functions) since
# the folder containing /radar_processing_scripts/ was already in the (global)
# sys.path variable (i.e. the project folder).
from radar_processing_scripts import utils
import warnings
warnings.filterwarnings('ignore')
# sys.path.insert(0, '../')
from radar_processing_scripts.radarmet import *
from radar_processing_scripts.read_cband_dwd import *

DATES = ["20170727", "20180923", "20170725", "20170719", "20170724", "20170726",
        "20190623", "20190720", "20180728", "20180809", "20190829", "20190830",
        "20181202", "20170720",
        ]
LOCATIONS = ['boo', 'eis', 'fld', 'mem', 'neu', 'ros', 'tur', 'umd', 'drs',
             'ess', 'fbg', 'hnr', 'isn', 'nhb', 'oft', 'pro'
             ]
# Info on the DWD scanning strategy:
# y https://www.dwd.de/EN/ourservices/radar_products/radar_products.html
# Scans 00-05 are the volume scans (5.5°, 4.5°, 3.5°, 2.5°, 1.5° and 0.5°),
# the rest are 8.0°, 12.0°, 17.0° and 25.0°
SCANS = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
# moments = ["DBZH", "DBZV", "TH", "TV", "ZDR", "UZDR",
#            "ZDR1", "UZDR1", "VRADH", "VRADV", "UVRADH", "UVRADV",
#            "FVRADH", "UFVRADH", "WRADH", "UWRADH", "UPHIDP", "KDP",
#            "RHOHV", "URHOHV", "SQIH", "SQIV", "SQI2H", "SQI2V",
#            "SQI3H", "SQI3V", "CPAH", "CPAV", "STDH", "STDV",
#            "CCORH", "CCORV", "SNRHC", "SNRVC",
#            # "DBSNRH", "DBSNRV",
#            "CMAP", "CFLAGS",
#            ]
moments = None

fmt = "%Y%m%d"

date = DATES[2]
loc = LOCATIONS[-1]
scan = SCANS[7]
case_studies = {cs: dt.datetime.strptime(cs, fmt) for cs in DATES}
mode = "vol5minng01"  # "vol5minng01" 'pcpng01'
path = header.dir_data_obs_realpep

date = "20210604"
loc = "pro"
case = "01"
scan = SCANS[7]
case_studies = {cs: dt.datetime.strptime(cs, fmt) for cs in [date]}
mode = "vol5minng10"
path = header.dir_data_obs + "OpHymet2-case" + case + "-" + date + "/"

# Start and End time
starttime = case_studies[date] + dt.timedelta(hours=0, minutes=0)
endtime = case_studies[date] + dt.timedelta(hours=24, minutes=0)

# Radardata filelist
file_list_gen = create_dwd_filelist(path=path,
                                    starttime=starttime,
                                    endtime=endtime,
                                    moments=moments,
                                    mode=mode,
                                    loc=loc,
                                    scan=scan,
                                    )
file_list = list(file_list_gen)

vol = wrl.io.open_odim(file_list, loader="h5py", chunks={})
swp = vol[0].data.pipe(wrl.georef.georeference_dataset)
swp.to_netcdf(header.dir_data_obs + '/test3.nc')
