#!/usr/bin/env python3.11
# #!/automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python3.11

# --------------------------------------------------------------------------- #
# Julian Steinheuer; 24.01.24                                                 #
# download_ERA5_data_DE.py                                                    #
#                                                                             #
# [...] Description here [...]                                                #
# --------------------------------------------------------------------------- #

import HEADER_RADAR_toolbox as header
import cdsapi
import os

# --------------------------------------------------------------------------- #
# from ECMWF: compute_ge
# --------------------------------------------------------------------------- #

DATES = ["20210604",  # case01
         "20210620", "20210621",  # case02
         "20210628", "20210629",  # case03
         "20220519", "20220520",  # case04
         "20220623", "20220624", "20220625",  # case05
         "20220626", "20220627", "20220628",  # case06+07
         "20220630", "20220701",  # case08
         "20210714",  # case09
         "20221222",  # case10
         ]
dir_out = header.dir_data_era5
overwrite = False
# date = DATES[0]
for date in DATES:

    year = date[0:4]
    mon = date[4:6]
    day = date[6:8]

    c = cdsapi.Client()
    file_out1 = dir_out + str(year) + str(mon) + str(day) + "-3D-T-q-ml.grib"
    if not overwrite and os.path.exists(file_out1):
        print('exists: ' + file_out1 + ' -> continue')
    elif not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out1.replace('grib', 'nc' + ' -> continue'))
    else:
        c.retrieve("reanalysis-era5-complete", {
            "class": "ea",
            "date": "%s-%s-%s" % (year, mon, day),
            "expver": "1",
            "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/"
                        "21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/"
                        "38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/"
                        "55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/"
                        "72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/"
                        "89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/"
                        "104/105/106/107/108/109/110/111/112/113/114/115/"
                        "116/117/118/119/120/121/122/123/124/125/126/127/"
                        "128/129/130/131/132/133/134/135/136/137",
            "levtype": "ml",
            "param": "130/133",  # Temperature (t) and specific humidity (q)
            "stream": "oper",
            "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/"
                    "06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/"
                    "12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/"
                    "18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type": "an",
            "grid": [0.25, 0.25],
            'area': [56, 2, 46.25, 17],
        }, file_out1)
    file_out2 = dir_out + str(year) + str(mon) + str(day) + \
                "-2D-z-lnsp-ml.grib"
    if not overwrite and os.path.exists(file_out2):
        print('exists: ' + file_out2 + ' -> continue')
    elif not overwrite and os.path.exists(file_out2.replace('grib', 'nc')):
        print('exists: ' + file_out2.replace('grib', 'nc' + ' -> continue'))
    else:
        c.retrieve("reanalysis-era5-complete", {
            "class": "ea",
            "date": "%s-%s-%s" % (year, mon, day),
            "expver": "1",
            "levelist": "1",
            "levtype": "ml",
            "param": "129/152",  # Geopotential (z) and ln of s. pres. (lnsp)
            "stream": "oper",
            "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/"
                    "06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/"
                    "12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/"
                    "18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type": "an",
            "grid": [0.25, 0.25],
            'area': [56, 2, 46.25, 17],
        }, file_out2)

    file_out3 = file_out1.replace('T-q-ml', 'z')
    if not overwrite and os.path.exists(file_out3):
        print('exists: ' + file_out3 + ' -> continue')
    elif not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out3.replace('grib', 'nc' + ' -> continue'))
    else:
        os.system('python compute_geopotential_on_ml.py ' + file_out1 +
                  ' ' + file_out2 + ' -o ' + file_out3)

    if not overwrite and os.path.exists(file_out1.replace('grib', 'nc')):
        print('exists: ' + file_out1.replace('grib', 'nc') + ' -> continue')
    else:
        os.system('cdo -f nc copy ' + file_out1 + ' ' +
                  file_out1.replace('grib', 'nc'))
        # os.system('rm ' + file_out1)

    # if not overwrite and os.path.exists(file_out2.replace('grib', 'nc')):
    #     print(
    #         'exists: ' + file_out2.replace('grib', 'nc') + ' -> continue')
    # else:
    #     os.system('cdo -f nc copy ' + file_out2 + ' ' +
    #               file_out2.replace('grib', 'nc'))
    #     # os.system('rm ' + file_out2)

    if not overwrite and os.path.exists(file_out3.replace('grib', 'nc')):
        print(
            'exists: ' + file_out3.replace('grib', 'nc') + ' -> continue')
    else:
        os.system('cdo -f nc copy ' + file_out3 + ' ' +
                  file_out3.replace('grib', 'nc'))
        # os.system('rm ' + file_out3)

    # os.system('cdo -f nc copy ' + file_out + ' ' +
    #           file_out.replace('grib', 'nc'))
    # os.system('rm ' + file_out)

    # file_out4 = dir_out + str(year) + str(mon) + str(day) + "-2D-SP-T2m.nc"
    # if not overwrite and os.path.exists(file_out4):
    #     print('exists: ' + file_out4 + ' -> continue')
    # else:
    #     c.retrieve('reanalysis-era5-single-levels', {
    #         'product_type': 'reanalysis',
    #         'variable': ['2m_temperature', 'surface_pressure',
    #                      'geopotential'],  # TODO: new
    #         # 'variable': ['2m_temperature', 'surface_pressure', ],  # old
    #         'year': year,
    #         'month': mon,
    #         'day': day,
    #         'time': ['00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
    #                  '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
    #                  '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
    #                  '18:00', '19:00', '20:00', '21:00', '22:00', '23:00', ],
    #         'area': [56, 2, 46.25, 17],
    #         'format': 'netcdf',
    #     }, file_out4)

