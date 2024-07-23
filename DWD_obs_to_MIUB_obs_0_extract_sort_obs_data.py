#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.
import sys
import warnings

warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", UserWarning)

import datetime as dt
import os
import glob
import argparse
import tarfile


def extract_file_name(fname, plist):
    parts = fname.split("-")
    sweep = parts[1].split("_")
    nparts = [parts[0]] + sweep + parts[2:]

    fdict = {}
    for i, key in enumerate(plist):
        fdict[key] = nparts[i]
    return fdict


def extract_file_name2(fname):
    fname=fname.split('/')[-1]  # TODO: new
    parts = fname.split("-")
    sweep = parts[1].split("_")
    nparts = [parts[0]] + sweep + parts[2:]
    return tuple(nparts)


def unpack_sort_dwd_obs(inpath, outpath, date_struct, log):
    tar_files = sorted(glob.glob(os.path.join(inpath, "*.tar*")))
    for tfile in tar_files:
        if tarfile.is_tarfile(tfile):
            th = tarfile.open(tfile)
            mem = th.getnames()
            for m in mem:
                (
                    dwd_type,
                    sweep_name,
                    sweep_type,
                    sweep_moments,
                    sweep_number,
                    dtime,
                    radar,
                    wmo,
                    ftype,
                ) = extract_file_name2(m)
                dtime = dt.datetime.strptime(dtime, "%Y%m%d%H%M%S00") # TODO: works for vol, pcp
                # dtime = dt.datetime.strptime(dtime, "%Y%m%d%H%M%S")  # TODO: for birdbath 240503
                date_path = dtime.strftime(date_struct)
                outp = os.path.join(
                    f"{date_path}",
                    f"{radar}",
                    f"{sweep_name}",
                    f"{sweep_number}",
                )
                outp = os.path.join(outpath, outp)
                if log:
                    print(m)
                if sys.version_info[1] == 12:
                    kwargs = dict(filter="data")
                else:
                    kwargs = {}
                th.extract(m, path=outp, **kwargs)
            th.close()

    # TODO: alle below is new
    tgz_files = sorted(glob.glob(os.path.join(inpath, "*.tgz*")))
    for tfile in tgz_files:
        th = tarfile.open(tfile, 'r:gz')
        mem = th.getnames()
        for m in mem:
            if len(m.split('/')[-1].split('_'))  < 4: # TODO: new
                continue

            (
                dwd_type,
                sweep_name,
                sweep_type,
                sweep_moments,
                sweep_number,
                dtime,
                radar,
                wmo,
                ftype,
            ) = extract_file_name2(m)
            dtime = dtime[::-1].zfill(16)[::-1]  # TODO: new!
            dtime = dt.datetime.strptime(dtime, "%Y%m%d%H%M%S00")
            date_path = dtime.strftime(date_struct)
            outp = os.path.join(
                f"{date_path}",
                f"{radar}",
                f"{sweep_name}",
                f"{sweep_number}",
            )
            outp = os.path.join(outpath, outp)
            if log:
                print(m)
            if sys.version_info[1] == 12:
                kwargs = dict(filter="data")
            else:
                kwargs = {}

            th.extract(m, path=outp, **kwargs)
            if m.split('/')[-1] != m.split('/')[:-1]:
                os.system('mv ' + outp + '/' + m + ' ' +
                          outp + '/' + m.split('/')[-1])
                os.system('rmdir ' + outp + '/' + m.split('/')[0])

        th.close()


def main():
    parser = argparse.ArgumentParser(
        description="Unpack and sort DWD obs data", add_help=True
    )

    d1group = parser.add_argument_group(description='specific date parameters')
    d1group.add_argument('-d', '--date-struct', required=True)

    iogroup = parser.add_argument_group(
        description="input and output folders, mandatory"
    )
    iogroup.add_argument("-i", "--input-folder", required=True)
    iogroup.add_argument("-o", "--output-folder", required=True)

    loggroup = parser.add_argument_group(
        description="logging"
    )

    loggroup.add_argument("-v", "--verbose", action="store_true")

    pargs = parser.parse_args()

    inpath = pargs.input_folder
    outpath = pargs.output_folder
    date_struct = pargs.date_struct
    log = pargs.verbose

    if not os.access(inpath, os.R_OK):
        raise IOError(f"No read access for input folder {inpath}.")

    if not os.access(outpath, os.W_OK):
        raise IOError(f"No write access for output folder {outpath}.")

    unpack_sort_dwd_obs(inpath, outpath, date_struct, log)


if __name__ == "__main__":
    main()


# 231126:
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210604 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case01-20210604 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210620 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case02-20210620 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210628 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case03-20210628 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220519 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220623 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case05-20220623 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220626 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case06+07-20220626 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220630 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case08-20220630 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /home/s6justei/mambaforge/envs/RADAR_toolbox/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20221222 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case10-20221222 -d "%Y/%Y-%m/%Y-%m-%d" -v

# 240503
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210604_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case01-20210604 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210620_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case02-20210620 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210628_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case03-20210628 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220519_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220623_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case05-20220623 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220626_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case06+07-20220626 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20220630_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case08-20220630 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20210714_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case09-20210714 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20221222_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case10-20221222 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/20170725_raw_birdbath -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-caseX-20170725 -d "%Y/%Y-%m/%Y-%m-%d" -v

# 240521
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20210604 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case01-20210604 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20210620 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case02-20210620 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20210628 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case03-20210628 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20220519 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20220623 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case05-20220623 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20220626 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case06+07-20220626 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20220630 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case08-20220630 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20210714 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case09-20210714 -d "%Y/%Y-%m/%Y-%m-%d" -v
# DONE:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/BB_20221222 -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-case10-20221222 -d "%Y/%Y-%m/%Y-%m-%d" -v

# 240722
# TODO:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/181223_BB -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-caseX-20181223 -d "%Y/%Y-%m/%Y-%m-%d" -v
# TODO:  /automount/agh/s6justei/mambaforge/envs/RADAR_toolbox_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/RADAR_toolbox/DWD_obs_to_MIUB_obs_0_extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/181223_PPI -o /automount/data02/agradar/operation_hydrometeors/data/obs/OpHymet2-caseX-20181223 -d "%Y/%Y-%m/%Y-%m-%d" -v

