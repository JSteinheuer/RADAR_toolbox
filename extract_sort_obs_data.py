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


# /automount/agh/s6justei/mambaforge/envs/hydrometeors_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/OP_Hydrometeors/extract_sort_obs_data.py -i /automount/realpep/upload/RealPEP-SPP/newdata/s6justei_2023/2021-06-28 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case03-20210628 -d "%Y%m%d" -v

# /automount/agh/s6justei/mambaforge/envs/hydrometeors_agh/bin/python /user/s6justei/PyCharm/PyCharmProjects/OP_Hydrometeors/extract_sort_obs_data.py -i /automount/realpep/upload/RealPEP-SPP/newdata/s6justei_2023/2022-05-19 -o /automount/agradar/operation_hydrometeors/data/obs/OpHymet2-case04-20220519 -d "%Y%m%d" -v

# /automount/agh/s6justei/mambaforge/envs/hydrometeors_agh/bin/python /automount/realpep/upload/RealPEP-SPP/scripts/extract_sort_obs_data.py -i /automount/ftp/wwwgast/spp-prom/OBS/OpHymet2_CASE06+07_2022-06-26/ -o /automount/agradar/operation_hydrometeors/data/obs -d "%Y%m%d" -v
