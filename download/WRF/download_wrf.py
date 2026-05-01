#!/usr/bin/env python
######################################################################
# Filename:    download_WRF.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to download Lader et al., 2020 SEAK WRF data
# https://registry.opendata.aws/wrf-se-ak-ar5/ (data link)
#
######################################################################
import os
import calendar
import argparse
import subprocess
from itertools import product
from multiprocessing import Pool

# -------------------------------
# CONFIG
# -------------------------------
BASE_OUT = "/cw3e/mead/projects/cwp140/data/downloads/SEAK-WRF"

MODEL_CONFIGS = {
    "cfsr": {
        "scenarios": ["historical"],
        "years": list(range(1981, 2020)),
    },
    "ccsm": {
        "scenarios": ["hist", "rcp85"],
        "years": {
            "hist": list(range(1981, 2011)),
            "rcp85": list(range(2031, 2061)),
        },
    },
    "gfdl": {
        "scenarios": ["hist", "rcp85"],
        "years": {
            "hist": list(range(1981, 2011)),
            "rcp85": list(range(2031, 2061)),
        },
    },
}

# -------------------------------
# BUILD JOB LIST
# -------------------------------
def build_job_list():
    jobs = []

    for model, cfg in MODEL_CONFIGS.items():
        for scenario in cfg["scenarios"]:

            if model == "cfsr":
                years = cfg["years"]
            else:
                years = cfg["years"][scenario]

            for year, month in product(years, range(1, 13)):
                jobs.append({
                    "model": model,
                    "scenario": scenario,
                    "year": year,
                    "month": month
                })

    return jobs


# -------------------------------
# PATH BUILDER
# -------------------------------
def get_s3_path(model, scenario, year):
    if model in ["ccsm", "gfdl"]:
        return f"s3://wrf-se-ak-ar5/{model}/{scenario}/daily/{year}/"
    elif model == "cfsr":
        return f"s3://wrf-se-ak-ar5/{model}/4km/daily/{year}/"
    else:
        raise ValueError(f"Unknown model: {model}")


# -------------------------------
# DOWNLOAD FUNCTION
# -------------------------------
def download_day(args):
    model, scenario, year, month, day = args

    day_str = str(day).zfill(2)
    month_str = str(month).zfill(2)

    s3_path = get_s3_path(model, scenario, year)
    fname = f"WRFDS_{year}-{month_str}-{day_str}.nc"

    input_path = f"{s3_path}{fname}"
    out_dir = os.path.join(BASE_OUT, model, scenario)
    os.makedirs(out_dir, exist_ok=True)

    output_path = os.path.join(out_dir, fname)

    # skip existing
    if os.path.exists(output_path):
        return f"SKIP {fname}"

    cmd = [
        "/cw3e/mead/projects/cwp140/aws/local/bin/aws", "s3", "cp",
        "--region", "us-west-2",
        input_path,
        output_path,
        "--no-sign-request"
    ]

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return f"DOWNLOADED {fname}"
    except subprocess.CalledProcessError as e:
        return f"FAILED {fname}: {e}"


# -------------------------------
# MAIN
# -------------------------------
def main(global_id, nproc):

    jobs = build_job_list()

    if global_id >= len(jobs):
        return
    
    job = jobs[global_id]

    model = job["model"]
    scenario = job["scenario"]
    year = job["year"]
    month = job["month"]

    print(f"Running: {model} {scenario} {year}-{month:02d}")

    # days in month
    ndays = calendar.monthrange(year, month)[1]

    args_list = [
        (model, scenario, year, month, day)
        for day in range(1, ndays + 1)
    ]

    # parallel download within node
    with Pool(nproc) as pool:
        results = pool.map(download_day, args_list)

    for r in results:
        print(r)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--nproc", type=int, default=4)

    args = parser.parse_args()

    global_id = args.task_id + args.offset

    main(global_id, args.nproc)