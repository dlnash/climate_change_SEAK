import pandas as pd
import yaml
from itertools import chain

# adjustable maximum jobs per config file
MAX_JOBS_PER_FILE = 999
conda_path = "/home/dnash/miniconda3/envs/SEAK-impacts/bin/python"
models = ["ccsm", "gfdl", "cfsr"]

for model in models:
    # set year ranges depending on model
    if model in ["ccsm", "gfdl"]:
        yr_start, yr_end = 2030, 2060
    elif model == "cfsr":
        yr_start, yr_end = 1980, 2019
    else:
        continue  # skip unknown models

    # create list of monthly dates
    start_date = pd.to_datetime(f"{yr_start}-01-01")
    end_date = pd.to_datetime(f"{yr_end}-12-31")
    date_lst = pd.date_range(start_date, end_date, freq="1M")

    jobcounter = 0
    filecounter = 0
    d_lst = []
    njob_lst = []

    # loop through dates to build config entries
    for date in date_lst:
        yr = date.strftime("%Y")
        month = date.strftime("%m")

        jobcounter += 1
        d = {
            f"job_{jobcounter}": {
                "year": yr,
                "month": month,
                "model": model,
            }
        }
        d_lst.append(d)

        # write to YAML once we reach MAX_JOBS_PER_FILE
        if jobcounter == MAX_JOBS_PER_FILE:
            filecounter += 1
            dest = dict(chain.from_iterable(map(dict.items, d_lst)))
            njob_lst.append(len(d_lst))
            with open(f"config_{model}_{filecounter}.yaml", "w") as file:
                yaml.dump(dest, file, allow_unicode=True, default_flow_style=None)

            # reset for next file
            jobcounter = 0
            d_lst = []

    # write remaining jobs if any left
    if d_lst:
        filecounter += 1
        dest = dict(chain.from_iterable(map(dict.items, d_lst)))
        njob_lst.append(len(d_lst))
        with open(f"config_{model}_{filecounter}.yaml", "w") as file:
            yaml.dump(dest, file, allow_unicode=True, default_flow_style=None)

    # now create calls_MODEL_N.txt for each config file
    for i, njobs in enumerate(njob_lst):
        call_str_lst = []
        for j in range(1, njobs + 1):
            call_string = f"{conda_path} -u getWRF_batch.py config_{model}_{i+1}.yaml 'job_{j}'"
            call_str_lst.append(call_string)

        with open(f"calls_{model}_{i+1}.txt", "w", encoding="utf-8") as f:
            for line in call_str_lst:
                f.write(line + "\n")
