######################################################################
# Filename:    create_job_configs.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to create .yaml configuration file to run job_array with slurm for preprocessing WRF Data
#
######################################################################

import yaml

MAX_JOBS_PER_FILE = 999
conda_path = "/home/dnash/miniconda3/envs/SEAK-impacts/bin/python"
models = ["ccsm", "gfdl", "cfsr"]
varnames = ['ivt', 'pcpt', 'freezing_level', 'uv925', 'snow']

for model in models:

    jobcounter = 0
    filecounter = 0
    d_lst = []
    njob_lst = []

    for varname in varnames:
        jobcounter += 1
        d_lst.append({
            f"job_{jobcounter}": {
                "model": model,
                "varname": varname
            }
        })

        if jobcounter == MAX_JOBS_PER_FILE:
            filecounter += 1
            dest = {k: v for dd in d_lst for k, v in dd.items()}
            njob_lst.append(len(d_lst))
            with open(f"config_{model}_{filecounter}.yaml", "w") as f:
                yaml.dump(dest, f, allow_unicode=True)
            jobcounter = 0
            d_lst = []

    # write remaining jobs
    if d_lst:
        filecounter += 1
        dest = {k: v for dd in d_lst for k, v in dd.items()}
        njob_lst.append(len(d_lst))
        with open(f"config_{model}_{filecounter}.yaml", "w") as f:
            yaml.dump(dest, f, allow_unicode=True)

    # create calls files
    for i, njobs in enumerate(njob_lst):
        call_lines = [
            f"{conda_path} -u compute_trends.py config_{model}_{i+1}.yaml 'job_{j}'"
            for j in range(1, njobs + 1)
        ]
        with open(f"calls_{model}_{i+1}.txt", "w", encoding="utf-8") as f:
            f.write("\n".join(call_lines))

