######################################################################
# Filename:    create_job_configs.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Create a single YAML configuration and call file
#              (plus grouped call files by model and variable)
#              to run job_array with Slurm for WRF data preprocessing.
######################################################################

import yaml
from collections import defaultdict

# ---------------- Configuration ----------------
conda_path = "/home/dnash/miniconda3/envs/SEAK-impacts/bin/python"
script_name = "compute_frequency_intensity_changes.py"  # update if needed

models = ["ccsm", "gfdl", "cfsr"]
varnames = ["ivt", "pcpt", "freezing_level", "uv925", "snow"]
seasons = ["DJF", "MAM", "JJA", "SON", "NDJFMA"]

# ---------------- Build job dictionary ----------------
config_dict = {}
job_counter = 0

# Group trackers for call files
calls_by_model = defaultdict(list)
calls_by_var = defaultdict(list)

for model in models:
    for varname in varnames:
        for ssn in seasons:
            job_counter += 1
            job_id = f"job_{job_counter}"
            config_dict[job_id] = {
                "model": model,
                "varname": varname,
                "season": ssn,
            }

            call_cmd = f"{conda_path} -u {script_name} config_all.yaml {job_id}"
            calls_by_model[model].append(call_cmd)
            calls_by_var[varname].append(call_cmd)

# ---------------- Write master config ----------------
config_filename = "config_all.yaml"
with open(config_filename, "w") as f:
    yaml.dump(config_dict, f, allow_unicode=True)
print(f"Wrote {config_filename} with {job_counter} jobs.")

# ---------------- Write master call file ----------------
call_filename = "calls_all.txt"
all_calls = [cmd for lst in calls_by_model.values() for cmd in lst]
with open(call_filename, "w", encoding="utf-8") as f:
    f.write("\n".join(all_calls))
print(f"Wrote {call_filename} with {job_counter} job calls.")

# ---------------- Write grouped call files ----------------
for model, cmds in calls_by_model.items():
    fname = f"calls_{model}.txt"
    with open(fname, "w", encoding="utf-8") as f:
        f.write("\n".join(cmds))
    print(f"  → Wrote {fname} ({len(cmds)} jobs)")

for varname, cmds in calls_by_var.items():
    fname = f"calls_{varname}.txt"
    with open(fname, "w", encoding="utf-8") as f:
        f.write("\n".join(cmds))
    print(f"  → Wrote {fname} ({len(cmds)} jobs)")

print("All grouped and master call files created successfully.")
