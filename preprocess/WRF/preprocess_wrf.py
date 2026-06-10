"""
Filename:    preprocess_wrf.py
Author:      Deanna Nash, dnash@ucsb.edu
Description: Preprocess daily variables from SEAK-WRF data and save as yearly NetCDF files.
"""

# --- Imports ---
import os
import sys
import glob
import xarray as xr
import argparse
import gc

# Add cwd for SLURM execution (script runs from a copied location)
sys.path.append(os.getcwd())

# Path to custom modules
sys.path.append('../../modules')
import globalvars
import wrf_preprocess as wrf_prep

# Path to custom modules
sys.path.append('../')
from build_jobs import build_job_list

# -------------------------------
# MAIN
# -------------------------------
def main(global_id):
    """Main preprocessing workflow."""

    jobs = build_job_list()

    if global_id >= len(jobs):
        return
    
    job = jobs[global_id]

    model = job["model"]
    scenario = job["scenario"]
    year = job["year"]

    # --- Map variables to preprocess functions ---
    varname_lst = ['ivt', 'pcpt', 'freezing_level', 'uv925', 'snow']

    preprocess_map = {
        "ivt": wrf_prep.preprocess_WRF_ivt,
        "freezing_level": wrf_prep.preprocess_WRF_freezing_level,
        "uv925": wrf_prep.preprocess_WRF_uv,
        "pcpt": wrf_prep.preprocess_WRF_pcpt,
        "snow": wrf_prep.preprocess_WRF_snow
    }

    # --- Define paths ---
    path_to_data = globalvars.path_to_data
    path_to_wrf = os.path.join(path_to_data, f"downloads/SEAK-WRF/{model}/{scenario}/")
    
    # --- Collect input files ---
    filenames = sorted(glob.glob(os.path.join(path_to_wrf, f"WRFDS_{year}*")))
    if not filenames:
        raise FileNotFoundError(f"No WRF files found for year {year} in {path_to_wrf}")

    # --- check to see if files exist before processing ---
    vars_to_process = []
    out_paths = {}
    
    for varname in varname_lst:
        path_to_out = os.path.join(
            path_to_data,
            f"preprocessed/SEAK-WRF/{model}/{scenario}/{varname}/"
        )
    
        out_file = os.path.join(
            path_to_out,
            f"WRFDS_{varname}_{year}.nc"
        )
    
        # store paths so we don’t recompute later
        out_paths[varname] = {
            "dir": path_to_out,
            "file": out_file
        }
    
        if not os.path.exists(out_file):
            vars_to_process.append(varname)
    
    # create directories ONCE for variables we will process
    for varname in vars_to_process:
        os.makedirs(out_paths[varname]["dir"], exist_ok=True)
    
    if not vars_to_process:
        print(f"All variables exist for {model} {scenario} {year}")
        return

    # --- Process datasets ---
    # initialize storage per variable
    ds_dict = {v: None for v in vars_to_process}

    for i, wrfin in enumerate(filenames):
        with xr.open_dataset(wrfin) as ds:
            for varname in vars_to_process:
                func = preprocess_map[varname]
                try:
                    ds_proc = func(ds, wrfin)
                    ds_proc = ds_proc.drop_vars("XTIME", errors="ignore")
                except Exception as e:
                    print(f"FAILED {varname} on {wrfin}: {e}")
                    continue

                if ds_dict[varname] is None:
                    ds_dict[varname] = ds_proc
                else:
                    old = ds_dict[varname]
                    ds_dict[varname] = xr.concat([old, ds_proc], dim="Time")
                    del old
    
                del ds_proc
    
        if i > 0 and i % 30 == 0:
            gc.collect()
    
    # --- write datasets to netCDF --- 
    for varname in vars_to_process:
    
        if ds_dict[varname] is None:
            print(f"Skipping {varname}, no data")
            continue
    
        out_file = out_paths[varname]["file"]
    
        new_ds = ds_dict[varname]
    
        encoding = {v: {"zlib": True, "complevel": 4} for v in new_ds.data_vars}
    
        print(f"{model} {scenario} {year} → {varname}")
        new_ds.to_netcdf(out_file, mode="w", format="NETCDF4", encoding=encoding)
    
        del new_ds
        del ds_dict[varname]
    
    gc.collect()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--print-jobs", action="store_true")
    parser.add_argument("--write-jobs", type=str, default=None)

    args = parser.parse_args()

    global_id = args.task_id + args.offset

    jobs = build_job_list()

    # --- debug: print job list ---
    if args.print_jobs:
        for i, job in enumerate(jobs):
            print(f"{i}: {job['model']} {job['scenario']} {job['year']}")
        sys.exit()
    
    # --- debug: write job list to file ---
    if args.write_jobs:
        with open(args.write_jobs, "w") as f:
            for i, job in enumerate(jobs):
                f.write(f"{i},{job['model']},{job['scenario']},{job['year']}\n")
        print(f"Wrote job list to {args.write_jobs}")
        sys.exit()

    main(global_id)
    