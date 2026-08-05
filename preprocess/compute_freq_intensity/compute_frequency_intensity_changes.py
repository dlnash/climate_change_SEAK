"""
Filename:    compute_frequency_intensity_changes.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Compute changes in the number of days where IVT >=250 kg m-1 s-1, Rainfall >= 5 mm day-1, and Snow > 3 mm day-1

"""
# --- Imports ---
import sys, os
import argparse
import xarray as xr
import glob

# Add cwd for SLURM execution (script runs from a copied location)
sys.path.append(os.getcwd())

# Path to modules
sys.path.append('../../modules/')
import globalvars
path_to_data = globalvars.path_to_data
from wrf_preprocess import preprocess_WRF_ros, compute_ros_frequency, compute_ros_intensity, save_ros_frequency, load_preprocessed_WRF_data
from utils import get_startmon_and_endmon, select_months_ds

# ============================================================
# Helper Functions
# ============================================================
def build_job_list():
    jobs = []

    models = {
        "cfsr": ["historical"],
        "ccsm": ["hist", "rcp85"],
        "gfdl": ["hist", "rcp85"],
    }

    varnames = [
        "ivt",
        "pcpt",
        "freezing_level",
        "uv925",
        "snow",
    ]

    seasons = ["ONDJFM"]

    for model, scenarios in models.items():
        for scenario in scenarios:
            for varname in varnames:
                for season in seasons:
                    jobs.append({
                        "model": model,
                        "scenario": scenario,
                        "varname": varname,
                        "season": season,
                    })

    return jobs
    

def main(global_id):
    """Main preprocessing workflow."""
    
    jobs = build_job_list()

    if global_id >= len(jobs):
        return

    job = jobs[global_id]

    model = job["model"]
    scenario = job["scenario"]
    varname = job["varname"]
    season = job["season"]
 
    # --- read the non-anomaly data ---
    ds = load_preprocessed_WRF_data(model, scenario, varname)

    # --- subset to specified season ---
    mon_s, mon_e = get_startmon_and_endmon(season)
    ds_ssn = select_months_ds(ds, mon_s, mon_e, time_varname='time')

    # --- compute 95th percentile for var for each year ---
    ds_95th = ds_ssn.groupby("time.year").quantile(0.95, dim="time").rename(year="time")

    # --- get clim ---
    ds_95th = ds_95th.mean('time')

    # --- add units to clim ---
    units_dict = {
                    'uv925': ('uv', 'm s$^{-1}$'),
                    'ivt': ('ivt', 'kg m$^{-1}$ s$^{-1}$'),
                    'pcpt': ('pcpt', 'mm d$^{-1}$'),
                    'freezing_level': ('freezing_level', 'm'),
                    'snow': ('snow', 'mm d$^{-1}$'),
                }
    varname, units = units_dict.get(varname, (varname, ''))
    ds_95th[varname].attrs['units'] = units

    # --- Save as netCDF ---
    save_netcdf(ds_95th, model, scenario, varname, season, option=None, filename_suffix="95th_percentile_clim")

    if varname == 'snow':
        # Work with clean copies so nothing carries over between runs
        ds_strict = ds.copy(deep=True)
        ds_flexible = ds.copy(deep=True)
    
        # Compute for 'strict'
        compute_ros_intensity(ds_strict, 'strict', season, model, scenario)
        save_ros_frequency(ds_strict, 'strict', season, model, scenario)
    
        # Compute for 'flexible'
        compute_ros_intensity(ds_flexible, 'flexible', season, model, scenario)
        save_ros_frequency(ds_flexible, 'flexible', season, model, scenario)

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
            print(f"{i}: {job['model']} {job['scenario']} {job['varname']} {job['season']}")
        sys.exit()
    
    # --- debug: write job list to file ---
    if args.write_jobs:
        with open(args.write_jobs, "w") as f:
            for i, job in enumerate(jobs):
                f.write(f"{i},{job['model']},{job['scenario']},{job['varname']},{job['season']}\n")
        print(f"Wrote job list to {args.write_jobs}")
        sys.exit()

    main(global_id)