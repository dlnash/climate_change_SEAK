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
from wrf_preprocess import preprocess_WRF_ros, compute_ros_frequency
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
    
def load_preprocessed_WRF_data(model, scenario, varname):

    datadir = os.path.join(
            path_to_data,
            f"preprocessed/SEAK-WRF/{model}/{scenario}/{varname}/"
        )
    
    fname_pattern = os.path.join(
        datadir,
        f"WRFDS_{varname}_*.nc"
    )

    if not glob.glob(fname_pattern):
        raise FileNotFoundError(
            f"No files found: {fname_pattern}"
        )
    
    ds = xr.open_mfdataset(fname_pattern,
                          engine='netcdf4',
                           combine='by_coords')
    
    ## rename coords
    ds = ds.rename({
        "Time": "time",
        "south_north": "y",
        "west_east": "x"
    })
    
    return ds
    
def save_netcdf(ds, model, scenario, varname, season, option, filename_suffix):
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(
    path_to_data,
    f"preprocessed/SEAK-WRF/{model}/{scenario}/trends/"
)
    os.makedirs(datadir, exist_ok=True)
    if option == None:
        outpath = os.path.join(datadir, f"{varname}_{model}_{season}_{filename_suffix}.nc")
    else:
        outpath = os.path.join(datadir, f"{varname}_{model}_{season}_{option}_{filename_suffix}.nc")
    ds.to_netcdf(path=outpath, mode='w', format='NETCDF4')
    print(f"Saved: {outpath}")

def compute_ros_intensity(ds, option, season, model, scenario):
    ds_ros_yearly = preprocess_WRF_ros(ds, temporal_resolution='yearly', option=option, season=season).mean('time')
    units_dict = {
                'ros': 'd yr$^{-1}$',
                'pcpt': 'mm d$^{-1}$',
                'snow': 'mm d$^{-1}$',
                'delsnowh': 'mm d$^{-1}$',
                'ros_intensity': 'mm d$^{-1}$',
            }

    for var, units in units_dict.items():
        if var in ds_ros_yearly:
            ds_ros_yearly[var].attrs['units'] = units

    save_netcdf(ds_ros_yearly, model, scenario, 'snow', season, option, "ros_intensity_clim")

def save_ros_frequency(ds, option, season, model, scenario):
    ds_ros_daily = preprocess_WRF_ros(ds, temporal_resolution='daily', option=option, season=season)
    ivt = load_preprocessed_WRF_data(model, scenario, 'ivt')
    mon_s, mon_e = get_startmon_and_endmon(season)
    ivt = select_months_ds(ivt, mon_s, mon_e, time_varname='time')
    ds_ros_daily = xr.merge([ds_ros_daily, ivt], compat="no_conflicts")
    ds_out = compute_ros_frequency(ds_ros_daily)
    save_netcdf(ds_out, model, scenario, 'snow', season, option, "ros_frequency_clim")


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