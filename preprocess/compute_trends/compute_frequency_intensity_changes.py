"""
Filename:    compute_frequency_intensity_changes.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Compute changes in the number of days where IVT >=250 kg m-1 s-1, Rainfall >= 5 mm day-1, and Snow > 3 mm day-1

"""
import sys, os
import pandas as pd
import xarray as xr
import yaml

# Add cwd for SLURM execution (script runs from a copied location)
sys.path.append(os.getcwd())

# Path to modules
sys.path.append('../../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data
from wrf_preprocess import preprocess_WRF_ros, compute_ros_frequency
from utils import get_startmon_and_endmon, select_months_ds

# ============================================================
# Helper Functions
# ============================================================
def save_netcdf(ds, model, varname, season, option, filename_suffix):
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/")
    os.makedirs(datadir, exist_ok=True)
    if option == None:
        outpath = os.path.join(datadir, f"{varname}_{model}_{season}_{filename_suffix}.nc")
    else:
        outpath = os.path.join(datadir, f"{varname}_{model}_{season}_{option}_{filename_suffix}.nc")
    ds.to_netcdf(path=outpath, mode='w', format='NETCDF4')
    print(f"Saved: {outpath}")

def compute_ros_intensity(ds, option, season, model):
    ds_ros_yearly = preprocess_WRF_ros(ds, temporal_resolution='yearly', option=option, season=season).mean('time')
    units_dict = {
                'ros': 'days yr$^{-1}$',
                'pcpt': 'mm day$^{-1}$',
                'snow': 'mm day$^{-1}$',
                'delsnowh': 'mm day$^{-1}$',
                'ros_intensity': 'mm day$^{-1}$',
            }

    for var, units in units_dict.items():
        if var in ds_ros_yearly:
            ds_ros_yearly[var].attrs['units'] = units

    save_netcdf(ds_ros_yearly, model, 'snow', season, option, "ros_intensity_clim")

def save_ros_frequency(ds, option, season, model):
    ds_ros_daily = preprocess_WRF_ros(ds, temporal_resolution='daily', option=option, season=season)
    ivt = load_preprocessed_WRF_data(model, 'ivt', anomaly=False)
    mon_s, mon_e = get_startmon_and_endmon(season)
    ivt = select_months_ds(ivt, mon_s, mon_e, time_varname='time')
    ds_ros_daily = xr.merge([ds_ros_daily, ivt], compat="no_conflicts")
    ds_out = compute_ros_frequency(ds_ros_daily)
    save_netcdf(ds_out, model, 'snow', season, option, "ros_frequency_clim")


def main(config_file: str, job_info: str):
    """Main preprocessing workflow."""
    
    # --- Load configuration ---
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ddict = config[job_info]
    varname = ddict["varname"]
    model = ddict["model"]
    season = ddict["season"]
 
    # --- read the non-anomaly data ---
    ds = load_preprocessed_WRF_data(model, varname, anomaly=False)

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
                    'pcpt': ('pcpt', 'mm day$^{-1}$'),
                    'freezing_level': ('freezing_level', 'm'),
                    'snow': ('snow', 'mm day$^{-1}$'),
                }
    varname, units = units_dict.get(varname, (varname, ''))
    ds_95th[varname].attrs['units'] = units

    # --- Save as netCDF ---
    save_netcdf(ds_95th, model, varname, season, option=None, filename_suffix="95th_percentile_clim")

    if varname == 'snow':
        option_lst = ['strict', 'flexible']
        for option in option_lst:
            # --- compute ROS intensity ---
            compute_ros_intensity(ds, option, season, model)
        
            # --- compute ROS frequency ---
            save_ros_frequency(ds, option, season, model)
            

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compute_frequency_changes.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2]
    main(config_file, job_info)