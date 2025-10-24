"""
Filename:    compute_frequency_changes.py
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
from wrf_preprocess import preprocess_WRF_ros
from utils import get_startmon_and_endmon, select_months_ds

# ============================================================
# Helper Functions
# ============================================================
def save_netcdf(ds, model, varname, season, filename_suffix):
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/")
    os.makedirs(datadir, exist_ok=True)
    outpath = os.path.join(datadir, f"{varname}_{model}_{season}_{filename_suffix}.nc")
    ds.to_netcdf(path=outpath, mode='w', format='NETCDF4')
    print(f"Saved: {outpath}")


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
    ds = select_months_ds(ds, mon_s, mon_e, time_varname='time')

    # --- compute 95th percentile for var for each year ---
    ds = ds.groupby("time.year").quantile(0.95, dim="time").rename(year="time")

    # --- get clim ---
    ds_clim = ds.mean('time')

    # --- add units to clim ---
    units_dict = {
                    'uv925': ('uv', 'm s$^{-1}$'),
                    'ivt': ('ivt', 'kg m$^{-1}$ s$^{-1}$'),
                    'pcpt': ('pcpt', 'mm day$^{-1}$'),
                    'freezing_level': ('freezing_level', 'm'),
                    'snow': ('snow', 'kg m$^{-2}$'),
                }
    varname, units = units_dict.get(varname, (varname, ''))
    ds_clim[varname].attrs['units'] = units

    # --- Save as netCDF ---
    save_netcdf(ds_clim, model, varname, season, "95th_percentile_clim")

    if varname == 'snow':

        # --- compute ROS intensity ---
        ds_ros_yearly = preprocess_WRF_ros(ds, temporal_resolution='yearly').mean('time')
        save_netcdf(ds_ros_yearly, model, 'snow', season, "ros_intensity_clim", path_to_data)
    
        # --- compute ROS frequency ---
        ds_ros_daily = preprocess_WRF_ros(ds, temporal_resolution='daily')
        ivt = load_preprocessed_WRF_data(model, 'ivt', anomaly=False)
        ivt = select_months_ds(ivt, mon_s, mon_e, time_varname='time')
        ds = xr.merge([ds_ros_daily, ivt], compat="no_conflicts")
        ds_out = compute_ros_frequency(ds)
        save_netcdf(ds_out, model, 'snow', season, "ros_frequency_clim")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compute_frequency_changes.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2], sys.argv[3]
    main(config_file, job_info)