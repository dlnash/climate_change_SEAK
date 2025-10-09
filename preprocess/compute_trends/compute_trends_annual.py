"""
Filename:    compute_trends_annual.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Compute annual trends for different variables and different models.

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
from xarrayMannKendall import compute_MK_trend_da

def main(config_file: str, job_info: str):
    """Main preprocessing workflow."""

    # --- Load configuration ---
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ddict = config[job_info]
    varname = ddict["varname"]
    model = ddict["model"]
 
    # --- read the non-anomaly data ---
    ds = load_preprocessed_WRF_data(model, varname, anomaly=False)

    # --- calculate rain (PCPT - SWE) ---
    ds['rain'] = ds['pcpt'] - ds['snow']
    
    # --- update snow depth to mm ---
    snowh_attrs = ds['snowh'].attrs
    ds['snowh'] = ds['snowh'] * 1000
    ds['snowh'].attrs = snowh_attrs
    ds['snowh'].attrs['units'] = 'mm'
    
    # --- calculate ROS intensity ---
    # ROS intensity is defined as the change in snow depth (snowh) for the ROS day
    ds['delsnow'] = ds['snow'].diff(dim='time') # n+1 - n
    ds['delsnowh'] = ds['snowh'].diff(dim='time') # n+1 - n
    
    # --- calculate ROS ---
    # A ROS event is defined as a day when rain (pcpt) is greater than 5 mm d−1, 
    # delta snowh (change in snow depth) decreases (< 0) and snow cover is true
    # Create the mask and add it as a new variable called "ros"
    ds['ros'] = ((ds['rain'] > 5) & (ds['delsnowh'] < 0) & (ds['snowc'] > 0)).astype(int)

    # --- calculate ROS intensity --- 
    # ROS intensity is the sum of rainfall (ds['rain']) and snowmelt (ds[delsnowh']*-1)
    ds['ros_intensity'] = ds['rain'].where(ds['ros'] == 1) + (ds['delsnowh'].where(ds['ros'] == 1))*-1

    # --- compute sum of rain-on-snow days per year ---
    da_ros = ds['ros'].groupby("time.year").sum('time').rename(year="time")
    
    # --- compute change in rainfall and snowmelt during ROS events
    rain_mask = ds['rain'].where(ds['ros'] == 1)
    # negative delsnowh indicates snow decreased between today and tomorrow
    # need to multiply by 1
    snowmelt_mask = (ds['delsnowh'].where(ds['ros'] == 1))*-1 

    da_rain = rain_mask.groupby("time.year").mean('time').rename(year="time")
    da_snowmelt = snowmelt_mask.groupby("time.year").mean('time').rename(year="time")
    da_ros_intensity = ds['ros_intensity'].groupby("time.year").mean('time').rename(year="time")

    # --- merge into single dataset ---
    ds = xr.merge([da_ros, da_rain, da_snowmelt, da_ros_intensity])

    # --- get clim ---
    ds_clim = ds.mean('time')

    # --- Save as netCDF ---
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/")
    os.makedirs(datadir, exist_ok=True) # ensure the directory exists
    outdir = os.path.join(datadir, f"{varname}_{model}_ros_clim.nc")
    ds_clim.to_netcdf(path=outdir, mode = 'w', format='NETCDF4')

    # --- get number of years in model - this will be added as an attribute ---
    n_years = ds.time.size  # since time now represents years
    
    # --- Compute trends ---
    vars_lst = [v for v in ds.data_vars if v not in ["lat", "lon"]]
    trend_lst = []
    for varname in vars_lst:
        ds[varname] = ds[varname].fillna(0)
        trend = compute_MK_trend_da(ds[varname])
        rename_map = {"trend": f"{varname}_trend", "signif": f"{varname}_signif", 
                      "p": f"{varname}_p", "std_error": f"{varname}_std_error"}
        trend = trend.rename(rename_map)
        trend_lst.append(trend)
    
    trend = xr.merge(trend_lst)
    trend.attrs['n_years'] = int(n_years)
    
    # --- Save as netCDF ---
    outdir = os.path.join(datadir, f"snow_{model}_ros_trends.nc")
    trend.to_netcdf(path=outdir, mode = 'w', format='NETCDF4')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compute_trends.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2]
    main(config_file, job_info)