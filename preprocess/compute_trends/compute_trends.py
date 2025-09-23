"""
Filename:    compute_trends.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Compute trends for different variables and different models.

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
    
    # --- compute 95th percentile for var for each year ---
    ds = ds.groupby("time.year").quantile(0.95, dim="time").rename(year="time")

    # --- get clim of 95th percentile ---
    ds_clim = ds.mean('time')

    # --- add units to clim ---
    if varname == 'uv925':
        units = 'm s$^{-1}$'
        varname = 'uv'
    elif varname == 'ivt':
        units = 'kg m$^{-1}$ s$^{-1}$'
    elif varname == 'pcpt':
        units = 'mm day$^{-1}$'
    elif varname == 'freezing_level':
        units = 'm'
    ds_clim[varname].attrs['units'] = units

    # --- Save as netCDF ---
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/")
    os.makedirs(datadir, exist_ok=True) # ensure the directory exists
    outdir = os.path.join(datadir, f"{varname}_{model}_clim.nc")
    ds_clim.to_netcdf(path=outdir, mode = 'w', format='NETCDF4')

    # --- get number of years in model - this will be added as an attribute ---
    n_years = ds.time.size  # since time now represents years
    
    # --- Compute trends ---
    vars_lst = [v for v in ds.data_vars if v not in ["lat", "lon"]]
    trend_lst = []
    for varname in vars_lst:
        if varname == 'freezing_level':
            ds[varname] = ds[varname].fillna(0)
        trend = compute_MK_trend_da(ds[varname])
        rename_map = {"trend": f"{varname}_trend", "signif": f"{varname}_signif", 
                      "p": f"{varname}_p", "std_error": f"{varname}_std_error"}
        trend = trend.rename(rename_map)
        trend_lst.append(trend)
    
    trend = xr.merge(trend_lst)
    trend.attrs['n_years'] = int(n_years)
    
    # --- Save as netCDF ---
    outdir = os.path.join(datadir, f"{varname}_{model}_trends.nc")
    trend.to_netcdf(path=outdir, mode = 'w', format='NETCDF4')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compute_trends.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2]
    main(config_file, job_info)