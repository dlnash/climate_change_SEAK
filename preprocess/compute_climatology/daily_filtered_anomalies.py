"""
Filename:    daily_filtered_anomalies.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions to filter annual climatology of daily gridded/time series data using harmonics.

"""
import sys, os
import numpy as np
import xarray as xr
import yaml
from pathlib import Path

# Add cwd for SLURM execution (script runs from a copied location)
sys.path.append(os.getcwd())

# Path to modules
sys.path.append('../../modules/')
import globalvars 
from harmonics import harmonic

def main(config_file: str, job_info: str):
    """Main preprocessing workflow."""

    # --- Load configuration ---
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ddict = config[job_info]
    varname = ddict["varname"]
    model = ddict["model"]

    # --- Define paths ---
    path_to_data = globalvars.path_to_data
    datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/{varname}/")

    # --- Read data ---
    print('Step 1: Reading data...')
    filename_pattern = os.path.join(datadir, f"WRFDS_{varname}_*.nc")
    print(filename_pattern)
    ds = xr.open_mfdataset(filename_pattern,
                           engine='netcdf4',
                           combine='by_coords')
    
    ds = ds.sortby('Time')
    
    # --- Calculate Annual Climatology ---
    print('Step 2: Calculating annual climatology...')
    clim_mean = ds.groupby('Time.dayofyear').mean('Time')
    clim_std = ds.groupby('Time.dayofyear').std('Time')
    
    # --- Save daily clim as netCDF ---
    clim_path = os.path.join(datadir, f"daily_mean_clim_{varname}.nc")
    clim_mean.load().to_netcdf(path=clim_path, mode = 'w', format='NETCDF4')
    
    # --- Save standard deviation as netCDF ---
    std_path = os.path.join(datadir, f"daily_std_clim_{varname}.nc")
    clim_std.load().to_netcdf(path=std_path, mode = 'w', format='NETCDF4')
    
    # --- Filter annual climatology using harmonics ---
    print('Step 3: Filtering annual climatology...')
    filtered_clim = harmonic(clim_mean)
    
    # --- Save filtered climatology as netCDF ---
    print('Step 4: Saving filtered climatology...')
    clim_path = os.path.join(datadir, f"filtered_daily_mean_clim_{varname}.nc")
    filtered_clim.load().to_netcdf(path=clim_path, mode = 'w', format='NETCDF4')
    
    # --- Calculate anomalies ---
    print('Step 5: Calculating anomalies...')
    anomalies = ds.unify_chunks().groupby('Time.dayofyear') - filtered_clim
    
    # --- Write anomalies to yearly file ---
    print('Step 6: Writing anomalies to yearly files...')
    years, datasets = zip(*anomalies.groupby('Time.year'))
    paths = [datadir +'anomalies/daily_filtered_anomalies_{0}_%s.nc'.format(varname) % y for y in years]
    # ensure the directory exists
    Path(os.path.join(datadir, "anomalies")).mkdir(parents=True, exist_ok=True)
    xr.save_mfdataset(datasets, paths)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python daily_filtered_anomalies.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2]
    main(config_file, job_info)