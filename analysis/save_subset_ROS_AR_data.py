import sys
import string
from pathlib import Path
import pandas as pd
import numpy as np
import gc

# --- Personal Modules ---
sys.path.append('../modules')
import globalvars
from utils import get_startmon_and_endmon, select_months_df, select_months_ds
from wrf_preprocess import load_preprocessed_WRF_data, preprocess_WRF_ros

study_start = "1981-01-01"
study_end   = "2019-12-31"

ds = load_preprocessed_WRF_data("cfsr", "historical", "snow")

# Compute for 'strict'
ds_ros = preprocess_WRF_ros(ds, 
                           temporal_resolution='daily', 
                           option="strict", 
                           season=None)

ds_ros = ds_ros['ros'].load()

outpath = '../out/cfsr_ros_daily.nc'
ds_ros.to_netcdf(path=outpath, mode='w', format='NETCDF4')
print(f"Saved: {outpath}")

del ds, ds_ros
gc.collect()

# --- Load AR dataset once ---
ar_ds = xr.open_dataset(globalvars.path_to_data + 'downloads/globalARcatalog_ERA5_1940-2024_v4.0.nc')
ar_ds = ar_ds.kidmap.squeeze()
ar_ds = ar_ds.assign_coords({"lon": (((ar_ds.lon + 180) % 360) - 180)}).sortby('lon')
ar_ds = ar_ds.sel(time=slice(study_start, study_end))
ar_ds = ar_ds.sel(lat=slice(62., 54.), lon=slice(-146., -127.))
ar_ds = ar_ds.resample(time="1D").max()
outpath = '../out/ar_ds_daily.nc'
ar_ds.to_netcdf(path=outpath, mode='w', format='NETCDF4')
print(f"Saved: {outpath}")