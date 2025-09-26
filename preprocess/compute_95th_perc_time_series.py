######################################################################
# Filename:    compute_95th_perc_time_series.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Reads the preprocessed WRF data (variables and models), subsets to the six communities, then saves as a single netCDF.
#
######################################################################


import sys, os
import numpy as np
import xarray as xr
import pandas as pd

# Path to modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data, subset_wrf_ds

models = ["cfsr", "gfdl", "era5"]
varnames = ["ivt", "pcpt", "freezing_level", "uv", "snow"]

model_ds_list = []  # list of datasets, one per model

for model in models:
    var_ds_list = []  # list of datasets for this model
    for varname in varnames:
        # handle uv special case
        load_varname = "uv925" if varname == "uv" else varname

        # --- load data ---
        ds = load_preprocessed_WRF_data(model, load_varname, anomaly=False)

        # --- subset to community locations ---
        ds = subset_wrf_ds(ds)

        # --- compute 95th percentile for var for each year ---
        ds = ds.groupby("time.year").quantile(0.95, dim="time").rename(year="time")

        var_ds_list.append(ds)

    # merge all variables for this model
    merged = xr.merge(var_ds_list)
    merged = merged.expand_dims(model=[model])  # add model dimension
    model_ds_list.append(merged)

# concat along model dimension
final_ds = xr.concat(model_ds_list, dim="model")

# save to netCDF
outpath = "../data/preprocessed_trends.nc"
final_ds.to_netcdf(outpath)

print(f"Saved preprocessed dataset to {outpath}")
