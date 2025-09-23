"""
Filename:    trend_plot.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Script to plot trends for each variable given a model name.
"""

import sys, os
import numpy as np
import xarray as xr
import pandas as pd

# matplotlib
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import cmocean
# cartopy
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar
import cartopy.feature as cfeature

# Path to modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data
from plotter import plot_trend_with_clim

# --- Settings ---
model = 'cfsr'
path_to_data = globalvars.path_to_data
varnames = ['ivt', 'pcpt', 'freezing_level', 'uv925']

# --- load dates list ---
fname = '../out/PCPT_95th_25perc-cov_dates.csv'
df = pd.read_csv(fname)
dates = pd.to_datetime(df['date'], format='%Y-%m-%d')

# --- Loop over all variables ---
for varname in varnames:
    print(f"Processing {varname}...")

    # Load full dataset (non-anomaly) for climatology
    ds_full = load_preprocessed_WRF_data(model, varname, anomaly=False)

    # Average climatology over selected dates
    ds_clim = ds_full.sel(time=dates.values).mean('time')

    # Get lons/lats from this dataset
    lons = ds_full.lon.values
    lats = ds_full.lat.values

    # Load precomputed trends file
    trend_fname = os.path.join(
        path_to_data,
        f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_trends.nc"
    )
    ds_trend = xr.open_dataset(trend_fname)

    # Make plot
    plot_trend_with_clim(
        ds_trend, ds_clim, varname, lons, lats, model,
        lonmin=-141., lonmax=-130., latmin=54.5, latmax=60.,
        sig_level=0.1
    )
