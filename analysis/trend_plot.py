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
varnames = ['ivt', 'pcpt', 'freezing_level', 'uv']

# --- Loop over all variables ---
for varname in varnames:
    print(f"Processing {varname}...")

    fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_clim.nc")
    ds_clim = xr.open_dataset(fname)

    # Get lons/lats from this dataset
    lons = ds_clim.lon.values
    lats = ds_clim.lat.values

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
