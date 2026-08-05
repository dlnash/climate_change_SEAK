#!/usr/bin/env python3
"""
Filename:    pull_numbers_for_manuscript.py
Author:      Deanna Nash (adapted for multi-model comparison)
Description: 
    Open preprocessed NetCDF files, subset each file to four defined 
    locations from a shapefile, compute min, max, and mean values 
    for each model, and calculate model–CFSR differences for each region.
"""
# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------

import os
import sys
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray  # for spatial clipping
from pyproj import CRS
from pathlib import Path

# Add path to custom modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_wrf_metrics
from plot_configs import PLOT_CONFIGS, MODELS, SCENARIOS
from wrf_rio import wrf_prepare_for_rio

# ---------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------
def compute_area_avg_max_min(
    path_to_data,
    plot_type,
    season,
    option=None,
):
    """
    Compute area-averaged statistics for each model and variable, including
    model minus CFSR differences.

    Parameters
    ----------
    models : list of str
        Model names (first should be 'cfsr' for baseline)
    season : str
        Season label (e.g., 'DJF')
    option: str
        choice in ROS category ('strict' or 'flexible')
    path_to_data : str
        Path to data directory
    fsuffix : str
        Filename suffix ('95th_percentile_clim', 'ros_intensity_clim', or 'ros_frequency_clim')

    Returns
    -------
    pandas.DataFrame
        DataFrame containing model, variable, polygon name, mean, min, max, and model–CFSR diff
    """

    cfg = PLOT_CONFIGS[plot_type]

    varnames = cfg["varnames"]
    var_labels = cfg["labels"]

    # --- load shapefile and subset ---
    fp = os.path.join(globalvars.path_to_data, 'downloads/AK_climate_divisions/AK_divisions_NAD83.shp')
    polys = gpd.read_file(fp)
    keep_names = ["Northeast Gulf", "North Panhandle", "Central Panhandle", "South Panhandle"]
    polys = polys[polys["Name"].isin(keep_names)]

    results = []

    # --- loop through variables ---
    for i, varname in enumerate(varnames):
        fields, titles, units = load_preprocessed_wrf_metrics(
            varname, 
            plot_type, 
            season, 
            option
        )

        # --- loop through MMM Historical, MMM Future, Model Difference --- 
        for k, field in enumerate(fields):

            # prepare for rioxarray
            # da_cfsr = ds_cfsr[varname]
            da = wrf_prepare_for_rio(field)
    
            # reproject shapefile to match WRF CRS
            polys = polys.to_crs(da.rio.crs)

            # print(f"\nVariable: {varname}, Model: {model}")
            # print(f"  Raster CRS: {da.rio.crs}")
            # print(f"  Polygon CRS: {polys.crs}")
            # print(f"  Raster bounds: {da.rio.bounds()}")
            # print(f"  Polygon bounds: {polys.geometry.bounds}")


            # --- loop through polygons ---
            for _, poly in polys.iterrows():
                # Clip both datasets to polygon
                da_clip = da.rio.clip([poly.geometry], polys.crs, drop=True)


                # Compute stats
                stats = lambda da: {
                    "mean": float(da.mean().values),
                    "min": float(da.min().values),
                    "max": float(da.max().values)
                }

                da_stats = stats(da_clip)

                results.append({
                    "season": season,
                    "plot_type": plot_type,
                    "variable": varname,
                    "region": poly["Name"],
                    "metric_name": titles[k],
                    "mean_model": da_stats["mean"],
                    "min_model": da_stats["min"],
                    "max_model": da_stats["max"],
                })

    df = pd.DataFrame(results)

    outfile = os.path.join(path_to_data, f"processed_summary_stats_with_diffs_{option}_{season}_{plot_type}.csv")
    df.to_csv(outfile, index=False)
    print(f"\n✅ Saved summary statistics with differences to {outfile}")


    return df


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    seasons = ["ONDJFM"]
    options = ["strict"]

    # --- 95th percentile climatology ---
    for season in seasons:

        print(
            f"Creating 95th percentile df: {season}"
        )

        compute_area_avg_max_min(
            path_to_data,
            "95th_percentile_clim",
            season
        )

    # --- ROS plots ---
    for option in options:
        for season in seasons:

            print(
                f"Creating ROS frequency df: "
                f"{option} {season}"
            )

            compute_area_avg_max_min(
                path_to_data,
                "ros_frequency_clim",
                season,
                option,
            )

            print(
                f"Creating ROS intensity df: "
                f"{option} {season}"
            )

            compute_area_avg_max_min(
                path_to_data,
                "ros_intensity_clim",
                season,
                option,
            )