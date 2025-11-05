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
from wrf_rio import wrf_prepare_for_rio

# ---------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------
def compute_area_avg_max_min(models, ssn, path_to_data, fsuffix):
    """
    Compute area-averaged statistics for each model and variable, including
    model minus CFSR differences.

    Parameters
    ----------
    models : list of str
        Model names (first should be 'cfsr' for baseline)
    ssn : str
        Season label (e.g., 'DJF')
    path_to_data : str
        Path to data directory
    fsuffix : str
        Filename suffix ('95th_percentile_clim', 'ros_intensity_clim', or 'ros_frequency_clim')

    Returns
    -------
    pandas.DataFrame
        DataFrame containing model, variable, polygon name, mean, min, max, and model–CFSR diff
    """

    # --- define variable names for each suffix ---
    if fsuffix == "flexible_ros_frequency_clim":
        varnames = ['ros', 'ivt', 'pcpt', 'delsnow', 'delsnowh']
    elif fsuffix == "flexible_ros_intensity_clim":
        varnames = ["ros", "pcpt", "delsnow", "delsnowh", "ros_intensity"]
    elif fsuffix == "95th_percentile_clim":
        varnames = ["ivt", "uv", "freezing_level", "pcpt", "snow"]
    else:
        raise ValueError(f"Unknown fsuffix: {fsuffix}")

    # --- load shapefile and subset ---
    fp = os.path.join(globalvars.path_to_data, 'downloads/AK_climate_divisions/AK_divisions_NAD83.shp')
    polys = gpd.read_file(fp)
    keep_names = ["Northeast Gulf", "North Panhandle", "Central Panhandle", "South Panhandle"]
    polys = polys[polys["Name"].isin(keep_names)]

    results = []

    # --- loop through variables ---
    for varname in varnames:
        # --- load CFSR reference ---
        if fsuffix == "95th_percentile_clim":
            cfsr_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/{varname}_cfsr_{ssn}_{fsuffix}.nc")
        else:
            cfsr_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/snow_cfsr_{ssn}_{fsuffix}.nc")

        if not os.path.exists(cfsr_path):
            print(f"⚠️ Missing CFSR file for {varname}: {cfsr_path}")
            continue

        ds_cfsr = xr.open_dataset(cfsr_path)
        if varname not in ds_cfsr:
            print(f"⚠️ Variable {varname} not in CFSR file")
            continue

        da_cfsr = ds_cfsr[varname]
        # prepare for rioxarray
        da_cfsr = wrf_prepare_for_rio(da_cfsr)

        # reproject shapefile to match WRF CRS
        polys = polys.to_crs(da_cfsr.rio.crs)

        # --- loop through each model (skip first since it's CFSR itself) ---
        for model in models:
            if model == "cfsr":
                continue

            if fsuffix == "95th_percentile_clim":
                model_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_{ssn}_{fsuffix}.nc")
            else:
                model_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/snow_{model}_{ssn}_{fsuffix}.nc")

            if not os.path.exists(model_path):
                print(f"⚠️ Missing file: {model_path}")
                continue

            ds_model = xr.open_dataset(model_path)
            if varname not in ds_model:
                print(f"⚠️ Variable {varname} not in {model_path}")
                continue

            da_model = ds_model[varname]
            # prepare for rioxarray
            da_model = wrf_prepare_for_rio(da_model)

            # --- compute model - CFSR difference before clipping ---
            # Align grids before subtracting
            da_model, da_cfsr_aligned = xr.align(da_model, da_cfsr, join="inner")
            da_diff = da_model - da_cfsr_aligned

            # print(f"\nVariable: {varname}, Model: {model}")
            # print(f"  Raster CRS: {da_model.rio.crs}")
            # print(f"  Polygon CRS: {polys.crs}")
            # print(f"  Raster bounds: {da_model.rio.bounds()}")
            # print(f"  Polygon bounds: {polys.geometry.bounds}")


            # --- loop through polygons ---
            for _, poly in polys.iterrows():
                # Clip both datasets to polygon
                model_clip = da_model.rio.clip([poly.geometry], polys.crs, drop=True)
                cfsr_clip  = da_cfsr_aligned.rio.clip([poly.geometry], polys.crs, drop=True)
                diff_clip  = da_diff.rio.clip([poly.geometry], polys.crs, drop=True)


                # Compute stats
                stats = lambda da: {
                    "mean": float(da.mean().values),
                    "min": float(da.min().values),
                    "max": float(da.max().values)
                }

                model_stats = stats(model_clip)
                cfsr_stats = stats(cfsr_clip)
                diff_stats = stats(diff_clip)

                results.append({
                    "season": ssn,
                    "fsuffix": fsuffix,
                    "variable": varname,
                    "region": poly["Name"],
                    "model": model,
                    "mean_model": model_stats["mean"],
                    "min_model": model_stats["min"],
                    "max_model": model_stats["max"],
                    "mean_cfsr": cfsr_stats["mean"],
                    "min_cfsr": cfsr_stats["min"],
                    "max_cfsr": cfsr_stats["max"],
                    "mean_diff": diff_stats["mean"],
                    "min_diff": diff_stats["min"],
                    "max_diff": diff_stats["max"]
                })

    df = pd.DataFrame(results)
    return df


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    ssn_lst = ["NDJFMA"]
    fsuffix_lst = ["95th_percentile_clim", "flexible_ros_intensity_clim", "flexible_ros_frequency_clim"]

    all_results = []
    for fsuffix in fsuffix_lst:
        for ssn in ssn_lst:
            print(f"Processing {fsuffix} for {ssn}")
            df = compute_area_avg_max_min(models, ssn, path_to_data, fsuffix)
            all_results.append(df)

    final_df = pd.concat(all_results, ignore_index=True)
    outfile = os.path.join(path_to_data, "processed_summary_stats_with_diffs.csv")
    final_df.to_csv(outfile, index=False)
    print(f"\n✅ Saved summary statistics with differences to {outfile}")
