#!/usr/bin/env python3
"""
Filename:    landslide_ros_ar_analysis_batch.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: 
    Loop through landslides in database and:
    1) Generate a time series plot for each.
    2) Record max ROS and AR values during the event window.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import os
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cmocean.cm as cmo

# Local modules
sys.path.append('../modules/')
import globalvars

# ---------------------------------------------------------------------
# Update ROS choice
# ---------------------------------------------------------------------
study_start = "1981-01-01"
study_end   = "2019-12-31"
# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def find_nearest_indices(ds, lat, lon):
    """Find nearest grid indices on 2D curvilinear grid."""
    dist = (ds['lat'] - lat) ** 2 + (ds['lon'] - lon) ** 2
    iy, ix = np.unravel_index(dist.argmin(), dist.shape)
    return iy, ix

def calculate_climatology_exposure(
    ros_ds,
    ar_ds,
    lat,
    lon,
):
    """
    Calculate the number of AR and ROS days at a landslide location
    over the full climatological study period.

    ROS:
        ROS is identified at the nearest WRF/CFSR grid cell to the
        landslide location.

    AR:
        AR is identified if an AR occurs anywhere within +/- 1 degree
        latitude and longitude of the landslide location.

    Parameters
    ----------
    ros_ds : xarray.DataArray
        Daily ROS dataset with 2D lat/lon coordinates.
    ar_ds : xarray.DataArray
        AR dataset with lat/lon coordinates.
    lat : float
        Landslide latitude.
    lon : float
        Landslide longitude.
        
    Returns
    -------
    ros_days : int
        Number of ROS days at the landslide location.
    ar_days : int
        Number of AR days within the +/- 1 degree box.
    """

    # ---------------------------------------------------------------
    # ROS exposure
    # ---------------------------------------------------------------

    # Find nearest WRF/CFSR grid cell
    iy, ix = find_nearest_indices(ros_ds, lat, lon)

    ros_point = ros_ds.isel(y=iy, x=ix)
    print(ros_point)

    # Count days with ROS
    ros_days = int((ros_point == 1).sum().item())

    # ---------------------------------------------------------------
    # AR exposure
    # ---------------------------------------------------------------

    ar_sub = ar_ds.sel(
        lat=slice(lat + 1, lat - 1),
        lon=slice(lon - 1, lon + 1)
    )

    # Collapse to daily maximum and across the spatial box
    ar_sub = ar_sub.max(["lat", "lon"])

    # Count days with an AR
    ar_days = int((ar_sub > 0).sum().item())

    return ros_days, ar_days

# ---------------------------------------------------------------------
# Main Loop
# ---------------------------------------------------------------------
if __name__ == "__main__":
    # --- Load landslide database ---
    fname = globalvars.path_to_data + 'downloads/SEAK_News_Reported_Landslides.csv'
    df = pd.read_csv(fname)
    df = df.set_index(pd.to_datetime(df['Day_min'], format='%m/%d/%y', errors='coerce'))
    df = df.loc[(df.index >= study_start) & (df.index <= study_end)]

    # --- Get ROS days --- 
    ## how many days in the study period was there a ROS day across SEAK?
    fname = '../out/cfsr_ros_daily.nc'
    ds_ros = xr.open_dataset(fname)
    ds_ros = ds_ros.sel(time=slice(study_start, study_end))
    print(ds_ros['ros'])

    # --- Load AR dataset once ---
    fname = '../out/ar_ds_daily.nc'
    ar_ds = xr.open_dataset(fname)
    print(ar_ds['kidmap'])

    # --- Prepare outputs ---
    summary = []

    # --- Iterate through landslides ---
    for i, row in df.iterrows():
        lat = row['Lat']
        lon = row['Lon']
        
        # --- Calculate climatological exposure ---
        ros_clim_days, ar_clim_days = calculate_climatology_exposure(
            ds_ros['ros'],
            ar_ds['kidmap'],
            lat,
            lon,
        )

        start_date = pd.to_datetime(row['Day_min'], format='%m/%d/%y', errors='coerce') - pd.Timedelta(days=1)
        end_date = pd.to_datetime(row['Day_max'], format='%m/%d/%y', errors='coerce') + pd.Timedelta(days=1)

        print(start_date, end_date, lat, lon)
        # Subset WRF data
        ds_sub = ds_ros.sel(time=slice(start_date, end_date))
        iy, ix = find_nearest_indices(ds_sub, lat, lon)
        ds_point = ds_sub.isel(y=iy, x=ix)
        print(ds_point)
        # Subset AR data
        ar = ar_ds.sel(time=slice(start_date, end_date))
        ar = ar.sel(lat=slice(lat+1, lat-1), lon=slice(lon-1, lon+1))
        print(ar)
        ar = ar['kidmap'].max(["lat", "lon"])
        
        ar_mask = (ar > 0).values
        print(ar_mask)
        # Compute summary stats
        max_ros = float(ds_point['ros'].max())
        max_ar = float(ar_mask.max())
        summary.append({
            'start_date': start_date,
            'end_date': end_date,
            'lat': lat,
            'lon': lon,
            'max_ros': max_ros,
            'max_ar': max_ar,
            'ros_climatology_days': ros_clim_days,
            'ar_climatology_days': ar_clim_days
        })

    # --- Save summary DataFrame ---
    outname = f'../out/landslide_summary_strict_cfsr.csv'
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(outname, index=False)
    print(f"✅ Saved summary to: {outname}")


    # --- Compute ROS / AR stats dataframe ---
    total_ros_exposure = summary_df['ros_climatology_days'].sum()
    total_ar_exposure = summary_df['ar_climatology_days'].sum()
    
    n_ls_ros = summary_df['max_ros'].sum()
    n_ls_ar = summary_df['max_ar'].sum()

    n_landslides = len(summary_df)
    
    ros_percent = n_ls_ros / n_landslides * 100
    ar_percent = n_ls_ar / n_landslides * 100
    
    F_ros = n_ls_ros / total_ros_exposure
    F_ar = n_ls_ar / total_ar_exposure
    
    relative_frequency = F_ros / F_ar
    relative_increase = (relative_frequency - 1) * 100


    # --- Build Dataframe --- 
    table_df = pd.DataFrame({
        "Statistic": [
            "Total landslides",
            "AR days",
            "ROS days",
            "AR-associated landslides",
            "ROS-associated landslides",
            "AR-associated landslides (%)",
            "ROS-associated landslides (%)",
            "Landslides per AR day",
            "Landslides per ROS day",
            "ROS-to-AR frequency ratio",
            "ROS-to-AR frequency increase (%)",
        ],
        "Value": [
            n_landslides,
            total_ar_exposure,
            total_ros_exposure,
            n_ls_ar,
            n_ls_ros,
            ar_percent,
            ros_percent,
            F_ar,
            F_ros,
            relative_frequency,
            relative_increase,
        ]
    })
    
    # table_df["Value"] = [
    #     f"{int(n_landslides):d}",
    #     f"{total_ar_exposure:d}",
    #     f"{total_ros_exposure:d}",
    #     f"{n_ls_ar:d}",
    #     f"{n_ls_ros:d}",
    #     f"{ar_percent:.1f}",
    #     f"{ros_percent:.1f}",
    #     f"{F_ar:.4f}",
    #     f"{F_ros:.4f}",
    #     f"{relative_frequency:.2f}",
    #     f"{relative_increase:.1f}",
    # ]
    
    table_df.to_csv(
        f"../out/landslide_normalized_frequency_statistics_{study_start}_{study_end}.csv",
        index=False
    )
