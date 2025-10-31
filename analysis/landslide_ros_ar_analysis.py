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
from wrf_utils import load_preprocessed_WRF_data
from wrf_preprocess import preprocess_WRF_ros
from plotter import set_font

# ---------------------------------------------------------------------
# Update ROS choice
# ---------------------------------------------------------------------
option = 'strict' # strict or flexible

# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def find_nearest_indices(ds, lat, lon):
    """Find nearest grid indices on 2D curvilinear grid."""
    dist = (ds['lat'] - lat) ** 2 + (ds['lon'] - lon) ** 2
    iy, ix = np.unravel_index(dist.argmin(), dist.shape)
    return iy, ix

def plot_ros_ar_timeseries(ds_point, ar, start_date, end_date, outpath):
    """Plot ROS, AR, IVT, Rain, Snow, and ΔSnowH time series."""
    print(start_date, end_date)
    current_dpi = 300
    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi) ** 0.3
    set_font(current_dpi, scaling_factor)

    fig = plt.figure(figsize=(12, 6), dpi=current_dpi)
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax1 = fig.add_subplot(gs[0, 0])

    # Extract variables
    prec = np.squeeze(ds_point['pcpt'].to_numpy())
    deltaswe = np.squeeze(ds_point['delsnow'].to_numpy())
    deltasnowh = np.squeeze(ds_point['delsnowh'].to_numpy())
    ivt_ts = np.squeeze(ds_point['ivt'].to_numpy())
    ros = np.squeeze(ds_point['ros'].to_numpy()).astype(bool)
    time = pd.to_datetime(ds_point['time'].values)

    # Plot bars
    bar_width = 0.4
    x = np.arange(len(time))  # or use time index if it's regular daily
    ax1.bar(x - bar_width/2, prec, color='#1f77b4', width=bar_width, label='Precip (mm)')
    ax1.bar(x + bar_width/2, deltaswe, color='#89CFF0', width=bar_width, label='ΔSWE (mm)')
    ax1.set_xticks(x)
    ax1.set_xticklabels([pd.to_datetime(t).strftime('%b %d') for t in time.values])

    ax1.set_ylabel('Precipitation / Snow (mm)')
    ax1.set_xticks(x)
    ax1.set_xticklabels([t.strftime('%b %d') for t in time])
    ax1.axhline(0, color='black', linestyle='--', linewidth=1)
    ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.4)
    # ax1.set_ylim(-40, 40)

    # IVT line
    ax2 = ax1.twinx()
    ax2.plot(x, ivt_ts, color='#62bfa4', linewidth=1.8, label='IVT (kg m$^{-1}$ s$^{-1}$)')
    ax2.set_ylabel('IVT (kg m$^{-1}$ s$^{-1}$)', color='#62bfa4')
    ax2.tick_params(axis='y', labelcolor='#62bfa4')
    # ax2.set_ylim(0, 500)

    # Shading for ROS and AR
    half_day = 0.5
    for t, r in zip(x, ros):
        if r:
            ax1.axvspan(t - half_day, t + half_day, color='lightgray', alpha=0.4, zorder=0, ymin=0.5, ymax=1)
    for t, a in zip(x, ar):
        if a:
            ax1.axvspan(t - half_day, t + half_day, color='tab:green', alpha=0.4, zorder=0, ymin=0, ymax=0.5)

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='lower center', bbox_to_anchor=(0.5, -0.12),
               ncol=4, frameon=False, fontsize=8, handlelength=1.2)

    # Save
    fig.savefig(outpath, bbox_inches='tight', dpi=current_dpi)
    plt.close(fig)
    print(f"✅ Saved: {outpath}")

# ---------------------------------------------------------------------
# Main Loop
# ---------------------------------------------------------------------
if __name__ == "__main__":
    # --- Load landslide database ---
    fname = globalvars.path_to_data + 'downloads/SEAK_News_Reported_Landslides.csv'
    df = pd.read_csv(fname)
    df = df.set_index(pd.to_datetime(df['Day_min'], format='%m/%d/%y', errors='coerce'))
    df = df.loc[(df.index >= '1980-01-01') & (df.index <= '2019-12-31')]

    # --- Load WRF data once ---
    ds = load_preprocessed_WRF_data('cfsr', 'snow', anomaly=False)
    ds = preprocess_WRF_ros(ds, temporal_resolution='daily', option=option)
    ivt = load_preprocessed_WRF_data('cfsr', 'ivt', anomaly=False)
    ds = xr.merge([ds, ivt], compat="no_conflicts")

    # --- Load AR dataset once ---
    ar_ds = xr.open_dataset(globalvars.path_to_data + 'downloads/globalARcatalog_ERA5_1940-2024_v4.0.nc')
    ar_ds = ar_ds.kidmap.squeeze()
    ar_ds = ar_ds.assign_coords({"lon": (((ar_ds.lon + 180) % 360) - 180)}).sortby('lon')

    # --- Prepare outputs ---
    figs_dir = Path(f'../figs/landslides_{option}/')
    figs_dir.mkdir(parents=True, exist_ok=True)
    summary = []

    # --- Iterate through landslides ---
    for i, row in df.iterrows():
        start_date = pd.to_datetime(row['Day_min'], format='%m/%d/%y', errors='coerce') - pd.Timedelta(days=1)
        end_date = pd.to_datetime(row['Day_max'], format='%m/%d/%y', errors='coerce') + pd.Timedelta(days=1)
        lat = row['Lat']
        lon = row['Lon']
        print(start_date, end_date, lat, lon)
        # Subset WRF data
        ds_sub = ds.sel(time=slice(start_date, end_date))
        iy, ix = find_nearest_indices(ds_sub, lat, lon)
        ds_point = ds_sub.isel(y=iy, x=ix)

        # Subset AR data
        ar = ar_ds.sel(time=slice(start_date, end_date))
        # ar = ar.sel(lat=lat, lon=lon, method='nearest')
        ar = ar.sel(lat=slice(lat+1, lat-1), lon=slice(lon-1, lon+1))
        ar = ar.resample(time="1D").max(["lat", "lon"])
        ar_mask = (ar > 0).values.squeeze()

        # Save figure
        outpath = figs_dir / f"landslide_{i.strftime('%Y%m%d')}_{lat:.2f}_{lon:.2f}.png"
        plot_ros_ar_timeseries(ds_point, ar_mask, start_date, end_date, outpath)

        # Compute summary stats
        max_ros = float(ds_point['ros'].max())
        max_ar = float(ar_mask.max())
        summary.append({
            'start_date': start_date,
            'end_date': end_date,
            'lat': lat,
            'lon': lon,
            'max_ros': max_ros,
            'max_ar': max_ar
        })

    # --- Save summary DataFrame ---
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f'../out/landslide_summary_{option}.csv', index=False)
    print(f"✅ Saved summary to: {figs_dir / 'landslide_summary.csv'}")
