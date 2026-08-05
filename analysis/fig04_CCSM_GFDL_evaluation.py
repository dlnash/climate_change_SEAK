#!/usr/bin/env python3
"""
Filename:    fig04_CCSM_GFDL_evaluation.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Compare annual 95th percentile precipitation from historical CCSM and GFDL simulations to four station observation locations
"""

import sys
import string
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt

# --- Personal Modules ---
sys.path.append('../modules')
import globalvars
from utils import get_startmon_and_endmon, select_months_df, select_months_ds
from wrf_preprocess import load_preprocessed_WRF_data

# --- HELPER FUNCTIONS ---
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # km

    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return R * c
    
def nearest_grid(ds, lat0, lon0, method="haversine"):
    if method == "haversine":
        dist = haversine(lat0, lon0, ds.lat.values, ds.lon.values)
    else:
        dist = (ds.lat.values - lat0) ** 2 + (ds.lon.values - lon0) ** 2

    iy, ix = np.unravel_index(np.argmin(dist), dist.shape)

    return iy, ix

# --- CONFIG ---
season = 'ONDJFM'
station_data = {
    'JUNEAU AIRPORT, AK US':
    {
        'latitude': 58.354,
        'longitude': -134.55606,
        'elevation': 5.9, # m
    },
    'YAKUTAT AIRPORT, AK US':
    {
            'latitude': 59.51211, 
            'longitude': -139.67103,
            'elevation': 12.6, # m
        },
    'AUKE BAY, AK US':
    {
            'latitude': 58.3814, 
            'longitude': -134.645, 
            'elevation': 13.4, # m
        },

    'ANNETTE AIRPORT, AK US':
    {
            'latitude': 55.0389,
            'longitude': -131.5787,
            'elevation': 33.2, # m
        },

}

# --- Read station data ---
fname = Path(globalvars.path_to_data) / 'downloads' / 'GHCN-Daily_SEAK.csv'
df = pd.read_csv(fname)

gauge_timeseries = {}
for name in df['NAME'].unique():
    idx = (df['NAME'] == name)
    station = df.loc[idx].copy()
    # --- set DATE to index
    station['DATE'] = pd.to_datetime(station['DATE'])
    station = station.set_index('DATE')

    # --- remove values below 0.01 ---
    # station = station[station["PRCP"] * 25.4 >= 0.1] 

    # --- subset to specified season ---
    mon_s, mon_e = get_startmon_and_endmon(season)
    df_ssn = select_months_df(station, mon_s, mon_e)
    
    # --- Compute annual 95th percentile of PRCP grouped by year ---
    df_95th = df_ssn.groupby(df_ssn.index.year)['PRCP'].quantile(0.95)
    gauge_timeseries[name] = df_95th * 25.4 # convert from inches to mm


# --- Read WRF data ---
model_lst = ['ccsm', 'gfdl']
records = []

for model in model_lst:

    ds = load_preprocessed_WRF_data(model, "hist", "pcpt")

    for station_name, info in station_data.items():

        lat = info["latitude"]
        lon = info["longitude"]
        print(lat, lon)

        iy, ix = nearest_grid(ds, lat, lon)

        ds_station = ds.pcpt.isel(
            y=iy,
            x=ix,
        )
        print(ds_station.lat.values, ds_station.lon.values)

        # --- subset to specified season ---
        mon_s, mon_e = get_startmon_and_endmon(season)
        ds_station = select_months_ds(ds_station, mon_s, mon_e, time_varname='time')

        # --- remove values below 0.01 ---
        ds_station = ds_station.compute()
        # ds_station = ds_station.where(ds_station >= 0.1, drop=True)
    
        # --- compute 95th percentile for var for each year ---
        station_95 = ds_station.groupby("time.year").quantile(0.95, dim="time")

        for year in station_95.year.values:
        
            records.append(
                {
                    "station": station_name,
                    "model": model,
                    "year": int(year),
                    "wrf95": float(
                        station_95.sel(year=year).values
                    ),
                    "obs95": float(
                        gauge_timeseries[station_name].loc[year]
                    ),
                }
            )

validation_df = pd.DataFrame(records)

labels = list(string.ascii_lowercase)  # ['a', 'b', 'c', ...]
bbox_dict = dict(facecolor='white', edgecolor='k', boxstyle='circle,pad=0.3', alpha=1.)
stations = {
    "JUNEAU AIRPORT, AK US": ("o", "tab:blue"),
    "YAKUTAT AIRPORT, AK US": ("s", "tab:orange"),
    "AUKE BAY, AK US": ("^", "tab:green"),
    "ANNETTE AIRPORT, AK US": ("D", "tab:red"),
}

stats_records = []

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10,5),
    sharex=True,
    sharey=True,
)

for i, (ax, model) in enumerate(zip(axes, ["ccsm", "gfdl"])):

    sub = validation_df[validation_df.model == model]

    # -----------------------------
    # Plot each station
    # -----------------------------
    for station, (marker, color) in stations.items():

        s = sub[sub.station == station]

        ax.scatter(
            s.obs95,
            s.wrf95,
            marker=marker,
            color=color,
            edgecolor="k",
            s=60,
            alpha=0.8,
            label=station.split(",")[0],
        )

    # -----------------------------
    # Overall statistics
    # -----------------------------
    bias = (sub.wrf95 - sub.obs95).mean()

    rmse = np.sqrt(
        mean_squared_error(
            sub.obs95,
            sub.wrf95,
        )
    )

    r = np.corrcoef(
        sub.obs95,
        sub.wrf95,
    )[0,1]

    stats_records.append(
        {
            "Model": model.upper(),
            "Bias (mm/day)": bias,
            "RMSE (mm/day)": rmse,
            "Correlation": r,
            "N": len(sub),
        }
    )

    stats_text = (
        f"$r$ = {r:.2f}\n"
        f"RMSE = {rmse:.1f}\n"
        f"Bias = {bias:.1f}"
    )

    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(
            facecolor="white",
            edgecolor="0.7",
            alpha=0.9,
        ),
    )

    # -----------------------------
    # 1:1 line
    # -----------------------------
    mn = min(sub.obs95.min(), sub.wrf95.min())
    mx = max(sub.obs95.max(), sub.wrf95.max())

    ax.plot(
        [mn, mx],
        [mn, mx],
        "k--",
        lw=1,
    )

    ax.set_xlim(mn-1.5, mx+1.5)
    ax.set_ylim(mn-1.5, mx+1.5)

    ax.set_title(f" Historical {model.upper()}")
    ax.set_xlabel("Observed annual 95th percentile (mm day$^{-1}$)")

    # a, b labels
    ax.text(
        0.05, 0.96, labels[i],
        transform=ax.transAxes,  # position relative to the axes
        va='top', ha='left',
        fontsize=12,
        bbox=bbox_dict
    )

axes[0].set_ylabel("WRF annual 95th percentile (mm day$^{-1}$)")

handles, labels = axes[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))

fig.legend(
    by_label.values(),
    by_label.keys(),
    loc="lower center",
    ncol=4,
)

plt.tight_layout(rect=[0,0.08,1,1])

plt.show()

output_path = Path('../figs/CCSM_GFDL_evaluation.png')
fig.savefig(output_path, bbox_inches='tight', dpi=fig.dpi)

print(f"Saved figure to: {output_path.resolve()}")

# =====================================================
# Create summary dataframe for supplement
# =====================================================

stats_df = pd.DataFrame(stats_records)

print(stats_df)

stats_df.to_latex(
    "validation_statistics.tex",
    index=False,
    float_format="%.2f",
)