import sys
import string
from pathlib import Path
import xarray as xr
import pandas as pd
import numpy as np

# --- Personal Modules ---
sys.path.append('../modules')
import globalvars
from utils import get_startmon_and_endmon, select_months_df, select_months_ds
import rioxarray  # for spatial clipping
from wrf_rio import wrf_prepare_for_rio
from wrf_crs import load_climate_division_shapefile

# --- Helper Functions --- 
def clip_ds_to_shapefile(ds, x_dim="x", y_dim="y"):

    # --- load shapefile and subset ---
    polys = load_climate_division_shapefile()
    # Dissolve into one multipolygon
    study_area = polys.dissolve()
    
    # --- prepare for rioxarray ---
    da = wrf_prepare_for_rio(ds, x_dim, y_dim)
    
    # --- Reproject once after reading the raster ---
    study_area = study_area.to_crs(da.rio.crs)

    # --- Clip both datasets to polygon ---
    da_clip = da.rio.clip(
        study_area.geometry,
        study_area.crs,
        drop=True
    )

    fill_value = -9223372036854775808
    da_clip = da_clip.astype(float)
    da_clip = da_clip.where(da_clip != fill_value, np.nan)

    return da_clip

# --- Config ---
season="ONDJFM"
study_start = "1981-01-01"
study_end   = "2019-12-31"

# --- Get ROS days --- 
## how many days in the study period was there a ROS day across SEAK?
fname = '../out/cfsr_ros_daily.nc'
ds_ros = xr.open_dataset(fname)
ds_ros = ds_ros.sel(time=slice(study_start, study_end))

print("Original ROS values:")
print(np.unique(ds_ros["ros"].values))

ros = clip_ds_to_shapefile(ds_ros["ros"])

print("Clipped ROS values:")
print(np.unique(ros.values))
ros_days = (ros == 1).any(dim=["y", "x"])

# --- debugging ---
ros_count_per_day = (ros == 1).sum(dim=["y", "x"])

print("ROS grid cells per day:")
print(f"  min:    {ros_count_per_day.min().item():.0f}")
print(f"  max:    {ros_count_per_day.max().item():.0f}")
print(f"  mean:   {ros_count_per_day.mean().item():.2f}")
print(f"  median: {ros_count_per_day.median().item():.0f}")
print(f"  total:  {ros_count_per_day.sum().item():.0f}")

print("\nNumber of days with ROS:")
print(f"  {int((ros_count_per_day > 0).sum().item()):,}")

print("\nROS days with number of affected grid cells:")
print(
    ros_count_per_day
    .where(ros_count_per_day > 0, drop=True)
    .to_series()
    .describe()
)

n_ros_gridcell_days = int((ros == 1).sum().item())
n_valid_gridcell_days = int(ros.notnull().sum().item())

print(f"ROS grid-cell days: {n_ros_gridcell_days:,}")
print(f"Valid grid-cell days: {n_valid_gridcell_days:,}")
print(f"Fraction of grid-cell days with ROS: "
      f"{n_ros_gridcell_days / n_valid_gridcell_days * 100:.2f}%")

print(f"ROS days: {int(ros_days.sum().item()):,}")
print(f"Fraction of study days with ROS: "
      f"{ros_days.sum().item() / ros.sizes['time'] * 100:.2f}%")

## how many days in the study period was there an AR across SEAK?
# --- Load AR dataset once ---
ar_ds = xr.open_dataset(globalvars.path_to_data + 'downloads/globalARcatalog_ERA5_1940-2024_v4.0.nc')
ar_ds = ar_ds.kidmap.squeeze()
ar_ds = ar_ds.assign_coords({"lon": (((ar_ds.lon + 180) % 360) - 180)}).sortby('lon')
ar_ds = ar_ds.sel(time=slice(study_start, study_end))
ar_ds = ar_ds.sel(lat=slice(62., 54.), lon=slice(-146., -127.))
ar_ds = ar_ds.resample(time="1D").max()

ar_ds = clip_ds_to_shapefile(ar_ds, x_dim="lon", y_dim="lat")
ar_days = ar_ds.notnull().any(dim=["lat", "lon"])

# --- Compute relative normalized frequency ---
fname = '../out/landslide_summary_strict.csv'
df = pd.read_csv(fname)

# --- Summary of Numbers ---
n_ros_days = int(ros_days.sum().item())
n_ar_days = int(ar_days.sum().item())
print(f"Number of ROS days: {n_ros_days}")
print(f"Number of AR days: {n_ar_days}")

n_ls_ros = int(df['max_ros'].sum())
n_ls_ar = int(df['max_ar'].sum())

n_landslides = len(df)

ros_percent = n_ls_ros / n_landslides * 100
ar_percent = n_ls_ar / n_landslides * 100

F_ros = n_ls_ros / n_ros_days
F_ar = n_ls_ar / n_ar_days

relative_frequency = F_ros / F_ar
relative_increase = (F_ros / F_ar - 1) * 100

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
        n_ar_days,
        n_ros_days,
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

table_df["Value"] = [
    f"{n_landslides:d}",
    f"{n_ar_days:d}",
    f"{n_ros_days:d}",
    f"{n_ls_ar:d}",
    f"{n_ls_ros:d}",
    f"{ar_percent:.1f}",
    f"{ros_percent:.1f}",
    f"{F_ar:.2f}",
    f"{F_ros:.2f}",
    f"{relative_frequency:.2f}",
    f"{relative_increase:.1f}",
]

table_df.to_csv(
    f"../out/landslide_normalized_frequency_statistics_{study_start}_{study_end}.csv",
    index=False
)