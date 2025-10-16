#!/usr/bin/env python3
"""
Filename:    fig4_ros_case.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: 
    Create a four-panel figure showing:
      (1) ROS time series,
      (2) Southeast Alaska ROS map,
      (3) North Pacific IVT map,
      (4) Landslide image.
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import os
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import cmocean.cm as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Local modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data
from wrf_preprocess import preprocess_WRF_ros
from plotter import set_font

# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def find_nearest_indices(ds, lat, lon):
    """Find nearest grid indices on 2D curvilinear grid."""
    dist = (ds['lat'] - lat) ** 2 + (ds['lon'] - lon) ** 2
    iy, ix = np.unravel_index(dist.argmin(), dist.shape)
    return iy, ix


# ---------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------
# --- WRF data ---
ds = load_preprocessed_WRF_data('cfsr', 'snow', anomaly=False)
ds = preprocess_WRF_ros(ds, temporal_resolution='daily')
ivt = load_preprocessed_WRF_data('cfsr', 'ivt', anomaly=False)
ds = xr.merge([ds, ivt])

# --- Case selection ---
start_date, end_date = '2018-12-30', '2019-01-04'
ds = ds.sel(time=slice(start_date, end_date))

# --- Location ---
lat_lst = [55.683099]
lon_lst = [-132.571308]
indices = [find_nearest_indices(ds, la, lo) for la, lo in zip(lat_lst, lon_lst)]
iy = [i[0] for i in indices]
ix = [i[1] for i in indices]
ds_point = ds.isel(y=xr.DataArray(iy, dims="location"),
                   x=xr.DataArray(ix, dims="location"))

# --- IVT Data ---
ivt_files = [
    '/cw3e/mead/projects/cwp162/data/preprocessed/ERA5/ivt/201812_IVT.nc',
    '/cw3e/mead/projects/cwp162/data/preprocessed/ERA5/ivt/201901_IVT.nc'
]
ivt = xr.open_mfdataset(ivt_files, engine='netcdf4')
ivt = ivt.sel(time=slice(start_date, end_date)).mean('time')


# ---------------------------------------------------------------------
# Setup Figure
# ---------------------------------------------------------------------
current_dpi = 300
base_dpi = 100
scaling_factor = (current_dpi / base_dpi) ** 0.3
set_font(current_dpi, scaling_factor)

fig = plt.figure(figsize=(10, 10), dpi=current_dpi)
gs = gridspec.GridSpec(
    nrows=3, ncols=2,
    width_ratios=[1, 1],
    height_ratios=[1, 1, 0.05],
    wspace=0.05, hspace=0.2
)
cax = fig.add_subplot(gs[-1, -1])  # colorbar axis

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.Mercator()


# ---------------------------------------------------------------------
# (1) ROS Time Series Panel
# ---------------------------------------------------------------------
ax1 = fig.add_subplot(gs[0, 0])

# Extract point time series
time = ds_point['time']
ivt_ts = ds_point['ivt']
rain = ds_point['rain']
snow = ds_point['snow']
deltasnowh = ds_point['delsnowh']
ros = ds_point['ros']

# Left axis: precip and snow variables
ax1.plot(time, rain, label='Rain (mm)', color='tab:blue', lw=1.5)
ax1.plot(time, snow, label='Snow (mm)', color='tab:cyan', lw=1.5)
ax1.plot(time, deltasnowh, label='ΔSnow Depth (mm)', color='tab:orange', lw=1.5)
ax1.set_ylabel('Precipitation / Snow (mm)', color='black')

# Right axis: IVT
ax2 = ax1.twinx()
ax2.plot(time, ivt_ts, label='IVT (kg m$^{-1}$ s$^{-1}$)',
         color='tab:purple', lw=1.8)
ax2.set_ylabel('IVT (kg m$^{-1}$ s$^{-1}$)', color='tab:purple')
ax2.tick_params(axis='y', labelcolor='tab:purple')

# Shade where ROS == 1 (centered)
ros_mask = ros.values.squeeze().astype(bool)
time_np = np.array(time.values)
half_day = np.timedelta64(12, 'h')
for t, r in zip(time_np, ros_mask):
    if r:
        ax1.axvspan(t - half_day, t + half_day, color='lightgray', alpha=0.5, zorder=0)

ax1.axhline(0, color='black', linestyle='--', lw=1)
ax1.set_xlabel('Time')

# Legend
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', frameon=False)


# ---------------------------------------------------------------------
# (2) ROS Map Panel (Southeast Alaska)
# ---------------------------------------------------------------------
ax_ros = fig.add_subplot(gs[1, 0], projection=mapcrs)
ax_ros.set_extent([-142, -130, 54, 61], crs=datacrs)

# Base map
for feat in [cfeature.COASTLINE, cfeature.BORDERS]:
    ax_ros.add_feature(feat, lw=0.8)
ax_ros.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax_ros.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)

# Contour plot
levels_ros = np.arange(0, 3.5, 0.5)
cf_ros = ax_ros.contourf(
    ds['lon'], ds['lat'], ds['ros'].sum('time'),
    levels=levels_ros, cmap='Blues',
    transform=datacrs, extend='neither'
)


# ---------------------------------------------------------------------
# (3) IVT Map Panel (North Pacific)
# ---------------------------------------------------------------------
ax_ivt = fig.add_subplot(gs[1, 1], projection=mapcrs)
ax_ivt.set_extent([-160, -120, 40, 60], crs=datacrs)

for feat in [cfeature.COASTLINE, cfeature.BORDERS]:
    ax_ivt.add_feature(feat, lw=0.8)
ax_ivt.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax_ivt.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)

levels_ivt = np.arange(250, 1250, 250)
cf_ivt = ax_ivt.contourf(
    ivt['lon'], ivt['lat'], ivt['IVT'],
    levels=levels_ivt, cmap=cmo.deep,
    transform=datacrs, extend='neither'
)

# IVT vectors
uvec = ivt['uIVT']
vvec = ivt['vIVT']
mask = ivt['IVT'] >= 250
Q = ax_ivt.quiver(
    ivt['lon'], ivt['lat'],
    uvec.where(mask), vvec.where(mask),
    transform=datacrs, color='k', regrid_shape=20,
    scale=250, scale_units='xy'
)
ax_ivt.quiverkey(Q, 0.02, -0.1, 250, '250 kg m$^{-1}$ s$^{-1}$', labelpos='E',
                 coordinates='axes', fontproperties={'size': 6})


# ---------------------------------------------------------------------
# (4) Landslide Image Panel
# ---------------------------------------------------------------------
ax_img = fig.add_subplot(gs[0, 1])
image = mpimg.imread('img_landslide.jpg')
ax_img.imshow(image)
ax_img.axis('off')


# ---------------------------------------------------------------------
# Colorbar & Save
# ---------------------------------------------------------------------
cb = fig.colorbar(cf_ivt, cax=cax, orientation='horizontal')
cb.set_label('IVT (kg m$^{-1}$ s$^{-1}$)', fontsize=11)
cb.ax.tick_params(labelsize=10)

output_path = Path('../figs/fig4_ros_case.png')
fig.savefig(output_path, bbox_inches='tight', dpi=fig.dpi)
print(f"✅ Saved figure to: {output_path.resolve()}")
