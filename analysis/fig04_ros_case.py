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
import string
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
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
ds = preprocess_WRF_ros(ds, temporal_resolution='daily', option="flexible")
ivt = load_preprocessed_WRF_data('cfsr', 'ivt', anomaly=False)
ds = xr.merge([ds, ivt], compat="no_conflicts")

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
ivt_file = '/cw3e/mead/projects/cwp140/data/downloads/ERA5/20190101_IVT.nc'
ivt = xr.open_dataset(ivt_file)

rename_dict = {'viwve': 'uIVT', 'viwvn': 'vIVT', 
               'longitude': 'lon', 'latitude':'lat',
              'valid_time': 'time'}
ivt = ivt.rename(rename_dict)
ivt['IVT'] = np.sqrt(ivt['uIVT']**2 + ivt['vIVT']**2)

ivt = ivt.mean('time')


# ---------------------------------------------------------------------
# Setup Figure
# ---------------------------------------------------------------------
current_dpi = 300
base_dpi = 100
scaling_factor = (current_dpi / base_dpi) ** 0.3
set_font(current_dpi, scaling_factor)

fig = plt.figure(figsize=(12, 12), dpi=current_dpi)
gs = gridspec.GridSpec(
    nrows=2, ncols=3,
    width_ratios=[1, 1, 0.05],
    height_ratios=[1, 1],
    wspace=0.05, hspace=0.08
)
cax = fig.add_subplot(gs[0, -1])  # colorbar axis

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.Mercator()

# ---------------------------------------------------------------------
# (1) Landslide Image Panel
# ---------------------------------------------------------------------
ax_img = fig.add_subplot(gs[0, 0])
image = mpimg.imread('../img_landslide.jpg')
ax_img.imshow(image)
ax_img.axis('off')

# ---------------------------------------------------------------------
# (2) IVT Map Panel (North Pacific)
# ---------------------------------------------------------------------
ax_ivt = fig.add_subplot(gs[0, 1], projection=mapcrs)
ax_ivt.set_extent([-160, -125, 37, 60], crs=datacrs)

for feat in [cfeature.COASTLINE, cfeature.BORDERS]:
    ax_ivt.add_feature(feat, lw=0.8)
ax_ivt.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax_ivt.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)

levels_ivt = np.arange(250, 950, 100)
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
    uvec.where(mask).values, vvec.where(mask).values,
    transform=datacrs,
    color='k',
    regrid_shape=15,
    scale=7000  # bigger = shorter arrows
)
ax_ivt.quiverkey(Q, X=0.8, Y=0.03, U=250, label='250 kg m$^{-1}$ s$^{-1}$',
                 labelpos='E', coordinates='axes', fontproperties={'size':6})


# --- add landslide point ---
ax_ivt.plot(lon_lst[0], lat_lst[0], 'ro', alpha=0.4,
            markersize=5, transform=datacrs, zorder=201)

# ---------------------------------------------------------------------
# (3) ROS Map Panel (Southeast Alaska)
# ---------------------------------------------------------------------
ax_ros = fig.add_subplot(gs[1, 0], projection=mapcrs)
ax_ros.set_extent([-141, -130, 54.5, 60], crs=datacrs)

# Base map
for feat in [cfeature.COASTLINE, cfeature.BORDERS]:
    ax_ros.add_feature(feat, lw=0.8)
ax_ros.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax_ros.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)

# Contour plot
levels_ros = np.arange(1, 6, 3)
cf_ros = ax_ros.contourf(
    ds['lon'], ds['lat'], ds['ros'].sum('time'),
    levels=levels_ros, cmap='Blues',
    transform=datacrs, extend='neither'
)

ax_ros.plot(lon_lst[0], lat_lst[0], 'ro', alpha=0.4,
            markersize=5, transform=datacrs, zorder=201)




# ---------------------------------------------------------------------
# (4) Plot ROS Time Series
# ---------------------------------------------------------------------
ax1 = fig.add_subplot(gs[1, 1:])

# Ensure arrays are 1D
prec = np.squeeze(ds_point['pcpt'].to_numpy())
swe = np.squeeze(ds_point['snow'].to_numpy())
time = pd.to_datetime(ds_point['time'].values)
ivt_ts = ds_point['ivt'].values
deltasnowh = ds_point['delsnowh'].values*-1
ros = ds_point['ros'].values

# Plot bars
bar_width = 0.4
x = np.arange(len(time))  # or use time index if it's regular daily
ax1.bar(x - bar_width/2, prec, color='#1f77b4', width=bar_width, label='Precip (mm)')
ax1.bar(x + bar_width/2, swe, color='#89CFF0', width=bar_width, label='SWE (mm)')
ax1.set_xticks(x)
ax1.set_xticklabels([pd.to_datetime(t).strftime('%b %d') for t in time.values])

# --- Overlay ΔSnowH line ---
ax1.plot(x, deltasnowh, label='ΔSnow Height (mm)', color='tab:purple', linewidth=1.8)

ax1.set_ylabel('Precipitation / Snow (mm)')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xlim(x.min() + .5, x.max() - 0.5)
ax1.set_ylim(-39, 39)
ax1.margins(x=0.05)

# --- IVT line on right axis ---
ax2 = ax1.twinx()
ax2.plot(x, ivt_ts, color='#62bfa4', label='IVT (kg m$^{-1}$ s$^{-1}$)', linewidth=1.8)
ax2.set_ylabel('IVT (kg m$^{-1}$ s$^{-1}$)', color='#62bfa4')
ax2.tick_params(axis='y', labelcolor='#62bfa4')
ax2.set_ylim(0, 500)

# --- Shade where ROS == 1 (centered on each day) ---
ros_mask = ds_point['ros'].values.squeeze().astype(bool)
time_np = x
half_day = 0.5

for t, r in zip(time_np, ros_mask):
    if r:
        ax1.axvspan(t - half_day, t + half_day, color='lightgray', alpha=0.4, zorder=0)

# --- Add horizontal reference line at 0 ---
ax1.axhline(0, color='black', linestyle='--', linewidth=1)

# --- Combine legends from both axes ---
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()

# --- Legend at bottom with 4 columns ---
ax1.legend(lines_1 + lines_2, labels_1 + labels_2,
    loc='lower center',      # Place below plot
    bbox_to_anchor=(0.5, -0.12),  # Center it horizontally, move down a bit
    ncol=4,                  # Four columns
    frameon=False,           # No box around legend
    fontsize=8,              # Slightly smaller font
    handlelength=1.2,        # Shorter legend lines
)
# --- Formatting tweaks ---
ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.4)
pos = ax1.get_position()
ax1.set_position([pos.x0 + 0.05, pos.y0 + 0.015, pos.width * 0.9, pos.height * 0.925])




# ---------------------------------------------------------------------
# Add a, b, c labels
# ---------------------------------------------------------------------
labels = list(string.ascii_lowercase)  # ['a', 'b', 'c', ...]
bbox_dict = dict(facecolor='white', edgecolor='k', boxstyle='circle,pad=0.3', alpha=1.)
axes = [ax_img, ax_ivt, ax_ros, ax1]  # your panel axes

# Loop over axes and labels
for ax, label in zip(axes, labels):
    ax.text(
        0.05, 0.96, label,
        transform=ax.transAxes,  # position relative to the axes
        va='top', ha='left',
        fontsize=12,
        bbox=bbox_dict
    )
    
# ---------------------------------------------------------------------
# Colorbar & Save
# ---------------------------------------------------------------------
cb = fig.colorbar(cf_ivt, cax=cax, orientation='vertical')
cb.set_label('IVT (kg m$^{-1}$ s$^{-1}$)', fontsize=11)
cb.ax.tick_params(labelsize=10)

output_path = Path('../figs/ros_case.png')
fig.savefig(output_path, bbox_inches='tight', dpi=fig.dpi)
print(f"✅ Saved figure to: {output_path.resolve()}")
