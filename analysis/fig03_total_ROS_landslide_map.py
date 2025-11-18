#!/usr/bin/env python3
"""
Filename:    fig03_total_ROS_landslide_map.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Plot a topographic map of Southeast Alaska with labeled terrain features and communities.
"""

# --- Standard Library Imports ---
import sys
from pathlib import Path
import textwrap

# --- Third-Party Imports ---
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo

# --- Personal Modules ---
sys.path.append('../modules')
from plotter import draw_basemap, set_font
import globalvars
from wrf_utils import load_preprocessed_WRF_data
from wrf_preprocess import preprocess_WRF_ros

# ---------------------------------------------------------------------
# Update ROS choice
# ---------------------------------------------------------------------
option = 'strict' # strict or flexible

# ---------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------
# --- read the non-anomaly snow data ---
ds = load_preprocessed_WRF_data('cfsr', 'snow', anomaly=False)

# --- compute ROS information ---
ds = preprocess_WRF_ros(ds, temporal_resolution='daily', option=option, season='ONDJFM')

# --- Sum over time: number of ROS days per grid cell ---
ros = ds['ros'].groupby('time.year').sum('time').rename({'year': 'time'})
ros.attrs['units'] = "ROS (d yr$^{{-1}}$)"
ros_sum = ros.mean('time', keep_attrs=True)
print(ros_sum.values.max())
print(ros_sum.quantile(q=0.98).values)
# print(ros_sum.mean(skipna=True).values)

# --- read summary of landslide information --- 
df = pd.read_csv(f'../out/landslide_summary_{option}.csv')
# Sort so that points with ar=1 and ros=1 are plotted last
# df = df.sort_values(by=['max_ros'], ascending=[False])

df2 = df.loc[df['max_ros']==1.]

df_lst = [df, df2]
# ---------------------------------------------------------------------
# Map and Label Definitions
# ---------------------------------------------------------------------

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.Mercator()
extent = [-141., -130., 54.5, 60.]
xticks = np.arange(extent[0], extent[1] + 2, 2)
yticks = np.arange(extent[2] - 0.5, extent[3] + 1, 1)


# ---------------------------------------------------------------------
# Create Figure and GridSpec Layout
# ---------------------------------------------------------------------
current_dpi=300
base_dpi=100
scaling_factor = (current_dpi / base_dpi)**0.3
set_font(current_dpi, scaling_factor)
fig = plt.figure(figsize=(8, 6.5), dpi=current_dpi)
gs = gridspec.GridSpec(
    nrows=1, ncols=2,  # one for map, one for colorbar
    width_ratios=[1, 0.04],  # colorbar narrower
    wspace=0.05
)

# Main map axis
ax = fig.add_subplot(gs[0, 0], projection=mapcrs)

# Colorbar axis
cax = fig.add_subplot(gs[0, 1])

# ---------------------------------------------------------------------
# Plot Data
# ---------------------------------------------------------------------

# Draw base map
ax = draw_basemap(
    ax, extent=extent, xticks=xticks, yticks=yticks,
    left_lats=True, right_lats=False,
    mask_ocean=False, coastline=True
)

# --- Contour plot ---
# Define contour levels for discrete ROS counts (0, 1, 2, 3)
levels = np.arange(1, 11, 1)  # finer levels between 0–3
cmap = cmo.dense

cf = ax.contourf(
    ros_sum['lon'], ros_sum['lat'], ros_sum.values,
    levels=levels,
    cmap=cmap,
    extend='max',
    transform=ccrs.PlateCarree()
)
for df in df_lst:
    # Plot landslide markers
    for i, row in df.iterrows():
        x = row['lon']
        y = row['lat']
        # Determine fill color
        marker_color = 'yellow' if row['max_ar'] == 1 else 'gray'
        
        # Determine outline color and width
        marker_outline = 'red' if row['max_ros'] == 1 else 'black'
        marker_edgewidth = .5 if row['max_ros'] == 1 else 0.3
        
        # Plot marker with both fill and edge colors
        ax.plot(
            x, y,
            marker='o',
            markersize=6,
            markerfacecolor=marker_color,
            markeredgecolor=marker_outline,
            markeredgewidth=marker_edgewidth,
            transform=datacrs,
            alpha=0.9,
            zorder=201
        )

# --- Create legend handles ---
ar_handle = mlines.Line2D([], [], color='yellow', marker='o', linestyle='None',
                          markersize=6, markerfacecolor='yellow', markeredgecolor='black',
                          label='AR-related')
non_ar_handle = mlines.Line2D([], [], color='gray', marker='o', linestyle='None',
                              markersize=6, markerfacecolor='gray', markeredgecolor='black',
                              label='Not AR-related')

ros_handle = mlines.Line2D([], [], color='red', marker='o', linestyle='None',
                           markersize=6, markerfacecolor='none', markeredgecolor='red',
                           label='ROS event')

# --- Add legend to plot ---
ax.legend(handles=[ar_handle, non_ar_handle, ros_handle],
          loc='lower left', frameon=True, fontsize=8, title='Landslide Type')
    
# ---------------------------------------------------------------------
# Colorbar and Save Figure
# ---------------------------------------------------------------------
cb = fig.colorbar(cf, cax=cax, orientation='vertical')
cb.set_label('Mean ROS Frequency (d yr$^{-1}$)', fontsize=11)
cb.ax.tick_params(labelsize=10)

output_path = Path(f'../figs/cfsr_ros_landslide_{option}.png')
fig.savefig(output_path, bbox_inches='tight', dpi=fig.dpi)

print(f"Saved figure to: {output_path.resolve()}")
