#!/usr/bin/env python3
"""
Filename:    fig1_topographic_map.py
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
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- Personal Modules ---
sys.path.append('../modules')
from plotter import draw_basemap, set_font
import globalvars
from wrf_utils import filter_vars


# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def make_land_cmap(vmin=0, vmax=3000):
    """
    Return a terrain-style colormap and normalization for land elevations only.
    """
    colors_land = plt.cm.terrain(np.linspace(0.25, 1, 256))
    cmap = mcolors.LinearSegmentedColormap.from_list('terrain_land', colors_land)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    return cmap, norm


# ---------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------
data_path = Path(globalvars.path_to_data) / 'downloads' / 'SEAK-WRF' / 'geo_southeast.nc'
elev = xr.open_dataset(data_path)
elev = filter_vars(elev.squeeze(), data_path, "hgt")

# ---------------------------------------------------------------------
# Map and Label Definitions
# ---------------------------------------------------------------------
communities = {
    'Hoonah': (-135.4519, 58.1122, 'center'),
    'Skagway': (-135.3277, 59.4538, 'left'),
    'Klukwan': (-135.8894, 59.3988, 'right'),
    'Yakutat': (-139.6710, 59.5121, 'center'),
    'Craig': (-133.1358, 55.4769, 'right'),
    'Kasaan': (-132.4009, 55.5400, 'left'),
}

feature_labels = {
    'Lynn Canal': (-135.1132, 58.7004),
    'Dixon Entrance': (-131.5954, 54.7378),
    'Gulf of Alaska': (-138.0, 57.5),
    'Glacier Bay National Park and Preserve': (-137.1026, 58.43),
}

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
cmap, norm = make_land_cmap(vmin=0, vmax=2500)

# Draw base map
ax = draw_basemap(
    ax, extent=extent, xticks=xticks, yticks=yticks,
    left_lats=True, right_lats=False,
    mask_ocean=True, coastline=False
)

# Plot topography
cf = ax.pcolormesh(
    elev.lon, elev.lat, elev.hgt.where(elev.hgt > 0),
    rasterized=False, cmap=cmap, norm=norm,
    transform=datacrs, alpha=0.6
)

# Plot community markers and labels
for name, (x, y, ha) in communities.items():
    ax.plot(x, y, 'ro', markersize=5, transform=datacrs, zorder=201)
    ax.annotate(
        name, (x, y), xytext=(0, 12),
        textcoords="offset points",
        ha=ha,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.5, alpha=0.8),
        xycoords=datacrs._as_mpl_transform(ax),
        zorder=200,
    )

# Plot geographic feature labels
style = {'color': 'black', 'fontweight': 'normal'}
for label, (x, y) in feature_labels.items():
    ax.annotate(
        textwrap.fill(label, 11),
        (x, y),
        textcoords="offset points",
        xytext=(0, 0),
        ha='center',
        xycoords=datacrs._as_mpl_transform(ax),
        path_effects=[pe.withStroke(linewidth=1.25, foreground="white")],
        zorder=200,
        **style
    )

    
# ---------------------------------------------------------------------
# Colorbar and Save Figure
# ---------------------------------------------------------------------
cb = fig.colorbar(cf, cax=cax, orientation='vertical')
cb.set_label('Elevation (m)', fontsize=11)
cb.ax.tick_params(labelsize=10)

output_path = Path('../figs/topographic_map.png')
fig.savefig(output_path, bbox_inches='tight', dpi=fig.dpi)

print(f"Saved figure to: {output_path.resolve()}")
