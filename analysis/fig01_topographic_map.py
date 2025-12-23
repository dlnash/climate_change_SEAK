#!/usr/bin/env python3
"""
Filename:    fig1_topographic_map.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Plot a topographic map of Southeast Alaska with labeled terrain features and communities.
TODO: Bold the text and add POW
"""

# --- Standard Library Imports ---
import sys, os
from pathlib import Path
import textwrap

# --- Third-Party Imports ---
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature

# --- Personal Modules ---
sys.path.append('../modules')
from plotter import draw_basemap, set_font
import globalvars
from wrf_utils import filter_vars

# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def make_land_cmap(vmin=0, vmax=2500, step=250):
    """
    Return a segmented terrain-style colormap and normalization for land elevations.
    Colors are discrete, with steps defined by 'step'.
    """
    # Define levels (edges of color bins)
    cflevs = np.arange(vmin, vmax + step, step)

    # Sample terrain colormap over a subset to focus on land colors
    colors_land = plt.cm.terrain(np.linspace(0.25, 1, len(cflevs) - 1))

    # Create a ListedColormap for discrete colors
    cmap = mcolors.ListedColormap(colors_land, name='terrain_land')

    # Create a BoundaryNorm to map data values into the discrete intervals
    norm = mcolors.BoundaryNorm(boundaries=cflevs, ncolors=cmap.N, clip=True)

    return cmap, norm

# ---------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------
data_path = Path(globalvars.path_to_data) / 'downloads' / 'SEAK-WRF' / 'geo_southeast.nc'
elev = xr.open_dataset(data_path)
elev = filter_vars(elev.squeeze(), data_path, "hgt")
elev = elev.hgt

# --- load shapefile and subset ---
fp = os.path.join(globalvars.path_to_data, 'downloads/AK_climate_divisions/AK_divisions_NAD83.shp')
polys = gpd.read_file(fp)
keep_names = ["Northeast Gulf", "North Panhandle", "Central Panhandle", "South Panhandle"]
polys = polys[polys["Name"].isin(keep_names)]
polys = polys.to_crs(epsg=4326)

# ---------------------------------------------------------------------
# Map and Label Definitions
# ---------------------------------------------------------------------

feature_labels = {
    'Chichagof Island': (-135.833242, 57.7),
    'Baranof Island': (-135.2, 56.971734),
    'Prince of Wales Island': (-132.684034, 55.3)
}

label_offsets = {
    "Northeast Gulf": {"dx": 110, "dy": -80},  # southwest of centroid
    "North Panhandle": {"dx": 0, "dy": 0},       # east of centroid
    "Central Panhandle": {"dx": 0, "dy": 0},     # east
    "South Panhandle": {"dx": 0, "dy": 0},       # east
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
cmap, norm = make_land_cmap()

# Draw base map
ax = draw_basemap(
    ax, extent=extent, xticks=xticks, yticks=yticks,
    left_lats=True, right_lats=False,
    mask_ocean=True, coastline=False
)

ax.set_extent(extent, crs=datacrs)  # extent in lon/lat degrees
ax.set_aspect('auto')  # or 'equal' if you want correct lat/lon ratio

# Plot topography
cf = ax.pcolormesh(
    elev.lon, elev.lat, elev.where(elev > 0),
    rasterized=False, cmap=cmap, norm=norm,
    transform=datacrs, alpha=0.6
)

# Plot Climate Divisions
polys.crs = 'epsg:4326'
ec_lst = ['#000000', '#DDAA33', '#AA3377', '#004488']
for idx, (i, poly) in enumerate(polys.iterrows()):
    feature = ShapelyFeature([poly.geometry], ccrs.PlateCarree(),
                             edgecolor=ec_lst[idx], facecolor='none', linewidth=1.)
    ax.add_feature(feature, zorder=200)

    geom = poly.geometry
    centroid = geom.centroid
    offsets = label_offsets.get(poly["Name"], {"dx": 0, "dy": 0})
    print(poly["Name"], offsets, centroid.x, centroid.y)
    ax.annotate(
        poly["Name"],
        xy=(centroid.x, centroid.y),
        xytext=(offsets["dx"], offsets["dy"]),  # shift in points
        textcoords="offset points",
        transform=datacrs,
        fontweight='bold',
        ha="center",
        va="center",
        color=ec_lst[idx],
        path_effects=[pe.withStroke(linewidth=2, foreground="white")],
        zorder=215
    )

# Plot geographic feature labels
style = {'color': 'black', 'fontweight': 'bold'}
for label, (x, y) in feature_labels.items():
    ax.annotate(
        textwrap.fill(label, 11),
        (x, y),
        textcoords="offset points",
        xytext=(0, 0),
        ha='center',
        fontsize=8,
        xycoords=datacrs._as_mpl_transform(ax),
        path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
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
