"""
Filename:    plotter.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for plotting
"""

# Import Python modules

import os, sys
import numpy as np
import itertools
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import colorsys
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm # Linear interpolation for color maps
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.projections import get_projection_class
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar # different way to handle colorbar
import seaborn as sns
import cmocean
import cmocean.cm as cmo
from datetime import timedelta
import textwrap
from matplotlib.gridspec import GridSpec
import itertools ## need this for the cbarticks
from matplotlib import font_manager as fm

## import personal modules
import cw3ecmaps as ccmaps
import globalvars

def set_font(current_dpi, scaling_factor):
    fm.fontManager.addfont(globalvars.path_to_repo+'modules/helvetica.ttc')

    plt.rcParams.update({
                    'font.family' : 'Helvetica',
                    'figure.dpi': current_dpi,
                    'font.size': 8 * scaling_factor, #changes axes tick label
                    'axes.labelsize': 8 * scaling_factor,
                    'axes.titlesize': 8 * scaling_factor,
                    'xtick.labelsize': 8 * scaling_factor,#do nothing
                    'ytick.labelsize': 8 * scaling_factor, #do nothing
                    'legend.fontsize': 5 * scaling_factor,
                    'lines.linewidth': 0.7 * scaling_factor,
                    'axes.linewidth': 0.2 * scaling_factor,
                    'legend.fontsize': 12 * scaling_factor,
                    'xtick.major.width': 0.8 * scaling_factor,
                    'ytick.major.width': 0.8 * scaling_factor,
                    'xtick.minor.width': 0.6 * scaling_factor,
                    'ytick.minor.width': 0.6 * scaling_factor,
                    'lines.markersize': 6 * scaling_factor
                })
    
def make_brgr_white_cmap(cflevs, white_range):
    """
    Create a 'BrGr'-style diverging colormap with white at the center.

    Parameters
    ----------
    cflevs : array-like
        Contour levels (e.g., np.arange(-10, 11, 2)).
    white_range : tuple (low, high)
        Range of values to make white (e.g., (-2, 2)).

    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        Custom colormap with white center.
    norm : matplotlib.colors.BoundaryNorm
        Normalization for use in contourf/pcolormesh.
    """

    # Use Brewer "BrBG" diverging colormap
    base_cmap = plt.get_cmap('BrBG', len(cflevs) - 1)
    colors = base_cmap(np.arange(len(cflevs) - 1))

    # Identify color bins that fall inside the white range
    mask = (cflevs[:-1] >= white_range[0]) & (cflevs[1:] <= white_range[1])

    # Set those bins to white
    colors[mask] = [1, 1, 1, 1]

    # Build new colormap and norm
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(cflevs, ncolors=cmap.N, clip=True)

    return cmap, norm


def plot_terrain(ax, ext):
    fname = '/expanse/nfs/cw3e/cwp140/downloads/ETOPO1_Bed_c_gmt4.grd'
    datacrs = ccrs.PlateCarree()
    grid = xr.open_dataset(fname)
    grid = grid.where(grid.z > 0) # mask below sea level
    grid = grid.sel(x=slice(ext[0], ext[1]), y=slice(ext[2], ext[3]))
    cs = ax.pcolormesh(grid.x, grid.y, grid.z,
                        cmap=cmo.gray_r, transform=datacrs, alpha=0.7)
    
    return ax
    
def draw_basemap(ax, datacrs=ccrs.PlateCarree(), extent=None, xticks=None, yticks=None, grid=False, left_lats=True, right_lats=False, bottom_lons=True, mask_ocean=False, coastline=True):
    """
    Creates and returns a background map on which to plot data. 
    
    Map features include continents and country borders.
    Option to set lat/lon tickmarks and draw gridlines.
    
    Parameters
    ----------
    ax : 
        plot Axes on which to draw the basemap
    
    datacrs : 
        crs that the data comes in (usually ccrs.PlateCarree())
        
    extent : float
        Set map extent to [lonmin, lonmax, latmin, latmax] 
        Default: None (uses global extent)
        
    grid : bool
        Whether to draw grid lines. Default: False
        
    xticks : float
        array of xtick locations (longitude tick marks)
    
    yticks : float
        array of ytick locations (latitude tick marks)
        
    left_lats : bool
        Whether to add latitude labels on the left side. Default: True
        
    right_lats : bool
        Whether to add latitude labels on the right side. Default: False
        
    Returns
    -------
    ax :
        plot Axes with Basemap
    
    Notes
    -----
    - Grayscale colors can be set using 0 (black) to 1 (white)
    - Alpha sets transparency (0 is transparent, 1 is solid)
    
    """
    ## some style dictionaries
    kw_ticklabels = {'size': 10, 'color': 'dimgray', 'weight': 'light'}
    kw_grid = {'linewidth': .5, 'color': 'k', 'linestyle': '--', 'alpha': 0.4}
    kw_ticks = {'length': 4, 'width': 0.5, 'pad': 2, 'color': 'black',
                         'labelsize': 10, 'labelcolor': 'dimgray'}

    # Use map projection (CRS) of the given Axes
    mapcrs = ax.projection    
    
    # Add map features (continents and country borders)
    ax.add_feature(cfeature.LAND, facecolor='0.9')      
    ax.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)
    if coastline == True:
        ax.add_feature(cfeature.COASTLINE, edgecolor='0.4', linewidth=0.8)
    if mask_ocean == True:
        ax.add_feature(cfeature.OCEAN, edgecolor='0.4', zorder=12, facecolor='white') # mask ocean
        
    ## Tickmarks/Labels
    ## Add in meridian and parallels
    if mapcrs == ccrs.NorthPolarStereo():
        gl = ax.gridlines(draw_labels=False,
                      linewidth=.5, color='black', alpha=0.5, linestyle='--')
    elif mapcrs == ccrs.SouthPolarStereo():
        gl = ax.gridlines(draw_labels=False,
                      linewidth=.5, color='black', alpha=0.5, linestyle='--')
        
    else:
        gl = ax.gridlines(crs=datacrs, draw_labels=True, **kw_grid)
        gl.top_labels = False
        gl.left_labels = left_lats
        gl.right_labels = right_lats
        gl.bottom_labels = bottom_lons
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = kw_ticklabels
        gl.ylabel_style = kw_ticklabels
    
    ## Gridlines
    # Draw gridlines if requested
    if (grid == True):
        gl.xlines = True
        gl.ylines = True
    if (grid == False):
        gl.xlines = False
        gl.ylines = False
            

    # apply tick parameters
    ax.set_xticks(xticks, crs=datacrs)
    ax.set_yticks(yticks, crs=datacrs)
    plt.yticks(color='w', size=1) # hack: make the ytick labels white so the ticks show up but not the labels
    plt.xticks(color='w', size=1) # hack: make the ytick labels white so the ticks show up but not the labels
    ax.ticklabel_format(axis='both', style='plain')

    ## Map Extent
    # If no extent is given, use global extent
    if extent is None:        
        ax.set_global()
        extent = [-180., 180., -90., 90.]
    # If extent is given, set map extent to lat/lon bounding box
    else:
        ax.set_extent(extent, crs=datacrs)
    
    return ax

def plot_trend_with_clim(ds, ds_clim, varname, lons, lats, model,
                         lonmin=-141., lonmax=-130., latmin=54.5, latmax=60.,
                         sig_level=0.1):
    """
    Generic plotting function for <varname>_trend and climatology.
    Handles scalar fields (pcpt, freezing level, etc.) and vector fields (IVT, winds).
    Produces a 2-column figure: climatology (left) and trend (right).
    """
    # --- Setup ---
    mapcrs = ccrs.PlateCarree()
    datacrs = ccrs.PlateCarree()
    dx = np.arange(lonmin, lonmax+3, 3)
    dy = np.arange(latmin, latmax+1, 1)

    fig = plt.figure(figsize=(9, 3.))  # wider for two columns
    fig.dpi = 300
    fname = f'../figs/{model}_{varname}_clim_trend'
    fmt = 'png'

    nrows, ncols = 2, 2
    gs = GridSpec(nrows, ncols,
                  height_ratios=[1, 0.05],  # map row + colorbar row
                  width_ratios=[1, 1],      # climatology + trend
                  wspace=0.05, hspace=0.05)

    # === Left column: Climatology ===
    if varname == 'ivt':
        cflevs = np.arange(150, 450, 25)
        cmap = cmo.deep
    elif varname == 'pcpt':
        cflevs = np.arange(0, 110, 10)
        cmap = cmo.rain
    elif varname == 'uv':
        cflevs = np.arange(0, 25, 2)
        cmap = cmo.dense
    elif varname == 'freezing_level':
        cflevs = np.arange(2500, 3700, 100)
        cmap = cmocean.tools.crop_by_percent(cmo.ice, 20, which='min', N=None)
    ax0 = fig.add_subplot(gs[0, 0], projection=mapcrs)
    ax0 = draw_basemap(ax0, extent=[lonmin, lonmax, latmin, latmax],
                       xticks=dx, yticks=dy, left_lats=True,
                       right_lats=False, bottom_lons=True)

    cfield = ds_clim[varname].values
    # cflevs = np.linspace(np.nanmin(cfield), np.nanmax(cfield), 21)
    
    cf0 = ax0.contourf(lons, lats, cfield, transform=datacrs,
                       levels=cflevs, cmap=cmap, extend='max')
    ax0.set_title(f"{varname.upper()} avg 95th percentile", loc='left')

    # --- Colorbar for climatology ---
    cbax0 = plt.subplot(gs[1, 0])
    cb0 = Colorbar(ax=cbax0, mappable=cf0, orientation='horizontal', ticklocation='bottom')
    cb0.set_label(fr'{varname.upper()} ({ds_clim[varname].attrs['units']})', fontsize=10)
    cb0.ax.tick_params(labelsize=10)

    # === Right column: Trend (with optional vectors) ===
    ax1 = fig.add_subplot(gs[0, 1], projection=mapcrs)
    ax1 = draw_basemap(ax1, extent=[lonmin, lonmax, latmin, latmax],
                       xticks=dx, yticks=dy, left_lats=False,
                       right_lats=False, bottom_lons=True)

    if varname in ["ivt", "uv"]:
        # IVT or winds (need u and v components)
        if varname == "ivt":
            ukey, vkey, pkey = "ivtu_trend", "ivtv_trend", "ivt_p"
            ckey = "ivt_trend"
            cflevs = np.arange(-10, 12, 2)
        elif varname == "uv":
            ukey, vkey, pkey = "u_trend", "v_trend", "uv_p"
            ckey = "u_trend"  # or another scalar field
            cflevs = np.arange(-10, 12, 2)
            
        uvec = ds[ukey].where(ds[pkey] <= sig_level).values
        vvec = ds[vkey].where(ds[pkey] <= sig_level).values

        # --- compute percent change based on number of years and the clim
        perc_change = (ds[ckey].values*ds.attrs['n_years']/ds_clim[varname].values)*100.
        print(f'Minimum % change: {np.nanmin(perc_change)}, Maximum % change: {np.nanmax(perc_change)}')
        # cflevs = np.linspace(np.nanmin(ds[ckey].values), np.nanmax(ds[ckey].values), 0.05)
        cf1 = ax1.contourf(lons, lats, perc_change, transform=datacrs,
                           levels=cflevs, cmap='BrBG', extend='both')

        Q = ax1.quiver(lons, lats, uvec, vvec, transform=datacrs,
                       color='k', regrid_shape=13, pivot='middle',
                       angles='xy', scale_units='xy', scale=1, units='xy')
        per_yr = 'yr$^{-1}$'
        ax1.quiverkey(Q, 0.65, 1.05, 0.25, f'0.25 {ds_clim[varname].attrs['units']} {per_yr}', labelpos='E',
                      coordinates='axes', fontproperties={'size': 8.0})

    else:
        ckey = f"{varname}_trend"
        pkey = f"{varname}_p"
        field = ds[ckey].where(ds[pkey] <= sig_level).values

        if varname == "pcpt":
            cflevs = np.arange(-40, 50, 10)
            cmap = 'BrBG'
        elif varname == 'freezing_level':
            cflevs = np.arange(-10, 12, 2)
            cmap = cmocean.tools.crop_by_percent(cmo.balance, 20, which='both', N=None)

            
        # --- compute percent change based on number of years and the clim
        perc_change = (field*ds.attrs['n_years']/ds_clim[varname].values)*100.
        print(f'Minimum % change: {np.nanmin(perc_change)}, Maximum % change: {np.nanmax(perc_change)}')
        
        cf1 = ax1.contourf(lons, lats, perc_change, transform=datacrs,
                           levels=cflevs, cmap=cmap, extend='both')

    ax1.set_title(f"{varname.upper()} Trend", loc='left')

    # --- Colorbar for trend ---
    cbax1 = plt.subplot(gs[1, 1])
    cb1 = Colorbar(ax=cbax1, mappable=cf1, orientation='horizontal', ticklocation='bottom')
    cb1.set_label(fr'$\Delta$ {varname.upper()} (%)', fontsize=10)
    cb1.ax.tick_params(labelsize=10)

    # Save
    fig.savefig(f"{fname}.{fmt}", bbox_inches='tight', dpi=fig.dpi)
    plt.show()
