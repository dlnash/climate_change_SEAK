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
    
    if extent is None:
        ax.set_global()
    else:
        ax.set_extent(extent, crs=datacrs)
    
    # Add map features (continents and country borders)
    ax.add_feature(cfeature.LAND, facecolor='0.9')      
    ax.add_feature(cfeature.BORDERS, edgecolor='0.5', linewidth=0.4, zorder=199)
    if coastline == True:
        ax.add_feature(cfeature.COASTLINE, edgecolor='0.4', linewidth=0.4)
    if mask_ocean == True:
        ocean = cfeature.NaturalEarthFeature('physical', 'ocean', \
        scale='50m', edgecolor='none', facecolor='#89C2D9')
        ax.add_feature(ocean)
        
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
    
    
    # Gridlines
    if grid:
        gl.xlines = True
        gl.ylines = True
    else:
        gl.xlines = False
        gl.ylines = False
    
    # Add tick marks (no labels)
    ax.set_xticks(xticks, crs=datacrs)
    ax.set_yticks(yticks, crs=datacrs)
    ax.tick_params(labelbottom=False, labelleft=False, length=3, width=0.3, color='k')

    return ax
