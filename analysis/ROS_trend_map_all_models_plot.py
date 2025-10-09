"""
Filename:    ROS_trend_map_plot.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Script to plot trends for each variable across all models,
             with shared colorbars per column.
"""

import sys, os
import numpy as np
import xarray as xr
import pandas as pd

# matplotlib
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import cmocean
# cartopy
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar
import cartopy.feature as cfeature

# Path to modules
sys.path.append('../modules/')
import globalvars
from plotter import make_brgr_white_cmap


# ============================================================
# Plot function
# ============================================================
def plot_trend_multi_models(models, varname, path_to_data,
                            lonmin=-141., lonmax=-130.,
                            latmin=54.5, latmax=60., sig_level=0.1):
    """
    Create a multi-row, 3-column plot:
        rows = models
        col 0 = climatology
        col 1 = trend
        col 2 = cfsr clim - future model clim
    Shared colorbars per column.
    """

    # --- Setup ---
    mapcrs = ccrs.PlateCarree()
    datacrs = ccrs.PlateCarree()
    dx = np.arange(lonmin, lonmax+3, 3)
    dy = np.arange(latmin, latmax+1, 1)

    nrows = len(models)
    ncols = 3
    fig = plt.figure(figsize=(12, 2*nrows))
    fig.dpi = 300
    gs = GridSpec(nrows+1, ncols, height_ratios=[1]*nrows + [0.05],
                  width_ratios=[1]*ncols, wspace=0.05, hspace=0.05)

    # Store colorbar handles
    cfs_clim = []
    cfs_diff = []
    cfs_trend = []

    for i, model in enumerate(models):
        print(f"Processing {varname} for {model}...")

        # Load climatology
        fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/snow_{model}_ros_clim.nc")
        ds_clim = xr.open_dataset(fname)
        # ds_clim = ds_clim.rename({'delsnowh': 'snowmelt'})
        if model == 'cfsr':
            cfsr_clim = ds_clim
            yr_range = "(1980-2019)"
        else:
            clim_diff = ds_clim - cfsr_clim
            yr_range = "(2030-2060)"

        ds_clim['rain'].attrs['units'] = 'mm day-1'
        ds_clim['delsnowh'].attrs['units'] = 'mm day-1'
        ds_clim['ros_intensity'].attrs['units'] = 'mm day-1'
        ds_clim['ros'].attrs['units'] = '#'

        # Load trend
        trend_fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/snow_{model}_ros_trends.nc")
        ds_trend = xr.open_dataset(trend_fname)
        # ds_trend = ds_trend.rename({'delsnowh': 'snowmelt'})

        lons = ds_clim.lon.values
        lats = ds_clim.lat.values

        # === Col 0: Climatology ===
        ax0 = fig.add_subplot(gs[i, 0], projection=mapcrs)
        ax0.set_extent([lonmin, lonmax, latmin, latmax])
        ax0.coastlines(resolution="50m")
        ax0.gridlines(draw_labels=False)
        ax0.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)

        # pick contour levels and cmap
        if varname == 'ros':
            cflevs = np.arange(0, 20, 1)
            cmap = cmo.rain
        else:
            cflevs = np.arange(0, 110, 10)
            cmap = cmo.rain
        # elif (varname == 'delsnow'):
        #     cflevs = np.arange(-25, 25, 5)
        #     white_range = (-5, 5)
        #     cmap, norm = make_brgr_white_cmap(cflevs, white_range) 
        # elif (varname == 'delsnowh'):
        #     cflevs = np.arange(-100, 120, 20)
        #     white_range = (-20, 20)
        #     cmap, norm = make_brgr_white_cmap(cflevs, white_range)
        # elif varname == 'rain':
        #     cflevs = np.arange(0, 110, 10)
        #     cmap = cmo.rain

        cfield = ds_clim[varname].values
        print(np.nanmin(cfield), np.nanmax(cfield))
        cf0 = ax0.contourf(lons, lats, cfield, transform=datacrs,
                           levels=cflevs, cmap=cmap, extend='both')
        cfs_clim.append(cf0)

        if i == 0:
            ax0.set_title(f"{varname.upper()} avg", loc="center")
        if i == nrows-1:
            ax0.set_xlabel("Longitude")
        if i == nrows//2:
            ax0.set_ylabel("Climatology")
        hlabel = f"{model.upper()} {yr_range}"
        ax0.text(-0.1, 0.5, hlabel, transform=ax0.transAxes,
                 rotation=90, ha="right", va="center", fontsize=10)

        # === Col 1: Trend ===
        ax1 = fig.add_subplot(gs[i, 1], projection=mapcrs)
        ax1.set_extent([lonmin, lonmax, latmin, latmax])
        ax1.coastlines(resolution="50m")
        ax1.gridlines(draw_labels=False)
        ax1.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)
        cflevs_tr = np.arange(-100, 110, 10)
        
        cmap_tr = 'BrBG'
        ckey = f"{varname}_trend"
        pkey = f"{varname}_p"
        field = ds_trend[ckey].where(ds_trend[pkey] <= sig_level).values
        
        perc_change = (field*ds_trend.attrs['n_years']/ds_clim[varname].values)*100.
        cf1 = ax1.contourf(lons, lats, perc_change, transform=datacrs,
                           levels=cflevs_tr, cmap=cmap_tr, extend="both")

        cfs_trend.append(cf1)

        if i == 0:
            ax1.set_title(f"{varname.upper()} Trend", loc="center")
        if i == nrows-1:
            ax1.set_xlabel("Longitude")

        # === Col 2: Climatology Difference ===
        if i > 0:
            ax2 = fig.add_subplot(gs[i, 2], projection=mapcrs)
            ax2.set_extent([lonmin, lonmax, latmin, latmax])
            ax2.coastlines(resolution="50m")
            ax2.gridlines(draw_labels=False)
            ax2.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)
            
            # pick contour levels and cmap
            if varname == 'ros':
                cflevs = np.arange(-5, 6, 1)
                white_range = (-1, 1)
            else:
                cflevs = np.arange(-20, 24, 4)
                white_range = (-4, 4)
            # elif (varname == 'delsnow'):
            #     cflevs = np.arange(-4, 5, 1)
            #     white_range = (-1, 1)
            # elif (varname == 'delsnowh'):
            #     cflevs = np.arange(-4, 5, 1)
            #     white_range = (-1, 1)
            # elif varname == 'rain':
            #     cflevs = np.arange(-20, 24, 4)
            #     white_range = (-4, 4)

            cmap, norm = make_brgr_white_cmap(cflevs, white_range)
    
            cfield = clim_diff[varname].values
            cf0 = ax2.contourf(lons, lats, cfield, transform=datacrs,
                               levels=cflevs, cmap=cmap, norm=norm, extend='both')
            cfs_diff.append(cf0)
    
            if i == 1:
                ax2.set_title(f"Future - CFSR", loc="center")
            if i == nrows-1:
                ax2.set_xlabel("Longitude")
            if i == nrows//2:
                ax2.set_ylabel("Climatology")

    # === Shared colorbars ===
    cbax0 = fig.add_subplot(gs[-1, 0])
    cb0 = Colorbar(ax=cbax0, mappable=cfs_clim[0], orientation="horizontal", ticklocation="bottom")
    cb0.set_label(fr'{varname.upper()} ({ds_clim[varname].attrs["units"]})', fontsize=10)

    cbax1 = fig.add_subplot(gs[-1, 1])
    cb1 = Colorbar(ax=cbax1, mappable=cfs_trend[0], orientation="horizontal", ticklocation="bottom")
    cb1.set_label(fr'$\Delta$ {varname.upper()} (%)', fontsize=10)

    cbax2 = fig.add_subplot(gs[-1, 2])
    cb2 = Colorbar(ax=cbax2, mappable=cfs_diff[0], orientation="horizontal", ticklocation="bottom")
    cb2.set_label(fr'$\Delta$ {varname.upper()} ({ds_clim[varname].attrs["units"]})', fontsize=10)

    # Save figure
    outname = f"../figs/ROS_MODEL_COMPARISON_{varname}_clim_trend.png"
    fig.savefig(outname, bbox_inches="tight", dpi=fig.dpi)
    plt.close(fig)
    print(f"Saved {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    varnames = ["ros", "ros_intensity", "delsnowh", "rain"]

    for varname in varnames:
        plot_trend_multi_models(models, varname, path_to_data,sig_level=1)
