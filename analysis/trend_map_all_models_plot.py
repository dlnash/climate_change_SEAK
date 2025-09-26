"""
Filename:    trend_map_plot.py
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


# ============================================================
# Plot function
# ============================================================
def plot_trend_multi_models(models, varname, path_to_data,
                            lonmin=-141., lonmax=-130.,
                            latmin=54.5, latmax=60., sig_level=0.1):
    """
    Create a multi-row, 2-column plot:
        rows = models
        col 0 = climatology
        col 1 = trend
    Shared colorbars per column.
    """

    # --- Setup ---
    mapcrs = ccrs.PlateCarree()
    datacrs = ccrs.PlateCarree()
    dx = np.arange(lonmin, lonmax+3, 3)
    dy = np.arange(latmin, latmax+1, 1)

    nrows = len(models)
    ncols = 2
    fig = plt.figure(figsize=(9, 3*nrows))
    fig.dpi = 300
    gs = GridSpec(nrows+1, ncols, height_ratios=[1]*nrows + [0.05],
                  width_ratios=[1, 1], wspace=0.05, hspace=0.05)

    # Store colorbar handles
    cfs_clim = []
    cfs_trend = []

    for i, model in enumerate(models):
        print(f"Processing {varname} for {model}...")

        # Load climatology
        fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_clim.nc")
        ds_clim = xr.open_dataset(fname)

        # Load trend
        trend_fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_trends.nc")
        ds_trend = xr.open_dataset(trend_fname)

        lons = ds_clim.lon.values
        lats = ds_clim.lat.values

        # === Left col: Climatology ===
        ax0 = fig.add_subplot(gs[i, 0], projection=mapcrs)
        ax0.set_extent([lonmin, lonmax, latmin, latmax])
        ax0.coastlines(resolution="50m")
        ax0.gridlines(draw_labels=False)
        ax0.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)

        # pick contour levels and cmap
        if varname == 'ivt':
            cflevs = np.arange(150, 750, 100); cmap = cmo.deep
        elif (varname == 'pcpt'):
            cflevs = np.arange(0, 110, 10); cmap = cmo.rain
        elif (varname == 'uv'):
            cflevs = np.arange(0, 55, 5); cmap = cmo.dense
        elif varname == 'freezing_level':
            cflevs = np.arange(2500, 4600, 100)
            cmap = cmocean.tools.crop_by_percent(cmo.ice, 20, which='min')
        elif varname == 'snow':
            cflevs = np.arange(0, 1100, 100); cmap = cmo.rain

        cfield = ds_clim[varname].values
        cf0 = ax0.contourf(lons, lats, cfield, transform=datacrs,
                           levels=cflevs, cmap=cmap, extend='max')
        cfs_clim.append(cf0)

        if i == 0:
            ax0.set_title(f"{varname.upper()} avg 95th percentile", loc="center")
        if i == nrows-1:
            ax0.set_xlabel("Longitude")
        if i == nrows//2:
            ax0.set_ylabel("Climatology")
        ax0.text(-0.1, 0.5, model.upper(), transform=ax0.transAxes,
                 rotation=90, ha="right", va="center", fontsize=10)

        # === Right col: Trend ===
        ax1 = fig.add_subplot(gs[i, 1], projection=mapcrs)
        ax1.set_extent([lonmin, lonmax, latmin, latmax])
        ax1.coastlines(resolution="50m")
        ax1.gridlines(draw_labels=False)
        ax1.add_feature(cfeature.BORDERS, edgecolor='0.4', linewidth=0.8)
        cflevs_tr = np.arange(-20, 22, 4)
        if varname in ["ivt", "uv"]:
            if varname == "ivt":
                ukey, vkey, pkey, ckey = "ivtu_trend", "ivtv_trend", "ivt_p", "ivt_trend"
                qscale=10
                # cflevs_tr = np.arange(-10, 12, 2)
            else:
                ukey, vkey, pkey, ckey = "u_trend", "v_trend", "uv_p", "u_trend"
                qscale = 0.25
                # cflevs_tr = np.arange(-10, 12, 2)

            uvec = ds_trend[ukey].where(ds_trend[pkey] <= sig_level).values
            vvec = ds_trend[vkey].where(ds_trend[pkey] <= sig_level).values
            perc_change = (ds_trend[ckey].values *
                           ds_trend.attrs['n_years'] /
                           ds_clim[varname].values) * 100.

            cf1 = ax1.contourf(lons, lats, perc_change, transform=datacrs,
                               levels=cflevs_tr, cmap="BrBG", extend="both")
            ax1.quiver(lons, lats, uvec, vvec, transform=datacrs,
                       color="k", regrid_shape=13, pivot="middle",
                       angles="xy", scale_units="xy", scale=qscale, units="xy")

        else:
            ckey = f"{varname}_trend"
            pkey = f"{varname}_p"
            field = ds_trend[ckey].where(ds_trend[pkey] <= sig_level).values

            if (varname == "pcpt") | (varname == "snow"):
                # cflevs_tr = np.arange(-40, 50, 10); 
                cmap_tr = "BrBG"
            elif varname == "freezing_level":
                # cflevs_tr = np.arange(-10, 12, 2)
                cmap_tr = cmocean.tools.crop_by_percent(cmo.balance, 20, which="both")

            perc_change = (field*ds_trend.attrs['n_years']/ds_clim[varname].values)*100.
            cf1 = ax1.contourf(lons, lats, perc_change, transform=datacrs,
                               levels=cflevs_tr, cmap=cmap_tr, extend="both")

        cfs_trend.append(cf1)

        if i == 0:
            ax1.set_title(f"{varname.upper()} Trend", loc="center")
        if i == nrows-1:
            ax1.set_xlabel("Longitude")

    # === Shared colorbars ===
    cbax0 = fig.add_subplot(gs[-1, 0])
    cb0 = Colorbar(ax=cbax0, mappable=cfs_clim[0], orientation="horizontal", ticklocation="bottom")
    cb0.set_label(fr'{varname.upper()} ({ds_clim[varname].attrs["units"]})', fontsize=10)

    cbax1 = fig.add_subplot(gs[-1, 1])
    cb1 = Colorbar(ax=cbax1, mappable=cfs_trend[0], orientation="horizontal", ticklocation="bottom")
    cb1.set_label(fr'$\Delta$ {varname.upper()} (%)', fontsize=10)

    # Save figure
    outname = f"../figs/MODEL_COMPARISON_{varname}_clim_trend.png"
    fig.savefig(outname, bbox_inches="tight", dpi=fig.dpi)
    plt.close(fig)
    print(f"Saved {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    varnames = ["ivt", "pcpt", "freezing_level", "uv", "snow"]

    for varname in varnames:
        plot_trend_multi_models(models, varname, path_to_data,sig_level=1)
