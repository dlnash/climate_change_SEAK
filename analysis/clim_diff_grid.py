"""
Filename:    clim_diff_grid.py
Author:      Deanna Nash (adapted for multi-model comparison)
Description: Plot 3-column x 5-row comparison of CFSR climatology and model differences
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar

# Path to modules
import sys
sys.path.append('../modules/')
import globalvars

# ============================================================
# Plot Function
# ============================================================
def plot_clim_diff_grid(models, varnames, path_to_data,
                        lonmin=-141., lonmax=-130.,
                        latmin=54.5, latmax=60.):
    """
    Create a 5-row x 3-column plot:
        Rows = variables
        Col 1 = CFSR climatology
        Col 2 = CCSM - CFSR
        Col 3 = GFDL - CFSR
    Each row has:
        - one colorbar for CFSR (col 1)
        - one shared colorbar for CCSM & GFDL (cols 2–3)
    """

    # --- Setup ---
    mapcrs = ccrs.Mercator()
    datacrs = ccrs.PlateCarree()
    fig = plt.figure(figsize=(13, 14))
    fig.dpi = 300
    nrows = len(varnames)
    ncols = len(models)

    gs = GridSpec(nrows, ncols+2, height_ratios=[1]*nrows, width_ratios=[1, 0.04, 1, 1, 0.04],
                  hspace=0.25, wspace=0.08)

    # --- Load reference CFSR climatology once ---
    cfsr_ds = {}
    for var in varnames:
        fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/{var}_cfsr_clim.nc")
        cfsr_ds[var] = xr.open_dataset(fname)

    # ============================================================
    # MAIN LOOP
    # ============================================================
    for i, varname in enumerate(varnames):
        print(f"Plotting {varname}...")
        lons = cfsr_ds[varname].lon.values
        lats = cfsr_ds[varname].lat.values

        # ----------------------------
        # Get base field and color levels
        # ----------------------------
        if varname == 'ivt':
            levs_clim = np.arange(150, 750, 100); cmap_clim = cmo.deep
            levs_diff = np.arange(-250, 300, 50)
        elif varname == 'pcpt':
            levs_clim = np.arange(0, 110, 10); cmap_clim = cmo.rain
            levs_diff = np.arange(-20, 24, 4)
        elif varname == 'uv':
            levs_clim = np.arange(0, 55, 5); cmap_clim = cmo.dense
            levs_diff = np.arange(-5, 6, 1)
        elif varname == 'freezing_level':
            levs_clim = np.arange(2500, 4600, 100)
            cmap_clim = cmocean.tools.crop_by_percent(cmo.ice, 20, which='min')
            levs_diff = np.arange(-2000, 2250, 500)
        elif varname == 'snow':
            levs_clim = np.arange(0, 1100, 100); cmap_clim = cmo.rain
            levs_diff = np.arange(-500, 525, 50)
        else:
            levs_clim = np.linspace(0, 1, 10)
            levs_diff = np.linspace(-1, 1, 10)
            cmap_clim = cmo.tempo

        # ============================================================
        # --- Column 1: CFSR Climatology ---
        # ============================================================
        ax0 = fig.add_subplot(gs[i, 0], projection=mapcrs)
        ax0.set_extent([lonmin, lonmax, latmin, latmax])
        ax0.coastlines(resolution="50m")
        ax0.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor='0.4')
        ax0.gridlines(draw_labels=False)
        cfield = cfsr_ds[varname][varname].values
        cf0 = ax0.contourf(lons, lats, cfield, levels=levs_clim, cmap=cmap_clim,
                           transform=datacrs, extend='max')
        
        if i == 0:
            ax0.set_title(f"CFSR", fontsize=11)
        if i == nrows - 1:
            ax0.set_xlabel("Longitude")
        if i == nrows // 2:
            ax0.set_ylabel("Variable")

        # Add variable name on left, centered horizontally
        ax0.text(-0.3, 0.5, varname.upper(), transform=ax0.transAxes,
                 rotation=90, ha='right', va='center', fontsize=11, fontweight='bold')

        # --- CFSR colorbar (vertical to the right) ---
        cbax_cfsr = fig.add_subplot(gs[i, 1])
        cb0 = Colorbar(ax=cbax_cfsr, mappable=cf0, orientation="vertical")
        cb0.set_label(f"{varname.upper()} ({cfsr_ds[varname][varname].attrs.get('units', '')})", fontsize=9)

        # ============================================================
        # --- Columns 2 & 3: CCSM & GFDL Differences ---
        # ============================================================
        diff_handles = []
        for j, model in enumerate(models[1:], start=2):
            ax = fig.add_subplot(gs[i, j], projection=mapcrs)
            ax.set_extent([lonmin, lonmax, latmin, latmax])
            ax.coastlines(resolution="50m")
            ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor='0.4')
            ax.gridlines(draw_labels=False)

            fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_clim.nc")
            ds_model = xr.open_dataset(fname)
            diff = ds_model[varname] - cfsr_ds[varname][varname]
            cf = ax.contourf(lons, lats, diff, levels=levs_diff, cmap=cmo.balance,
                             transform=datacrs, extend='both')
            diff_handles.append(cf)
            if i == 0:
                ax.set_title(f"{model.upper()} - CFSR", fontsize=11)
            if i == nrows - 1:
                ax.set_xlabel("Longitude")

        # --- Shared colorbar for both diffs (cols 2 & 3) ---
        cbax_diff = fig.add_subplot(gs[i, -1])
        cb1 = Colorbar(ax=cbax_diff, mappable=diff_handles[0], orientation="vertical")
        cb1.set_label(f"Δ{varname.upper()} ({cfsr_ds[varname][varname].attrs.get('units', '')})", fontsize=9)

    # ============================================================
    # Save Figure
    # ============================================================
    outname = f"../figs/MODEL_DIFF_GRID.png"
    fig.savefig(outname, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved figure: {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    varnames = ["ivt", "pcpt", "freezing_level", "uv", "snow"]
    plot_clim_diff_grid(models, varnames, path_to_data)
