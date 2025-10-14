"""
Filename:    clim_trend_grid.py
Author:      Deanna Nash (adapted for multi-model comparison)
Description: Plot 3-column x 5-row of model and variable trends
"""

import os
import string
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
from plotter import set_font

# ============================================================
# Plot Function
# ============================================================
def plot_clim_trend_grid(models, varnames, path_to_data,
                        lonmin=-141., lonmax=-130.,
                        latmin=54.5, latmax=60., sig_level = 0.1):
    """
    Create a 5-row x 3-column plot:
        Rows = variables
        Col 1 = CFSR
        Col 2 = CCSM
        Col 3 = GFDL
    Each row has:
        - one shared colorbar
    """

    # --- Setup ---
    mapcrs = ccrs.Mercator()
    datacrs = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9.5, 14))
    fig.dpi = 300
    current_dpi=300
    base_dpi=100
    scaling_factor = (current_dpi / base_dpi)**0.3
    set_font(current_dpi, scaling_factor)
    
    nrows = len(varnames)
    ncols = len(models)

    gs = GridSpec(nrows, ncols+1, height_ratios=[1]*nrows, width_ratios=[1, 1, 1, 0.05],
                  hspace=0.05, wspace=0.05)

    # ============================================================
    # MAIN LOOP
    # ============================================================
    varname_lbl = ['IVT', 'UV', 'Freezing Level', 'Precipitation', 'SWE']
    labels = list(string.ascii_lowercase)  # ['a', 'b', 'c', ...]
    label_idx = 0
    bbox_dict = dict(facecolor='white', edgecolor='k', boxstyle='circle,pad=0.3', alpha=1.)
    
    for i, varname in enumerate(varnames):
        for j, model in enumerate(models):
            print(f"Plotting {varname}...")
            # Load climatology
            fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_clim.nc")
            ds_clim = xr.open_dataset(fname)
            
            # Load trend
            trend_fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_trends.nc")
            ds_trend = xr.open_dataset(trend_fname)
            lons = ds_clim.lon.values
            lats = ds_clim.lat.values
            
            # ----------------------------
            # Get base field and color levels
            # ----------------------------
            cflevs_tr = np.arange(-20, 22, 4)
            if (varname == "pcpt") | (varname == "snow"):
                # cflevs_tr = np.arange(-40, 50, 10); 
                cmap_tr = "BrBG"
            elif varname == "freezing_level":
                # cflevs_tr = np.arange(-10, 12, 2)
                cmap_tr = cmocean.tools.crop_by_percent(cmo.balance, 20, which="both")
            else:
                cmap_tr = 'BrBG'
    
            # ============================================================
            # --- Column: Trends ---
            # ============================================================
            ax0 = fig.add_subplot(gs[i, j], projection=mapcrs)
            ax0.set_extent([lonmin, lonmax, latmin, latmax])
            ax0.coastlines(resolution="50m")
            ax0.add_feature(cfeature.BORDERS, linewidth=0.75, edgecolor='k')
                
            ckey = f"{varname}_trend"
            pkey = f"{varname}_p"
            
            field = ds_trend[ckey]
            perc_change = (field * ds_trend.attrs['n_years'] / ds_clim[varname].values) * 100.0
            
            # --- Plot all values ---
            cf = ax0.contourf(lons, lats, perc_change, transform=datacrs,
                              levels=cflevs_tr, cmap=cmap_tr, extend="both")
            
            # --- Overlay transparency where NOT significant ---
            if pkey in ds_trend:
                sig_mask = ds_trend[pkey] > sig_level
                # Add semi-transparent gray overlay
                ax0.contourf(lons, lats, sig_mask, levels=[0.5, 1.5], colors='none',
                             hatches=['////'], alpha=0, transform=datacrs)
            # if pkey in ds_trend:
            #     nonsig_mask = ds_trend[pkey] > sig_level
            #     ax0.contourf(lons, lats, nonsig_mask, levels=[0.5, 1.5],
            #                  colors='gray', alpha=0.3, transform=datacrs)
            # --- Overlay stippling where NOT significant ---
            # if pkey in ds_trend:
            #     nonsig_mask = (ds_trend[pkey] > sig_level)
            #     # Plot stippling as black dots at those gridpoints
            #     y, x = np.meshgrid(lats, lons, indexing="ij")
            #     ax0.scatter(x[nonsig_mask], y[nonsig_mask],
            #                 s=2, color='k', alpha=0.3, transform=datacrs)


            
            if i == 0:
                ax0.set_title(model.upper())
    
            # --- add a, b, c labels --- 
            ax0.text(0.05, 0.96, f"{labels[label_idx]}", transform=ax0.transAxes,
                     va='top', ha='left', bbox=bbox_dict)
            label_idx += 1

        # --- colorbar (vertical to the right) ---
        cbax = fig.add_subplot(gs[i, -1])
        cb0 = Colorbar(ax=cbax, mappable=cf, orientation="vertical")
        cb0.set_label(fr'$\Delta$ {varname_lbl[i]} (%)')

    # ============================================================
    # Save Figure
    # ============================================================
    outname = f"../figs/MODEL_TREND_GRID.png"
    fig.savefig(outname, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved figure: {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    varnames = ["ivt", "uv", "freezing_level", "pcpt", "snow"]
    plot_clim_trend_grid(models, varnames, path_to_data)
