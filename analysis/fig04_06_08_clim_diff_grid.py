"""
Filename:    plot_change_grids.py
Author:      Deanna Nash (refactored)
Description:
Generalized plotting script for:
    - 95th percentile climatology changes
    - ROS frequency changes
    - ROS intensity changes

Columns:
    1. CCSM future - historical
    2. GFDL future - historical
    3. Multi-model mean Δ
"""

# ============================================================
# Imports
# ============================================================
import os
import string
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar

# Path to modules
import sys
sys.path.append('../modules/')

import globalvars
from wrf_utils import load_preprocessed_wrf_metrics
from plotter import set_font
from colormaps import get_colormap_and_levels
from plot_configs import PLOT_CONFIGS, MODELS, SCENARIOS


# ============================================================
# Plot Function
# ============================================================
def plot_change_grid(
    path_to_data,
    plot_type,
    season,
    option=None,
    lonmin=-141.,
    lonmax=-130.,
    latmin=54.5,
    latmax=60.
):
    """
    Create a 5-row x 3-column plot:

    Columns:
        1 = Multi-model mean (historical)
        2 = Multi model mean (future - historical)
        3 = CCSM - GFDL (future)
    """

    cfg = PLOT_CONFIGS[plot_type]

    varnames = cfg["varnames"]
    var_labels = cfg["labels"]

    # ============================================================
    # Figure setup
    # ============================================================
    mapcrs = ccrs.Mercator()
    datacrs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 14))
    fig.dpi = 300

    current_dpi = 300
    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi)**0.3

    set_font(current_dpi, scaling_factor)

    nrows = len(varnames)
    ncols = 3

    gs = GridSpec(
        nrows, 
        ncols+3, 
        height_ratios=[1]*nrows, 
        width_ratios=[1, 0.05, 0.25, 1, 1, 0.05], 
        hspace=0.05, 
        wspace=0.05
    )

    # -------------------------------------------------------------
    # Labels
    # -------------------------------------------------------------

    labels = list(string.ascii_lowercase)
    label_idx = 0

    bbox_dict = dict(
        facecolor='white',
        edgecolor='k',
        boxstyle='circle,pad=0.3',
        alpha=1.
    )

    # ============================================================
    # Main Loop
    # ============================================================
    for i, varname in enumerate(varnames):

        fields, titles, units = load_preprocessed_wrf_metrics(
            varname, 
            plot_type, 
            season, 
            option
        )
        
        lons = fields[0].lon.values
        lats = fields[0].lat.values
    
        (
            levs_clim,
            cmap_clim,
            norm_clim,
            levs_diff,
            cmap_diff,
            norm_diff,
        ) = get_colormap_and_levels(
            cfg["cmap_key"],
            varname
        )
        print(f"Plotting {varname}")

        # --------------------------------------------------------
        # Plot columns
        # --------------------------------------------------------
        grid_cols = [0, 3, 4]

        for k, gcol in enumerate(grid_cols):
        
            ax = fig.add_subplot(
                gs[i, gcol],
                projection=mapcrs,
            )
        
            ax.set_extent([lonmin, lonmax, latmin, latmax])
        
            ax.coastlines(resolution="50m")
        
            ax.add_feature(
                cfeature.BORDERS,
                linewidth=0.75,
                edgecolor="k",
            )
        
            # -----------------------------
            # Historical climatology
            # -----------------------------
            if k == 0:
        
                cf = ax.contourf(
                    lons,
                    lats,
                    fields[k],
                    levels=levs_clim,
                    cmap=cmap_clim,
                    norm=norm_clim,
                    transform=datacrs,
                    extend=cfg["extend"],
                )
        
            # -----------------------------
            # Difference panels
            # -----------------------------
            else:
        
                cf = ax.contourf(
                    lons,
                    lats,
                    fields[k],
                    levels=levs_diff,
                    cmap=cmap_diff,
                    norm=norm_diff,
                    transform=datacrs,
                    extend="both",
                )
        
            if i == 0:
                ax.set_title(titles[k])
        
            ax.text(
                0.05,
                0.96,
                labels[label_idx],
                transform=ax.transAxes,
                va="top",
                ha="left",
                bbox=bbox_dict,
            )
        
            label_idx += 1
    
            # --------------------------------------------------------
            # Colorbar
            # --------------------------------------------------------
            
            if k == 0:

                cbax = fig.add_subplot(gs[i, 1])
            
                cb = Colorbar(
                    ax=cbax,
                    mappable=cf,
                    orientation="vertical",
                )
            
                if plot_type == "ros_frequency_clim":
                    if varname == "ros":
                        cb.set_label("ROS (days yr$^{-1}$)")
                    else:
                        cb.set_label(units)
                else:
                    cb.set_label(
                        f"{var_labels[i]} ({units})"
                    )

            elif k == 2:

                cbax = fig.add_subplot(gs[i, -1])
            
                cb = Colorbar(
                    ax=cbax,
                    mappable=cf,
                    orientation="vertical",
                )
            
                if plot_type == "ros_frequency_clim":
                    if varname == "ros":
                        cb.set_label(r"$\Delta$ ROS (days yr$^{-1}$)")
                    else:
                        cb.set_label(rf"$\Delta$ {units}")
                else:
                    cb.set_label(
                        rf"$\Delta$ {var_labels[i]} ({units})"
                    )

    # ============================================================
    # Save Figure
    # ============================================================
    outname = cfg["outfile"].format(
        season=season,
        option=option,
        plot_name="map",
    )

    os.makedirs(
        os.path.dirname(outname),
        exist_ok=True
    )

    fig.savefig(
        outname,
        bbox_inches="tight",
        dpi=300
    )

    plt.close(fig)

    print(f"Saved figure: {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    seasons = ["ONDJFM"]
    options = ["strict"]

    # --- 95th percentile climatology ---
    for season in seasons:

        print(
            f"Creating 95th percentile plot: {season}"
        )

        plot_change_grid(
            path_to_data,
            "95th_percentile_clim",
            season
        )

    # --- ROS plots ---
    for option in options:
        for season in seasons:

            print(
                f"Creating ROS frequency plot: "
                f"{option} {season}"
            )

            plot_change_grid(
                path_to_data,
                "ros_frequency_clim",
                season,
                option,
            )

            print(
                f"Creating ROS intensity plot: "
                f"{option} {season}"
            )

            plot_change_grid(
                path_to_data,
                "ros_intensity_clim",
                season,
                option,
            )