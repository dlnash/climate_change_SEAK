#!/usr/bin/env python3
"""
Filename:    scatterplots.py
Author:      Deanna Nash (dnash@ucsd.edu)

Description:
    Scatter/hexbin plots of projected changes:
        Column 1 = CCSM (Future - Historical)
        Column 2 = GFDL (Future - Historical)
        Column 3 = Multi-model mean Δ

    x-axis = Elevation
    y-axis = Latitude
    color = projected change
"""

# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import os
import sys
import string
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar

# Add path to custom modules
sys.path.append('../modules/')
import globalvars
from plotter import set_font
from wrf_utils import load_preprocessed_wrf_metrics, read_wrf_terrain
from colormaps import get_colormap_and_levels
from plot_configs import PLOT_CONFIGS, MODELS, SCENARIOS

import rioxarray  # for spatial clipping
from wrf_rio import wrf_prepare_for_rio
from wrf_crs import load_climate_division_shapefile

# ---------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------
def plot_clim_diff_scatter(
                           path_to_data,
                           plot_type,
                           ssn, 
                           option=None, 
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

    # --- load shapefile and subset ---
    polys = load_climate_division_shapefile()
    # Dissolve into one multipolygon
    study_area = polys.dissolve()

    # -------------------------------------------------------------
    # Figure setup
    # -------------------------------------------------------------
    fig = plt.figure(figsize=(9.75, 14))
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

    mask, x_flat = read_wrf_terrain()

    # =============================================================
    # Main Loop
    # =============================================================
    for i, varname in enumerate(varnames):

        fields, titles, units = load_preprocessed_wrf_metrics(
            varname, 
            plot_type, 
            season, 
            option
        )

        # ---------------------------------------------------------
        # Colormap
        # ---------------------------------------------------------
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

        diff_handles = []

        # ---------------------------------------------------------
        # Plot columns
        # ---------------------------------------------------------
        grid_cols = [0, 3, 4]

        for k, gcol in enumerate(grid_cols):

            # prepare for rioxarray
            da = wrf_prepare_for_rio(fields[k])
            
            # Reproject once after reading the raster
            study_area = study_area.to_crs(da.rio.crs)

            # Clip both datasets to polygon
            da_clip = da.rio.clip(
                study_area.geometry,
                study_area.crs,
                drop=True
            )

            lat2d = da_clip.lat.values
            y_flat = lat2d.flatten()[mask]
            cfield = da_clip.values.flatten()[mask]
        
            ax = fig.add_subplot(
                gs[i, gcol],
            )

            ax.set_xlim(0, 1900)
            ax.set_ylim(54.5, 60)

            # -----------------------------
            # Historical climatology
            # -----------------------------
            if k == 0:
                cf = ax.hexbin(
                    x_flat,
                    y_flat,
                    C=cfield,
                    gridsize=50,
                    cmap=cmap_clim,
                    norm=norm_clim,
                    reduce_C_function=np.mean,
                    mincnt=1
                )
                
            # -----------------------------
            # Difference panels
            # -----------------------------
            else:

                cf = ax.hexbin(
                    x_flat,
                    y_flat,
                    C=cfield,
                    gridsize=50,
                    cmap=cmap_diff,
                    norm=norm_diff,
                    reduce_C_function=np.mean,
                    mincnt=1
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

            if i == nrows - 1:
                ax.set_xlabel("Elevation (m)")
            else:
                ax.set_xticklabels([])
            if k == 0:
                ax.set_ylabel("Latitude (°N)")
            else:
                ax.set_yticklabels([])

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
        plot_name="hexbin",
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


# =============================================================
# Driver
# =============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    seasons = ["ONDJFM"]
    options = ["strict"]

    # --- 95th percentile climatology ---
    for season in seasons:

        print(
            f"Creating 95th percentile plot: {season}"
        )

        plot_clim_diff_scatter(
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

            plot_clim_diff_scatter(
                path_to_data,
                "ros_frequency_clim",
                season,
                option,
            )

            print(
                f"Creating ROS intensity plot: "
                f"{option} {season}"
            )

            plot_clim_diff_scatter(
                path_to_data,
                "ros_intensity_clim",
                season,
                option,
            )
