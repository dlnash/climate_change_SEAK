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
from plotter import set_font
from colormaps import get_colormap_and_levels


# ============================================================
# CONFIG
# ============================================================
PLOT_CONFIGS = {

    "95th_percentile_clim": {
        "varnames": [
            "ivt",
            "uv",
            "freezing_level",
            "pcpt",
            "snow"
        ],

        "labels": [
            "IVT",
            r"$\overline{\mathrm{UV}}$",
            "Freezing Level",
            "Precipitation",
            "SWE"
        ],

        "filename_template":
            "{var}_{model}_{season}_95th_percentile_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/clim/{season}_95th_percentile_clim_change.png",

        "cmap_key":
            "95th_percentile_clim",

        "extend":
            "both",
    },

    "ros_frequency_clim": {
        "varnames": [
            "ros",
            "ivt",
            "pcpt",
            "delsnow",
            "delsnowh"
        ],

        "labels": [
            "ROS",
            "IVT",
            "Precip",
            r"$\Delta$ SWE",
            "Snowmelt"
        ],

        "filename_template":
            "snow_{model}_{season}_{option}_ros_frequency_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/ros_{option}/{season}/{season}_ROS_FREQ_CHANGE.png",

        "cmap_key":
            "ros_frequency_clim",

        "extend":
            "max",
    },

    "ros_intensity_clim": {
        "varnames": [
            "ros",
            "pcpt",
            "snow",
            "delsnowh",
            "ros_intensity"
        ],

        "labels": [
            "ROS",
            "Precip",
            "SWE",
            "Snowmelt",
            "Intensity"
        ],

        "filename_template":
            "snow_{model}_{season}_{option}_ros_intensity_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/ros_{option}/{season}/{season}_ROS_INTENSITY_CHANGE.png",

        "cmap_key":
            "ros_intensity_clim",

        "extend":
            "both",
    },
}


MODELS = ["ccsm", "gfdl"]
SCENARIOS = {
    "historical": "hist",
    "future": "rcp85"
}


# ============================================================
# Helper Functions
# ============================================================
def load_dataset(
    path_to_data,
    model,
    scenario,
    plot_type,
    varname,
    season,
    option=None
):
    """Load dataset for a given plot type."""

    cfg = PLOT_CONFIGS[plot_type]

    datadir = os.path.join(
        path_to_data,
        "preprocessed",
        "SEAK-WRF",
        cfg["subdir"].format(
            model=model,
            scenario=scenario
        )
    )

    fname = cfg["filename_template"].format(
        var=varname,
        model=model,
        season=season,
        option=option
    )

    path = os.path.join(datadir, fname)

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    return xr.open_dataset(path)


# ============================================================
# Plot Function
# ============================================================
def plot_change_grid(
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
        1 = CCSM future - historical
        2 = GFDL future - historical
        3 = Multi-model mean Δ
    """

    cfg = PLOT_CONFIGS[plot_type]

    varnames = cfg["varnames"]
    var_labels = cfg["labels"]

    path_to_data = globalvars.path_to_data

    # ============================================================
    # Setup
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

    gs = GridSpec(
        nrows,
        4,
        width_ratios=[1, 1, 1, 0.05],
        hspace=0.05,
        wspace=0.05
    )

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

        print(f"Plotting {varname}")

        # --------------------------------------------------------
        # Load datasets
        # --------------------------------------------------------
        hist = {}
        future = {}

        for model in MODELS:

            hist[model] = load_dataset(
                path_to_data,
                model,
                "hist",
                plot_type,
                varname,
                season,
                option
            )

            future[model] = load_dataset(
                path_to_data,
                model,
                "rcp85",
                plot_type,
                varname,
                season,
                option
            )

        # --------------------------------------------------------
        # Compute differences
        # --------------------------------------------------------
        ccsm_diff = (
            future["ccsm"][varname]
            - hist["ccsm"][varname]
        )

        gfdl_diff = (
            future["gfdl"][varname]
            - hist["gfdl"][varname]
        )

        mmm_diff = xr.concat(
            [ccsm_diff, gfdl_diff],
            dim="model"
        ).mean("model")

        fields = [
            ccsm_diff,
            gfdl_diff,
            mmm_diff
        ]

        titles = [
            "ΔCCSM (Future − Historical)",
            "ΔGFDL (Future − Historical)",
            "Multi-model Mean Δ"
        ]

        lons = ccsm_diff.lon.values
        lats = ccsm_diff.lat.values

        (
            _,
            _,
            _,
            levs_diff,
            cmap_diff,
            norm_diff
        ) = get_colormap_and_levels(
            cfg["cmap_key"],
            varname
        )

        # --------------------------------------------------------
        # Plot columns
        # --------------------------------------------------------
        diff_handles = []

        for j in range(3):

            ax = fig.add_subplot(
                gs[i, j],
                projection=mapcrs
            )

            ax.set_extent([
                lonmin,
                lonmax,
                latmin,
                latmax
            ])

            ax.coastlines(resolution="50m")

            ax.add_feature(
                cfeature.BORDERS,
                linewidth=0.75,
                edgecolor='k'
            )

            cf = ax.contourf(
                lons,
                lats,
                fields[j],
                levels=levs_diff,
                cmap=cmap_diff,
                norm=norm_diff,
                transform=datacrs,
                extend='both'
            )

            diff_handles.append(cf)

            if i == 0:
                ax.set_title(titles[j])

            # panel labels
            ax.text(
                0.05,
                0.96,
                labels[label_idx],
                transform=ax.transAxes,
                va='top',
                ha='left',
                bbox=bbox_dict
            )

            label_idx += 1

        # --------------------------------------------------------
        # Colorbar
        # --------------------------------------------------------
        cbax = fig.add_subplot(gs[i, -1])

        cb = Colorbar(
            ax=cbax,
            mappable=diff_handles[0],
            orientation="vertical"
        )

        units = (
            "days yr$^{-1}$"
            if varname == "ros"
            else hist["ccsm"][varname]
            .attrs.get("units", "")
        )
        if plot_type == "ros_frequency_clim":
            cb.set_label(
                f"Δ {units}"
            )
        else:
            cb.set_label(
                f"Δ{var_labels[i]} ({units})"
            )

    # ============================================================
    # Save Figure
    # ============================================================
    outname = cfg["outfile"].format(
        season=season,
        option=option
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

    seasons = ["ONDJFM"]
    options = ["strict", "flexible"]

    # --- 95th percentile climatology ---
    for season in seasons:

        print(
            f"Creating 95th percentile plot: {season}"
        )

        plot_change_grid(
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
                "ros_frequency_clim",
                season,
                option
            )

            print(
                f"Creating ROS intensity plot: "
                f"{option} {season}"
            )

            plot_change_grid(
                "ros_intensity_clim",
                season,
                option
            )