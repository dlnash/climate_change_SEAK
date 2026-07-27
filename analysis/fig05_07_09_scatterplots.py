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
from wrf_utils import filter_vars
from colormaps import get_colormap_and_levels


# ---------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------
def plot_clim_diff_scatter(ssn, option, path_to_data,
                           fsuffix, plot_type):

    # -------------------------------------------------------------
    # Variable configuration
    # -------------------------------------------------------------
    if fsuffix == "ros_frequency_clim":
        varnames = ['ros', 'ivt', 'pcpt', 'delsnow', 'delsnowh']
        varname_lbl = ['ROS', 'IVT', 'Precip', 'ΔSWE', 'Snowmelt']

    elif fsuffix == "ros_intensity_clim":
        varnames = ["ros", "pcpt", "snow",
                    "delsnowh", "ros_intensity"]
        varname_lbl = ['ROS', 'Precip',
                       'SWE', 'Snowmelt', 'Intensity']

    elif fsuffix == "95th_percentile_clim":
        varnames = [
            "ivt",
            "uv",
            "freezing_level",
            "pcpt",
            "snow"
        ]

        varname_lbl = [
            'IVT',
            r'$\overline{\mathrm{UV}}$',
            'Freezing Level',
            'Precipitation',
            'SWE'
        ]

    else:
        raise ValueError(f"Unknown fsuffix: {fsuffix}")

    # -------------------------------------------------------------
    # Figure setup
    # -------------------------------------------------------------
    fig = plt.figure(figsize=(9, 14))
    fig.dpi = 300

    current_dpi = 300
    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi)**0.3
    set_font(current_dpi, scaling_factor)

    nrows = len(varnames)

    gs = GridSpec(
        nrows, 4,
        width_ratios=[1, 1, 1, 0.05],
        hspace=0.15,
        wspace=0.08
    )

    # -------------------------------------------------------------
    # Read terrain
    # -------------------------------------------------------------
    elev_fname = os.path.join(
        path_to_data,
        "downloads/SEAK-WRF/geo_southeast.nc"
    )

    elev_ds = xr.open_dataset(elev_fname)
    elev_ds = filter_vars(
        elev_ds.squeeze(),
        elev_fname,
        "hgt"
    )

    terrain = elev_ds["hgt"].values
    landmask = elev_ds["landmask"].values
    land_mask_bool = landmask == 1

    mask = land_mask_bool.flatten()
    x_flat = terrain.flatten()[mask]

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

    column_titles = [
        "ΔCCSM (Future − Historical)",
        "ΔGFDL (Future − Historical)",
        "Multi-model Mean Δ"
    ]

    # =============================================================
    # Main Loop
    # =============================================================
    for i, varname in enumerate(varnames):

        print(f"Plotting {varname}...")

        # ---------------------------------------------------------
        # Load datasets
        # ---------------------------------------------------------
        def load_model(model, scenario):

            if fsuffix == "95th_percentile_clim":
                fname = os.path.join(
                    path_to_data,
                    f"preprocessed/SEAK-WRF/"
                    f"{model}/{scenario}/trends/"
                    f"{varname}_{model}_{ssn}_{fsuffix}.nc"
                )

            else:
                fname = os.path.join(
                    path_to_data,
                    f"preprocessed/SEAK-WRF/"
                    f"{model}/{scenario}/trends/"
                    f"snow_{model}_{ssn}_{option}_{fsuffix}.nc"
                )

            return xr.open_dataset(fname)

        ccsm_hist = load_model("ccsm", "hist")
        ccsm_future = load_model("ccsm", "rcp85")

        gfdl_hist = load_model("gfdl", "hist")
        gfdl_future = load_model("gfdl", "rcp85")

        # ---------------------------------------------------------
        # Compute differences
        # ---------------------------------------------------------
        diff_ccsm = (
            ccsm_future[varname]
            - ccsm_hist[varname]
        )

        diff_gfdl = (
            gfdl_future[varname]
            - gfdl_hist[varname]
        )

        diff_mean = xr.concat(
            [diff_ccsm, diff_gfdl],
            dim="model"
        ).mean("model")

        diffs = [
            diff_ccsm,
            diff_gfdl,
            diff_mean
        ]

        lat2d = diff_ccsm.lat.values
        y_flat = lat2d.flatten()[mask]

        # ---------------------------------------------------------
        # Colormap
        # ---------------------------------------------------------
        (
            _,
            _,
            _,
            levs_diff,
            cmap_diff,
            norm_diff
        ) = get_colormap_and_levels(
            fsuffix,
            varname
        )

        diff_handles = []

        # ---------------------------------------------------------
        # Columns
        # ---------------------------------------------------------
        for j, diff in enumerate(diffs):

            ax = fig.add_subplot(gs[i, j])

            ax.set_xlim(0, 1900)
            ax.set_ylim(54.5, 60)

            cfield = diff.values.flatten()[mask]

            if plot_type == "hexbin":

                sc = ax.hexbin(
                    x_flat,
                    y_flat,
                    C=cfield,
                    gridsize=50,
                    cmap=cmap_diff,
                    norm=norm_diff,
                    reduce_C_function=np.mean,
                    mincnt=1
                )

            else:

                sc = ax.scatter(
                    x_flat,
                    y_flat,
                    c=cfield,
                    cmap=cmap_diff,
                    vmin=levs_diff.min(),
                    vmax=levs_diff.max(),
                    s=0.5
                )

            diff_handles.append(sc)

            if i == 0:
                ax.set_title(column_titles[j])

            if i == nrows - 1:
                ax.set_xlabel("Elevation (m)")

            if j == 0:
                ax.set_ylabel("Latitude (°N)")
            else:
                ax.set_yticklabels([])

            ax.text(
                0.05,
                0.95,
                labels[label_idx],
                transform=ax.transAxes,
                va='top',
                ha='left',
                bbox=bbox_dict
            )

            label_idx += 1

        # ---------------------------------------------------------
        # Shared colorbar
        # ---------------------------------------------------------
        cbax = fig.add_subplot(gs[i, -1])

        cb = Colorbar(
            ax=cbax,
            mappable=diff_handles[0],
            orientation="vertical",
            extend='both'
        )

        units = ccsm_hist[varname].attrs.get(
            "units",
            ""
        )

        if fsuffix == "ros_frequency_clim":
            cb.set_label(
                f"Δ {units}"
            )
        else:
            cb.set_label(
                f"Δ{varname_lbl[i]} ({units})"
            )

    # =============================================================
    # Save Figure
    # =============================================================
    if fsuffix == "95th_percentile_clim":
        outname = (
            f"../figs/clim/"
            f"{ssn}_{fsuffix}_{plot_type}.png"
        )
    else:
        outname = (
            f"../figs/ros_{option}/{ssn}/"
            f"{ssn}_{fsuffix}_{plot_type}.png"
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

    ssn_lst = ["ONDJFM"]

    options = ["strict", "flexible"]

    fsuffix_lst = [
        "95th_percentile_clim",
        "ros_intensity_clim",
        "ros_frequency_clim"
    ]

    for option in options:
        for ssn in ssn_lst:
            for fsuffix in fsuffix_lst:
                for plot_type in ["hexbin"]:

                    print(
                        f"Creating plot for "
                        f"{option}, "
                        f"{ssn}, "
                        f"{fsuffix}, "
                        f"{plot_type}"
                    )

                    plot_clim_diff_scatter(
                        ssn,
                        option,
                        path_to_data,
                        fsuffix,
                        plot_type
                    )