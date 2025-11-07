#!/usr/bin/env python3
"""
Filename:    scatterplots.py
Author:      Deanna Nash (dnash@ucsd.edu)
Description: 
    Open preprocessed NetCDF files, open terrain file, and combine.
    Mask ocean points.
    For each variable, create scatterplots of values, sorted by elevation (x-axis) and latitude (y-axis).
    Colors of each point match the colors of the maps.
"""
# ---------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------
import os
import sys
import string
import pandas as pd
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
def plot_clim_diff_scatter(models, ssn, option, path_to_data, fsuffix, plot_type):
    """
    Create scatter plots for climatology and differences:
        Rows = variables
        Col 1 = CFSR climatology
        Col 2 = CCSM - CFSR
        Col 3 = GFDL - CFSR
    Scatter x-axis = Elevation, y-axis = latitude, color = variable value or difference
    """

    # --- define variable names for each suffix ---
    if fsuffix == "ros_frequency_clim":
        varnames = ['ros', 'ivt', 'pcpt', 'snow', 'delsnowh']
        varname_lbl = ['ROS', 'IVT', 'Precip', 'SWE', 'Snowmelt']
    elif fsuffix == "ros_intensity_clim":
        varnames = ["ros", "pcpt", "snow", "delsnowh", "ros_intensity"]
        varname_lbl = ['ROS', 'Precip', 'SWE', 'Snowmelt', 'Intensity']
    elif fsuffix == "95th_percentile_clim":
        varnames = ["ivt", "uv", "freezing_level", "pcpt", "snow"]
        varname_lbl = ['IVT', 'UV', 'Freezing Level', 'Precipitation', 'SWE']
    else:
        raise ValueError(f"Unknown fsuffix: {fsuffix}")

    fig = plt.figure(figsize=(10, 14))
    fig.dpi = 300
    current_dpi = 300
    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi)**0.3
    set_font(current_dpi, scaling_factor)

    nrows = len(varnames)
    ncols = len(models)

    gs = GridSpec(
        nrows, ncols + 3,
        height_ratios=[1]*nrows,
        width_ratios=[1, 0.05, 0.25, 1, 1, 0.05],
        hspace=0.15, wspace=0.05
    )

    # Load CFSR datasets
    cfsr_ds = {}
    for var in varnames:
        if fsuffix == "95th_percentile_clim":
            cfsr_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/{var}_cfsr_{ssn}_{fsuffix}.nc")
        else:
            cfsr_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/snow_cfsr_{ssn}_{option}_{fsuffix}.nc")

        if not os.path.exists(cfsr_path):
            print(f"⚠️ Missing CFSR file for {varname}: {cfsr_path}")
            continue
        cfsr_ds[var] = xr.open_dataset(cfsr_path)

    # Read elevation
    elev_fname = os.path.join(path_to_data, "downloads/SEAK-WRF/geo_southeast.nc")
    elev_ds = xr.open_dataset(elev_fname)
    elev_ds = filter_vars(elev_ds.squeeze(), elev_fname, "hgt")
    terrain = elev_ds["hgt"].values
    landmask = elev_ds["landmask"].values
    # --- Create mask for land points only ---
    land_mask_bool = landmask == 1

    
    labels = list(string.ascii_lowercase)
    label_idx = 0
    bbox_dict = dict(facecolor='white', edgecolor='k', boxstyle='circle,pad=0.3', alpha=1.)

    for i, varname in enumerate(varnames):
        print(f"Plotting scatter for {varname}...")

        lat2d = cfsr_ds[varname].lat.values
        lon2d = cfsr_ds[varname].lon.values

        # --- Apply land mask before flattening ---
        mask = land_mask_bool.flatten()
        x_flat = terrain.flatten()[mask]
        y_flat = lat2d.flatten()[mask]
        cfield = cfsr_ds[varname][varname].values.flatten()[mask]

        # Get base field and color levels
        levs_clim, cmap_clim, levs_diff, cmap_diff = get_colormap_and_levels(fsuffix, varname)

        # --- Column 1: CFSR scatter ---
        ax0 = fig.add_subplot(gs[i, 0])

        if plot_type == "hexbin":
            sc0 = ax0.hexbin(
                x_flat,
                y_flat,
                C=cfield,
                gridsize=35,        # adjust for resolution
                cmap=cmap_clim,
                vmin=levs_clim.min(), vmax=levs_clim.max(),
                reduce_C_function=np.mean,  # how to aggregate the data in each bin
                mincnt=1
            )
        else:
            sc0 = ax0.scatter(x_flat, y_flat, c=cfield, cmap=cmap_clim,
                          vmin=levs_clim.min(), vmax=levs_clim.max(), s=0.5)

        if i == 0:
            ax0.set_title("CFSR")
        if i == nrows - 1:
            ax0.set_xlabel("Elevation (m)")
        ax0.set_ylabel("Latitude (°N)")
        ax0.set_ylim(54.5, 60)
        ax0.set_xlim(0, 3300)
        ax0.text(0.05, 0.95, labels[label_idx], transform=ax0.transAxes,
                 va='top', ha='left', bbox=bbox_dict)
        label_idx += 1

        # --- Colorbar for CFSR Plots ---
        cbax_cfsr = fig.add_subplot(gs[i, 1])
        cb0 = Colorbar(ax=cbax_cfsr, mappable=sc0, orientation="vertical")
        if (fsuffix == "95th_percentile_clim") | (fsuffix == "ros_intensity_clim"):
            cb0.set_label(f"{varname_lbl[i]} ({cfsr_ds[varname][varname].attrs.get('units', '')})")
        else:
            cb0.set_label(f"{cfsr_ds[varname][varname].attrs.get('units', '')}")

        # --- Columns 2 & 3: CCSM and GFDL diffs ---
        diff_handles = []
        for j, model in enumerate(models[1:], start=3):
            ax = fig.add_subplot(gs[i, j])
            
            if fsuffix == "95th_percentile_clim":
                model_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/{varname}_{model}_{ssn}_{fsuffix}.nc")
            else:
                model_path = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/snow_{model}_{ssn}_{option}_{fsuffix}.nc")

            if not os.path.exists(model_path):
                print(f"⚠️ Missing file: {model_path}")
                continue

            ds_model = xr.open_dataset(model_path)
            if varname not in ds_model:
                print(f"⚠️ Variable {varname} not in {model_path}")
                continue
            
            diff = (ds_model[varname] - cfsr_ds[varname][varname]).values.flatten()[mask]

            if plot_type == "hexbin":
                sc = ax.hexbin(
                    x_flat,
                    y_flat,
                    C=diff,
                    gridsize=35,        # adjust for resolution
                    cmap=cmap_diff,
                    reduce_C_function=np.mean,  # how to aggregate the data in each bin
                    vmin=levs_diff.min(), vmax=levs_diff.max(),
                    mincnt=1
                )
            else:
                sc = ax.scatter(x_flat, y_flat, c=diff, cmap=cmap_diff,
                            vmin=levs_diff.min(), vmax=levs_diff.max(), s=0.5)
            
            diff_handles.append(sc)
            if i == 0:
                ax.set_title(f"{model.upper()} - CFSR")
            if i == nrows - 1:
                ax.set_xlabel("Elevation (m)")
            
            ax.set_ylim(54.5, 60)
            ax.set_xlim(0, 3300)
            ax.set_yticklabels([])
            ax.text(0.05, 0.95, labels[label_idx], transform=ax.transAxes,
                    va='top', ha='left', bbox=bbox_dict)
            label_idx += 1

        # --- Colorbar for Difference Plots ---
        cbax_diff = fig.add_subplot(gs[i, -1])
        cb1 = Colorbar(ax=cbax_diff, mappable=diff_handles[0], orientation="vertical")
        if (fsuffix == "95th_percentile_clim") | (fsuffix == "ros_intensity_clim"):
            cb1.set_label(f"Δ{varname_lbl[i]} ({cfsr_ds[varname][varname].attrs.get('units', '')})")
        else:
            cb1.set_label(f"Δ{cfsr_ds[varname][varname].attrs.get('units', '')}")
        
    if fsuffix == "95th_percentile_clim":
        outname = f"../figs/clim/{ssn}_{fsuffix}_{plot_type}.png"
    else:
        outname = f"../figs/ros_{option}/{ssn}/{ssn}_{fsuffix}_{plot_type}.png"
    fig.savefig(outname, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved figure: {outname}")



# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    ssn_lst = ["NDJFMA"]
    options = ["strict", "flexible"]
    fsuffix_lst = ["95th_percentile_clim", "ros_intensity_clim", "ros_frequency_clim"]
    
    for option in options:
        for ssn in ssn_lst:
            for fsuffix in fsuffix_lst:
                for plot_type in ["hexbin", "scatter"]:
                    print(f"Creating plot for {option}, {ssn}, {fsuffix}, {plot_type}")
                    plot_clim_diff_scatter(models, ssn, option, path_to_data, fsuffix, plot_type)
