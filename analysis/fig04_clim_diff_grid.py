"""
Filename:    clim_diff_grid.py
Author:      Deanna Nash
Description: Plot future-historical differences for CCSM and GFDL
               Col 1 = CCSM Δ (future - historical)
               Col 2 = GFDL Δ (future - historical)
               Col 3 = Multi-model mean Δ
"""

import os
import string
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar

import sys
sys.path.append('../modules/')

import globalvars
from plotter import set_font
from colormaps import get_colormap_and_levels


# ============================================================
# Helper
# ============================================================
def load_climatology(path_to_data, model, scenario, varname, season):

    fname = os.path.join(
        path_to_data,
        f"preprocessed/SEAK-WRF/"
        f"{model}/{scenario}/trends/"
        f"{varname}_{model}_{season}_95th_percentile_clim.nc"
    )

    if not os.path.exists(fname):
        raise FileNotFoundError(fname)

    return xr.open_dataset(fname)


# ============================================================
# Plot Function
# ============================================================
def plot_clim_diff_grid(
    varnames,
    ssn,
    path_to_data,
    lonmin=-141.,
    lonmax=-130.,
    latmin=54.5,
    latmax=60.,
):

    mapcrs = ccrs.Mercator()
    datacrs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 14), dpi=300)

    current_dpi = 300
    base_dpi = 100
    scaling_factor = (current_dpi / base_dpi) ** 0.3
    set_font(current_dpi, scaling_factor)

    nrows = len(varnames)
    ncols = 3

    gs = GridSpec(
        nrows,
        ncols + 1,
        width_ratios=[1, 1, 1, 0.05],
        hspace=0.05,
        wspace=0.05
    )

    model_titles = [
        "CCSM Future − Historical",
        "GFDL Future − Historical",
        "Multi-model Mean Δ",
    ]

    varname_lbl = [
        "IVT",
        r"$\overline{UV}$",
        "Freezing Level",
        "Precipitation",
        "SWE",
    ]

    labels = list(string.ascii_lowercase)
    label_idx = 0

    bbox_dict = dict(
        facecolor='white',
        edgecolor='k',
        boxstyle='circle,pad=0.3',
        alpha=1.
    )

    # ============================================================
    # MAIN LOOP
    # ============================================================
    for i, varname in enumerate(varnames):

        print(f"Plotting {varname}...")

        # --------------------------------------------------------
        # Load historical + future climatologies
        # --------------------------------------------------------
        ds_ccsm_hist = load_climatology(
            path_to_data, "ccsm", "hist", varname, ssn
        )

        ds_ccsm_future = load_climatology(
            path_to_data, "ccsm", "rcp85", varname, ssn
        )

        ds_gfdl_hist = load_climatology(
            path_to_data, "gfdl", "hist", varname, ssn
        )

        ds_gfdl_future = load_climatology(
            path_to_data, "gfdl", "rcp85", varname, ssn
        )

        # --------------------------------------------------------
        # Compute deltas
        # --------------------------------------------------------
        ccsm_diff = (
            ds_ccsm_future[varname]
            - ds_ccsm_hist[varname]
        )

        gfdl_diff = (
            ds_gfdl_future[varname]
            - ds_gfdl_hist[varname]
        )

        multi_model_mean = (
            ccsm_diff + gfdl_diff
        ) / 2

        fields = [
            ccsm_diff,
            gfdl_diff,
            multi_model_mean,
        ]

        lons = ds_ccsm_hist.lon.values
        lats = ds_ccsm_hist.lat.values

        # --------------------------------------------------------
        # Colormap
        # --------------------------------------------------------
        (
            _,
            _,
            _,
            levs_diff,
            cmap_diff,
            norm_diff
        ) = get_colormap_and_levels(
            "95th_percentile_clim",
            varname
        )

        # --------------------------------------------------------
        # Plot columns
        # --------------------------------------------------------
        cf = None

        for j in range(ncols):

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

            if i == 0:
                ax.set_title(model_titles[j])

            if i == nrows - 1:
                ax.set_xlabel("Longitude")

            # subplot label
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
        # Shared row colorbar
        # --------------------------------------------------------
        cbax = fig.add_subplot(gs[i, -1])

        cb = Colorbar(
            ax=cbax,
            mappable=cf,
            orientation="vertical"
        )

        units = ds_ccsm_hist[varname].attrs.get(
            "units",
            ""
        )

        cb.set_label(
            f"Δ{varname_lbl[i]} ({units})"
        )

    # ============================================================
    # Save
    # ============================================================
    outname = os.path.join(
        "../figs/clim",
        f"{ssn}_95th_percentile_clim_delta.png"
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

    varnames = [
        "ivt",
        "uv",
        "freezing_level",
        "pcpt",
        "snow",
    ]

    ssn_lst = ["ONDJFM"]

    for ssn in ssn_lst:
        print(f"Creating plot for {ssn}")

        plot_clim_diff_grid(
            varnames,
            ssn,
            path_to_data
        )