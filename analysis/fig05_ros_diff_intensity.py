"""
Filename:    ros_diff_intensity.py
Author:      Deanna Nash (adapted for multi-model comparison)
Description: Plot 3-column x 5-row comparison of CFSR climatology and model differences for ROS frequency and intensity
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
from plotter import set_font, make_brgr_white_cmap

# ============================================================
# Plot Function
# ============================================================
def plot_ros_diff_intensity(models, varnames, ssn, option, path_to_data,
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
    fig = plt.figure(figsize=(10, 14))
    fig.dpi = 300
    current_dpi=300
    base_dpi=100
    scaling_factor = (current_dpi / base_dpi)**0.3
    set_font(current_dpi, scaling_factor)
    
    nrows = len(varnames)
    ncols = len(models)

    gs = GridSpec(nrows, ncols+3, height_ratios=[1]*nrows, width_ratios=[1, 0.05, 0.25, 1, 1, 0.05],
                  hspace=0.05, wspace=0.05)

    # --- Load reference CFSR climatology once ---
    fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/cfsr/trends/snow_cfsr_{ssn}_{option}_ros_intensity_clim.nc")
    cfsr_ds = xr.open_dataset(fname)

    cfsr_ds['pcpt'].attrs['units'] = 'mm day-1'
    cfsr_ds['delsnow'].attrs['units'] = 'mm day-1'
    cfsr_ds['delsnowh'].attrs['units'] = 'mm day-1'
    cfsr_ds['ros_intensity'].attrs['units'] = 'mm day-1'
    cfsr_ds['ros'].attrs['units'] = '#'

    # ============================================================
    # MAIN LOOP
    # ============================================================
    varname_lbl = ['ROS', 'Precip', 'SWE', 'Snowmelt', 'Intensity']
    labels = list(string.ascii_lowercase)  # ['a', 'b', 'c', ...]
    label_idx = 0
    bbox_dict = dict(facecolor='white', edgecolor='k', boxstyle='circle,pad=0.3', alpha=1.)

    for i, varname in enumerate(varnames):
        print(f"Plotting {varname}...")
        lons = cfsr_ds.lon.values
        lats = cfsr_ds.lat.values

        # ----------------------------
        # Get base field and color levels
        # ----------------------------
        if (varname == 'ros'):
            levs_clim = np.arange(0, 11, 1)
            cmap_clim = cmo.rain
            levs_diff = np.arange(-5, 6, 1)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-1, 1))
        elif (varname == 'pcpt'):
            levs_clim = np.arange(0, 110, 10)
            cmap_clim = cmo.rain
            levs_diff = np.arange(-20, 24, 4)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-4, 4))
        elif (varname == 'delsnow'):
            levs_clim = np.arange(0, 44, 4)
            cmap_clim = cmo.rain
            levs_diff = np.arange(-20, 24, 4)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-4, 4))
        else:
            levs_clim = np.arange(0, 220, 20)
            cmap_clim = cmo.rain
            levs_diff = np.arange(-40, 48, 8)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-8, 8))

        # ============================================================
        # --- Column 1: CFSR Climatology ---
        # ============================================================
        ax0 = fig.add_subplot(gs[i, 0], projection=mapcrs)
        ax0.set_extent([lonmin, lonmax, latmin, latmax])
        ax0.coastlines(resolution="50m")
        ax0.add_feature(cfeature.BORDERS, linewidth=0.75, edgecolor='k')
        # ax0.gridlines(draw_labels=False)
        
        
        cfield = cfsr_ds[varname].values
        cf0 = ax0.contourf(lons, lats, cfield, levels=levs_clim, cmap=cmap_clim,
                           transform=datacrs, extend='max')
        
        if i == 0:
            ax0.set_title(f"CFSR")

        # --- add a, b, c labels --- 
        ax0.text(0.05, 0.96, f"{labels[label_idx]}", transform=ax0.transAxes,
                 va='top', ha='left', bbox=bbox_dict)
        label_idx += 1

        # --- CFSR colorbar (vertical to the right) ---
        cbax_cfsr = fig.add_subplot(gs[i, 1])
        cb0 = Colorbar(ax=cbax_cfsr, mappable=cf0, orientation="vertical")
        cb0.set_label(f"{varname_lbl[i]} ({cfsr_ds[varname].attrs.get('units', '')})")

        # ============================================================
        # --- Columns 2 & 3: CCSM & GFDL Differences ---
        # ============================================================
        diff_handles = []
        for j, model in enumerate(models[1:], start=3):
            ax = fig.add_subplot(gs[i, j], projection=mapcrs)
            ax.set_extent([lonmin, lonmax, latmin, latmax])
            ax.coastlines(resolution="50m")
            ax.add_feature(cfeature.BORDERS, linewidth=0.75, edgecolor='k')
            # ax.gridlines(draw_labels=False)

            fname = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/trends/snow_{model}_{ssn}_{option}_ros_intensity_clim.nc")
            ds_model = xr.open_dataset(fname)
            diff = ds_model[varname] - cfsr_ds[varname]
            cf = ax.contourf(lons, lats, diff, levels=levs_diff, cmap=cmap_diff,
                             norm=norm, transform=datacrs, extend='both')
            diff_handles.append(cf)
            if i == 0:
                ax.set_title(f"{model.upper()} - CFSR")
            
            # --- add a, b, c labels --- 
            ax.text(0.05, 0.96, f"{labels[label_idx]}", transform=ax.transAxes,
                    va='top', ha='left', bbox=bbox_dict)
            label_idx += 1


        # --- Shared colorbar for both diffs (cols 2 & 3) ---
        cbax_diff = fig.add_subplot(gs[i, -1])
        cb1 = Colorbar(ax=cbax_diff, mappable=diff_handles[0], orientation="vertical")
        cb1.set_label(f"Δ{varname_lbl[i]} ({cfsr_ds[varname].attrs.get('units', '')})")

    # ============================================================
    # Save Figure
    # ============================================================
    outname = f"../figs/ros_{option}/{ssn}/{ssn}_ROS_INTENSITY_DIFF.png"
    fig.savefig(outname, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved figure: {outname}")


# ============================================================
# Driver
# ============================================================
if __name__ == "__main__":
    path_to_data = globalvars.path_to_data
    models = ["cfsr", "ccsm", "gfdl"]
    varnames = ["ros", "pcpt", "delsnow", "delsnowh", "ros_intensity"]
    ssn_lst = ["DJF", "MAM", "JJA", "SON", "NDJFMA"]
    ssn_lst = ["NDJFMA"]
    option_lst = ['strict', 'flexible']
    for option in option_lst:
        for ssn in ssn_lst:
            print(f"Creating plot for {option} {ssn}")
            plot_ros_diff_intensity(models, varnames, ssn, option, path_to_data)
                