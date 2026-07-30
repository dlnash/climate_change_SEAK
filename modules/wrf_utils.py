"""
Filename:    wrf_utils.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Shared constants and helper functions for WRF preprocessing.

"""
import sys, os
import xarray as xr
import numpy as np
from typing import Dict, List, Optional
import globalvars
path_to_data = globalvars.path_to_data
from plot_configs import PLOT_CONFIGS, MODELS, SCENARIOS

# Variables needed for each preprocessing target
KEEP_VARS: Dict[str, List[str]] = {
    "ivt": ["QVAPOR", "U", "V", "PSFC", "lat", "lon"],
    "uv925": ["U", "V", "interp_levels", "lat", "lon"],
    "freezing_level": ["T", "GHT", "lat", "lon"],
    "pcpt": ["PCPT", "lat", "lon"],
    "snow": ["SNOW", "SNOWH", "SNOWC", "PCPT", "lat", "lon"],
    "hgt": ["HGT_M", "LANDMASK", "XLAT_M", "XLONG_M"]
}

# Optional renaming rules for consistency
RENAME_MAP: Dict[str, Dict[str, str]] = {
    "uv925": {"interp_levels": "isobaricInhPa", "U": "u", "V": "v"},
    "pcpt": {"PCPT": "pcpt"},
    "ivt": {'interp_levels': 'isobaricInhPa', 'QVAPOR': 'q', 'V': 'v', 'U': 'u', 'PSFC': 'sp'},
    "freezing_level": {'interp_levels': 'isobaricInhPa'},
    "snow": {"SNOW": "snow", "PCPT": "pcpt", "SNOWH": "snowh", "SNOWC": "snowc"},
    "hgt": {"HGT_M": "hgt", "LANDMASK": "landmask", "XLAT_M": "lat", "XLONG_M": "lon"},
    # other variables could also go here if needed
}


def filter_vars(
    ds: xr.Dataset,
    fname: str,
    varname: str,
    rename: Optional[Dict[str, str]] = None,
    find_date_func=None,
) -> xr.Dataset:
    """
    Keep only required variables, rename them, and add Time coordinate.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset.
    fname : str
        Filename (used to extract valid date/time).
    varname : str
        Variable key (must exist in KEEP_VARS).
    rename : dict, optional
        Mapping of {old_name: new_name} for consistency.
    find_date_func : callable, optional
        Function that extracts datetime from filename.

    Returns
    -------
    xr.Dataset
        Filtered dataset with standardized variable names and Time coordinate.
    """
    if varname not in KEEP_VARS:
        raise ValueError(f"Unknown variable type: {varname}")

    # Filter to keep list
    keep = KEEP_VARS[varname]
    ds = ds[[v for v in keep if v in ds]]

    if not ds.data_vars:
        raise KeyError(f"No variables found for {varname} in dataset.")

    # Apply renaming rules
    rename_map = rename or RENAME_MAP.get(varname, {})
    ds = ds.rename({k: v for k, v in rename_map.items() if k in ds})

    # Add time coordinate
    if find_date_func:
        ds = ds.assign_coords(Time=find_date_func(fname))

    # Drop XTIME coord safely (some of the files don't have it)
    ds = ds.drop_vars(['XTIME'], errors="ignore")
    # Make lat and lon variables coordinates and replace XLAT and XLONG
    ds = ds.assign_coords(
        lat=(("south_north", "west_east"), ds["lat"].values),
        lon=(("south_north", "west_east"), ds["lon"].values),
    )

    ds = ds.drop_vars(["XLAT", "XLONG"], errors="ignore")

    return ds


def find_nearest_indices(ds, lat, lon):
    # Function to find nearest grid indices on 2D curvilinear grid
    dist = (ds['lat'] - lat)**2 + (ds['lon'] - lon)**2
    iy, ix = np.unravel_index(dist.argmin(), dist.shape)
    return iy, ix

def subset_wrf_ds(ds):
    # input lists
    lon_lst = [-135.4519, -135.3277, -135.8894, -139.671, -133.1358, -132.4009]
    lat_lst = [58.1122, 59.4538, 59.3988, 59.5121, 55.4769, 55.5400]
    lbl_lst = ['Hoonah', 'Skagway', 'Klukwan', 'Yakutat', 'Craig', 'Kasaan']
    
    # Get indices for each requested point
    indices = [find_nearest_indices(ds, la, lo) for la, lo in zip(lat_lst, lon_lst)]
    iy = [i[0] for i in indices]
    ix = [i[1] for i in indices]
    
    # Subset the dataset
    subset = ds.isel(
        y=xr.DataArray(iy, dims="location"),
        x=xr.DataArray(ix, dims="location")
    )
    
    # Attach the requested coordinates and labels to the new "location" dim
    subset = subset.assign_coords(
        location=lbl_lst,
        lat=("location", lat_lst),
        lon=("location", lon_lst)
    )

    return subset


def load_preprocessed_dataset(
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

def load_preprocessed_wrf_metrics(varname, plot_type, season, option):
    # --------------------------------------------------------
    # Load datasets
    # --------------------------------------------------------
    hist = {}
    future = {}

    for model in MODELS:

        hist[model] = load_preprocessed_dataset(
            path_to_data,
            model,
            "hist",
            plot_type,
            varname,
            season,
            option
        )

        future[model] = load_preprocessed_dataset(
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
    # Historical MMM
    hist_mmm = xr.concat(
        [
            hist["ccsm"][varname],
            hist["gfdl"][varname]
        ],
        dim="model"
    ).mean("model")
    
    # Future MMM
    future_mmm = xr.concat(
        [
            future["ccsm"][varname],
            future["gfdl"][varname]
        ],
        dim="model"
    ).mean("model")
    
    # Climate change signal
    mmm_change = future_mmm - hist_mmm
    
    # Future model disagreement
    future_difference = (
        future["ccsm"][varname]
        - future["gfdl"][varname]
    )

    fields = [
        hist_mmm,
        mmm_change,
        future_difference
    ]
    
    titles = [
        "Historical\nMulti-model Mean",
        "Future Change\nMulti-model Mean",
        "Future Projection\n(CCSM - GFDL)"
    ]

    units = (
                "days yr$^{-1}$"
                if varname == "ros"
                else hist["ccsm"][varname]
                .attrs.get("units", "")
            )

    return fields, titles, units

def read_wrf_terrain():
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

    return mask, x_flat