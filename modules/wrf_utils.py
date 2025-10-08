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

# Variables needed for each preprocessing target
KEEP_VARS: Dict[str, List[str]] = {
    "ivt": ["QVAPOR", "U", "V", "PSFC", "lat", "lon"],
    "uv925": ["U", "V", "interp_levels", "lat", "lon"],
    "freezing_level": ["T", "GHT", "lat", "lon"],
    "pcpt": ["PCPT", "lat", "lon"],
    "snow": ["SNOW", "SNOWH", "SNOWC", "PCPT", "lat", "lon"]
}

# Optional renaming rules for consistency
RENAME_MAP: Dict[str, Dict[str, str]] = {
    "uv925": {"interp_levels": "isobaricInhPa", "U": "u", "V": "v"},
    "pcpt": {"PCPT": "pcpt"},
    "ivt": {'interp_levels': 'isobaricInhPa', 'QVAPOR': 'q', 'V': 'v', 'U': 'u', 'PSFC': 'sp'},
    "freezing_level": {'interp_levels': 'isobaricInhPa'},
    "snow": {"SNOW": "snow", "PCPT": "pcpt", "SNOWH": "snowh", "SNOWC": "snowc"}
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

def load_preprocessed_WRF_data(model, varname, anomaly=False):
    
    if anomaly == False:
        ## read the non-anomaly data
        datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/{varname}/")
        fname_pattern = f'WRFDS_{varname}_*.nc'
    else:
        ## read the anomaly data
        datadir = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/{varname}/anomalies/")
        fname_pattern = f'daily_filtered_anomalies_{varname}_*.nc'
    
    ds = xr.open_mfdataset(datadir+fname_pattern,
                          engine='netcdf4',
                           combine='by_coords')
    
    ## rename time variable
    ds = ds.rename({'Time': 'time'})
    
    ## rename dims from south_north to lat and west_east to lon
    ds = ds.rename_dims({'south_north': 'y', 'west_east': 'x'})
    
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


