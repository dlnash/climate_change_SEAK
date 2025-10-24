#!/usr/bin/python3
"""
Filename:    wrf_preprocess.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: preprocess functions for AK 4 km WRF simulations downloaded from https://registry.opendata.aws/wrf-se-alaska-snap/
"""

## Imports
import numpy as np
import xarray as xr
import pandas as pd
from wrf import interplevel

from wrf_utils import filter_vars
from time_helpers import find_date_based_on_filename

def calc_IVT_manual(ds):
    '''
    Calculate IVT manually (not using scipy.integrate)
    This is in case you need to remove values below the surface
     '''
    lats = ds.lat
    lons = ds.lon
    pressure = ds.isobaricInhPa.values*100 # convert from hPa to Pa
    dp = np.diff(pressure) # delta pressure
    g = 9.81 # gravity constant
    
    qu_lst = []
    qv_lst = []
    # enumerate through pressure levels so we select the layers
    for i, pres in enumerate(ds.isobaricInhPa.values[:-1]):
        pres2 = ds.isobaricInhPa.values[i+1]
        tmp = ds.sel(isobaricInhPa=[pres, pres2]) # select layer
        tmp = tmp.mean(dim='isobaricInhPa', skipna=True) # average q, u, v in layer
        # calculate ivtu in layer
        qu = ((tmp.q*tmp.u*dp[i])/g)*-1
        qu_lst.append(qu)
        # calculate ivtv in layer
        qv = ((tmp.q*tmp.v*dp[i])/g)*-1
        qv_lst.append(qv)
    
    ## add up u component of ivt from each layer
    qu = xr.concat(qu_lst, pd.Index(pressure[:-1], name="pres"))
    qu = qu.sum('pres')
    qu.name = 'ivtu'
    
    # ## add up v component of ivt from each layer
    qv = xr.concat(qv_lst, pd.Index(pressure[:-1], name="pres"))
    qv = qv.sum('pres')
    qv.name = 'ivtv'
    
    ## calculate IVT magnitude
    ivt = np.sqrt(qu**2 + qv**2)
    ivt.name = 'ivt'

    ds = xr.merge([qu, qv, ivt], compat='no_conflicts')

    return ds
    
def preprocess_WRF_ivt(ds: xr.Dataset, fname: str) -> xr.Dataset:
    """
    Preprocess WRF dataset to compute IVTu and IVTv.

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset.
    fname : str
        Filename (used to extract valid date/time).

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with only IVTu and IVTv.
    """
    
    ds = filter_vars(ds, fname, "ivt", find_date_func=find_date_based_on_filename)

    # subset to the levels we need for interpolation
    ds = ds.sel(isobaricInhPa=slice(300, 1000)) 
    ds = ds.reindex(isobaricInhPa=ds.isobaricInhPa[::-1]) ## flip pressure levels
    ## mask values below surface pressure
    print('Masking values below surface ....')
    varlst = ['q', 'u', 'v']
    for i, varname in enumerate(varlst):
        ds[varname] = ds[varname].where(ds[varname].isobaricInhPa < ds.sp/100., drop=False)
    ds = calc_IVT_manual(ds)
    
    return ds

def preprocess_WRF_uv(ds: xr.Dataset, fname: str, lev: float = 925.0) -> xr.Dataset:
    """
    Preprocess WRF dataset to extract 850-hPa winds (u and v).

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset.
    fname : str
        Filename (used to extract valid date/time).

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with only U and V at the requesting.
    """
    
    ds = filter_vars(ds, fname, "uv925", find_date_func=find_date_based_on_filename)

    # Subset to requested pressure level
    if "isobaricInhPa" not in ds.dims:
        raise KeyError("'isobaricInhPa' dimension not found in dataset")
    if lev not in ds["isobaricInhPa"].values:
        raise ValueError(f"Level {lev} hPa not available in dataset")

    ds = ds.sel(isobaricInhPa=lev)
    
    ## calculate UV magnitude
    ds['uv'] = np.sqrt(ds['u']**2 + ds['v']**2)

    return ds

def preprocess_WRF_pcpt(ds: xr.Dataset, fname: str) -> xr.Dataset:
    """
    Preprocess WRF dataset to extract precipitation (PCPT).

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset.
    fname : str
        Filename (used to extract valid date/time).

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with precipitation only.
    """
    return filter_vars(ds, fname, "pcpt", find_date_func=find_date_based_on_filename)

def preprocess_WRF_snow(ds: xr.Dataset, fname: str) -> xr.Dataset:
    """
    Preprocess WRF dataset to extract snow related variables (i.e., SWE, snow depth, snow cover).

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset.
    fname : str
        Filename (used to extract valid date/time).

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with snow only.
    """

    ds = filter_vars(ds, fname, "snow", find_date_func=find_date_based_on_filename)

    return ds

def preprocess_WRF_freezing_level(ds: xr.Dataset, fname: str) -> xr.Dataset:
    """
    Preprocess WRF dataset to compute freezing level or the height in m at which temperature was 0*C.

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset.
    fname : str
        Filename (used to extract valid date/time).

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with freezing level in m.
    """

    ds = filter_vars(ds, fname, "freezing_level", find_date_func=find_date_based_on_filename)
    ds = ds.sel(isobaricInhPa=slice(200, 1000)) # only interested in freezing level below 200 hPa
    ## need 2 3D arrays for input
    ## reshape output arrays
    ntime, nlev, nlat, nlon = ds.GHT.shape
    gh = ds.GHT.values # height in meters
    t = ds.T.values # temperature in *C
    
    # interpolate gh to temperature = 0
    interp_var = interplevel(gh, t, [0])
    ds['freezing_level'] = (("Time", "south_north", "west_east"), interp_var.data[np.newaxis, :, :])
    keep = ["freezing_level", "lat", "lon", "time"]
    ds = ds[[v for v in keep if v in ds]]
    
    return ds

def preprocess_WRF_ros(ds: xr.Dataset, temporal_resolution: str = 'daily') -> xr.Dataset:
    """
    Preprocess WRF dataset to compute rain-on-snow (ROS) diagnostics.

    Parameters
    ----------
    ds : xr.Dataset
        Input WRF dataset (already preprocessed snow dataset).
    temporal_resolution : str, optional
        Temporal resolution of output. Options:
        - 'daily' (default): returns daily time steps.
        - 'yearly': aggregates ROS counts (sum) and other variables (mean) per year.

    Returns
    -------
    xr.Dataset
        Preprocessed dataset with variables:
        - rain: rainfall (mm)
        - delsnowh: daily snow depth change (mm)
        - ros: binary indicator (1 = ROS event)
        - ros_intensity: combined rainfall + snowmelt (mm)
        If yearly, variables are aggregated accordingly.
    """

    # --- Convert snow depth to mm ---
    ds['snowh'] = ds['snowh'] * 1000
    ds['snowh'].attrs.update({'units': 'mm', 'long_name': 'Snow depth'})

    # --- Snow change per timestep ---
    ds['delsnow'] = (ds['snow'].diff('time'))* -1 # positive values indicate melting, sublimation, or compaction
    ds['delsnowh'] = (ds['snowh'].diff('time'))* -1 # positive values mean melt
    ds['delsnowh'].attrs.update({'units': 'mm', 'long_name': 'Change in snow depth'})

    # --- Rain-on-snow (ROS) indicator ---
    ds['ros'] = ((ds['pcpt'] > 25.4) & (ds['snowc'] > 0) & (ds['delsnow'] > 6.35)).astype(int)
    ds['ros'].attrs.update({'long_name': 'Rain-on-snow event', 'description': '1 = ROS event'})

    # --- ROS intensity (rain + snowmelt) ---
    ds['ros_intensity'] = ds['pcpt'].where(ds['ros'] == 1) + (ds['delsnowh'].where(ds['ros'] == 1)) 
    ds['ros_intensity'].attrs.update({'units': 'mm', 'long_name': 'ROS intensity (prec + snowmelt)'})

    # --- Temporal aggregation ---
    if temporal_resolution.lower() == 'yearly':
        # --- 1. Sum total ROS events per year ---
        da_ros = ds['ros'].groupby('time.year').sum('time')
    
        # --- 2. Identify which variables to average ---
        vars_to_process = ['pcpt', 'delsnow', 'delsnowh', 'ros_intensity']
    
        yearly_means = []
        for var in vars_to_process:
            da = ds[var].where(ds['ros'] == 1)                # keep only ROS days
            da_yearly = da.groupby('time.year').mean('time')  # average per year
            yearly_means.append(da_yearly)
    
        # --- 3. Merge results into single dataset ---
        ds_out = xr.merge([da_ros] + yearly_means)
        ds_out = ds_out.rename({'year': 'time'})

    elif temporal_resolution.lower() == 'daily':
        ds_out = ds  # retain daily data

    else:
        raise ValueError("temporal_resolution must be either 'daily' or 'yearly'.")

    return ds_out

def compute_ros_frequency(ds):
    vars_info = [
        ('pcpt', 25.4, 'Precip'),
        ('delsnow', 6.35, 'ΔSWE'),
        ('ivt', 250., 'IVT'),
        ('delsnowh', 25.4, 'Snowmelt'),
    ]

    freq_lst = [ds['ros'].groupby('time.year').sum('time').rename({'year': 'time'})]

    for var, thres, label in vars_info:
        exceed = ds[var] > thres
        freq = exceed.groupby('time.year').sum('time').rename({'year': 'time'})
        freq.attrs['units'] = f"days yr$^{{-1}}$ {label} > {int(thres)}"
        freq_lst.append(freq)

    ds_out = xr.merge(freq_lst).mean('time', keep_attrs=True)

    return ds_out



