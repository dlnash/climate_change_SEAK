"""
Filename:    utils.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: utility functions
"""

import math
import glob
import re
import numpy as np

def get_startmon_and_endmon(ssn):
    ## set start_mon and end_mon based on ssn
    if ssn == 'DJF':
        start_mon, end_mon = (12, 2)
    elif ssn == 'MAM':
        start_mon, end_mon = (3, 5)
    elif ssn == 'JJA':
        start_mon, end_mon = (6, 8)
    elif ssn == 'SON':
        start_mon, end_mon = (9, 11)
    elif ssn == 'NDJFMA':
        start_mon, end_mon = (11, 4)
    elif ssn == 'MJJASO':
        start_mon, end_mon = (5, 10)

    return start_mon, end_mon

def roundPartial(value, resolution):
    return np.round(value / resolution) * resolution

def round_latlon_degree(df, res):

    df['lon-round'] = roundPartial(df['longitude'], res)
    df['lat-round'] = roundPartial(df['latitude'], res)
    
    return df

def select_months_ds(ds, mon_s, mon_e, time_varname='time'):    
    # Select months from xarray dataset
    if mon_s > mon_e:
        idx = (ds[time_varname].dt.month >= mon_s) | (ds[time_varname].dt.month <= mon_e)
    else:
        idx = (ds[time_varname].dt.month >= mon_s) & (ds[time_varname].dt.month <= mon_e)
    
    if time_varname == 'time':
        ds = ds.sel(time=idx)
    elif time_varname == 'start_date':
        ds = ds.sel(start_date=idx)
    elif time_varname == 'date':
        ds = ds.sel(date=idx)
    return ds

def select_months_df(df, mon_s, mon_e):
    # Select months from pandas dataframe
    if mon_s > mon_e:
        idx = (df.index.month >= mon_s) | (df.index.month <= mon_e)
    else:
        idx = (df.index.month >= mon_s) & (df.index.month <= mon_e)

    df = df.loc[idx]
    
    return df 