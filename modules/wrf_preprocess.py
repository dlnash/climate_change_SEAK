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
import re

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

    ds = xr.merge([qu, qv, ivt, lats, lons])

    return ds
    
def preprocess_WRF_ivt(ds, fname):
    drop_varlst = ['SLP', 'U10', 'V10', 'T2', 'TSLB',
                   'Q2', 'TSK', 'SNOW', 'SNOWH', 'SNOWC',
                   'PW', 'LH', 'HFX', 'ALBEDO', 'LWDNB',
                   'LWDNBC', 'LWUPB', 'LWUPBC', 'SWDNB', 'SWDNBC',
                   'SWUPB', 'SWUPBC', 'PCPT', 'ACSNOW', 'TMAX',
                   'TMIN', 'T', 'CLDFRA', 'GHT', 'SMOIS', 'SH2O']
    ds = ds.drop_vars(drop_varlst)
    ds = ds.rename({'interp_levels': 'isobaricInhPa', 'QVAPOR': 'q', 'V': 'v', 'U': 'u', 'PSFC': 'sp'})
    ds = ds.sel(isobaricInhPa=slice(300, 1000)) # subset to the levels we need for interpolation
    ds = ds.reindex(isobaricInhPa=ds.isobaricInhPa[::-1]) ## flip pressure levels
    ## mask values below surface pressure
    print('Masking values below surface ....')
    varlst = ['q', 'u', 'v']
    for i, varname in enumerate(varlst):
        ds[varname] = ds[varname].where(ds[varname].isobaricInhPa < ds.sp/100., drop=False)
    ds = calc_IVT_manual(ds)

    numbers = re.findall(r'\d+', fname)
    df = pd.DataFrame([numbers], columns=['Year', 'Month', 'Day'])
    date = pd.to_datetime(df)
    date = date.values

    ds['Time'] = date
    
    return ds
