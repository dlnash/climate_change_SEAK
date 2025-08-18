"""
Filename:    daily_filtered_anomalies.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions to filter annual climatology of daily gridded/time series data using harmonics.

"""
import numpy as np
import xarray as xr

def calc_harmonics_tseries(tseries):
    ''' calculates sum of first two harmonics to approximate annual cycle'''
    mtot = tseries.size
    a0 = tseries.mean().item() # calculate mean of series
    t = np.arange(0, mtot, 1, dtype='float')

    omegas = []
    coef_a = []
    coef_b = []
    
    # compute coefficients
    for p in np.arange(1, 4):
        omegas.append((2.*np.pi*np.float(p))/np.float(mtot))
        coef_a.append((2./mtot)*(tseries*np.cos((2*np.pi*p*t)/mtot)).sum())
        coef_b.append((2./mtot)*(tseries*np.sin((2*np.pi*p*t)/mtot)).sum())
        
    ## first harmonic
    tmp1 = coef_a[0].item()*np.cos(omegas[0]*t) + coef_b[0].item()*np.sin(omegas[0]*t)
    tsrecSC1 = a0 + tmp1
    ## second harmonic
    tmp2 = coef_a[1].item()*np.cos(omegas[1]*t) + coef_b[1].item()*np.sin(omegas[1]*t)
    tsrecSC2 = a0 + tmp2
    ## first plus second harmonic
    tsrecSC = a0 + tmp1 + tmp2
    
    return(tsrecSC1, tsrecSC2, tsrecSC)


def fourier(ds,p):
    time = ds.dayofyear
    n = len(time) 
 
    ds,time = xr.broadcast(ds,time)

    f = 2.*np.pi*p/n
    ft = f*time

    sum_a = ds*np.cos(ft - 1.)
    sum_b = ds*np.sin(ft - 1.)
    coef_a = (2./n)*sum_a.sum('dayofyear',skipna=True)
    coef_b = (2./n)*sum_b.sum('dayofyear',skipna=True)

    return ft,coef_a, coef_b

def harmonic(ds):
    ''' use this function to use the first two harmonics to approximate the annual cycle for an xarray ds'''
    a0 = ds.mean('dayofyear',skipna=True)

    #-First Harmonic
    p     = 1
    ft,coef_a,coef_b = fourier(ds,p)
    harm1 = a0 + coef_a*np.cos(ft-1.) + coef_b*np.sin(ft-1.)

    #-Second Harmonic
    p     = 2
    ft,coef_a,coef_b = fourier(ds,p)
    harm2 = a0 + coef_a*np.cos(ft-1.) + coef_b*np.sin(ft-1.)

    #-First plus second
    combo = harm1 + coef_a*np.cos(ft-1.) + coef_b*np.sin(ft-1.)

    return combo