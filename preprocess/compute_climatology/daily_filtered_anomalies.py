"""
Filename:    daily_filtered_anomalies.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions to filter annual climatology of daily gridded/time series data using harmonics.

"""
import sys
import numpy as np
import xarray as xr

# Path to modules
sys.path.append('../../modules/')

from harmonics import harmonic

#######################
# Variables to Update #
#######################

model = 'CCSM'
varname = 'ivt' ## 'ivt' or '700z'

### Constants ###
outdir = '/expanse/nfs/cw3e/cwp140/preprocessed/SEAK-WRF/{0}/{1}/'.format(model, varname)
datadir = '/expanse/nfs/cw3e/cwp140/preprocessed/SEAK-WRF/{0}/{1}/'.format(model, varname)
fmt = '.nc'

print('Step 1: Reading data...')
filename_pattern = datadir + 'WRFDS_{0}_*.nc'.format(varname)
print(filename_pattern)
ds = xr.open_mfdataset(filename_pattern,
                       engine='netcdf4',
                       combine='by_coords')

print('ds size in GB {:0.2f}\n'.format(ds.nbytes / 1e9))
ds = ds.sortby('Time')

## calculate annual climatology
print('Step 2: Calculating annual climatology...')
clim_mean = ds.groupby('Time.dayofyear').mean('Time')
clim_std = ds.groupby('Time.dayofyear').std('Time')

# # ### Save Daily Climatology as netcdf

## Save Daily Climatology as netcdf
clim_path = outdir + 'daily_mean_clim_' + varname + fmt
clim_mean.load().to_netcdf(path=clim_path, mode = 'w', format='NETCDF4')

## Save Standard Deviation as netcdf
std_path = outdir + 'daily_std_clim_' + varname + fmt
clim_std.load().to_netcdf(path=std_path, mode = 'w', format='NETCDF4')

## filter annual climatology
print('Step 3: Filtering annual climatology...')
filtered_clim = harmonic(clim_mean)

## Save Filtered Climatology as netcdf
print('Step 4: Saving filtered climatology...')
clim_path = outdir + 'filtered_daily_mean_clim_' + varname + fmt
filtered_clim.load().to_netcdf(path=clim_path, mode = 'w', format='NETCDF4')

## Calculate Anomalies
print('Step 5: Calculating anomalies...')
anomalies = ds.unify_chunks().groupby('Time.dayofyear') - filtered_clim
# long_term_mean = anomalies.mean('time', skipna=True)

## Write anomalies to yearly file
print('Step 6: Writing anomalies to yearly files...')
years, datasets = zip(*anomalies.groupby('Time.year'))
paths = [outdir +'anomalies/daily_filtered_anomalies_{0}_%s.nc'.format(varname) % y for y in years]
xr.save_mfdataset(datasets, paths)