######################################################################
# Filename:    compute_prec_thres.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Reads 6-hr SEAK WRF precipitation data, resamples to daily, computes the 95th percentile, 
# then finds the dates where at least 25% of the land grid cells in the SEAK data had precipitation >= 95th percentile
# saves those dates as a csv for creating composites/trend plots etc.
#
######################################################################

import os, sys
import xarray as xr
import pandas as pd
import gc
sys.path.append('../modules/')
import globalvars

path_to_data = globalvars.path_to_data

## open all of the data
fname_pattern = path_to_data + f"preprocessed/SEAK-WRF/cfsr/pcpt/WRFDS_PCPT_*.nc"
ds = xr.open_mfdataset(fname_pattern)
ds = ds.resample(time='1d').sum()

## compute quantiles/percentiles
## need to rechunk so time is a single chunk
ds = ds.chunk(dict(time=-1))

## Calculate the percentiles
quantile = ds.quantile([0.95], dim=['time'], skipna=True)
quantile = quantile.compute()

# write to netCDF
fname = os.path.join(path_to_data, 'preprocessed/SEAK-WRF/cfsr/pcpt/pcpt_quantiles_daily.nc')
quantile.load().to_netcdf(path=fname, mode = 'w', format='NETCDF4')

lm = xr.open_dataset('/cw3e/mead/projects/cwp140/data/downloads/SEAK-WRF/geo_southeast.nc')
lm = lm['LANDMASK'].isel(Time=0)
## count how many grid cells total
total_grid_ct = lm.sum()
# Mask PCPT over ocean
wrf["PCPT"] = wrf["PCPT"].where(wrf["LANDMASK"] == 1)

# ## open all of the data
# fname_pattern = path_to_data + f"preprocessed/SEAK-WRF/cfsr/pcpt/WRFDS_PCPT_*.nc"
# wrf = xr.open_mfdataset(fname_pattern)
# wrf = wrf.resample(time='1d').sum()
# wrf = wrf.assign({"LANDMASK": (("y", "x"), lm.values)})

# ## open the precalculated quantiles
# fname = path_to_data +'preprocessed/SEAK-WRF/cfsr/pcpt/pcpt_quantiles_daily.nc'
# quantile = xr.open_dataset(fname)

## use where statement to get dates where precip >95th percentile
tmp = wrf.where(wrf['PCPT'] >= quantile['PCPT'])
tmp = tmp.compute()

## get number of grids where prec >=95th percentile
grid_ct = tmp['PCPT'].count(['y', 'x'])

## get list of dates when precipitation exceeded this threshold
grid_perc = (grid_ct/total_grid_ct)*100
final_dates = grid_perc.where(grid_perc >= 25, drop=True).time.values

d = {'date': final_dates}
df = pd.DataFrame(data=d)
df.to_csv(f"../out/PCPT_95th_25perc-cov_dates.csv")