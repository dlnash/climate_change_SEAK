"""
Filename:    preprocess_WRF.py
Author:      Deanna Nash, dnash@ucsb.edu
Description: preprocess daily variables from SEAK-WRF data and save as yearly nc files
"""

## Imports
import os, sys
import yaml
import glob
import numpy as np
import xarray as xr
import shutil

# Necessary to add cwd to path when script run
# by SLURM (since it executes a copy)
sys.path.append(os.getcwd())

# Path to modules
sys.path.append('../../modules')
# Import my modules
from wrf_preprocess import preprocess_WRF_ivt

path_to_wrf = '/expanse/lustre/scratch/dnash/temp_project/downloaded/WRF/'
path_to_out = '/expanse/lustre/scratch/dnash/temp_project/preprocessed/'

# command for selecting variables **include XTIME,XLONG,XLAT in var list
# ncks -v varname1,varname2 input_file output_file

### Input ###
### Variable to Process ###
config_file = str(sys.argv[1]) # this is the config file name
job_info = str(sys.argv[2]) # this is the job name

config = yaml.load(open(config_file), Loader=yaml.SafeLoader) # read the file
ddict = config[job_info] # pull the job info from the dict

year = ddict['year']
output_varname = ddict['varname']
model = 'CCSM'

# get list of filenames that contain data from that year from current year folder
filenames = []
for name in glob.glob(path_to_wrf + 'WRFDS_{0}*'.format(str(year))):
    filenames.append(name)
# sort filenames so they are in chronological order
filenames = sorted(filenames)

ds_lst = []
for i, wrfin in enumerate(filenames):
    ds = xr.open_dataset(wrfin)

    if output_varname == 'ivt':
        ds = preprocess_WRF_ivt(ds, wrfin)
        ds_lst.append(ds)

new_ds = xr.concat(ds_lst, dim='Time')

## if XTIME is a variable still, drop it
try:
    new_ds = new_ds.drop_vars(["XTIME"])
except:
    pass

# write to netCDF
print('Writing', output_varname, ' to netCDF')
fname = os.path.join(path_to_out, 'SEAK-WRF/{2}/{0}/WRFDS_{0}_{1}.nc').format(output_varname, str(year), model)
new_ds.to_netcdf(path=fname, mode = 'w', format='NETCDF4')

### COPY FROM LUSTRE TO MAIN SPACE
outname = '/expanse/nfs/cw3e/cwp140/preprocessed/SEAK-WRF/{0}/{1}/WRFDS_{1}_{2}.nc'.format(model, output_varname, str(year))
print('Copying preprocessed data...')
print('... {0} to {1}'.format(fname, outname))
shutil.copy(fname, outname) # copy file over to data folder