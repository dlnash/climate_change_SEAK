"""
Filename:    getWRF_batch.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Download Lader et al., 2020 SEAK WRF data based on input configuration dictionary.
"""
import sys
import yaml
import subprocess
import calendar

### Imports config name from argument when submit
yaml_doc = sys.argv[1]
config_name = sys.argv[2]

# import configuration file for dictionary choice
config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
ddict = config[config_name]

year = ddict['year']
month = ddict['month']

## loop through days here
def get_days_in_month(year, month):
    '''Returns a list of days in the given month and year.'''
    num_days = calendar.monthrange(int(year), int(month))[1]
    return list(range(1, num_days + 1))

day_lst = get_days_in_month(year, month)

for i, day in enumerate(day_lst):
    day = str(day).zfill(2)
    ## run download_WRF.sh to download data 
    bash_script = "/home/dnash/repos/climate_change_SEAK/download/WRF/download_WRF.sh"
    print(subprocess.run([bash_script, year, month, day]))
