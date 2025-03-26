"""
Filename:    getWRF_batch.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Download Lader et al., 2020 SEAK WRF data based on input configuration dictionary.
"""
import sys
import yaml
import subprocess

### Imports config name from argument when submit
yaml_doc = sys.argv[1]
config_name = sys.argv[2]

# import configuration file for dictionary choice
config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
ddict = config[config_name]

year = ddict['year']
month = ddict['month']
day = ddict['day']

## run download_WRF.sh to download data 
bash_script = "/home/dnash/repos/SEAK_AR_impacts/downloads/WRF/download_WRF.sh"
print(subprocess.run([bash_script, year, month, day]))
