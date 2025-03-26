######################################################################
# Filename:    create_job_configs.py
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to create .yaml configuration file to run job_array with slurm for downloading GEFSv12 Reforecast Data
#
######################################################################

## import libraries
import pandas as pd
import numpy as np
import yaml
from itertools import chain

yr_start = 2030
yr_end = 2060

## create list of dates to download in parallel
start_date = pd.to_datetime('{0}-01-01'.format(yr_start))
end_date = pd.to_datetime('{0}-12-31'.format(yr_start))

## make a list of dates between start_date and end_date
date_lst = pd.date_range(start_date, end_date, freq='1D')
			
jobcounter = 0
filecounter = 0
## loop through to create dictionary for each job
d_lst = []
dest_lst = []
njob_lst = []
for i, date in enumerate(date_lst):
    yr = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")
    
    jobcounter += 1
    d = {'job_{0}'.format(jobcounter):
         {'year': yr,
          'month': month,
          'day': day
          }}
    d_lst.append(d)
    
    if (jobcounter == 999):
        filecounter += 1
        ## merge all the dictionaries to one
        dest = dict(chain.from_iterable(map(dict.items, d_lst)))
        njob_lst.append(len(d_lst))
        ## write to .yaml file and close
        file=open("config_{0}.yaml".format(str(filecounter)),"w")
        yaml.dump(dest,file, allow_unicode=True, default_flow_style=None)
        file.close()
        
        ## reset jobcounter and d_lst
        jobcounter = 0
        d_lst = []
        
## now save the final config
filecounter += 1
## merge all the dictionaries to one
dest = dict(chain.from_iterable(map(dict.items, d_lst)))
njob_lst.append(len(d_lst))
## write to .yaml file and close
file=open("config_{0}.yaml".format(str(filecounter)),"w")
yaml.dump(dest,file, allow_unicode=True, default_flow_style=None)
file.close()

## create calls.txt for config_1(-8)

for i, njobs in enumerate(njob_lst):
    call_str_lst = []
    for j, job in enumerate(range(1, njobs+1, 1)):
        call_string = "python getWRF_batch.py config_{0}.yaml 'job_{1}'".format(i+1, j+1)
        call_str_lst.append(call_string)
        
    ## now write those lines to a text file
    with open('calls_{0}.txt'.format(i+1), 'w',encoding='utf-8') as f:
        for line in call_str_lst:
            f.write(line)
            f.write('\n')
        f.close()