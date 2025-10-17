#!/bin/bash
######################################################################
# Filename:    download_ERA5.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to download all the necessary ERA5 files
# To run, activate conda env with "conda activate cds", then run the script with "bash download_ERA5.sh"
######################################################################

# names of configuration dictionaries to loop through
array=(
ivt_hourly
)

# now loop through each configuration dictionary to download the ERA5 data
for i in ${!array[*]}
do 
    inconfig="${array[$i]}"
    python getERA5_batch.py ${inconfig}
done