#!/bin/bash
######################################################################
# Filename:    download_WRF.sh
# Author:      Deanna Nash dnash@ucsd.edu
# Description: Script to download Lader et al., 2020 SEAK WRF data
# https://registry.opendata.aws/wrf-se-ak-ar5/ (data link)
#
######################################################################

### Inputs
YEAR=$1
MONTH=$2
DAY=$3

# ### FOR CFSR REANALYSIS
# ### Set up paths
# PATH_TO_OUT="/expanse/lustre/scratch/dnash/temp_project/downloaded/WRF/"
# PATH_TO_DATA="s3://wrf-se-ak-ar5/cfsr/4km/daily/${YEAR}/"

# INPUT_FNAME="${PATH_TO_DATA}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
# OUTPUT_FNAME="${PATH_TO_OUT}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
# # echo ${OUTPUT_FNAME}
# /expanse/nfs/cw3e/cwp140/aws/local/bin/aws s3 cp --region us-west-2 ${INPUT_FNAME} ${OUTPUT_FNAME} --no-sign-request

### FOR CCSM RCP85 REANALYSIS
### Set up paths
PATH_TO_OUT="/expanse/lustre/scratch/dnash/temp_project/downloaded/WRF/"
PATH_TO_DATA="s3://wrf-se-ak-ar5/ccsm/rcp85/daily/${YEAR}/"

INPUT_FNAME="${PATH_TO_DATA}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
OUTPUT_FNAME="${PATH_TO_OUT}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
# echo ${OUTPUT_FNAME}
/expanse/nfs/cw3e/cwp140/aws/local/bin/aws s3 cp --region us-west-2 ${INPUT_FNAME} ${OUTPUT_FNAME} --no-sign-request