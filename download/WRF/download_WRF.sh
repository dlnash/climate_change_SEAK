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
MODEL=$4

### Directory to save data in
PATH_TO_OUT="/cw3e/mead/projects/cwp140/data/downloads/SEAK-WRF/${MODEL}/"

case "$MODEL" in
  ccsm|gfdl)
    ### FOR CCSM OR GFDL RCP85 REANALYSIS
    PATH_TO_DATA="s3://wrf-se-ak-ar5/${MODEL}/rcp85/daily/${YEAR}/"
    ;;
  cfsr)
    ### FOR CFSR 4 km
    PATH_TO_DATA="s3://wrf-se-ak-ar5/${MODEL}/4km/daily/${YEAR}/"
    ;;
  *)
    echo "Unknown model: $MODEL"
    exit 1
    ;;
esac

INPUT_FNAME="${PATH_TO_DATA}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
OUTPUT_FNAME="${PATH_TO_OUT}WRFDS_${YEAR}-${MONTH}-${DAY}.nc"
# echo ${OUTPUT_FNAME}
echo ${INPUT_FNAME}
/cw3e/mead/projects/cwp140/aws/local/bin/aws s3 cp --region us-west-2 ${INPUT_FNAME} ${OUTPUT_FNAME} --no-sign-request