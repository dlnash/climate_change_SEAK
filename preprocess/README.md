## Preprocess Data

These scripts can be run after downloading the data required as outlined in `../downloads/` .
For each of the models (cfsr, ccsm, gfdl) preprocess the necessary variables:
    - ivt
    - 850 hPa wind
    - freezing level
    - precipitation

1. Preprocess WRF data to extract variables needed for analysis with the sbatch script: `preprocess/WRF/run_preprocess_WRF.slurm`
2. Compute anomalies by removing the annual cycle via harmonics with the script: `preprocess/compute_climatology/daily_filtered_anomalies.py`
3. Find the dates of extreme precipitation with the script `preprocess/compute_prec_thres.py`. These will be the dates for the trend analysis.
4. Compute the trends of the variables with the sbatch script: `preprocess/compute_trends/compute_trends.py`