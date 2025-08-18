## Preprocess Data

These scripts can be run after downloading the data required as outlined in `../downloads/` .

1. Preprocess WRF data to extract variables needed for analysis with the sbatch script: `preprocess/WRF/run_preprocess_WRF.slurm`
    - ivt
    - 925 hPa wind
    - freezing level
2. Compute anomalies by removing the annual cycle via harmonics with the script: `preprocess/compute_climatology/daily_filtered_anomalies.py`
3. Compute the trends of the variables.