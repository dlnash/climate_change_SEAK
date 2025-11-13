## Preprocess Data

These scripts can be run after downloading the data required as outlined in `../downloads/` .
For each of the models (cfsr, ccsm, gfdl) preprocess the necessary variables:
    - ivt
    - 925 hPa wind
    - freezing level
    - precipitation
    - snow

1. Preprocess WRF data to extract variables needed for analysis with the sbatch script: `preprocess/WRF/run_preprocess_WRF.slurm`
2. Compute the differences between historical and future of the variables with the sbatch script: `preprocess/compute_freq_intensity/run_compute_freq_intensity_changes.slurm`