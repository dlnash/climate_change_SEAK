## Download Data

1. Download daily SEAK-WRF data (Lader et al., 2020) from AWS using the sbatch scripts in `/download/WRF/`
   To create the config and calls files, run `python create_job_config.py`. This creates six files.
    - `calls_ccsm_1.txt` and `config_ccsm_1.yaml` downloads all of the daily data from the 2031-2060 CCSM RCP 8.5 WRF run.
    - `calls_cfsr_1.txt` and `config_cfsr_1.yaml` downloads all of the daily data from the 1981-2019 CFSR historical run.
    - `calls_gfdl_1.txt` and `config_gfdl_1.yaml` downloads all of the daily data from the 2031-2060 GFDL RCP 8.5 WRF run.
   You can then use the SLURM script `run_download_WRF.slurm` to complete the downloads.