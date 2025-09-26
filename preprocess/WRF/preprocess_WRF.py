"""
Filename:    preprocess_WRF.py
Author:      Deanna Nash, dnash@ucsb.edu
Description: Preprocess daily variables from SEAK-WRF data and save as yearly NetCDF files.
"""

# --- Imports ---
import os
import sys
import glob
import yaml
import xarray as xr

# Add cwd for SLURM execution (script runs from a copied location)
sys.path.append(os.getcwd())

# Path to custom modules
sys.path.append('../../modules')
import globalvars
import wrf_preprocess as wrf_prep

def main(config_file: str, job_info: str):
    """Main preprocessing workflow."""

    # --- Load configuration ---
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    ddict = config[job_info]
    year = ddict["year"]
    output_varname = ddict["varname"]
    model = ddict["model"]

    # --- Define paths ---
    path_to_data = globalvars.path_to_data
    path_to_wrf = os.path.join(path_to_data, f"downloads/SEAK-WRF/{model}/")
    path_to_out = os.path.join(path_to_data, f"preprocessed/SEAK-WRF/{model}/{output_varname}/")
    os.makedirs(path_to_out, exist_ok=True)

    # --- Collect input files ---
    filenames = sorted(glob.glob(os.path.join(path_to_wrf, f"WRFDS_{year}*")))
    if not filenames:
        raise FileNotFoundError(f"No WRF files found for year {year} in {path_to_wrf}")

    # --- Map variables to preprocess functions ---
    preprocess_map = {
        "ivt": wrf_prep.preprocess_WRF_ivt,
        "freezing_level": wrf_prep.preprocess_WRF_freezing_level,
        "uv925": wrf_prep.preprocess_WRF_uv,
        "pcpt": wrf_prep.preprocess_WRF_pcpt,
        "snow": wrf_prep.preprocess_WRF_snow
    }

    if output_varname not in preprocess_map:
        raise ValueError(f"Unsupported variable: {output_varname}")

    preprocess_func = preprocess_map[output_varname]

    # --- Process datasets ---
    ds_list = []
    for wrfin in filenames:
        with xr.open_dataset(wrfin) as ds:
            ds_proc = preprocess_func(ds, wrfin)
            ds_list.append(ds_proc)

    # --- Concatenate along time dimension ---
    new_ds = xr.concat(ds_list, dim="Time")

    # Drop XTIME if present
    new_ds = new_ds.drop_vars("XTIME", errors="ignore")

    # --- Save to NetCDF ---
    out_file = os.path.join(path_to_out, f"WRFDS_{output_varname}_{year}.nc")
    print(f"Writing {output_varname} for {year} to {out_file}")
    new_ds.to_netcdf(out_file, mode="w", format="NETCDF4")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python preprocess_WRF.py <config_file> <job_info>")
        sys.exit(1)

    config_file, job_info = sys.argv[1], sys.argv[2]
    main(config_file, job_info)
