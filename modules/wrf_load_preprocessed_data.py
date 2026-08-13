## Imports
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd

import globalvars
path_to_data = globalvars.path_to_data
from plot_configs import PLOT_CONFIGS, MODELS, SCENARIOS

def load_preprocessed_dataset(
    path_to_data,
    model,
    scenario,
    plot_type,
    varname,
    season,
    option=None
):
    """Load dataset for a given plot type."""

    cfg = PLOT_CONFIGS[plot_type]

    datadir = os.path.join(
        path_to_data,
        "preprocessed",
        "SEAK-WRF",
        cfg["subdir"].format(
            model=model,
            scenario=scenario
        )
    )

    fname = cfg["filename_template"].format(
        var=varname,
        model=model,
        season=season,
        option=option
    )

    path = os.path.join(datadir, fname)

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    return xr.open_dataset(path)

def load_preprocessed_wrf_metrics(varname, plot_type, season, option):
    # --------------------------------------------------------
    # Load datasets
    # --------------------------------------------------------
    hist = {}
    future = {}

    for model in MODELS:

        hist[model] = load_preprocessed_dataset(
            path_to_data,
            model,
            "hist",
            plot_type,
            varname,
            season,
            option
        )

        future[model] = load_preprocessed_dataset(
            path_to_data,
            model,
            "rcp85",
            plot_type,
            varname,
            season,
            option
        )

    # --------------------------------------------------------
    # Compute differences
    # --------------------------------------------------------
    # Historical MMM
    hist_mmm = xr.concat(
        [
            hist["ccsm"][varname],
            hist["gfdl"][varname]
        ],
        dim="model"
    ).mean("model")
    
    # Future MMM
    future_mmm = xr.concat(
        [
            future["ccsm"][varname],
            future["gfdl"][varname]
        ],
        dim="model"
    ).mean("model")
    
    # Climate change signal
    mmm_change = future_mmm - hist_mmm
    
    # Future model disagreement
    future_difference = (
        future["ccsm"][varname]
        - future["gfdl"][varname]
    )

    fields = [
        hist_mmm,
        mmm_change,
        future_difference
    ]
    
    titles = [
        "Historical\nMulti-model Mean",
        "Future Change\nMulti-model Mean",
        "Future Projection\n(CCSM - GFDL)"
    ]

    units = (
                "days yr$^{-1}$"
                if varname == "ros"
                else hist["ccsm"][varname]
                .attrs.get("units", "")
            )

    return fields, titles, units

