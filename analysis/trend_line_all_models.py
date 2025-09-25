"""
Filename:    trend_line_plot.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Script to plot trend line plots for each variable and model,
             and compute R² values saved to CSV.
"""

import sys, os
import numpy as np
import xarray as xr
import pandas as pd

import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt

# Path to modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data, subset_wrf_ds

# ----------------------
# Settings
# ----------------------
models = ["cfsr", "ccsm", "gfdl"]       # your three models
varnames = ["ivt", "pcpt", "freezing_level", "uv"]
path_to_data = globalvars.path_to_data

# ----------------------
# Function to plot one panel
# ----------------------
def add_regplot(ax, df, varname, show_legend=False):
    """Add regression plots for each location, return dict of R² values."""
    xvar, yvar = "time", varname
    r2_dict = {}

    for loc, sub in df.groupby("location"):
        slope, intercept, r_value, p_value, std_err = linregress(sub[xvar], sub[yvar])
        r2_dict[loc] = round(r_value**2, 2)  # two decimal places

        sns.regplot(
            data=sub, x=xvar, y=yvar,
            ax=ax, label=loc,  # just location name, legend handled globally
            scatter=True, ci=None
        )

    ax.set_title(varname if varname != "uv925" else "uv")
    if not show_legend:
        ax.get_legend().remove() if ax.get_legend() else None
    return r2_dict

# ----------------------
# Prepare figure
# ----------------------
fig, axes = plt.subplots(
    nrows=len(models), ncols=len(varnames),
    figsize=(16, 12), dpi=300,
    sharex=False, sharey=False
)

# Store R² results
r2_results = {model: {} for model in models}

# Loop over models and variables
for i, model in enumerate(models):
    for j, varname in enumerate(varnames):
        ax = axes[i, j]

        # --- load data ---
        if varname == "uv":
            load_varname = "uv925"
        else:
            load_varname = varname
        ds = load_preprocessed_WRF_data(model, load_varname, anomaly=False)

        # --- subset to community locations ---
        ds = subset_wrf_ds(ds)

        # --- compute 95th percentile for var for each year ---
        ds = ds.groupby("time.year").quantile(0.95, dim="time").rename(year="time")

        # --- change to tidy dataframe format ---
        df = ds.to_dataframe().reset_index()

        # Plot + collect R²
        r2_dict = add_regplot(ax, df, varname, show_legend=(i == 0 and j == 0))
        r2_results[model][varname] = r2_dict  # store all communities

        # Label rows/cols
        if j == 0:
            ax.set_ylabel(model.upper())
        if i == len(models) - 1:
            ax.set_xlabel("Year")

# ----------------------
# Global legend
# ----------------------
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, title="Community", loc="lower center", ncol=6)

fig.tight_layout(rect=[0, 0.05, 1, 1])  # leave space for legend

# Save figure
fig.savefig("../figs/trend_lineplot_allmodels.png", bbox_inches="tight")
plt.close(fig)

# ----------------------
# Save R² to CSV (all communities separate)
# ----------------------
rows = []
for model in models:
    for var in varnames:
        for loc, r2 in r2_results[model][var].items():
            rows.append({"model": model, "variable": var, "location": loc, "R2": r2})

r2_df = pd.DataFrame(rows)
r2_df.to_csv("../figs/trend_lineplot_r2_full.csv", index=False)

