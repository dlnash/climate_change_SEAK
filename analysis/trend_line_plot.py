"""
Filename:    trend_line_plot.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Script to plot trend line plots for each variable given a model name.
"""

import sys, os
import numpy as np
import xarray as xr
import pandas as pd

import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Path to modules
sys.path.append('../modules/')
import globalvars
from wrf_utils import load_preprocessed_WRF_data, subset_wrf_ds

model = 'gfdl'
path_to_data = globalvars.path_to_data
varnames = ['ivt', 'pcpt', 'freezing_level', 'uv925']

df_lst = []
for varname in varnames:
    print(varname)
    # --- load data ---
    ds = load_preprocessed_WRF_data(model, varname, anomaly=False)
    
    # --- subset to community locations ---
    ds = subset_wrf_ds(ds)
    
    # --- compute 95th percentile for var for each year ---
    ds = ds.groupby("time.year").quantile(0.95, dim="time").rename(year="time")
    
    # --- change to tidy dataframe format ---
    df = ds.to_dataframe().reset_index()
    df_lst.append(df)

def add_regplot(ax, df, varname):
    xvar, yvar = "time", varname
    
    # Compute R² per location
    r2_dict = {}
    for loc, sub in df.groupby("location"):
        slope, intercept, r_value, p_value, std_err = linregress(sub[xvar], sub[yvar])
        r2_dict[loc] = r_value**2
        # plot regression line
        sns.regplot(
            data=sub, x=xvar, y=yvar,
            ax=ax, label=f"{loc} (R²={r2_dict[loc]:.2f})"
        )
    
    ax.legend(title="Location")
    ax.set_title(varname)
    return ax

# --- create figure with GridSpec ---
fig = plt.figure(figsize=(9, 9), dpi=300)
gs = GridSpec(2, 2, wspace=0.3, hspace=0.3)

for i, (varname, df) in enumerate(zip(varnames, df_lst)):
    row, col = divmod(i, 2)
    ax = fig.add_subplot(gs[row, col])
    if varname == "uv925":
        varname = "uv"
    add_regplot(ax, df, varname)

fig.savefig(f"../figs/{model}_lineplot_trend.png", bbox_inches="tight")
plt.show()


