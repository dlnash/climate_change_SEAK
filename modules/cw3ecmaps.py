"""
Filename:    cw3ecmaps.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Functions for cmaps from the CW3E website adapted from from https://github.com/samwisehawkins/nclcmaps
"""

import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter

__all__ = ['cw3ecmaps']


cw3e_cmaps =    {
                "ivt" :{ 
                        "colors":[[255, 255, 3],
                                [255, 229, 3],
                                [255, 201, 2],
                                [255, 175, 2],
                                [255, 131, 1],
                                [255, 79, 1], 
                                [255, 24, 1],  
                                [235, 1, 7],  
                                [185, 0, 55], 
                                [134, 0, 99],
                                [86, 0, 137]],
                        "bounds":[250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600, 5000],
                        "ticks":[250, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600],
                        "label": "IVT (kg m$^{-1}$ s$^{-1}$)",
                        },
                "iwv" : { 
                        "colors":[[5, 0, 206], 
                                [8, 25, 255], 
                                [9, 113, 255], 
                                [10, 196, 255],
                                [11, 249, 237],
                                [6, 226, 160], 
                                [3, 201, 81], 
                                [9, 184, 6], 
                                [97, 208, 3],  
                                [177, 232, 2], 
                                [255, 255, 3], 
                                [255, 225, 3], 
                                [255, 197, 2], 
                                [255, 166, 2], 
                                [255, 112, 1], 
                                [255, 56, 1], 
                                [255, 1, 0], 
                                [209, 0, 33], 
                                [155, 0, 81], 
                                [100, 0, 126], 
                                [128, 59, 167]],
                        "bounds":[20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 200],
                        "ticks":[20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60],
                        "label": "IWV (mm)"
                        },
}

def cmap(cbarname):
        data = np.array(cw3e_cmaps[cbarname]["colors"])
        data = data / np.max(data)
        cmap = ListedColormap(data, name=cbarname)
        bnds = cw3e_cmaps[cbarname]["bounds"]
        norm = mcolors.BoundaryNorm(bnds, cmap.N)
        cbarticks = cw3e_cmaps[cbarname]["ticks"]
        cbarlbl = cw3e_cmaps[cbarname]["label"]
        return cmap, norm, bnds, cbarticks, cbarlbl