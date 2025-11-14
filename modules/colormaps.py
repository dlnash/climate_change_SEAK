import numpy as np
import cmocean
import cmocean.cm as cmo
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt

# Optional: your custom diverging cmap function
from plotter import make_brgr_white_cmap  


def get_colormap_and_levels(clim_type, varname):
    """
    Return color levels, colormaps, and norms for a given climatology type and variable.

    Returns
    -------
    levs_clim : ndarray
    cmap_clim : Colormap
    norm_clim : BoundaryNorm
    levs_diff : ndarray
    cmap_diff : Colormap
    norm_diff : BoundaryNorm or None
    """

    norm_diff = None  # default
    norm_clim = None  # default

    # === 95th_percentile_clim ===
    if clim_type == "95th_percentile_clim":
        if varname == 'ivt':
            levs_clim = np.arange(150, 425, 25); cmap_clim = cmo.deep
            levs_diff = np.arange(0, 165, 15);  cmap_diff = cmo.deep
        elif varname == 'uv':
            levs_clim = np.arange(0, 22, 2);    cmap_clim = cmo.dense
            levs_diff = np.arange(-5, 6, 1);    cmap_diff = cmo.balance
        elif varname == 'freezing_level':
            levs_clim = np.arange(1500, 2600, 100)
            cmap_clim = cmocean.tools.crop_by_percent(cmo.ice, 20, which='min')
            levs_diff = np.arange(0, 550, 50);  cmap_diff = 'Reds'
        elif varname == 'pcpt':
            levs_clim = np.arange(0, 66, 6);    cmap_clim = cmo.rain
            levs_diff = np.arange(-15, 18, 3);  cmap_diff = 'BrBG'
        elif varname == 'snow':
            levs_clim = np.arange(0, 2200, 200); cmap_clim = cmo.rain
            levs_diff = np.arange(-500, 525, 100); cmap_diff = 'BrBG'
        else:
            levs_clim = np.linspace(0, 1, 10); cmap_clim = cmo.tempo
            levs_diff = np.linspace(-1, 1, 10); cmap_diff = cmo.balance

        # Add BoundaryNorm for clim and diff
        norm_clim = BoundaryNorm(levs_clim, ncolors=cmap_clim.N, clip=True)
        if isinstance(cmap_diff, str) or hasattr(cmap_diff, 'N'):
            norm_diff = BoundaryNorm(levs_diff, ncolors=plt.get_cmap(cmap_diff).N if isinstance(cmap_diff, str) else cmap_diff.N, clip=True)

    # === ros_intensity_clim ===
    elif clim_type == "ros_intensity_clim":
        if varname == 'ros':
            levs_clim = np.arange(1, 11, 1); cmap_clim = cmo.dense
            levs_diff = np.arange(-3, 3.5, 0.5)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-0.5, 0.5))
        elif varname == 'pcpt':
            levs_clim = np.arange(25, 80, 5); cmap_clim = cmo.rain
            levs_diff = np.arange(-25, 30, 5)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-5, 5))
        elif varname == 'snow':
            levs_clim = np.arange(0, 550, 50); cmap_clim = cmo.rain
            levs_diff = np.arange(-250, 300, 50)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-50, 50))
        elif varname == 'delsnowh':
            levs_clim = np.arange(0, 165, 15); cmap_clim = cmo.rain
            levs_diff = np.arange(-50, 60, 10)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-10, 10))
        elif varname == 'ros_intensity':
            levs_clim = np.arange(60, 180+12, 12); cmap_clim = cmo.rain
            levs_diff = np.arange(-100, 120, 20)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-10, 10))

        norm_clim = BoundaryNorm(levs_clim, ncolors=cmap_clim.N, clip=True)

    # === ros_frequency_clim ===
    elif clim_type == "ros_frequency_clim":
        if varname == 'ros':
            levs_clim = np.arange(1, 11, 1); cmap_clim = cmo.dense
            levs_diff = np.arange(-3, 3.5, 0.5)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-0.5, 0.5))
        elif varname == 'ivt':
            levs_clim, cmap_clim = np.arange(0, 33, 3), cmo.deep
            levs_diff = np.arange(-20, 24, 4)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-4, 4))
        elif varname == 'pcpt':
            levs_clim, cmap_clim = np.arange(0, 44, 4), cmo.rain
            levs_diff = np.arange(-10, 12, 2)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-2, 2))
        elif varname == 'delsnow':
            levs_clim, cmap_clim = np.arange(0, 20, 2), cmo.rain
            levs_diff = np.arange(-15, 18, 3)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-3, 3))
        elif varname == 'delsnowh':
            levs_clim, cmap_clim = np.arange(0, 165, 15), cmo.rain
            levs_diff = np.arange(-25, 30, 5)
            cmap_diff, norm_diff = make_brgr_white_cmap(levs_diff, (-5, 5))

        norm_clim = BoundaryNorm(levs_clim, ncolors=cmap_clim.N, clip=True)

    else:
        raise ValueError(f"Unknown climatology type: {clim_type}")

    return levs_clim, cmap_clim, norm_clim, levs_diff, cmap_diff, norm_diff
