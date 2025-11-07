# colormaps.py
import numpy as np
import cmocean
import cmocean.cm as cmo

# Optional: your custom function (imported in your scripts)
from plotter import make_brgr_white_cmap  


def get_colormap_and_levels(clim_type, varname):
    """
    Return color levels and colormaps for a given climatology type and variable.

    Parameters
    ----------
    clim_type : str
        One of ['95th_percentile_clim', 'ros_intensity_clim', 'ros_frequency_clim'].
    varname : str
        Variable name (e.g., 'ivt', 'uv', 'freezing_level', 'pcpt', 'snow', etc.)

    Returns
    -------
    levs_clim : ndarray
    cmap_clim : Colormap
    levs_diff : ndarray
    cmap_diff : Colormap or (Colormap, norm)
    """

    # === 95th_percentile_clim ===
    if clim_type == "95th_percentile_clim":
        if varname == 'ivt':
            levs_clim = np.arange(150, 425, 25); cmap_clim = cmo.deep
            levs_diff = np.arange(0, 132, 12); cmap_diff = cmo.deep
        elif varname == 'uv':
            levs_clim = np.arange(0, 33, 3); cmap_clim = cmo.dense
            levs_diff = np.arange(-5, 6, 1); cmap_diff = cmo.balance
        elif varname == 'freezing_level':
            levs_clim = np.arange(1500, 2600, 100)
            cmap_clim = cmocean.tools.crop_by_percent(cmo.ice, 20, which='min')
            levs_diff = np.arange(0, 660, 60); cmap_diff = 'Reds'
        elif varname == 'pcpt':
            levs_clim = np.arange(0, 110, 10); cmap_clim = cmo.rain
            levs_diff = np.arange(-20, 24, 4); cmap_diff = 'BrBG'
        elif varname == 'snow':
            levs_clim = np.arange(0, 2200, 200); cmap_clim = cmo.rain
            levs_diff = np.arange(-500, 525, 100); cmap_diff = 'BrBG'
        else:
            levs_clim = np.linspace(0, 1, 10)
            levs_diff = np.linspace(-1, 1, 10)
            cmap_clim = cmo.tempo
            cmap_diff = cmo.balance
        return levs_clim, cmap_clim, levs_diff, cmap_diff

    # === ros_intensity_clim ===
    elif clim_type == "ros_intensity_clim":
        if varname == 'ros':
            levs_clim = np.arange(0, 11, 1); cmap_clim = cmo.rain
            levs_diff = np.arange(-5, 6, 1)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-1, 1))
        elif varname == 'pcpt':
            levs_clim = np.arange(0, 66, 6); cmap_clim = cmo.rain
            levs_diff = np.arange(-20, 24, 4)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-4, 4))
        elif varname == 'snow':
            levs_clim = np.arange(0, 220, 20); cmap_clim = cmo.rain
            levs_diff = np.arange(-2000, 2400, 400)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-4, 4))
        else:
            levs_clim = np.arange(0, 220, 20); cmap_clim = cmo.rain
            levs_diff = np.arange(-40, 48, 8)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-8, 8))
        return levs_clim, cmap_clim, levs_diff, cmap_diff

    # === ros_frequency_clim ===
    elif clim_type == "ros_frequency_clim":
        if varname == 'ros':
            levs_clim, cmap_clim = np.arange(0, 11, 1), cmo.rain
            levs_diff = np.arange(-5, 6, 1)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-1, 1))
        elif varname == 'ivt':
            levs_clim, cmap_clim = np.arange(0, 44, 4), cmo.deep
            levs_diff = np.arange(-20, 24, 4)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-4, 4))
        elif varname == 'pcpt':
            levs_clim, cmap_clim = np.arange(0, 88, 8), cmo.rain
            levs_diff = np.arange(-10, 12, 2)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-2, 2))
        elif varname == 'snow':
            levs_clim, cmap_clim = np.arange(0, 275, 25), cmo.rain
            levs_diff = np.arange(-40, 48, 8)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-2, 2))
        elif varname == 'delsnowh':
            levs_clim, cmap_clim = np.arange(0, 185, 18), cmo.rain
            levs_diff = np.arange(-40, 48, 8)
            cmap_diff, norm = make_brgr_white_cmap(levs_diff, (-8, 8))
        else:
            levs_clim = np.linspace(0, 1, 10)
            levs_diff = np.linspace(-1, 1, 10)
            cmap_clim = cmo.tempo
            cmap_diff = cmo.balance
        return levs_clim, cmap_clim, levs_diff, cmap_diff

    else:
        raise ValueError(f"Unknown climatology type: {clim_type}")
