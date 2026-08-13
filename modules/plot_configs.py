# ============================================================
# CONFIG
# ============================================================
MODELS = ["ccsm", "gfdl"]
SCENARIOS = {
    "historical": "hist",
    "future": "rcp85"
}

PLOT_CONFIGS = {

    "95th_percentile_clim": {
        "varnames": [
            "ivt",
            "uv",
            "freezing_level",
            "pcpt",
            "snow"
        ],

        "labels": [
            "IVT",
            r"$\overline{\mathrm{UV}}$",
            "Freezing Level",
            "Precipitation",
            "SWE"
        ],

        "filename_template":
            "{var}_{model}_{season}_95th_percentile_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/clim/{season}_95th_percentile_clim_change_{plot_name}.png",

        "cmap_key":
            "95th_percentile_clim",

        "extend":
            "both",
    },

    "ros_frequency_clim": {
        "varnames": [
            "ros",
            "ivt",
            "pcpt",
            "delsnow",
            "delsnowh"
        ],

        "labels": [
            "ROS",
            "IVT",
            "Precip",
            r"$\Delta$ SWE",
            "Snowmelt"
        ],

        "filename_template":
            "snow_{model}_{season}_{option}_ros_frequency_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/ros_{option}/{season}/{season}_ROS_FREQ_CHANGE_{plot_name}.png",

        "cmap_key":
            "ros_frequency_clim",

        "extend":
            "both",
    },

    "ros_intensity_clim": {
        "varnames": [
            "ros",
            "pcpt",
            "snow",
            "delsnowh",
            "ros_intensity"
        ],

        "labels": [
            "ROS",
            "Precip",
            "SWE",
            "Snowmelt",
            "Intensity"
        ],

        "filename_template":
            "snow_{model}_{season}_{option}_ros_intensity_clim.nc",

        "subdir":
            "{model}/{scenario}/trends",

        "outfile":
            "../figs/ros_{option}/{season}/{season}_ROS_INTENSITY_CHANGE_{plot_name}.png",

        "cmap_key":
            "ros_intensity_clim",

        "extend":
            "both",
    },
}