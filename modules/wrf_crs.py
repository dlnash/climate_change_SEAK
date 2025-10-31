"""
Filename:    wrf_crs.py
Author:      Deanna Nash, dnash@ucsd.edu
Description: Create a map crs for the WRF data

"""
import xarray as xr
from pyproj import CRS
from pathlib import Path
import globalvars

def create_wrf_crs():
    data_path = Path(globalvars.path_to_data) / 'downloads' / 'SEAK-WRF' / 'geo_southeast.nc'
    wrf = xr.open_dataset(data_path)

    attrs = wrf.attrs
    map_proj = attrs.get("MAP_PROJ")
    
    if map_proj == 1:  # Lambert Conformal
        crs = CRS.from_proj4(
            f"+proj=lcc +lat_1={attrs['TRUELAT1']} +lat_2={attrs['TRUELAT2']} "
            f"+lat_0={attrs['CEN_LAT']} +lon_0={attrs['STAND_LON']} "
            f"+a=6370000 +b=6370000"
        )
    elif map_proj == 2:  # Polar Stereographic
        crs = CRS.from_proj4(
            f"+proj=stere +lat_0=90 +lat_ts={attrs['TRUELAT1']} "
            f"+lon_0={attrs['STAND_LON']} +a=6370000 +b=6370000"
        )
    elif map_proj == 3:  # Mercator
        crs = CRS.from_proj4(
            f"+proj=merc +lat_ts={attrs['TRUELAT1']} "
            f"+lon_0={attrs['STAND_LON']} +a=6370000 +b=6370000"
        )
    elif map_proj == 6:  # Lat/Lon
        crs = CRS.from_epsg(4326)
    else:
        raise ValueError(f"Unknown MAP_PROJ: {map_proj}")

    return crs
