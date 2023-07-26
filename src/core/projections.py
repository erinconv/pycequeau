# from __future__ import annotations

from pyproj import Transformer, Proj
import numpy as np
import xarray as xr
from osgeo import gdal, osr


def get_proj_code(ds: gdal.Dataset):
    datum = osr.SpatialReference(wkt=ds.GetProjection())
    EPSG_code = datum.GetAttrValue('AUTHORITY', 1)
    EPSG = f"EPSG:{EPSG_code}"
    return EPSG


def get_transformer(source_crs: str,
                    target_crs: str):
    proj_source = Proj(source_crs, preserve_units=False)
    proj_target = Proj(target_crs, preserve_units=False)
    return Transformer.from_proj(proj_source, proj_target)


def latlon_to_utm(lon: list,
                  lat: list,
                  target: str) -> tuple:
    transformer = get_transformer("EPSG:4326", target)
    lat_utm, lon_utm = transformer.transform(lat, lon)
    return lat_utm, lon_utm


def utm_to_latlon(x: list,
                  y: list,
                  source: str) -> tuple:
    transformer = get_transformer(source, "EPSG:4326")
    lat_utm, lon_utm = transformer.transform(y, x)
    return lat_utm, lon_utm



