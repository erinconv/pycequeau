# from __future__ import annotations

from pyproj import Transformer, Proj
import numpy as np
import xarray as xr
from osgeo import gdal, osr


def get_proj_code(ds: gdal.Dataset):
    """_summary_

    Args:
        ds (gdal.Dataset): _description_

    Returns:
        _type_: _description_
    """
    datum = osr.SpatialReference(wkt=ds.GetProjection())
    EPSG_code = datum.GetAttrValue('AUTHORITY', 1)
    EPSG = f"EPSG:{EPSG_code}"
    return EPSG


def get_transformer(source_crs: str,
                    target_crs: str):
    """_summary_

    Args:
        source_crs (str): _description_
        target_crs (str): _description_

    Returns:
        _type_: _description_
    """
    proj_source = Proj(source_crs, preserve_units=False)
    proj_target = Proj(target_crs, preserve_units=False)
    return Transformer.from_proj(proj_source, proj_target)


def latlon_to_utm(lon: list,
                  lat: list,
                  target: str) -> tuple:
    """_summary_

    Args:
        lon (list): _description_
        lat (list): _description_
        target (str): _description_

    Returns:
        tuple: _description_
    """
    transformer = get_transformer("EPSG:4326", target)
    lat_utm, lon_utm = transformer.transform(lat, lon)
    return lat_utm, lon_utm


def utm_to_latlon(x: list,
                  y: list,
                  source: str) -> tuple:
    """_summary_

    Args:
        x (list): _description_
        y (list): _description_
        source (str): _description_

    Returns:
        tuple: _description_
    """
    transformer = get_transformer(source, "EPSG:4326")
    lat_utm, lon_utm = transformer.transform(y, x)
    return lat_utm, lon_utm
