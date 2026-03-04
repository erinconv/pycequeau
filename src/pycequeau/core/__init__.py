"""Core module."""
from __future__ import annotations

from .cop_processor import CopernicusDEMProcessor
from .utils import (
    polygonize_raster,
    fix_geometry,
    rasterize_shp,
    rasterize_shp_as_byte,
    get_altitude_point,
    GetExtent,
    regrid_CE,
    rasterize_feature,
    get_index_list
)

__all__ = [
    "CopernicusDEMProcessor",
    "polygonize_raster",
    "fix_geometry", 
    "rasterize_shp",
    "rasterize_shp_as_byte",
    "get_altitude_point",
    "GetExtent",
    "regrid_CE",
    "rasterize_feature",
    "get_index_list"
]

# from . import missing