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
from .matlab_utils import (
    to_mat_compatible,
    dataframe_to_struct_array,
    mat_field_value,
    mat_to_py,
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
    "get_index_list",
    "to_mat_compatible",
    "dataframe_to_struct_array",
    "mat_field_value",
    "mat_to_py",
]

# from . import missing
