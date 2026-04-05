from typing import Any


import os
import pandas as pd
import numpy as np
import xarray as xr
from osgeo import gdal, ogr, osr

# GDAL sidecars appended to the full raster filename (e.g. merged_WBM.tif.aux.xml)
_RASTER_SIDECAR_SUFFIXES = (
    ".aux.xml",
    ".ovr",
    ".vat.dbf",
    ".aux",
)
# World file extensions use the stem only (e.g. merged_WBM.tfw)
_RASTER_WORLD_SUFFIXES = (".tfw", ".tifw")


def delete_raster_and_sidecars(tif_path: str) -> list[str]:
    """Delete a raster TIFF and all associated sidecar files (e.g. .aux.xml, .ovr).

    Removes the main file and common GDAL sidecars: .aux.xml, .ovr, .vat.dbf,
    .tfw, .tifw, .aux. Missing files are skipped without error.

    Args:
        tif_path: Path to the raster file (e.g. .../geographic/merged_WBM.tif).
            Can use forward or backslashes; .tif or .tiff.

    Returns:
        List of paths that were successfully deleted.
    """
    if not tif_path or not tif_path.strip():
        return []
    tif_path = os.path.normpath(tif_path.strip())
    if not (tif_path.lower().endswith(".tif") or tif_path.lower().endswith(".tiff")):
        return []
    directory = os.path.dirname(tif_path)
    base_name = os.path.basename(tif_path)  # e.g. merged_WBM.tif
    base_stem = os.path.splitext(base_name)[0]  # merged_WBM
    base_path = os.path.join(directory, base_stem)
    # Main raster: exact path (support both .tif and .tiff)
    candidates = [tif_path]
    # Sidecars that use the full raster filename as prefix (e.g. merged_WBM.tif.aux.xml)
    candidates.extend(
        os.path.join(directory, base_name + suffix)
        for suffix in _RASTER_SIDECAR_SUFFIXES
    )
    # World files use the stem only (e.g. merged_WBM.tfw)
    candidates.extend(
        os.path.join(directory, base_stem + suffix)
        for suffix in _RASTER_WORLD_SUFFIXES
    )
    deleted = []
    for path in candidates:
        if os.path.isfile(path):
            try:
                os.remove(path)
                deleted.append(path)
            except OSError:
                pass
    return deleted
