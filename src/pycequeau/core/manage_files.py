import os
import pandas as pd
import numpy as np
import xarray as xr
from osgeo import gdal, ogr, osr
from .netcdf import fix_calendar
from .units import (
    units_ERA,
    units_CORDEX
)

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


def dict_netCDF(files_path: str) -> dict:
    """_summary_

    Args:
        files_path (str): _description_

    Returns:
        dict: _description_
    """
    # List files whitin the folder
    files_list = os.listdir(files_path)
    filtered_list = [os.path.join(files_path, i)
                     for i in files_list if i.endswith(".nc")]
    vars_dict = {}
    for i, file in enumerate(filtered_list):
        var_name = "var"+str(i)
        vars_dict.update({var_name: xr.open_dataset(file,
                                                    engine="netcdf4")})
    return vars_dict


def get_CORDEX_Dataset(vars_dict: dict) -> xr.Dataset:
    """_summary_

    Args:
        vars_dict (dict): _description_

    Returns:
        xr.Dataset: _description_
    """
    keys = list(vars_dict.keys())
    # Create the dataset
    # dr1 = vars_dict.get(keys[0])
    # dr2 = vars_dict.get(keys[1])
    ds = xr.Dataset()
    for var in keys:
        # Get the var name

        # Check the data type to convert it to a commont data form
        # if vars_dict.get(var)[nc_file_var].dtype == "float64":
        #     # vars_dict.get(var)['lat'] = vars_dict.get(var)['lat'].astype(np.float32)
        #     # vars_dict.get(var)['lon'] = vars_dict.get(var)['lon'].astype(np.float32)
        #     ds = ds.merge(vars_dict.get(var).astype(np.float32))
        # else:
        #   ds = ds.merge(vars_dict.get(var))
        # Fix units
        dr = units_CORDEX(vars_dict.get(var))
        # Drop non needed variables
        variables = list(dr.variables)
        if "height" in variables:
            dr = dr.drop_vars("height")
        # Fix calendar
        dr = fix_calendar(dr)
        experiment = dr.attrs["experiment"]
        print(experiment)
        if dr.attrs["experiment"] == "historical":
            dr = dr.sel(time=slice("1980","2005"))
        elif experiment == "rcp85" or experiment== "RCP8.5":
            dr = dr.sel(time=slice("2040","2095"))
        
        nc_file_var = list(dr.keys())[0]
        ds = ds.merge(dr)
        # ds = xr.concat([ds,dr[nc_file_var]],dim="time")
        # Put attributes
        ds[nc_file_var].attrs = dr[nc_file_var].attrs
    return ds


def get_ERA_Dataset(vars_dict: dict) -> xr.Dataset:
    """_summary_

    Args:
        vars_dict (dict): _description_

    Returns:
        xr.Dataset: _description_
    """
    keys = list(vars_dict.keys())
    # Create the dataset
    # dr1 = vars_dict.get(keys[0])
    # dr2 = vars_dict.get(keys[1])
    ds = xr.Dataset()
    for var in keys:
        # Get the var name

        # Check the data type to convert it to a commont data form
        # if vars_dict.get(var)[nc_file_var].dtype == "float64":
        #     # vars_dict.get(var)['lat'] = vars_dict.get(var)['lat'].astype(np.float32)
        #     # vars_dict.get(var)['lon'] = vars_dict.get(var)['lon'].astype(np.float32)
        #     ds = ds.merge(vars_dict.get(var).astype(np.float32))
        # else:
        #   ds = ds.merge(vars_dict.get(var))
        # Fix units
        dr = units_ERA(vars_dict.get(var))
        nc_file_var = list(dr.keys())[0]
        # print(dr)
        ds = ds.merge(dr)
        # Put attributes
        ds[nc_file_var].attrs = dr[nc_file_var].attrs
    return ds
