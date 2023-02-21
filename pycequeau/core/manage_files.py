import os
import pandas as pd
import numpy as np
import xarray as xr
from osgeo import gdal, ogr, osr
from .units import (
    units_ERA,
    units_CORDEX
    )
import sys


def dict_netCDF(files_path: str) -> dict:
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
        nc_file_var =list(dr.keys())[0]
        ds = ds.merge(dr)
        # Put attributes
        ds[nc_file_var].attrs = dr[nc_file_var].attrs
    return ds

    
    pass


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
        nc_file_var =list(dr.keys())[0]
        # print(dr)
        ds = ds.merge(dr)
        # Put attributes
        ds[nc_file_var].attrs = dr[nc_file_var].attrs
    return ds
