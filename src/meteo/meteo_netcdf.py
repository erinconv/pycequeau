from __future__ import annotations

import os
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point
from osgeo import gdal, ogr
from src.core import projections
from src.core import manage_files
from src.physiographic.base import Basin
from src.core import utils as u
from .base import Meteo
from ._stations import (
    interpolation_netCDF,
    create_station_table,
    get_netCDF_grids,
    create_grid_var
)

__methods__ = [
    "Bilineal",
    "Voroni polygons",
    "Nearest neighbourgh"
]
# from ._retrieve_netcdf import ()


class StationNetCDF(Meteo):
    def __init__(self, basin_struct: Basin, ds: xr.Dataset) -> None:
        """_summary_

        Args:
            bassinVersant (Basin): _description_
            ds (xr.Dataset): _description_
        """
        super().__init__(basin_struct)
        self.ds = ds
        # Set the basin dimensions based on the information contained in the basssinVersant structure
        # Convert from km2 to m2
        CE_area = self.basin_struct.bassinVersant["superficieCE"]*1e6
        self.basin_struct.set_dimensions(np.sqrt(CE_area), np.sqrt(CE_area))
        self.lon_utm = None
        self.lat_utm = None

    @classmethod
    def charge_CORDEX_Meteo(cls, bassinVersant: Basin, vars_path: str) -> StationNetCDF:
        """_summary_

        Args:
            vars_path (str): _description_

        Returns:
            StationNetCDF: _description_
        """
        vars_dict = manage_files.dict_netCDF(vars_path)
        ds = manage_files.get_CORDEX_Dataset(vars_dict)
        # Construct object
        obj = cls(bassinVersant, ds)
        return obj

    @ classmethod
    def charge_ERA_Meteo(cls, bassinVersant: Basin, vars_path: str) -> StationNetCDF:
        """_summary_

        Args:
            vars_path (str): _description_

        Returns:
            StationNetCDF: _description_
        """
        vars_dict = manage_files.dict_netCDF(vars_path)
        ds = manage_files.get_ERA_Dataset(vars_dict)
        # Construct object
        obj = cls(bassinVersant, ds)
        return obj

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        """_summary_

        Args:
            ds (xr.Dataset): _description_

        Returns:
            xr.Dataset: _description_
        """
        # Get the var list
        var_list = list(ds.keys())
        CEs_name = os.path.join(basin_struct.project_path,"results","carreauxEntiers.csv")
        CEs_df = pd.read_csv(CEs_name,index_col=0)
        CEs_df.index = CEs_df["CEid"].values
        # Get the CEs and the i,j values from the bassinVersant structure.

        # Remove CE variable from the main dataset
        var_list.remove("CE")
        # Convert date to datenum
        datenum = np.array(list(pd.to_datetime(ds.time.values).map(
            lambda x: 366.0 + x.toordinal())), dtype=np.float32)
        # Create the dataset to store the variables in the CEQUEAU format
        for var_num, var_name in enumerate(var_list):
            if var_num == 0:
                dr = create_grid_var(ds, CEs_df, var_name, datenum)
            else:
                dr = dr.merge(create_grid_var(
                    ds, CEs_df, var_name, datenum))
            dr[var_name].attrs = ds[var_name].attrs
            # break
        return dr

    def interpolation(self, method: str,
                      name_meteo_grid_file="meteo_grid_points.shp") -> xr.DataArray:
        """_summary_

        Args:
            method (str): _description_
            name_meteo_grid_file (str): _description_

        Returns:
            xr.DataArray: _description_
        """
        # Get stations table
        self.stations_table(name_meteo_grid_file)
        # TODO: Allow to choose between different options
        dsi = interpolation_netCDF(self.ds,
                                   self.basin_struct.get_CEgrid,
                                   self.table,
                                   method)
        return dsi

    def stations_table(self, name_meteo_grid_file: str):
        """_summary_
        """
        # Open DEM to use it below
        DEM = gdal.Open(self.basin_struct.DEM, gdal.GA_ReadOnly)
        epsg_dem = projections.get_proj_code(DEM)
        # Create the CEgrid file
        ce_grid = gdal.Open(self.basin_struct.get_CEgrid, gdal.GA_ReadOnly)
        # ce_array = np.array(ce_grid.GetRasterBand(1).ReadAsArray())
        # CE_array = self.basin_struct.create_CEgrid()
        # CE_array = np.flipud(CE_array)
        # Open the shp file from the basin structure
        watershed = ogr.Open(self.basin_struct.Basin.replace(
            ".tif", ".shp"), gdal.GA_ReadOnly)
        # Get x-y pairs
        xy_pair = get_netCDF_grids(self.ds,
                                   DEM,
                                   watershed)
        # Conveert coordinates
        self.lon_utm, self.lat_utm = projections.latlon_to_utm(
            xy_pair[:, 1],  # -> Lon
            xy_pair[:, 0],  # -> Lat
            epsg_dem)

        # Create stations_table
        self.table = create_station_table(ce_grid,
                                          DEM,
                                          self.lon_utm,
                                          self.lat_utm,
                                          xy_pair)
        point = [Point(x, y)
                 for y, x in zip(self.table["lat_utm"], self.table["lon_utm"])]
        # Export stations table grid points as shp
        gpdf_table = gpd.GeoDataFrame(self.table,
                                      geometry=point,
                                      crs=epsg_dem)
        file_name = os.path.join(
            self.basin_struct.project_path, "geographic", name_meteo_grid_file)
        gpdf_table.to_file(file_name)
