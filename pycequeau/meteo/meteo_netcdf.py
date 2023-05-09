from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
import os
from osgeo import gdal, ogr
from .base import Meteo
from pycequeau.core import projections
from pycequeau.core import manage_files
from pycequeau.physiographic.base import Basin
from pycequeau.core import utils as u
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
        self.basin_struct.set_dimenssions(np.sqrt(CE_area), np.sqrt(CE_area))

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
        pass

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
        # Get the CEs and the i,j values from the bassinVersant structure.
        CEs = np.array(basin_struct.bassinVersant["carreauxEntiers"]["CEid"])
        # Here we substract 10 to use this vector as index in the meteo dataset
        i = np.array(
            basin_struct.bassinVersant["carreauxEntiers"]["i"]) - 10  # columns
        j = np.array(
            basin_struct.bassinVersant["carreauxEntiers"]["j"]) - 10  # rows

        # Remove CE variable from the main dataset
        var_list.remove("CE")
        # Convert date to datenum
        datenum = np.array(list(pd.to_datetime(ds.time.values).map(
            lambda x: 366.0 + x.toordinal())), dtype=np.float32)
        # Create the dataset to store the variables in the CEQUEAU format
        for var_num, var_name in enumerate(var_list):
            if var_num == 0:
                dr = create_grid_var(ds, j, i, CEs, var_name, datenum)
            else:
                dr = dr.merge(create_grid_var(
                    ds, j, i, CEs, var_name, datenum))
            dr[var_name].attrs = ds[var_name].attrs
            # break
        return dr

    def interpolation(self, method: str) -> xr.DataArray:
        """_summary_

        Args:
            method (str): _description_

        Returns:
            xr.DataArray: _description_
        """
        # Get stations table
        self.stations_table()
        # TODO: Allow to choose btween different options
        dsi = interpolation_netCDF(self.ds,
                                   self.basin_struct._CEgrid,
                                   self.table,
                                   method)
        return dsi

    def stations_table(self):
        """_summary_
        """
        # Create the CEgrid file
        CE_array = self.basin_struct.create_CEgrid()
        # Open the shp file from the basin structure
        watershed = ogr.Open(self.basin_struct._Basin, gdal.GA_ReadOnly)
        # Get x-y pairs
        xy_pair = get_netCDF_grids(self.ds,
                                   self.basin_struct._CEgrid,
                                   watershed)
        # Conveert coordinates
        epsg_dem = projections.get_proj_code(self.basin_struct._CEgrid)
        self.lat_utm, self.lon_utm = projections.latlon_to_utm(
            xy_pair[:, 0],
            xy_pair[:, 1],
            epsg_dem)
        # Open DEM to use it below
        DEM = gdal.Open(self.basin_struct._DEM, gdal.GA_ReadOnly)
        # Create stations_table
        self.table = create_station_table(self.basin_struct._CEgrid,
                                          DEM,
                                          self.lat_utm,
                                          self.lon_utm,
                                          xy_pair)
        # Export stations table grid points as csv
        # TODO: Export in the future as shp file
        file_name = os.path.join(
            self.basin_struct._project_path, "geographic", "meteo_stations.csv")
        self.table.to_csv(file_name, index=False)
