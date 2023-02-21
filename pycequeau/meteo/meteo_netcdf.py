from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from osgeo import gdal, ogr
from .base import MeteoStation
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


class StationNetCDF(MeteoStation):
    def __init__(self,
                 ds: xr.Dataset) -> None:
        """_summary_

        Args:
            ds (xr.Dataset): _description_
        """
        self.ds = ds

    @classmethod
    def charge_CORDEX_Meteo(cls,vars_path: str) -> StationNetCDF:
        """_summary_

        Args:
            vars_path (str): _description_

        Returns:
            StationNetCDF: _description_
        """
        vars_dict = manage_files.dict_netCDF(vars_path)
        ds = manage_files.get_CORDEX_Dataset(vars_dict)
        # Construct object
        obj = cls(ds)
        return obj
        pass
    
    @ classmethod
    def charge_ERA_Meteo(cls, vars_path: str) -> StationNetCDF:
        """_summary_

        Args:
            vars_path (str): _description_

        Returns:
            StationNetCDF: _description_
        """
        vars_dict = manage_files.dict_netCDF(vars_path)
        ds = manage_files.get_ERA_Dataset(vars_dict)
        # Construct object
        obj = cls(ds)
        return obj

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset) -> xr.Dataset:
        """_summary_

        Args:
            ds (xr.Dataset): _description_

        Returns:
            xr.Dataset: _description_
        """
        # Get unique values for CE
        CEs = np.unique(ds.CE.values)
        CEs = CEs[~np.isnan(CEs)]
        idx = []
        idy = []
        # Create empty dataset
        var_list = list(ds.keys())
        # Get the index where the CEs fall in
        for i in CEs:
            x, y = np.where(ds.CE.values == i)
            idx.append(x[0])
            idy.append(y[0])
        # Remove CE
        float()
        var_list.remove("CE")
        # Convert date to datenum
        datenum = np.array(list(pd.to_datetime(ds.time.values).map(
            lambda x: 366.0 + x.toordinal())), dtype=np.float32)
        for var_num, var_name in enumerate(var_list):
            if var_num == 0:
                dr = create_grid_var(ds, idx, idy, CEs, var_name, datenum)
            else:
                dr = dr.merge(create_grid_var(
                    ds, idx, idy, CEs, var_name, datenum))
            dr[var_name].attrs = ds[var_name].attrs
            # break
        return dr

    def charge_Basin(self, basin: Basin):
        self.CEgrid = basin.CEgrid
        self.watershed = basin.watershed
        self.DEM = basin.DEM
        self.fishnet = basin.fishnet
        self.grid_size = basin.grid_size

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
                                   self.CEregrid,
                                   self.table,
                                   method)
        return dsi

    def stations_table(self):
        """_summary_
        """
        # Get x-y pairs
        xy_pair = get_netCDF_grids(self.ds,
                                   self.CEgrid,
                                   self.watershed)
        # Conveert coordinates
        epsg_dem = projections.get_proj_code(self.CEgrid)
        self.lat_utm, self.lon_utm = projections.latlon_to_utm(
            xy_pair[:, 0],
            xy_pair[:, 1],
            epsg_dem)
        # Get the CE_Array regridded
        self.CEregrid = u.regrid_CE(self.CEgrid,
                                    self.fishnet,
                                    self.grid_size)
        # Create stations_table
        self.table = create_station_table(self.CEregrid,
                                          self.DEM,
                                          self.lat_utm,
                                          self.lon_utm,
                                          xy_pair)
