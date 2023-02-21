from __future__ import annotations

from osgeo import gdal, ogr
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Proj
from pycequeau.physiographic.base import Basin
from pycequeau.core import utils as u
# Base class for the meteo module


class MeteoStation:
    def __init__(self) -> None:
        pass

    # def __init__(self, *args, **kwargs) -> None:
    #     if self._check_inputs(*args, **kwargs):
    #         # print("Passed")
    #         # TODO: Validation code must be addedt
    #         pass

    @classmethod
    def cequeau_grid(cls, ds: xr.DataArray) -> xr.Dataset:
        
        grid = cls._cequeau_grid(ds)
        
        return grid

    @classmethod
    def Construct(csl, *args, **kwargs):
        pass

    @classmethod
    def _check_inputs(cls, *args, **kwargs) -> bool:
        """
        The order must follow the next structure:
        arg[0]: Basin name -> str
        arg[2]: Grid_size -> int
        ---------
        kwarg[0]: DEM -> gdal.Dataset 
        kwarg[1]: Cegrid -> gdal.Dataset
        kwarg[2]: Fishnet -> ogr.DataSource
        kwarg[3]: Watershed -> ogr.DataSource
        kwarg[4]: MeteoFile -> pd.DataFrame | np.array | xr.DataArray
        ---------
        Returns:  bool
        """
        if len(args) < 2 or len(args) > 2 or len(kwargs) < 5 or len(kwargs) > 5:
            raise ValueError("Incomplete number of arguments")
        kwargs_names = list(kwargs.keys())
        # print(isinstance(args[1],int))
        if not isinstance(args[0], str):
            raise TypeError(f"The expected object in {args[0]} is a str")
        if not isinstance(args[1], int):
            raise TypeError(f"The expected object in {args[1]} is an int")
        if not isinstance(kwargs[kwargs_names[0]], gdal.Dataset):
            raise TypeError(
                f"The expected object in {kwargs_names[0]} is a gdal.Dataset")
        if not isinstance(kwargs[kwargs_names[1]], gdal.Dataset):
            raise TypeError(
                f"The expected object in {kwargs_names[1]} is an ogr.DataSource")
        if not isinstance(kwargs[kwargs_names[2]], ogr.DataSource):
            raise TypeError(
                f"The expected object in {kwargs_names[2]} is an ogr.DataSource")
        if not isinstance(kwargs[kwargs_names[3]], ogr.DataSource):
            raise TypeError(
                f"The expected object in {kwargs_names[3]} is an ogr.DataSource")
        if not any(isinstance(kwargs[kwargs_names[4]], t) for t in [pd.DataFrame, np.ndarray, xr.DataArray, xr.Dataset]):
            raise TypeError(
                f"The expected object in {kwargs_names[4]} is a pd.DataFrame | np.ndarray | xr.DataArray | xr.Dataset")

        return True

        # if approve_score == 7:
        #     return True
        # # else:
        #     raise

    @classmethod
    def stattion_table(cls,
                       CEgrid: gdal.Dataset,
                       watershed: ogr.DataSource,
                       fishnet: ogr.DataSource,
                       grid_size: int,
                       **kwargs):
        cls._stations_table(CEgrid,
                            watershed,
                            fishnet,
                            grid_size,
                            **kwargs)
        pass

    @classmethod
    def interpolate(cls,
                    ds: xr.DataArray,
                    CEgrid: gdal.Dataset,
                    stations_table: pd.DataFrame):
        cls._interpolation(CEgrid, CEgrid, stations_table)
