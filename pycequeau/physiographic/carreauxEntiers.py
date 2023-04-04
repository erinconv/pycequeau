import numpy as np
import pandas as pd
from pycequeau.core import utils as u
import geopandas as gpd
from osgeo import gdal
import os



def find_grid_coordinates(CE_array: np.ndarray) -> pd.DataFrame:
    # Flip the array, since gdal stores it upside down
    CE_array = np.flipud(CE_array)
    i_res = CE_array.shape[0]
    j_res = CE_array.shape[1]
    i = np.linspace(0, i_res-1, i_res, dtype=int).astype("uint8")
    j = np.linspace(0, j_res-1, j_res, dtype=int).astype("uint8")
    # Create mesh grid with the i,j values
    jm, im = np.meshgrid(j,i)
    # Mask array based on the nondata value
    masked_CE = np.ma.masked_where(CE_array==0,CE_array)
    jm = jm[~masked_CE.mask]
    im = im[~masked_CE.mask]
    data = np.c_[CE_array[im,jm].astype("uint16"),im+10,jm+10]
    # Create table t store the coordinates
    coordinates = pd.DataFrame(columns=["CEid","i","j"],data=data)
    coordinates["CEid"] = pd.to_numeric(coordinates["CEid"])
    coordinates = coordinates.sort_values(by=["CEid"])
    return coordinates


