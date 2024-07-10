import numpy as np
import pandas as pd
import rasterstats as rs
from src.core import utils as u
from src.core import projections
import geopandas as gpd
from osgeo import gdal
import os


def find_grid_coordinates(CE_array: np.ndarray) -> pd.DataFrame:
    # Flip the array, since gdal stores it upside down
    # The array was flipped because in this step we use the arrays indexes
    # rather than geographical indexing
    CE_array = np.flipud(CE_array)
    i_res = CE_array.shape[1]  # columns
    j_res = CE_array.shape[0]  # rows
    i = np.arange(0, i_res, 1, dtype=int).astype("uint16")
    j = np.arange(0, j_res, 1, dtype=int).astype("uint16")
    # j = np.flip(j)
    # Create mesh grid with the i,j values
    im, jm = np.meshgrid(i, j)
    # idx, idy = np.where(CE_array != 0)
    # Mask array based on the nondata value
    masked_CE = np.ma.masked_where(CE_array == 0, CE_array)
    jm = jm[~masked_CE.mask]
    im = im[~masked_CE.mask]
    data = np.c_[CE_array[jm, im].astype("uint16"), im+10, jm+10]
    # Create table t store the coordinates
    coordinates = pd.DataFrame(columns=["CEid", "i", "j"], data=data)
    coordinates["CEid"] = pd.to_numeric(coordinates["CEid"])
    coordinates = coordinates.sort_values(by=["CEid"])
    return coordinates


def get_lat_lon_CE(CE_fishnet_name: str) -> np.ndarray:
    gdf = gpd.read_file(CE_fishnet_name)
    centroids = gdf.centroid
    x_coords = []
    y_coords = []
    for pp in centroids.values:
        x_coords.append(pp.x)
        y_coords.append(pp.y)
    epsg_code = gdf.crs.srs
    x, y = projections.utm_to_latlon(x_coords, y_coords,
                                     epsg_code)
    array_latlon = np.array([x, y]).T
    # print(centroids)
    # Convert utm to lat - lon

    return array_latlon


def ComputeMeanSlope(gdf: gpd.GeoDataFrame,
                     Slope_name: str,
                     attr: str) -> np.ndarray:
    gdf = gdf.sort_values(by=attr)
    zonal_stats = rs.zonal_stats(gdf.geometry,
                                 Slope_name,
                                 stats=["mean"])
    stat_array = np.zeros(len(zonal_stats))
    for i in range(len(zonal_stats)):
        stat_array[i] = zonal_stats[i]["mean"]

    # gdf["meanSlope"] = 
    return stat_array
