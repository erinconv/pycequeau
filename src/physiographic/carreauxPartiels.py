import numpy as np
import pandas as pd
from src.core import utils as u
import geopandas as gpd
from osgeo import gdal
import os


def get_CP_coordinates(carreuxEntiers: pd.DataFrame,
                       CPfishnet: gpd.GeoDataFrame) -> pd.DataFrame:
    # Create dataset to store the coordinates
    coordinates = pd.DataFrame(columns=["CPid", "i", "j"],
                               index=CPfishnet.index)
    coordinates["CPid"] = CPfishnet["newCPid"]
    columns_carreux = coordinates.columns.tolist()
    # Loop into the carreuxEntiers dataset
    for index, carreu in carreuxEntiers.iterrows():
        # Find the CEid in the CP fisnhet dataframe
        idx_CEs, = np.where(CPfishnet["newCEid"] == carreu["CEid"])
        # Place this values in the coordinates dataframe
        coordinates.iloc[idx_CEs, columns_carreux.index("i")] = carreu["i"]
        coordinates.iloc[idx_CEs, columns_carreux.index("j")] = carreu["j"]
    # Drop nan if it exist
    coordinates = coordinates.dropna(axis="index")
    return coordinates


def get_codes(CPfishnet: gpd.GeoDataFrame) -> pd.DataFrame:
    # Get unique CEids
    CEids = np.unique(CPfishnet["newCEid"])
    codes = np.zeros(len(CPfishnet))
    for i in range(len(CEids)):
        idx, = np.where(CPfishnet["newCEid"] == CEids[i])
        codes[idx] = list(range(65, 65+len(idx)))
    return codes


def cumulate_variables(outlet_routes: np.ndarray,
                       pctForet: np.ndarray,
                       pctLacRiviere: list,
                       pctMarais: list) -> pd.DataFrame:
    # "cumulPctSuperficieLacsAmont",
    #  "cumulPctSuperficieMaraisAmont","cumulPctSuperficieForetAmont"
    # Create the dataset to store the cumulates values
    cumulates = pd.DataFrame(columns=["cumulPctSuperficieLacsAmont",
                                      "cumulPctSuperficieMaraisAmont",
                                      "cumulPctSuperficieForetAmont"],
                             data=np.c_[pctLacRiviere,
                                        pctMarais,
                                        pctForet],
                             index=range(1, len(outlet_routes)+1))
    # Track the upstream cps
    upstreamCPs = []
    for i in range(1, len(outlet_routes)+1):
        # Create a copy of the main dataframe
        temp_df = outlet_routes.copy()
        # find the row of the current CP
        idx_row, _ = np.where(temp_df == i)
        # Mask the downstreams values
        temp_df = temp_df[idx_row, :]
        mask_df = temp_df < i
        # Drop the downstream values
        temp_df = temp_df[~mask_df]
        # Get the unique values
        temp_df = np.unique(temp_df).astype("uint16")
        # Cumulate all the variables
        cumulates.loc[i, "cumulPctSuperficieLacsAmont"] = np.sum(
            np.array(pctLacRiviere)[temp_df-1])
        cumulates.loc[i, "cumulPctSuperficieForetAmont"] = np.sum(
            pctForet[temp_df-1])
        cumulates.loc[i, "cumulPctSuperficieMaraisAmont"] = np.sum(
            np.array(pctMarais)[temp_df-1])

    return cumulates


def get_river_geometry(CPfishnet: gpd.GeoDataFrame,
                       rtable: pd.DataFrame) -> pd.DataFrame:
    rtable.index = CPfishnet.index.values
    # Create the dataframe to store the data
    river_geometry = pd.DataFrame(columns=["profondeurMin", "longueurCoursEauPrincipal",
                                           "largeurCoursEauPrincipal", "penteRiviere"],
                                  index=rtable.index)
    # calculated as function of upstream area (units = cm)
    river_geometry["profondeurMin"] = (
        0.0198*(np.power(CPfishnet["cumulArea"].values, 0.53)))*100.0
    # calculated as function of current CP (units = 1/100 km)
    river_geometry["longueurCoursEauPrincipal"] = np.power(
        CPfishnet["Area"]*1.0e-6, 0.5)*10.0
    # calculated as function of upstream area (units = 1/10m)
    river_geometry["largeurCoursEauPrincipal"] = (
        0.49*np.power(CPfishnet["cumulArea"].values, 0.6))*10.0
    # (units = 1/1000 metres/km)
    river_geometry["penteRiviere"] = 1000.0

    return river_geometry

## Get long_lat of CP for Canopy module in CEQUEAU
def get_lat_lon_CP(CP_fishnet_name: str) -> np.ndarray:
    gdf = gpd.read_file(CP_fishnet_name)
    centroids = gdf.centroid
    x_coords = []
    y_coords = []
    for pp in centroids.values:
        x_coords.append(pp.x)
        y_coords.append(pp.y)
    epsg_code = gdf.crs.srs
    x, y = projections.utm_to_latlon(x_coords, y_coords,
                                     epsg_code)
    array_latlonCP = np.array([x, y]).T
    # print(centroids)
    # Convert utm to lat - lon

    return array_latlonCP