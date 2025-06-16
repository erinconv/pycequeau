import os
import numpy as np
import pandas as pd
import rasterstats as rs
from src.core import utils as u
from src.core import projections
import geopandas as gpd
from osgeo import gdal
from shapely.geometry import Point, LineString


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


def create_cequeau_stream_network(project_path: str,
                                  CP_fishnet: gpd.GeoDataFrame,
                                  rtable: pd.DataFrame, 
                                  area_th=0.01) -> gpd.GeoDataFrame:
    """_summary_

    Args:
        area_th (float, optional): _description_. Defaults to 0.01.

    Returns:
        float: _description_
    """
    # Open the
    # CP_fishnet_name = os.path.join(project_path,
    #                                "geographic",
    #                                "CP_fishnet.shp")
    # CP_fishnet = gpd.read_file(CP_fishnet_name)
    # cp_struct_name = os.path.join(project_path,
    #                               "results",
    #                               "carreauxPartiels.csv")
    # carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
    CP_fishnet["x_c"] = CP_fishnet.centroid.x
    CP_fishnet["y_c"] = CP_fishnet.centroid.y

    # Array to store the line shape
    lines = np.zeros(len(CP_fishnet), dtype=LineString)
    # Mask small streams
    max_area = CP_fishnet["cumulArea"].max()
    idx_small_area = CP_fishnet["cumulArea"] > area_th*max_area
    CP_fishnet.index = CP_fishnet["newCPid"].values
    # carreaux_partiels.index = CP_fishnet["newCPid"].values
    df_line = pd.DataFrame(columns=["names", "cumulArea", "slope", "azimutCoursEau"],
                           data=np.empty(shape=[len(CP_fishnet), 4]),
                           index=CP_fishnet.index)
    df_line["cumulArea"] = CP_fishnet["cumulArea"].values
    # Save the outlet point of the river
    x_outlet, y_outlet = save_outlet_point(project_path, CP_fishnet)
    # Start creating the stream files
    for i, _ in CP_fishnet.iterrows():
        # cp_aval = carreaux_partiels.loc[i, "idCPAval"]
        cp_aval = rtable.loc[i, "downstreamCPs"]
        if cp_aval == 0:
            p1 = Point(CP_fishnet.loc[i, "x_c"],
                       CP_fishnet.loc[i, "y_c"])
            p2 = Point(x_outlet, y_outlet)
            lines[i-1] = LineString([p2, p1])
            # carreaux_partiels["altitudeMoy"]
            df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"
            slope = 0
            df_line.loc[i, "slope"] = slope
            # Compute river azimuth
            df_line.loc[i, "azimutCoursEau"] = compute_river_azimuth(p1, p2)
        else:
            p1 = Point(CP_fishnet.loc[i, "x_c"],
                       CP_fishnet.loc[i, "y_c"])
            p2 = Point(CP_fishnet.loc[cp_aval, "x_c"],
                       CP_fishnet.loc[cp_aval, "y_c"])
            lines[i-1] = LineString([p2, p1])
            df_line.loc[i, "names"] = f"CP{i} to CP{cp_aval}"
            dh = CP_fishnet.loc[i, "altitude"] - \
                CP_fishnet.loc[cp_aval, "altitude"]
            df_line.loc[i, "slope"] = abs(dh)/lines[i-1].length
            # Compute river azimuth
            df_line.loc[i, "azimutCoursEau"] = compute_river_azimuth(p1, p2)

    # df_line = df_line.loc[idx_small_area.values]
    gpdf_line = gpd.GeoDataFrame(df_line,
                                 geometry=lines,
                                 crs=CP_fishnet.geometry.crs)
    gpdf_line = gpdf_line.loc[idx_small_area.values]
    # Open the outlet routes
    outlet_file = os.path.join(project_path,
                               "results",
                               "outlet_routes.csv")
    outlet_routes = pd.read_csv(outlet_file, header=None)
    # find all the complete routes
    idy, = np.where(outlet_routes.iloc[:, -1] != 0)
    outlet_routes_complete = outlet_routes.iloc[idy, :]
    # Compute the length of all the routes
    lengths = []
    # Find the largest CP with stream
    idx_max = gpdf_line.index.max()
    for i in range(len(outlet_routes_complete)):
        routes = np.array(outlet_routes_complete.iloc[i, :])
        routes = routes[routes < idx_max]
        sub_gpdf = gpdf_line.loc[routes, :]
        lengths.append(sub_gpdf.geometry.length.sum())
    # Find the maximun length value
    idx_longest_path = lengths.index(max(lengths))
    main_route = np.array(outlet_routes_complete.iloc[idx_longest_path, :])
    main_route = main_route[main_route < idx_max]
    gpdf_line["main_path"] = 0
    gpdf_line.loc[main_route, "main_path"] = 1
    # Open the outlet routes
    streams_file = os.path.join(project_path,
                                "geographic",
                                "streams_cequeau.shp")
    gpdf_line.to_file(streams_file)
    # main_channel_length = max(lengths)
    # slope_channel = gpdf_line["slope"].where(
    #     gpdf_line["slope"] > 0).mean()
    return gpdf_line


def save_outlet_point(project_path: str, CP_fishnet: gpd.GeoDataFrame) -> tuple:
    FAC_path = os.path.join(project_path,
                            "geographic",
                            "FAC.tif")
    x_outlet, y_outlet = u.get_outlet_point(FAC_path)
    point_geom = np.array([Point(x_outlet, y_outlet)], dtype=Point)
    outet_df = pd.DataFrame(data=np.array([[x_outlet, y_outlet]]),
                            columns=["x", "y"])
    outlet_gdf = gpd.GeoDataFrame(outet_df,
                                  geometry=point_geom,
                                  crs=CP_fishnet.geometry.crs)
    # Save the outlet point
    point_file = os.path.join(project_path,
                              "geographic",
                              "outlet_point.shp")
    outlet_gdf.to_file(point_file)
    return x_outlet, y_outlet


def compute_river_azimuth(p1: Point, p2: Point):
    angle = np.arctan2((p2.y - p1.y), (p2.x - p1.x))
    return np.degrees(angle)

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
    array_latlon = np.array([x, y]).T
    # print(centroids)
    # Convert utm to lat - lon

    return array_latlon

def ComputeMeanCanopy(gdf: gpd.GeoDataFrame,
                      Canopy_name: str,
                      attr: str) -> np.ndarray:
    gdf = gdf.sort_values(by=attr)
    zonal_stats = rs.zonal_stats(gdf.geometry,
                                 Canopy_name,
                                 stats=["mean"])
    stat_canopy_array_cp = np.zeros(len(zonal_stats))
    for i in range(len(zonal_stats)):
        stat_canopy_array_cp[i] = zonal_stats[i]["mean"]
    # gdf["meanCanopy"]
    return  stat_canopy_array_cp
