import itertools
import sys
from math import ceil
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from osgeo import gdal
import geopandas as gpd
import rasterstats as rs
import warnings
# from src.core import utils as u


def convert_coords_to_index(df: gpd.GeoDataFrame,
                            dataset: gdal.Dataset) -> gpd.GeoDataFrame:
    """_summary_

    Args:
        df (gpd.GeoDataFrame): _description_
        dataset (gdal.Dataset): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    # band = dataset.GetRasterBand(1)
    # cols = dataset.RasterXSize
    # rows = dataset.RasterYSize
    df2 = pd.concat([df, pd.DataFrame(columns=["col_min", "row_min", "col_max", "row_max"],
                                      index=df.index.values)], axis=1)
    transform = dataset.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    for i in range(len(df)):
        df2.at[i, "col_min"] = ceil((df["minx"].iloc[i] - xOrigin)/pixelWidth)
        df2.at[i, "row_max"] = ceil((yOrigin - df["miny"].iloc[i])/pixelHeight)
        df2.at[i, "col_max"] = ceil((df["maxx"].iloc[i]-xOrigin)/pixelWidth)
        df2.at[i, "row_min"] = ceil((yOrigin - df["maxy"].iloc[i])/pixelHeight)
        # px_1 = int((df["minx"].iloc[i] - xOrigin)/pixelWidth)
        # py_1 = int((yOrigin - df["miny"].iloc[i])/pixelHeight)
        # px_2 = int((df["maxx"].iloc[i]-xOrigin)/pixelWidth)
        # py_2 = int((yOrigin - df["maxy"].iloc[i])/pixelHeight)

        # with rio.open(file) as dataset:
        #     py1, px1 = dataset.index(df["minx"].iloc[i], df["miny"].iloc[i])
        #     py2, px2 = dataset.index(df["maxx"].iloc[i], df["maxy"].iloc[i])
        #     a = 1
        # print(px1-px_1)
        # print(px2-px_2)
        # print(py1-py_1)
        # print(py2-py_2)
    return df2


def find_neighbors(gdf: gpd.GeoDataFrame,
                   id: str) -> gpd.GeoDataFrame:
    # https://gis.stackexchange.com/questions/281652/finding-all-neighbors-using-geopandas
    # Drop the column if it exist
    if 'NEIGHBORS' in gdf.columns:
        gdf = gdf.drop(columns=["NEIGHBORS", "KEEP"])
    # add NEIGHBORS column
    gdf = gdf.reindex(columns=gdf.columns.tolist() + ['NEIGHBORS', 'KEEP'])
    gdf["NEIGHBORS"] = ''
    gdf["KEEP"] = ''
    columns = gdf.columns.tolist()
    count = 0
    for index, CP in gdf.iterrows():
        # get 'not disjoint' countries
        neighbors = gdf[~gdf.geometry.disjoint(CP.geometry)][id].tolist()
        keep = gdf[~gdf.geometry.disjoint(CP.geometry)].Dissolve.tolist()
        # remove own name of the country from the list
        keep = [bool(numb) for numb, name in zip(
            keep, neighbors) if CP[id] != name]
        neighbors = [name for name in neighbors if CP[id] != name]
        # add names of neighbors as NEIGHBORS value
        # Catch an exception here
        try:
            gdf.at[index, 'NEIGHBORS'] = neighbors
            gdf.at[index, 'KEEP'] = keep
        except ValueError:
            if isinstance(neighbors, list):
                gdf["NEIGHBORS"].iloc[count] = neighbors
                gdf['KEEP'].iloc[count] = keep
            else:
                gdf["NEIGHBORS"].iloc[count] = list([int(neighbors)])
                gdf['KEEP'].iloc[count] = list([int(keep)])
        count += 1
    return gdf


def identify_small_CPs(CE_fishnet: gpd.GeoDataFrame,
                       CP_fishnet: gpd.GeoDataFrame,
                       thereshold: float):
    # Get the area of the CE grid
    CE_area = CE_fishnet.area[0]
    # Get area for each CP feature
    CP_fishnet = CP_fishnet.dropna(subset=['CEid'])
    CP_fishnet = CP_fishnet.explode()
    # CP_fishnet["CPid"] = range(1,len(CP_fishnet)+1)
    CP_fishnet["Area"] = CP_fishnet.area
    CP_fishnet = CP_fishnet.dropna(subset=['Area'])
    # Mask to drop values with extremly tiny areas
    mask_area = CP_fishnet["Area"] > 1.0
    CP_fishnet = CP_fishnet[mask_area]
    CP_fishnet["CPid"] = range(1, len(CP_fishnet)+1)
    mask_CP = CP_fishnet["Area"] < thereshold*CE_area
    CP_fishnet["Dissolve"] = 0
    CP_fishnet.at[mask_CP, "Dissolve"] = 1
    # Get the CP ousite the subbasin
    # This need to be changed for anohter
    mask_SUB = np.isnan(CP_fishnet["CATid"])
    index_drop = CP_fishnet.index[mask_SUB]
    # Drop this value from the main dataframe
    CP_fishnet = CP_fishnet.drop(index=index_drop)
    # CP_fishnet.at[mask_SUB, "Dissolve"] = 0
    CP_fishnet.index = CP_fishnet["CPid"].values
    # CP_fishnet["CPid"] = CP_fishnet.index
    CP_fishnet.index = CP_fishnet.index.rename("index")
    return CP_fishnet


def remove_border_CPs(CE_fishnet: gpd.GeoDataFrame,
                      CP_fishnet: gpd.GeoDataFrame,
                      FAC: str) -> list:
    """_summary_

    Args:
        CE_fishnet (gpd.GeoDataFrame): _description_
        CP_fishnet (gpd.GeoDataFrame): _description_
        FAC (str): _description_

    Returns:
        list: _description_
    """
    # Get the bounds
    # bounds = CE_fishnet["geometry"].bounds
    # CE_fishnet = pd.concat([CE_fishnet, bounds], axis=1)
    # FAC_dataset = gdal.Open(FAC)
    # CE_fishnet = convert_coords_to_index(CE_fishnet, FAC_dataset)

    CP_fishnet = pd.concat([CP_fishnet, pd.DataFrame(columns=["maxFAC"],
                                                     index=CP_fishnet.index)], axis=1)
    columnsCP = CP_fishnet.columns.tolist()
    for _, CE in CE_fishnet.iterrows():
        # Get the index for all the subbasin features
        idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
        # Get all features inside the CE
        CE_features = CP_fishnet.iloc[idx]
        # Check if the CE is empty or not
        if CE_features.empty:
            continue
        # If there is not features to dissolve, get rid off
        stats = rs.zonal_stats(CE_features, FAC, stats=['max'])
        CP_fishnet.iloc[idx, columnsCP.index("maxFAC")] = [
            s['max'] for s in stats]
        # Update values
        CE_features = CP_fishnet.iloc[idx]
        # Find neighbors
        CE_features = find_neighbors(CE_features, "CPid")
        for i, CP in CE_features.iterrows():
            # Take only the CP labeled with dissolve
            if CP["Dissolve"]:
                # Here I delete the Isolated CP in the border
                if CP["maxFAC"] is None:
                    CP_fishnet.at[i, "CPid"] = 0.0
                    CP_fishnet.at[i, "Dissolve"] = 0
                    CP_fishnet.at[i, "CEid"] = 0
                # Delete the CP less than 400 m2. DEM 20x20
                if CP["maxFAC"] is None and CP["Area"] <= 400:
                    CP_fishnet.at[i, "CPid"] = 0.0
                    CP_fishnet.at[i, "Dissolve"] = 0
                    CP_fishnet.at[i, "CEid"] = 0
                # Check if there are no neighbors
                if not CP["NEIGHBORS"]:
                    CP_fishnet.at[i, "CPid"] = 0.0
                    CP_fishnet.at[i, "maxFAC"] = None
                    CP_fishnet.at[i, "Dissolve"] = 0
    # Drop the dissolved values in this iterations
    CP_fishnet = CP_fishnet[CP_fishnet["CPid"] != 0]
    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    # Save file
    # CP_fishnet.index = range(len(CP_fishnet))
    CP_fishnet.index = CP_fishnet.index.rename("ind")
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    CP_fishnet.at[:, "Area"] = CP_fishnet.area

    return CP_fishnet, CE_fishnet


def remove_smallCP(CE_fishnet: gpd.GeoDataFrame,
                   CP_fishnet: gpd.GeoDataFrame,
                   thereshold: float) -> gpd.GeoDataFrame:
    """_summary_

    Args:
        CE_fishnet (gpd.GeoDataFrame): _description_
        CP_fishnet (gpd.GeoDataFrame): _description_
        thereshold (float): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    CP_fishnet.index = CP_fishnet["CPid"].values
    CE_area = CE_fishnet.area[0]
    mask_CP = CP_fishnet["Area"] < thereshold*CE_area
    CP_fishnet["Dissolve"] = 0
    CP_fishnet.loc[mask_CP, "Dissolve"] = 1
    count = 0
    while CP_fishnet["Dissolve"].max() == 1:
        # dissolve_vals = CP_fishnet["Dissolve"].values
        for _, CE in CE_fishnet.iterrows():
            # Track the already processed CPs
            tracked_cps = []
            # Get the index for all the subbasin features
            idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
            # Get all features inside the CE
            CE_features = CP_fishnet.iloc[idx]
            # Check if there are not CPs to dissolve to avoid computing the
            # neighbors and save time
            if CE_features["Dissolve"].max() == 0:
                continue
            CE_features = find_neighbors(CE_features, "CPid")
            for i, CP in CE_features.iterrows():
                # Check if the CP was already processed
                if CP["CPid"] in tracked_cps:
                    continue
                # Check if need to dissolve
                if CP["Dissolve"]:
                    # Get the neighbors and the condition
                    neighbors = CE_features.loc[i, "NEIGHBORS"]
                    dissolve = CE_features.loc[i, "KEEP"]
                    # Check if there are neighbors to dissolve also
                    if True in dissolve:
                        # Filter by dissolve or not dissolve
                        to_dissolve = list(
                            itertools.compress(neighbors, dissolve))
                        not_dissolve = list(itertools.compress(
                            neighbors, np.logical_not(dissolve)))
                        # 1- There are only CPs to dissolve.
                        if not not_dissolve:
                            to_dissolve.append(i)
                            idx_max = CE_features.loc[to_dissolve,
                                                      "maxFAC"].idxmax()
                            CP_fishnet.at[to_dissolve, "CPid"] = idx_max
                            # Track this CPs
                            tracked_cps.extend(to_dissolve)
                        # 2- There are both, neighbors to dissolve and not to
                        elif to_dissolve and not_dissolve:
                            to_dissolve.append(i)
                            # idx_max = CE_features.loc[not_dissolve, "maxFAC"].idxmax()
                            idx_max = CE_features.loc[to_dissolve, "maxFAC"].idxmax(
                            )
                            CP_fishnet.at[to_dissolve, "CPid"] = idx_max
                            # Track this CP
                            tracked_cps.extend(to_dissolve)
                    else:
                        if not CP["NEIGHBORS"]:
                            # This means the CPs does not have neighbors
                            CP_fishnet.at[i, "CPid"] = 0
                            continue
                        # Here the CPs with no neighbors to dissolve are processed
                        idx_max = CE_features.loc[CP["NEIGHBORS"], "maxFAC"].idxmax(
                        )
                        CP_fishnet.at[i, "CPid"] = idx_max

        CP_fishnet = CP_fishnet[CP_fishnet["CPid"] != 0]
        CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
        # Re-do the dissolve labeling
        CP_fishnet["Area"] = CP_fishnet.area
        mask_CP = CP_fishnet["Area"] < thereshold*CE_area
        CP_fishnet["Dissolve"] = 0
        CP_fishnet.loc[mask_CP, "Dissolve"] = 1
        # Save file
        CP_fishnet.index = CP_fishnet.index.rename("index")
        CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
        count += 1
        # Break the loop if a maximun number of iterations was exeeced
        if count > 15:
            break
    return CP_fishnet


def dissolve_pixels(CE_fishnet: gpd.GeoDataFrame,
                    CP_fishnet: gpd.GeoDataFrame,
                    pixel_size_area) -> gpd.GeoDataFrame:
    """
    The ressult from the previous process might lead to have multipolygon features
    We need to  make sure this is going to be well dissolved.
    Also, all the pixel-size features need to be merged with the bigger neighbors
    in order to avoid errors while computing zonal area in further steps.

    Args:
        CE_fishnet (gpd.GeoDataFrame): _description_
        CP_fishnet (gpd.GeoDataFrame): _description_
        area_th (_type_): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    #
    CP_fishnet = CP_fishnet.explode()
    # mask_pixel_feats = CP_fishnet["Area"] > 2*abs(pixel_size_area)
    CP_fishnet["Area"] = CP_fishnet.area
    # Drop values less than three pixels size
    mask_pixels = CP_fishnet["Area"] <= 3*abs(pixel_size_area)
    CP_fishnet["Dissolve"] = 0
    CP_fishnet.loc[mask_pixels, "Dissolve"] = 1
    CP_fishnet.loc[:, "CPid"] = range(1, len(CP_fishnet)+1)
    CP_fishnet.index = range(1, len(CP_fishnet)+1)
    for _, CE in CE_fishnet.iterrows():
        # Get the index for all the subbasin features
        idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
        # Get all features inside the CE
        CE_features = CP_fishnet.iloc[idx]
        # Check if the CE is empty
        if CE_features.empty or CE_features["Dissolve"].max() == 0:
            continue
        tracked_cps = []
        CE_features = find_neighbors(CE_features, "CPid")
        for i, CP in CE_features.iterrows():
            # Track the already processed CPs
            # Check if the CP was already processed
            if CP["CPid"] in tracked_cps:
                continue
            # Check if need to dissolve
            if CP["Dissolve"]:
                # Get the neighbors and the condition
                neighbors = CE_features.loc[i, "NEIGHBORS"]
                idx_max_area = CE_features.loc[neighbors, "maxFAC"].idxmax()
                CP_fishnet.loc[CP["CPid"], "CPid"] = idx_max_area
                tracked_cps.append(CP["CPid"])
            continue

    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    # Re-do the dissolve labeling
    # Save file
    CP_fishnet.index = CP_fishnet.index.rename("index")
    CP_fishnet.index = range(1, len(CP_fishnet)+1)
    CP_fishnet["CPid"] = CP_fishnet.index.values
    CP_fishnet["Area"] = CP_fishnet.area

    # Still find some pixels after dissolving them
    mask_pixels = CP_fishnet["Area"] <= 3*abs(pixel_size_area)
    CP_fishnet["Dissolve"] = 0
    CP_fishnet.loc[mask_pixels, "Dissolve"] = 1
    df = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    CP_fishnet = CP_fishnet[~mask_pixels]
    df = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    CP_fishnet.loc[:, "CPid"] = range(1, len(CP_fishnet)+1)
    CP_fishnet.index = range(1, len(CP_fishnet)+1)
    # df = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    return CP_fishnet


def force_4CP(CE_fishnet: gpd.GeoDataFrame,
              CP_fishnet: gpd.GeoDataFrame,
              area_th: float) -> gpd.GeoDataFrame:
    # Explode the fishnet
    CP_fishnet = CP_fishnet.explode()
    CP_fishnet.index = CP_fishnet["CPid"].values
    CE_area = CE_fishnet.area.max()
    mask_CP = CP_fishnet["Area"] < area_th*CE_area
    CP_fishnet.loc[mask_CP, "Dissolve"] = 1
    CP_fishnet.index = range(1, len(CP_fishnet)+1)
    CP_fishnet["CPid"] = CP_fishnet.index.values
    for _, CE in CE_fishnet.iterrows():

        # Get the index for all the subbasin features
        idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
        # Get all features inside the CE
        CE_features = CP_fishnet.iloc[idx]
        # Check if the CE is empty
        if CE_features.empty:
            continue
        # Get the CPid
        while len(CE_features) > 4:
            # Find neighbors
            CE_features = find_neighbors(CE_features, "CPid")
            # CE_features.loc[:, "Area"] = CE_features.area
            columns = CE_features.columns.tolist()
            # Find the smallest CP
            idx_small, = np.where(
                CE_features.loc[:, "Area"].values == np.amin(CE_features["Area"]))
            # Find the neighbor with the highest flow accumulation
            neighbors = CE_features.iloc[idx_small,
                                         columns.index("NEIGHBORS")].values[0]
            # Check that there are no more neighbours
            # if not neighbors:
            #     continue
            maxFAC = np.amax(CE_features.loc[neighbors, "maxFAC"])
            idx_maxFAC, = np.where(
                CE_features.loc[neighbors, "maxFAC"] == maxFAC)
            # Get the cpid of the replaced value
            idx_replaced = CE_features.iloc[idx_small, columns.index("CPid")]
            # Replace the value in the main dataframe
            CP_fishnet.loc[idx_replaced, "CPid"] = neighbors[idx_maxFAC[0]]
            # Now dissolve all CE_features and rearange the things
            CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
            CP_fishnet.index = CP_fishnet.index.rename("index")
            CP_fishnet.index = range(1, len(CP_fishnet)+1)
            CP_fishnet["CPid"] = CP_fishnet.index.values
            CP_fishnet["Area"] = CP_fishnet.area
            # Replace the CPid value into the small CP
            idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
            CE_features = CP_fishnet.iloc[idx]
            # df = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
            # CE_features.iloc[idx_small, columns.index(
            #     "CPid")] = neighbors[idx_maxFAC[0]]

    # Save file
    # CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    # CP_fishnet.index = CP_fishnet.index.rename("index")
    # CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    # CP_fishnet.at[:, "Area"] = CP_fishnet.area

    # Given that this is the final dissolving step, here I will check for
    # duplicate geometries and drop them all if that's the case
    CP_fishnet['geometry_str'] = CP_fishnet['geometry'].apply(str)
    CP_fishnet = CP_fishnet.explode()
    CP_fishnet['normalized_geometry'] = CP_fishnet['geometry'].apply(
        lambda geom: str(Polygon(geom.exterior.coords)))
    CP_fishnet = CP_fishnet.dissolve(by='normalized_geometry')
    CP_fishnet.index = CP_fishnet.index.rename("index")
    CP_fishnet.index = range(1, len(CP_fishnet)+1)
    CP_fishnet["CPid"] = CP_fishnet.index.values
    CP_fishnet = CP_fishnet.drop(columns=["geometry_str"])
    # df1 = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    # indexes_to_skip,processed_indexes = u.drop_duplicated_geometries(CP_fishnet["geometry"])
    return CP_fishnet


# def compute_flow_path(flow_dir: np.ndarray,
#                       flow_acc: np.ndarray,
#                       flow_th: float) -> np.ndarray:
#     # Mask the flow accumulation based on the threshold
#     flow_accu_mask = (flow_acc > flow_th)
#     rows, cols = np.indices(flow_dir.shape)
#     stream_network = np.zeros_like(flow_dir)
#     counter = 1
#     for row in range(flow_dir.shape[0]):
#         for col in range(flow_dir.shape[1]):
#             if flow_accu_mask[row, col]:
#                 next_row = row
#                 next_col = col
#                 # counter = 1
#                 while True:
#                     direction = flow_accu_mask[next_row, next_col]
#                     if direction == 0:
#                         break
#                     next_row += [-1, -1, 0, 1, 1, 1, 0, -1][direction - 1]
#                     next_col += [0, 1, 1, 1, 0, -1, -1, -1][direction - 1]
#                     if stream_network[next_row, next_col] > 0:
#                         break
#                 stream_network[row, col] = counter
#         # counter +=1
#     return stream_network

# Compute the mean altitude within each CE and CP
def mean_altitudes(CE_fishnet: gpd.GeoDataFrame,
                   CP_fishnet: gpd.GeoDataFrame,
                   DEM: str):
    # Add altitude column to each dataset
    CE_fishnet = CE_fishnet.reindex(
        columns=CE_fishnet.columns.tolist() + ['altitude'])
    CE_fishnet["altitude"] = None
    CP_fishnet = CP_fishnet.reindex(
        columns=CP_fishnet.columns.tolist() + ['altitude'])
    CP_fishnet["altitude"] = None

    # Compute the zonal statistics
    stats_CE = rs.zonal_stats(CE_fishnet, DEM, stats=['mean'])
    CE_fishnet.loc[:, "altitude"] = [s['mean'] for s in stats_CE]
    stats_CP = rs.zonal_stats(CP_fishnet, DEM, stats=['mean'])
    CP_fishnet.loc[:, "altitude"] = [s['mean'] for s in stats_CP]
    return CP_fishnet, CE_fishnet

# def main_path(flow_dir: np.ndarray,
#               flow_acc: np.ndarray,
#               flow_th: float) -> np.ndarray:
#     # Create a mask to extract only the cells with flow accumulation greater than a certain threshold
#     flow_accumulation_mask = (flow_acc > flow_th)

#     # Find the outlet cell of the catchment (i.e., the cell with the minimum flow accumulation)
#     outlet_row, outlet_col = np.unravel_index(
#         np.argmin(flow_acc), flow_acc.shape)

#     # Create an empty list to store the cells in the main stream
#     main_stream_cells = []

#     # Trace the main stream from the outlet cell to the start of the main stream
#     next_row, next_col = outlet_row, outlet_col
#     while True:
#         main_stream_cells.append((next_row, next_col))
#         direction = flow_dir[next_row, next_col]
#         if direction == 0:  # Check if the current cell is a sink cell and exit the loop if it is
#             break
#         # Calculate the row and column indices of the next cell in the main stream
#         next_row += [-1, -1, 0, 1, 1, 1, 0, -1][direction - 1]
#         next_col += [0, 1, 1, 1, 0, -1, -1, -1][direction - 1]
#         # Check if the next row or column index is out of bounds and exit the loop if it is
#         if next_row < 0 or next_row >= flow_dir.shape[0] or next_col < 0 or next_col >= flow_dir.shape[1]:
#             break
#         # Check if the current cell has more than one upstream cell and exit the loop if it does
#         upstream_count = 0
#         for i, j in [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]:
#             row = next_row + i
#             col = next_col + j
#             if row >= 0 and row < flow_dir.shape[0] and col >= 0 and col < flow_dir.shape[1]:
#                 if flow_dir[row, col] == (8 - direction):
#                     upstream_count += 1
#                     if upstream_count > 1:
#                         break
#         if upstream_count > 1:
#             break

#     # Create a new raster to represent the main stream cells
#     main_stream_raster = np.zeros_like(flow_dir)
#     for cell in main_stream_cells:
#         main_stream_raster[cell[0], cell[1]] = 1


def routing_table(CP_fishnet: gpd.GeoDataFrame,
                  CE_fishnet: gpd.GeoDataFrame,
                  FAC: str,
                  CP_array: np.ndarray,
                  CE_array: np.ndarray,
                  ncols: float,
                  nrows: float) -> pd.DataFrame:
    """_summary_

    Args:
        CP_fishnet (gpd.GeoDataFrame): _description_
        CE_fishnet (gpd.GeoDataFrame): _description_
        FAC (str): _description_
        CP_array (np.ndarray): _description_
        CE_array (np.ndarray): _description_
        ncols (float): _description_
        nrows (float): _description_

    Returns:
        pd.DataFrame: _description_
    """
    # Renaming index to make sure it starts from 1
    CP_fishnet["CPid"] = range(1, len(CP_fishnet)+1)
    CP_fishnet = CP_fishnet.sort_values(by=["CPid"])
    CP_fishnet.index = CP_fishnet["CPid"].values
    # Create dataframe to store the routing data
    routing = pd.DataFrame(columns=["CPid", "inCPid", "oldCP", "FAC"],
                           index=CP_fishnet.index)
    # Assign values to the routing dataframe
    # Compute the boundaries for each square within the CE and CP fishnets
    CE_fishnet = pd.concat([CE_fishnet,
                            CE_fishnet["geometry"].bounds], axis=1)
    routing["CPid"] = CP_fishnet["CPid"].values
    # Get the FAC array
    FAC_dataset = gdal.Open(FAC, gdal.GA_ReadOnly)
    band = FAC_dataset.GetRasterBand(1)
    FAC_array = band.ReadAsArray()
    FAC_array[FAC_array < 0] = 0
    # Get the CE square indexes
    CE_fishnet = convert_coords_to_index(CE_fishnet, FAC_dataset)
    CE_columns = CE_fishnet.columns.tolist()
    # df1 = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    # df2 = pd.DataFrame(CE_fishnet.drop(columns='geometry'))
    # Loop into each CE

    for index, feat in CP_fishnet.iterrows():
        # Find the rows and cols where the CP value is stored
        # Find the index which correspond to the CEid in the main dataframe
        idx_CE, = np.where(CE_fishnet["CEid"] == feat["CEid"])
        # Maybe the CE does not exist?
        if len(idx_CE) == 0:
            continue
        # Slice the CE array using the corners into the main CE dataset
        row_min = CE_fishnet.iloc[idx_CE[0], CE_columns.index("row_min")]
        row_max = CE_fishnet.iloc[idx_CE[0], CE_columns.index("row_max")]
        col_min = CE_fishnet.iloc[idx_CE[0], CE_columns.index("col_min")]
        col_max = CE_fishnet.iloc[idx_CE[0], CE_columns.index("col_max")]
        CE = CE_array[row_min:row_max,
                      col_min:col_max]
        # Check the size of the sliced subset.
        # Sometimes the sliced subset is smaller than the actual subset size
        # because the indexes are derived from float values and some rounding
        # artifacts can affect the precision of the index
        if CE.shape[0] != nrows:
            # row_min = row_min
            row_max = row_min + nrows
            CE = CE_array[row_min:row_max,
                          col_min:col_max]
        elif CE.shape[1] != ncols:
            # col_min = col_min
            col_max = col_min + ncols
            CE = CE_array[row_min:row_max,
                          col_min:col_max]

        # Get the unique rows and cols that exist in the sliced array
        unique_rows, _ = np.unique(CE, return_counts=True, axis=0)
        if unique_rows.shape[0] > 1:
            diff = np.diff(CE, axis=0)
            idx_diff_zero, _ = np.where(diff != 0)
            # Compare the value with respect to the number of cols
            if idx_diff_zero[0] > 0:
                # if counts_rows[0] < counts_rows[1]:
                val = -1
            else:
                val = 1
            row_min = row_min + val
            row_max = row_min + nrows
            # row_max = row_min-1+nrows
            CE = CE_array[row_min:row_max,
                          col_min:col_max]
        unique_cols, _ = np.unique(CE, return_counts=True, axis=1)
        if unique_cols.shape[1] > 1:
            # Find the differences along the column
            diff = np.diff(CE, axis=1)
            _, idy_diff_zero = np.where(diff != 0)
            # Compare the value with respect to the number of cols
            if idy_diff_zero[0] > 0:
                # if counts_cols[0] > counts_cols[1]:
                val = -1
            else:
                val = 1
            col_min = col_min + val
            col_max = col_min + ncols
            # col_max = col_min-1+ncols
            CE = CE_array[row_min:row_max,
                          col_min:col_max]
        # Find the CP values whitin the CE subset
        CE = CE_array[row_min-1:row_max+1,
                      col_min-1:col_max+1]
        CP = CP_array[row_min-1:row_max+1,
                      col_min-1:col_max+1]
        subFAC = FAC_array[row_min-1:row_max+1,
                           col_min-1:col_max+1]
        # Mask the FAC based on this CP
        mask_CP = (CP == feat['CPid']).astype('uint8')
        # Apply mask
        # Get the outlet indexes for the CP
        outlet_row, outlet_col = np.unravel_index(
            np.argmax(subFAC*mask_CP), subFAC.shape)
        if isinstance(outlet_row, np.ndarray):
            sys.exit("There is more than one outlet point")

        # Create mask for subFAC
        mask_FAC = np.zeros(subFAC.shape).astype("uint8")
        mask_FAC[outlet_row - 1:outlet_row + 2,
                 outlet_col - 1:outlet_col + 2] = 1
        # Get the FAC values on the mask
        # Find location of next pixel into which outlet flows
        inlet_row, inlet_col = np.unravel_index(
            np.argmax(subFAC*mask_FAC), subFAC.shape)
        # Add the CP where it discharges
        routing.loc[index, "inCPid"] = CP[inlet_row, inlet_col]
        routing.loc[index, "FAC"] = subFAC[outlet_row, outlet_col]
        if routing.loc[index, "CPid"] == routing.loc[index, "inCPid"]:
            outlet_point = routing.loc[index, "CPid"]
        if routing.loc[index, "inCPid"] == 0:
            continue
    return routing


def get_rtable(CP_fishnet: gpd.GeoDataFrame,
               routing: pd.DataFrame):
    """_summary_

    Args:
        CP_fishnet (gpd.GeoDataFrame): _description_
        routing (pd.DataFrame): _description_

    Returns:
        _type_: _description_
    """
    # *There are probably cases where a given CP drains into a non existence CP.
    # *So, here we make sure that we drop all the CP where this happens into the main data frame.
    # *This is because the CPs in the border can be so tiny that they do not account for the
    # *area threshold that we defined.
    # Get the CPs where the its downstream is zero
    idx_zero_inCP = routing["inCPid"] == 0
    index_zero_flow_in = routing.index[idx_zero_inCP]
    # Create the rouring table here
    rtable = pd.DataFrame(columns=["oldCPid", "newCPid",
                                   "upstreamCPs", "oldupstreams"],
                          index=routing.index)
    routing["diff"] = routing["CPid"] - routing["inCPid"]
    # Find where the difference is zero
    # This will help us to find the outlet CP
    idx_outlet = routing.index[routing["diff"] == 0].tolist()
    if len(idx_outlet) > 1:
        assert False, "There are more than 1 outlet point."
    elif len(idx_outlet) == 0:
        assert False, "The outlet point was not found"
    # assert len(idx_outlet) == 1, "There are more than 1 outlet point."
    # idx_outlet = routing.mask(routing["diff"] != 0)
    # Add this value in the rtable as the first value
    rtable.at[1, "oldCPid"] = routing.loc[idx_outlet, "CPid"].values[0]
    rtable.at[1, "newCPid"] = 1
    # Get column list
    columns = rtable.columns.tolist()
    # Set to zero the already tracked CP
    routing.loc[idx_outlet, "inCPid"] = -99999
    new_id_counter = 1
    for i, _ in routing.iterrows():
        # Find the upstream CPs
        index_routing = rtable.loc[i, "oldCPid"] == routing["inCPid"]
        idx_outlet = routing.index[index_routing]
        if not idx_outlet.empty:
            # Find the upstream CPs
            upstreams_cps = routing.loc[idx_outlet, "CPid"].values
            # Rename the upstream CPs
            upstreams_cps_newid = list(
                range(new_id_counter+1, new_id_counter+len(idx_outlet)+1))
            # Append the CPs vertically. +1 to step out the current index
            rtable.iloc[new_id_counter:new_id_counter +
                        len(idx_outlet), columns.index("newCPid")] = upstreams_cps_newid
            rtable.iloc[new_id_counter:new_id_counter +
                        len(idx_outlet), columns.index("oldCPid")] = upstreams_cps
            rtable.iloc[i-1,
                        columns.index("upstreamCPs")] = upstreams_cps_newid
            rtable.iloc[i-1, columns.index("oldupstreams")] = upstreams_cps
            new_id_counter += len(idx_outlet)
            routing.loc[idx_outlet, "inCPid"] = -99999
    # Check if there is NAN values in the rtable
    idx_nan_inCP = rtable["oldCPid"].isnull()
    index_nan = rtable.index[idx_nan_inCP]
    if len(index_nan) > 0:
        # warnings.warn(
        #     f"The CPs {index_nan} drop water into non existent CP \n")
        # warnings.warn(
        #     "This might no be a problem. But it may be a source of error in the next steps \n")
        # assert False, "Error found in the definition of routing function"
        # Fill this nan values in the newCPid with a flag to keep tracking it in the next steps
        # CP_fishnet.loc[index_zero_flow_in, "CPid"] = -9999
        # rtable.loc[index_nan, "newCPid"] = -9999
        rtable.loc[index_nan, "newCPid"] = index_nan
        rtable.loc[index_nan, "oldCPid"] = index_zero_flow_in
    rtable['newCPid'] = rtable['newCPid'].astype('int')
    rtable['oldCPid'] = rtable['oldCPid'].astype('int')
    # idx = np.where()
    # rtable.index = CP_fishnet.index
    # CP_fishnet = CP_fishnet.join(rtable[["oldCPid"]])
    # CP_fishnet = CP_fishnet.dropna(subset=['oldCPid'])

    # CP_fishnet = CP_fishnet.sort_values(by=["CPid"])
    # rtable = rtable.dropna(subset=['oldCPid'])
    # idx_nan_inCP = rtable["oldCPid"].isnull()
    # index_drop = rtable.index[idx_nan_inCP]
    # idx_fishnet_drop, = np.where(CP_fishnet["CPid"] == index_drop)
    # index_fishnet_drop = CP_fishnet.index[idx_fishnet_drop]
    # CP_fishnet = CP_fishnet.drop(index=rtable[""])
    # rtable.index = range(1, len(rtable)+1)
    # CP_fishnet.index = rtable.index
    # CP_fishnet["CPid"] = CP_fishnet.index
    # CP_fishnet = CP_fishnet.drop(columns=["minx", "miny",
    #                                       "maxx", "maxy",
    #                                       "col_min", "row_min",
    #                                       "col_max", "row_max"])
    # rtable = rtable.sort_values(by=["newCPid"])
    return rtable, CP_fishnet


def renumber_fishnets(CP_fishnet: gpd.GeoDataFrame,
                      CE_fishnet: gpd.GeoDataFrame,
                      rtable: pd.DataFrame,
                      routing: pd.DataFrame) -> gpd.GeoDataFrame:
    # Renumbering the CPs
    CP_fishnet["newCPid"] = 0
    CP_fishnet["newCEid"] = 0
    CE_fishnet["newCEid"] = 0
    CE_track_list = []
    CP_fishnet = CP_fishnet.sort_values(by=["CPid"])
    CP_fishnet.index = CP_fishnet["CPid"].values
    rtable = rtable.sort_values(by=["oldCPid"])
    df = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    rtable.index = rtable["oldCPid"].values
    # CP_fishnet["CPid"] = range(1, len(CP_fishnet)+1)
    # CP_fishnet["newCPid"] = pd.to_numeric(,dtype_backend="Int64")
    # CP_fishnet["newCPid"] = rtable["newCPid"].values
    # get the columns to use them as indexers
    # columns_rtable = rtable.columns.tolist()
    columns_CPfishnet = CP_fishnet.columns.tolist()
    columns_CEfishnet = CE_fishnet.columns.tolist()
    idx_zero_inCP = routing["inCPid"] == 0
    index_zero_flow_in = routing.index[idx_zero_inCP]
    # Iterate over the dataframe to rename CP fishnet
    for i, _ in CP_fishnet.iterrows():
        # Check if the current CP drops into a non existent
        # if i in index_zero_flow_in:
        #     CP_fishnet.loc[idx_old, "newCPid"] = 0
        #     continue
        # Find the index of the old CPid in the main dataframe
        index_routing = CP_fishnet["CPid"] == rtable.loc[i, "oldCPid"]
        idx_old = CP_fishnet.index[index_routing]
        # idx_old, = np.where(
        #     CP_fishnet["CPid"].values == rtable.loc[i, "oldCPid"])
        # Replace the value with the new CPid value
        CP_fishnet.loc[idx_old, "newCPid"] = rtable.loc[i, "newCPid"]
    # Change the index in the main dataframe
    CP_fishnet.index = CP_fishnet["newCPid"].values
    # df1 = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    # Sort values
    CP_fishnet = CP_fishnet.sort_values(by=["newCPid"])
    CP_fishnet["newCEid"] = np.array(CP_fishnet["newCEid"]).astype(int)
    # Renumbering the CEs. This is possible since the values are
    # already sorted in the main dataframe based on the new CPids
    CE_id = 0
    for i, CP in CP_fishnet.iterrows():
        # Check if the CE has already been tracked
        if CP["CEid"] in CE_track_list:
            continue
        # Track the CE in the CP fishnet
        CE_track_list.append(CP["CEid"])
        # Find the value in the dataset
        idx_new, = np.where(CE_fishnet["CEid"] == CE_track_list[CE_id])
        # Add the value in the column
        CE_fishnet.iloc[idx_new, columns_CEfishnet.index("newCEid")] = CE_id+1
        # Find the value in the CPfishnet dataset
        idx_new2, = np.where(CP_fishnet["CEid"] == CE_track_list[CE_id])
        CP_fishnet.iloc[idx_new2, columns_CPfishnet.index("newCEid")] = CE_id+1
        CE_id += 1
    CE_fishnet = CE_fishnet.sort_values(by=["newCEid"])
    CE_fishnet.index = CE_fishnet["newCEid"].values
    # Find where the CE is zero
    idx_drop_CE, = np.where(CE_fishnet["newCEid"] == 0)
    # Find index
    indexes_drop = CE_fishnet.index[idx_drop_CE]
    # Drop by  index
    CE_fishnet = CE_fishnet.drop(index=indexes_drop)
    # # Re-do the indexing to fix it
    # CP_fishnet = CP_fishnet.sort_values(by=["newCPid"])
    # CP_fishnet["newCPid"] = range(1, len(CP_fishnet)+1)
    # CP_fishnet.index = range(1, len(CP_fishnet)+1)
    # df1 = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    return CP_fishnet, CE_fishnet


def get_downstream_CP(rtable: pd.DataFrame) -> pd.DataFrame:
    # Create an empty table to store the values
    rtable.index = rtable["newCPid"].values
    downstreamCPs = pd.DataFrame(columns=["downstreamCPs"],
                                 index=rtable.index)
    # Convert lists into an array.
    # Get the lenght of each list to get the maximum value
    lists_len = [len(i)
                 for i in rtable["upstreamCPs"].values if isinstance(i, list)]
    # Create a zero array to store the values
    list_upstream_array = np.zeros(
        [len(rtable), max(lists_len)]).astype("uint16")
    # Fill up the array using the lists in the main dataframe
    list_to_upstream_append = []
    for i, _ in rtable.iterrows():
        if isinstance(rtable.loc[i, "upstreamCPs"], list):
            list_up = rtable.loc[i, "upstreamCPs"]
            list_upstream_array[i-1, 0:len(list_up)] = list_up
            list_to_upstream_append.append(
                list_upstream_array[i-1, :].tolist())
        else:
            list_to_upstream_append.append(
                list_upstream_array[i-1, :].tolist())
            # rtable.loc[i, "upstreamCPs"] = list_upstream_array[i-1, :].tolist()
    rtable["upstreamCPs"] = list_to_upstream_append
    # Now identify the Downstream CPs
    for i, table in rtable.iterrows():
        # Find the index of the CP in wwhich the current CP drains
        idx_down, _ = np.where(list_upstream_array == table["newCPid"])
        # Check if the list is empty
        if len(idx_down) == 0:
            downstreamCPs.loc[i, "downstreamCPs"] = 0
        else:
            downstreamCPs.loc[i,
                              "downstreamCPs"] = rtable.loc[idx_down[0]+1, "newCPid"]
    # Concat to the rtable
    rtable = pd.concat([rtable, downstreamCPs], axis=1)
    return rtable


def outlet_routes(rtable: pd.DataFrame) -> pd.DataFrame:

    # Create the array to store the lists
    allroute_lists = pd.DataFrame(columns=["outletRoutes"],
                                  index=rtable.index)
    # Start looping the array
    columns_rtable = rtable.columns.tolist()
    for i, _ in rtable.iterrows():
        route_list = []
        count_ups = 0
        # Set the upstream CP to
        up = rtable.loc[i, "newCPid"]
        # Set the number of the upstream CPs to 1
        # nu = up
        while up > 1:
            # append Nth CP to list
            route_list.append(up)
            up = rtable.loc[up, "downstreamCPs"]
            count_ups += 1
            route_list = np.unique(route_list)
            route_list = route_list.tolist()
            if count_ups > 10*len(rtable):
                # TODO: Check what would rise this issue. Although this does not affect the correct functioning of the process
                # raise ValueError
                break
            # go down from Nth position until end
            # up = rtable.iloc[up-1, "downstreamCPs"]
            # up = rtable.iloc[up-1, 4]
            # nu -= 1
        route_list.append(1)
        allroute_lists.loc[i, "outletRoutes"] = route_list
    # Convert it into a numpy array
    lists_len = [
        len(i) for i in allroute_lists["outletRoutes"].values if isinstance(i, list)]
    # Create a zero array to store the values
    allroute_lists_array = np.zeros([len(allroute_lists), max(lists_len)])
    # Fill up the array using the lists in the main dataframe
    for i in range(len(allroute_lists)):
        if isinstance(allroute_lists.loc[i+1, "outletRoutes"], list):
            list_up = allroute_lists.loc[i+1, "outletRoutes"]
            allroute_lists_array[i, 0:len(list_up)] = list_up
    return allroute_lists_array


def cumulative_areas(CP_fishnet: gpd.GeoDataFrame,
                     CE_fishnet: gpd.GeoDataFrame,
                     out_routes: np.ndarray) -> pd.DataFrame:
    """_summary_

    Args:
        CP_fishnet (gpd.GeoDataFrame): _description_
        CE_fishnet (gpd.GeoDataFrame): _description_
        out_routes (np.ndarray): _description_

    Returns:
        pd.DataFrame: _description_
    """
    # Update areas of the CP and get the CE area
    CE_area = CE_fishnet.area[1]
    CP_fishnet.indeẋ = range(1, len(CP_fishnet)+1)
    df1 = pd.DataFrame(CP_fishnet.drop(columns='geometry'))
    CP_fishnet["Area"] = CP_fishnet.area
    # Get the percentage of that area
    CP_fishnet["pctSurface"] = (CP_fishnet["Area"]/CE_area)*100
    # Cumulative areas
    CP_fishnet["cumulPctSurf"] = 0.0
    upstreamCPs = []
    for i in range(len(out_routes)):
        if i == 0:
            upstreamCPs.append(0)
            continue
        # Create a copy of the main dataframe
        temp_df = out_routes.copy()
        # find the row of the current CP
        idx_row, _ = np.where(temp_df == i)
        # Mask the downstreams values
        temp_df = temp_df[idx_row, :]
        mask_df = temp_df < i
        # Drop the downstream values
        temp_df = temp_df[~mask_df]
        # Get the unique values
        temp_df = np.unique(temp_df).astype("uint16")
        # party = CP_fishnet.loc[temp_df, "pctSurface"]
        # sum = CP_fishnet.loc[temp_df, "pctSurface"].sum()
        upstreamCPs.append(temp_df.tolist())
        CP_fishnet.loc[i, "cumulPctSurf"] = CP_fishnet.loc[temp_df,
                                                           "pctSurface"].sum()
        # Compute cumul areas in km2
        CP_fishnet.loc[i, "cumulArea"] = CP_fishnet.loc[temp_df,
                                                        "Area"].sum()*1e-6
        # CP_fishnet.loc[i, "cumulPctSurf"] = CP_fishnet.iloc[temp_df, 10].sum()

    # Now, fix the last upstream CP with the area
    CP_fishnet.loc[i+1, "cumulPctSurf"] = CP_fishnet.loc[i+1,
                                                         "pctSurface"]
    # Compute cumul areas in km2
    CP_fishnet.loc[i+1, "cumulArea"] = CP_fishnet.loc[i+1,
                                                      "Area"]*1e-6
    # Convert it into a numpy array
    lists_len = [
        len(i) for i in upstreamCPs if isinstance(i, list)]
    # Create a zero array to store the values
    upstreamcps_array = np.zeros([len(upstreamCPs), max(lists_len)])
    # Fill up the array using the lists in the main dataframe
    for i in range(len(upstreamCPs)):
        if isinstance(upstreamCPs[i], list):
            list_up = upstreamCPs[i]
            upstreamcps_array[i, 0:len(list_up)] = list_up
    upstreamcps_array = np.delete(upstreamcps_array, 0, 0)
    upstreamcps_array = np.pad(upstreamcps_array, ((0, 1), (0, 0)))
    # Assign the value to the last entry
    upstreamcps_array[-1, 0] = len(upstreamCPs)
    # return allroute_lists_array
    return CP_fishnet, upstreamcps_array
