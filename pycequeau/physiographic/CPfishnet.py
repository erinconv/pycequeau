import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal
from pycequeau.core import utils as u
import rasterstats as rs
from math import ceil
import itertools
import sys


def convert_coords_to_index(df: gpd.GeoDataFrame,
                            dataset: gdal.Dataset) -> gpd.GeoDataFrame:
    band = dataset.GetRasterBand(1)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    data = band.ReadAsArray(0, 0, cols, rows)
    df2 = pd.concat([df, pd.DataFrame(columns=["col_min", "row_min", "col_max", "row_max"],
                                      index=df.index)], axis=1)
    transform = dataset.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    for i in range(len(df)):
        df2.at[i, "col_min"] = ceil((df["minx"].iloc[i] - xOrigin)/pixelWidth)
        df2.at[i, "row_max"] = int((yOrigin - df["miny"].iloc[i])/pixelHeight)
        df2.at[i, "col_max"] = ceil((df["maxx"].iloc[i]-xOrigin)/pixelWidth)
        df2.at[i, "row_min"] = int((yOrigin - df["maxy"].iloc[i])/pixelHeight)
    return df2


def find_neighbors(gdf: gpd.GeoDataFrame,
                   id: str) -> gpd.GeoDataFrame:
    # https://gis.stackexchange.com/questions/281652/finding-all-neighbors-using-geopandas
    # Drop the column if it exist
    if 'NEIGHBORS' in gdf.columns:
        gdf = gdf.drop(columns=["NEIGHBORS", "KEEP"])
    # add NEIGHBORS column
    gdf = gdf.reindex(columns=gdf.columns.tolist() + ['NEIGHBORS', 'KEEP'])
    gdf["NEIGHBORS"] = None
    gdf["KEEP"] = None
    for index, CP in gdf.iterrows():
        # get 'not disjoint' countries
        neighbors = gdf[~gdf.geometry.disjoint(CP.geometry)][id].tolist()
        keep = gdf[~gdf.geometry.disjoint(CP.geometry)].Dissolve.tolist()
        # remove own name of the country from the list
        keep = [bool(numb) for numb, name in zip(
            keep, neighbors) if CP[id] != name]
        neighbors = [name for name in neighbors if CP[id] != name]
        # add names of neighbors as NEIGHBORS value
        gdf.at[index, "NEIGHBORS"] = neighbors
        gdf.at[index, "KEEP"] = keep
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
    CP_fishnet.at[mask_SUB, "Dissolve"] = 0
    CP_fishnet.index = CP_fishnet["CPid"].values
    # CP_fishnet["CPid"] = CP_fishnet.index
    CP_fishnet.index = CP_fishnet.index.rename("ind")
    return CP_fishnet


def remove_border_CPs(CE_fishnet: gpd.GeoDataFrame,
                      CP_fishnet: gpd.GeoDataFrame,
                      FAC: str) -> list:
    # Get the area of the CE grid
    CE_area = CE_fishnet.area[0]
    # Add the bounds for each polygon

    bounds = CE_fishnet["geometry"].bounds
    CE_fishnet = pd.concat([CE_fishnet, bounds], axis=1)
    FAC_dataset = gdal.Open(FAC)
    CE_fishnet = convert_coords_to_index(CE_fishnet, FAC_dataset)

    CP_fishnet = pd.concat([CP_fishnet, pd.DataFrame(columns=["maxFAC"],
                                                     index=CP_fishnet.index)], axis=1)
    columnsCP = CP_fishnet.columns.tolist()
    for index, CE in CE_fishnet.iterrows():
        # Get the index for all the subbasin features
        idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
        # Get all features inside the CE
        CE_features = CP_fishnet.iloc[idx]
        # If there is not features to dissolve, get rid off
        stats = rs.zonal_stats(CE_features, FAC, stats=['max'])
        CP_fishnet.iloc[idx, columnsCP.index("maxFAC")] = [
            s['max'] for s in stats]
        # Update values
        CE_features = CP_fishnet.iloc[idx]
        # Find neighbors
        CE_features = find_neighbors(CE_features, "CPid")
        # print(CE_features)
        for i, CP in CE_features.iterrows():
            # Take only the CP labeled with dissolve
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
            if CE_features.loc[CP["NEIGHBORS"], "maxFAC"].isnull().values.any() and CP["Dissolve"] == 1:
                CP_fishnet.at[i, "CPid"] = 0.0
                CP_fishnet.at[i, "Dissolve"] = 0
                CP_fishnet.at[i, "CEid"] = 0
            # Get the maximum index value of the neigbourhs in the subset
            # idx_max = CE_features["maxFAC"][CP["NEIGHBORS"]].isnull().values.any()

    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    # Save file
    # CP_fishnet.index = range(len(CP_fishnet))
    CP_fishnet.index = CP_fishnet.index.rename("ind")
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    CP_fishnet.at[:, "Area"] = CP_fishnet.area

    return CP_fishnet, CE_fishnet


def remove_smallCP(CE_fishnet: gpd.GeoDataFrame,
                   CP_fishnet: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    CP_fishnet.index = CP_fishnet["CPid"].values

    for index, CE in CE_fishnet.iterrows():
        # Get the index for all the subbasin features
        idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
        # Get all features inside the CE
        CE_features = CP_fishnet.iloc[idx]
        CE_features = find_neighbors(CE_features, "CPid")
        for i, CP in CE_features.iterrows():
            # Do not ehck the already labaled CPs
            if CP["CPid"] == 0:
                continue

            # Check if need to dissolve
            if CP["Dissolve"]:
                # Check if there are no neighbors
                if not CP["NEIGHBORS"]:
                    CP_fishnet.at[i, "CPid"] = 0.0
                    CP_fishnet.at[i, "maxFAC"] = None
                    CP_fishnet.at[i, "Dissolve"] = 0
                else:
                    # There are specific cases for dissolving the CP
                    # Here those cases are depicte
                    # Columns in dataframe
                    # columns = CE_features.columns.tolist()
                    # Find all the CP to dissolve
                    # idx_dissolve, = np.where(tuple(CE_features["Dissolve"]==1))
                    # Get the neighbors and the condition
                    neighbors = CE_features.loc[i, "NEIGHBORS"]
                    dissolve = CE_features.loc[i, "KEEP"]
                    # Check if there are neighbors to dissolve also
                    if True in dissolve:
                        # Filter by dissole or not dissolve
                        to_dissolve = list(
                            itertools.compress(neighbors, dissolve))
                        not_dissolve = list(itertools.compress(
                            neighbors, np.logical_not(dissolve)))
                        # 1- There are only neighbors to dissolve:
                        if to_dissolve and not not_dissolve:
                            # idx_max = CE_features.loc[to_dissolve, "maxFAC"].idxmax()
                            idx_max = CE_features.loc[to_dissolve, "Area"].idxmax(
                            )
                            CP_fishnet.at[to_dissolve, "CPid"] = idx_max
                            CP_fishnet.at[to_dissolve, "Dissolve"] = 0
                            CE_features.at[to_dissolve, "Dissolve"] = 0
                            # CE_features.at[to_dissolve, "CPid"] = idx_max
                            # CE_features = find_neighbors(CE_features, "CPid")
                        # CE 545 - interesting
                        # 2- There are both, neighbors to dissolve and not to
                        elif to_dissolve and not_dissolve:
                            # idx_max = CE_features.loc[not_dissolve, "maxFAC"].idxmax()
                            idx_max = CE_features.loc[not_dissolve, "Area"].idxmin(
                            )
                            CE_features.at[to_dissolve, "Dissolve"] = 0
                            CP_fishnet.at[to_dissolve, "CPid"] = idx_max
                            CP_fishnet.at[to_dissolve, "Dissolve"] = 0
                            # CE_features.at[to_dissolve, "CPid"] = idx_max
                            # CE_features = find_neighbors(CE_features, "CPid")
                    else:
                        # Here the CPs with no neighbors to dissolve are processed
                        idx_max = CE_features.loc[CP["NEIGHBORS"], "maxFAC"].idxmax(
                        )
                        CP_fishnet.at[i, "CPid"] = idx_max
                        CP_fishnet.at[i, "CEid"] = CE["CEid"]
                        CP_fishnet.at[i, "Dissolve"] = 0

    # Save file
    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    CP_fishnet.index = CP_fishnet.index.rename("ind")
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    CP_fishnet.at[:, "Area"] = CP_fishnet.area
    return CP_fishnet


def dissolve_pixels(CE_fishnet: gpd.GeoDataFrame,
                    CP_fishnet: gpd.GeoDataFrame,
                    area_th) -> gpd.GeoDataFrame:

    # Here to dissolve the CE
    CEarea = CE_fishnet.area[0]
    # This part drops the CEs where there exist multipolygons
    idx_multi, = np.where(
        CP_fishnet["geometry"].geom_type.values == "MultiPolygon")
    # Detele the zero entry, which is non data
    idx_multi = np.delete(idx_multi, 0)
    SubCP_fishnet = CP_fishnet.iloc[idx_multi]
    CEsDrop = np.unique(SubCP_fishnet["CEid"].values)
    SubCP_fishnet = CP_fishnet.loc[CP_fishnet["CEid"].isin(CEsDrop)]
    CP_fishnet = CP_fishnet.loc[~CP_fishnet["CEid"].isin(CEsDrop)]

    # Start looping in the CEs that need to be dissolved
    for CE in CEsDrop:
        # Get the index for all the subbasin features
        idx, = np.where(SubCP_fishnet["CEid"] == CE)
        CE_features = SubCP_fishnet.iloc[idx].explode()
        # Replace the index values. Select large values to avoid coinciding with the
        # CPs already storaged in the main dataset
        CP_vals = np.random.randint(99999, 999999, len(CE_features), dtype=int)
        # Compute the features of each CP to find the neighbors
        CE_features.index = CP_vals
        CE_features.at[:, "CPid"] = CP_vals
        CE_features.at[:, "Area"] = CE_features.area
        mask_CP = CE_features["Area"] < area_th*CEarea
        CE_features.loc[mask_CP.values, "Dissolve"] = 1
        CE_features = find_neighbors(CE_features, "CPid")
        columns = CE_features.columns.tolist()
        # Find features with an area of 400m
        idx_400km, = np.where(CE_features.loc[:, "Area"] <= 400.0)
        CP_400km_neighs = CE_features.iloc[idx_400km, columns.index(
            "NEIGHBORS")].values
        if len(idx_400km) == 1:
            # CE_features.iloc[idx_400km,columns.index("CPid")] = CP_400km_neighs[0][0]
            CE_features.iloc[idx_400km, columns.index("Dissolve")] = 0
            # Loop to check multi polygons
            flag_neig = True
            count = 0
            while flag_neig:
                # NEIGHBORS
                CE_features.iloc[idx_400km, columns.index(
                    "CPid")] = CP_400km_neighs[0][count]
                dissolve = CE_features.dissolve(by="CPid", aggfunc="max")
                if "MultiPolygon" not in dissolve.geom_type.values:
                    flag_neig = False
                count += 1

            if 'NEIGHBORS' in CE_features.columns:
                CE_features = CE_features.drop(columns=["NEIGHBORS", "KEEP"])
            CP_fishnet = pd.concat(
                [CP_fishnet, CE_features.iloc[idx_400km, :]], axis=0)
            # Drop the values
            CE_features = CE_features.drop(
                index=CE_features.iloc[idx_400km].index)
        elif len(idx_400km) > 1:
            for ind in range(len(idx_400km)):
                # CE_features.iloc[idx_400km[ind],columns.index("CPid")] = CP_400km_neighs[ind][0]
                CE_features.iloc[idx_400km[ind], columns.index("Dissolve")] = 0
                # Loop to check multi polygons
                flag_neig = True
                count = 0
                while flag_neig:
                    # NEIGHBORS
                    CE_features.iloc[idx_400km[ind], columns.index(
                        "CPid")] = CP_400km_neighs[ind][count]
                    dissolve = CE_features.dissolve(by="CPid", aggfunc="max")
                    if "MultiPolygon" not in dissolve.geom_type.values:
                        flag_neig = False
                    count += 1
            # Add only this CPS
            if 'NEIGHBORS' in CE_features.columns:
                CE_features = CE_features.drop(columns=["NEIGHBORS", "KEEP"])
            CP_fishnet = pd.concat(
                [CP_fishnet, CE_features.iloc[idx_400km[:], :]], axis=0)
            # Drop the values
            CE_features = CE_features.drop(
                index=CE_features.iloc[idx_400km].index)
        CP_fishnet = pd.concat([CP_fishnet, CE_features], axis=0)

    # Dissolve to make sure everything is restarted
    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    CP_fishnet.index = CP_fishnet.index.rename("index")
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    CP_fishnet = pd.concat([CP_fishnet, CE_features], axis=0)
    if 'NEIGHBORS' in CP_fishnet.columns:
        CP_fishnet = CP_fishnet.drop(columns=["NEIGHBORS", "KEEP"])
    # Update the mask of the already merged values
    mask_CP = CP_fishnet.loc[:, "Area"] < area_th*CEarea
    CP_fishnet.loc[mask_CP, "Dissolve"] = 1
    CP_fishnet.loc[np.logical_not(mask_CP), "Dissolve"] = 0
    # Treat the left-over cases
    if CP_fishnet["Dissolve"].any() == 1:
        idx_CES = np.where(CP_fishnet["Dissolve"] == 1)
        CEsDrop = np.unique(CP_fishnet.iloc[idx_CES]["CEid"].values)
        left_overs_CEs = CP_fishnet.loc[CP_fishnet["CEid"].isin(CEsDrop)]
        CP_fishnet = CP_fishnet.loc[~CP_fishnet["CEid"].isin(CEsDrop)]
        # Loop
        for CE in CEsDrop:
            idx, = np.where(left_overs_CEs["CEid"] == CE)
            CE_features = left_overs_CEs.iloc[idx]
            CE_features = find_neighbors(CE_features, "CPid")
            columns = CE_features.columns.tolist()
            # Check the cases.
            unique, counts = np.unique(
                CE_features.loc[:, "Dissolve"].values, return_counts=True)
            # 1- There is only one CP left to merge:
            if counts[1] == 1:
                # Find the CP to dissolve
                idx_dissolve, = np.where(
                    CE_features.loc[:, "Dissolve"].values == 1)
                # Find the neighbor with the maximum FAC
                neig_list = CE_features.iloc[idx_dissolve, columns.index(
                    "NEIGHBORS")].values[0]
                idx_FAC = CE_features.loc[neig_list, "maxFAC"].idxmax()
                # Replace the CP values within the CE
                CE_features.iloc[idx_dissolve,
                                 columns.index("CPid")] = idx_FAC.real
                CE_features.iloc[idx_dissolve, columns.index("Dissolve")] = 0
                # Merge it with the main dataset
                if 'NEIGHBORS' in CE_features.columns:
                    CE_features = CE_features.drop(
                        columns=["NEIGHBORS", "KEEP"])
                CP_fishnet = pd.concat([CP_fishnet, CE_features], axis=0)
            else:
                # Sort by dissolve or not
                CE_features = CE_features.sort_values(
                    by='Dissolve', ascending=False)
                # Loop into each CP
                # for idx_CP in range(len(CE_features)):
                for index, CP in CE_features.iterrows():
                    # Create an exception here since at each iteration, the dataframe
                    # is being reduced. Just to make sure it does not crash in run time
                    try:
                        CE_features.loc[index, "Dissolve"] == 1
                    except:
                        # Continue to jump to te next CP
                        continue
                    # Check if this needs to be dissolve
                    if CE_features.loc[index, "Dissolve"] == 1:
                        # Check the neigbors cases also
                        neig_list = CE_features.loc[index, "NEIGHBORS"]
                        # Only one neigh.
                        if len(neig_list) == 1:
                            # Check if this neighboor need to be dissolved also
                            if CE_features.loc[neig_list, "KEEP"].values:
                                # Set the two dissolve values to zero
                                CE_features.loc[index, "Dissolve"] = 0
                                CE_features.loc[neig_list, "Dissolve"] = 0
                                # Now merge the two CPs and drop them
                                CE_features.loc[index, "CPid"] = neig_list[0]
                                CP_fishnet = pd.concat(
                                    [CP_fishnet, CE_features.loc[[index]]], axis=0)
                                CP_fishnet = pd.concat(
                                    [CP_fishnet, CE_features.loc[neig_list]], axis=0)
                                # Drop the values
                                CE_features = CE_features.drop(index=index)
                                CE_features = CE_features.drop(
                                    index=CE_features.loc[neig_list].index)
                    else:
                        CP_fishnet = pd.concat(
                            [CP_fishnet, CE_features], axis=0)
                        continue
    # Drop the unnecesary columns
    if 'NEIGHBORS' in CP_fishnet.columns:
        CP_fishnet = CP_fishnet.drop(columns=["NEIGHBORS", "KEEP"])
    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    CP_fishnet.index = CP_fishnet.index.rename("index")
    CP_fishnet.index = range(len(CP_fishnet))
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values

    return CP_fishnet


def force_4CP(CE_fishnet: gpd.GeoDataFrame,
              CP_fishnet: gpd.GeoDataFrame,
              area_th: float) -> gpd.GeoDataFrame:

    mask_CP = CP_fishnet["Area"] < area_th*CE_fishnet.area[0]
    CP_fishnet.at[mask_CP, "Dissolve"] = 1

    for index, CE in CE_fishnet.iterrows():

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
            CE_features.at[:, "Area"] = CE_features.area
            columns = CE_features.columns.tolist()
            # Find the smallest CP
            idx_small, = np.where(
                CE_features.loc[:, "Area"].values == np.amin(CE_features["Area"]))
            # Find the neighbor with the highest flow accumulation
            neighbors = CE_features.iloc[idx_small,
                                         columns.index("NEIGHBORS")].values[0]
            maxFAC = np.amax(CE_features.loc[neighbors, "maxFAC"])
            idx_maxFAC, = np.where(
                CE_features.loc[neighbors, "maxFAC"] == maxFAC)
            # Get the cpid of the replaced value
            idx_replaced = CE_features.iloc[idx_small, columns.index("CPid")]
            # Replace the value in the main dataframe
            CP_fishnet.loc[idx_replaced, "CPid"] = neighbors[idx_maxFAC[0]]
            # Replace the CPid value into the small CP
            CE_features.iloc[idx_small, columns.index(
                "CPid")] = neighbors[idx_maxFAC[0]]
            # Now dissolve all CE_features and rearange the things
            CE_features = CE_features.dissolve(by="CPid", aggfunc="max")
            CE_features.index = CE_features.index.rename("index")
            CE_features["CPid"] = CE_features.index.values

    # Save file
    CP_fishnet = CP_fishnet.dissolve(by="CPid", aggfunc="max")
    CP_fishnet.index = CP_fishnet.index.rename("index")
    CP_fishnet.loc[:, "CPid"] = CP_fishnet.index.values
    CP_fishnet.at[:, "Area"] = CP_fishnet.area
    return CP_fishnet


def compute_flow_path(flow_dir: np.ndarray,
                      flow_acc: np.ndarray,
                      flow_th: float) -> np.ndarray:
    # Mask the flow accumulation based on the threshold
    flow_accu_mask = (flow_acc > flow_th)
    rows, cols = np.indices(flow_dir.shape)
    stream_network = np.zeros_like(flow_dir)
    counter = 1
    for row in range(flow_dir.shape[0]):
        for col in range(flow_dir.shape[1]):
            if flow_accu_mask[row, col]:
                next_row = row
                next_col = col
                # counter = 1
                while True:
                    direction = flow_accu_mask[next_row, next_col]
                    if direction == 0:
                        break
                    next_row += [-1, -1, 0, 1, 1, 1, 0, -1][direction - 1]
                    next_col += [0, 1, 1, 1, 0, -1, -1, -1][direction - 1]
                    if stream_network[next_row, next_col] > 0:
                        break
                stream_network[row, col] = counter
        # counter +=1
    return stream_network


def main_path(flow_dir: np.ndarray,
              flow_acc: np.ndarray,
              flow_th: float) -> np.ndarray:
    # Create a mask to extract only the cells with flow accumulation greater than a certain threshold
    flow_accumulation_mask = (flow_acc > flow_th)

    # Find the outlet cell of the catchment (i.e., the cell with the minimum flow accumulation)
    outlet_row, outlet_col = np.unravel_index(
        np.argmin(flow_acc), flow_acc.shape)

    # Create an empty list to store the cells in the main stream
    main_stream_cells = []

    # Trace the main stream from the outlet cell to the start of the main stream
    next_row, next_col = outlet_row, outlet_col
    while True:
        main_stream_cells.append((next_row, next_col))
        direction = flow_dir[next_row, next_col]
        if direction == 0:  # Check if the current cell is a sink cell and exit the loop if it is
            break
        # Calculate the row and column indices of the next cell in the main stream
        next_row += [-1, -1, 0, 1, 1, 1, 0, -1][direction - 1]
        next_col += [0, 1, 1, 1, 0, -1, -1, -1][direction - 1]
        # Check if the next row or column index is out of bounds and exit the loop if it is
        if next_row < 0 or next_row >= flow_dir.shape[0] or next_col < 0 or next_col >= flow_dir.shape[1]:
            break
        # Check if the current cell has more than one upstream cell and exit the loop if it does
        upstream_count = 0
        for i, j in [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]:
            row = next_row + i
            col = next_col + j
            if row >= 0 and row < flow_dir.shape[0] and col >= 0 and col < flow_dir.shape[1]:
                if flow_dir[row, col] == (8 - direction):
                    upstream_count += 1
                    if upstream_count > 1:
                        break
        if upstream_count > 1:
            break

    # Create a new raster to represent the main stream cells
    main_stream_raster = np.zeros_like(flow_dir)
    for cell in main_stream_cells:
        main_stream_raster[cell[0], cell[1]] = 1


def routing_table(CE_fishnet: gpd.GeoDataFrame,
                  CP_fishnet: gpd.GeoDataFrame,
                  FAC: str,
                  DIR: str,
                  CP_array: np.ndarray,
                  flow_th: float) -> pd.DataFrame:

    # CP_fishnet.index = range(len(CP_fishnet))
    # CP_fishnet["CPid"] = range(len(CP_fishnet))
    
    
    routing = pd.DataFrame(columns=["CPid", "inCPid", "outlet_row",
                                    "outlet_col", "inlet_row", "inlet_col"],
                           index=CP_fishnet.index)
    routing["CPid"] = CP_fishnet["CPid"]

    CP_fishnet = pd.concat([CP_fishnet,
                            CP_fishnet["geometry"].bounds], axis=1)
    # Get the FAC array
    FAC_dataset = gdal.Open(FAC, gdal.GA_ReadOnly)
    band = FAC_dataset.GetRasterBand(1)
    FAC_array = band.ReadAsArray()
    FAC_array[FAC_array < 0] = 0
    # Get the DIR array
    CP_fishnet = convert_coords_to_index(CP_fishnet, FAC_dataset)
    # Loop into each CE
    for index, feat in CP_fishnet.iterrows():
        # Do not work on the non data value
        if feat["CPid"] == 0:
            continue
        # Select the CP and the FAC based on the boundaries
        CP = CP_array[feat["row_min"]-1:feat["row_max"]+1,
                      feat["col_min"]-1:feat["col_max"]+1]
        subFAC = FAC_array[feat["row_min"]-1:feat["row_max"]+1,
                           feat["col_min"]-1:feat["col_max"]+1]
        # Mask the FAC based on this CP
        mask_CP = (CP == feat['CPid']).astype('uint8')
        # Apply mask
        # mask_subFAC = subFAC*mask_CP
        # Get the outlet indexes for the CP
        outlet_row, outlet_col = np.unravel_index(
            np.argmax(subFAC*mask_CP), subFAC.shape)
        if isinstance(outlet_row, np.ndarray):
            sys.exit("There is more than one outlet point")

        # Create mask for subFAC
        mask_FAC = np.zeros(subFAC.shape).astype("uint8")
        mask_FAC[outlet_row - 1:outlet_row + 2,
                 outlet_col - 1:outlet_col + 2] = 1
        # masked_FAC = subFAC*mask_FAC
        # Get the FAC values on the mask
        # Find location of next pixel into which outlet flows
        inlet_row, inlet_col = np.unravel_index(
            np.argmax(subFAC*mask_FAC), subFAC.shape)
        # Add the CP where it discharges
        routing.at[index, "inCPid"] = CP[inlet_row, inlet_col]
        # Add the coordinates into the dataframe
        routing.at[index, "outlet_col"] = outlet_col-1
        routing.at[index, "outlet_row"] = outlet_row-1
        routing.at[index, "inlet_row"] = inlet_row-1
        routing.at[index, "inlet_col"] = inlet_col-1
    # Create the rouring table here
    rtable = pd.DataFrame(columns=["oldCPid", "newCPid", "upstreamCPs"],
                          index=range(len(CP_fishnet)))
    routing["diff"] = routing["CPid"] - routing["inCPid"]
    # Find where the difference is zero
    # This will help us to find the outlet CP
    idx_outlet = routing.index[routing["diff"] == 0].values
    # Add this value in the rtable as the first value
    rtable.loc[[1], "oldCPid"] = routing.loc[idx_outlet, "CPid"].values[0]
    rtable.loc[[1], "newCPid"] = 1
    # Get column list
    columns = rtable.columns.tolist()
    # Set to zero the already tracked CP
    routing.loc[idx_outlet, :"inCPid"] = 0
    new_id_counter = 2
    for i, _ in routing.iterrows():
        if i == 0:
            rtable.loc[i, "newCPid"] = 0
            continue
        # Find the upstream CPs
        idx_outlet = routing.index[routing["inCPid"]
                                   == rtable.at[i, "oldCPid"]]
        # Set to nan the already tracked CP
        routing.loc[idx_outlet, "inCPid"] = -99999
        # print(idx_outlet)
        # print(rtable.at[i,"oldCPid"],i)
        if not idx_outlet.empty:
            # Find the upstream CPs
            upstreams_cps = routing.loc[idx_outlet, "CPid"].values
            # Rename the upstream CPs
            upstreams_cps_newid = list(
                range(new_id_counter, new_id_counter+len(idx_outlet)))
            # Append the CPs vertically. +1 to step out the current index
            rtable.iloc[new_id_counter:new_id_counter +
                        len(idx_outlet), columns.index("newCPid")] = upstreams_cps_newid
            rtable.iloc[new_id_counter:new_id_counter +
                        len(idx_outlet), columns.index("oldCPid")] = upstreams_cps
            rtable.iloc[i, columns.index("upstreamCPs")] = upstreams_cps_newid
            new_id_counter += len(idx_outlet)
            pass
    return rtable


def get_downstream_CP(rtable: pd.DataFrame) -> pd.DataFrame:
    # Create an empty table to store the values
    downstreamCPs = pd.DataFrame(columns=["downstreamCPs"],
                                 index=rtable.index)
    # Convert lists into an array.
    # Get the lenght of each list to get the maximum value
    lists_len = [len(i)
                 for i in rtable["upstreamCPs"].values if isinstance(i, list)]
    # Create a zero array to store the values
    list_upstream_array = np.zeros([len(rtable), max(lists_len)])
    # Fill up the array using the lists in the main dataframe
    for i in range(len(rtable)):
        if isinstance(rtable.loc[i, "upstreamCPs"], list):
            list_up = rtable.loc[i, "upstreamCPs"]
            list_upstream_array[i, 0:len(list_up)] = list_up
        rtable.loc[i, "upstreamCPs"] = list_upstream_array[i,:].tolist()
    # Now identify the Downstream CPs
    for i, table in rtable.iterrows():
        # Skip the first since it is non data value
        if i == 0:
            downstreamCPs.loc[i, "downstreamCPs"] = 0
            continue
        # Find the index of the CP in wwhich the current CP drains
        idx_down, _ = np.where(list_upstream_array == table["newCPid"])
        # Check if the list is empty
        if len(idx_down) == 0:
            downstreamCPs.loc[i, "downstreamCPs"] = 0
        else:
            downstreamCPs.loc[i,
                              "downstreamCPs"] = rtable.loc[idx_down[0], "newCPid"]
    # Concat to the rtable
    rtable = pd.concat([rtable, downstreamCPs], axis=1)
    return rtable


def outlet_routes(rtable: pd.DataFrame) -> pd.DataFrame:

    # Create the array to store the lists
    allroute_lists = pd.DataFrame(columns=["outletRoutes"],
                                  index=rtable.index)
    # Start looping the array
    for i, _ in allroute_lists.iterrows():
        # Skip the zero index
        if i == 0:
            continue
        route_list = []
        # Set the number of the upstream CPs to 1
        nu = 1
        # Set the upstream CP to
        up = rtable.loc[i, "newCPid"]

        while up != 1:
            # append Nth CP to list
            route_list.append(up)
            # go down from Nth position until end
            up = rtable.loc[up, "downstreamCPs"]
        route_list.append(1)
        allroute_lists.loc[i, "outletRoutes"] = route_list
    # Convert it into a numpy array
    lists_len = [
        len(i) for i in allroute_lists["outletRoutes"].values if isinstance(i, list)]
    # Create a zero array to store the values
    allroute_lists_array = np.zeros([len(allroute_lists), max(lists_len)])
    # Fill up the array using the lists in the main dataframe
    for i in range(len(allroute_lists)):
        if isinstance(allroute_lists.loc[i, "outletRoutes"], list):
            list_up = allroute_lists.loc[i, "outletRoutes"]
            allroute_lists_array[i, 0:len(list_up)] = list_up
    return allroute_lists_array


def cumulative_areas(CP_fishnet: gpd.GeoDataFrame,
                     CE_fishnet: gpd.GeoDataFrame,
                     outlet_routes: np.ndarray) -> pd.DataFrame:
    # Update areas of the CP and get the CE area
    CE_area = CE_fishnet.area[0]
    CP_fishnet["Area"] = CP_fishnet.area
    # Get the percentage of that area
    CP_fishnet["pctSurface"] = (CP_fishnet["Area"]/CE_area)*100
    # Cumulative areas
    CP_fishnet["cumulPctSurf"] = 0.0
    upstreamCPs = []
    for i in range(len(outlet_routes)):
        if i == 0:
            upstreamCPs.append(0)
            continue
        # Create a copy of the main dataframe
        temp_df = outlet_routes.copy()
        # find the row of the current CP
        idx_row , _ = np.where(temp_df == i)
        # Mask the downstreams values
        temp_df = temp_df[idx_row,:] 
        mask_df = temp_df < i
        # Drop the downstream values
        temp_df = temp_df[~mask_df]
        # Get the unique values
        temp_df = np.unique(temp_df).astype("uint16")
        party = CP_fishnet.loc[temp_df,"pctSurface"]
        sum = CP_fishnet.loc[temp_df,"pctSurface"].sum()
        upstreamCPs.append(temp_df.tolist())
        CP_fishnet.loc[i,"cumulPctSurf"] = CP_fishnet.loc[temp_df,"pctSurface"].sum()
    return CP_fishnet, upstreamCPs

def renumber_fishnets(CP_fishnet: gpd.GeoDataFrame,
                     CE_fishnet: gpd.GeoDataFrame,
                     rtable: pd.DataFrame) -> gpd.GeoDataFrame:
    # Renumbering the CPs
    CP_fishnet["newCPid"] = 0
    CP_fishnet["newCEid"] = 0
    CE_fishnet["newCEid"] = 0
    CE_track_list = []
    # get the columns to use them as indexers
    columns_rtable = rtable.columns.tolist()
    columns_CPfishnet = CP_fishnet.columns.tolist()
    columns_CEfishnet = CE_fishnet.columns.tolist()
    # Iterate over the dataframe to rename CP fishnet
    for i in range(len(CP_fishnet)):
        # Skip the non data first value
        if i == 0:
            CP_fishnet.iloc[0,columns_CPfishnet.index("newCPid")] = 0
            continue
        # Find the index of the old CPid in the main dataframe
        idx_old, = np.where(CP_fishnet["CPid"].values == rtable.iloc[i,columns_rtable.index("oldCPid")])
        # Replace the value with the new CPid value
        CP_fishnet.iloc[idx_old,columns_CPfishnet.index("newCPid")] = int(rtable.iloc[i,columns_rtable.index("newCPid")])
    # Change the index in the main dataframe
    CP_fishnet.index = CP_fishnet["newCPid"].values
    # Sort values
    CP_fishnet = CP_fishnet.sort_values(by=["newCPid"])
    # Renumbering the CEs. This is possible since the values are
    # already sorted in the main dataframe based on the new CPids
    CE_id = 1
    for i in range(len(CP_fishnet)):
        # Skip the first
        if i == 0:
            CE_track_list.append(0)
            continue
        # Check if the CE has already been tracked
        if CP_fishnet.loc[i,"CEid"] in CE_track_list:
            pass
        else:
            # Track the CE in the CP fishnet
            CE_track_list.append(CP_fishnet.loc[i,"CEid"])
            # Find the value in the dataset
            idx_new, = np.where(CE_fishnet["CEid"] == CE_track_list[CE_id])
            # Add the value in the column
            CE_fishnet.iloc[idx_new,columns_CEfishnet.index("newCEid")] = CE_id
            # Find the value in the CPfishnet dataset
            idx_new2, = np.where(CP_fishnet["CEid"] == CE_track_list[CE_id])
            CP_fishnet.iloc[idx_new2,columns_CPfishnet.index("newCEid")] = CE_id
            CE_id += 1
    CE_fishnet = CE_fishnet.sort_values(by=["newCEid"])
    return CP_fishnet, CE_fishnet


def get_atitudes():
    pass

# Add the data to the routing table

# subDIR = subDIR*mask_FAC
# Get the stream network
# stream = compute_flow_path(subDIR,subFAC,np.nanmean(subFAC))

# outlet_row, outlet_col = np.where(subFAC == np.amax(subFAC))
# if isinstance(outlet_row,np.ndarray) > 1:
#     print("AAA")
#     break
# else:
#     print(feat['CPid'],outlet_row, outlet_col)
#     print(subFAC[outlet_row, outlet_col])
# print(np.unique(CP))
# print(CP)
# print(subDIR)
# print(stream.astype("uint8"))
# print(subFAC)
# print(mask_FAC)
#     # Get the index for all the subbasin features
#     idx, = np.where(CP_fishnet["CEid"] == CE["CEid"])
#     # Get all features inside the CE
#     CE_features = CP_fishnet.iloc[idx]
#     # Check if the CE is empty
#     if CE_features.empty:
#         continue

# Get the ids for each subbasin
# id_sub = np.unique(CP_fishnet["OBJECTID"][mask_CP])
# id_CE = np.unique(CP_fishnet["id"][mask_CP])
# CP_fishnet = CP_fishnet.dissolve(by='CPid', aggfunc='first')
# Loop into each subbasin
# for CE in id_CE:
#     # Get the index for all the subbasin features
#     idx, = np.where(CP_fishnet["id"] == 527)
#     CE_features = CP_fishnet.iloc[idx]
#     for i in range(len(CE_features)):
#         # print(CE_features["Dissolve"])
#         if CE_features["Dissolve"].iloc[i]==1:
#             print(CE_features["geometry"].touches(CE_features["geometry"],align=True))
#     # print(CE_features["Dissolve"].iloc[1])
#     break
# for
# Mask nan
# mask_NaN = CP_fishnet[mask_CP]
# print(CP_fishnet[mask_CP])
# print(CP_area[mask_CP])
# print(len(CP_fishnet))
# print(len(CP_fishnet.explode()))
# mask = CP_fishnet['CPid']==759
# print(CP_fishnet[mask])
# return
