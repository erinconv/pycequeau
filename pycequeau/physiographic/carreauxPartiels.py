import numpy as np
import pandas as pd
from pycequeau.core import utils as u
import geopandas as gpd
from osgeo import gdal
import os


def get_CP_coordinates(carreuxEntiers: pd.DataFrame,
                       CPfishnet: gpd.GeoDataFrame)->pd.DataFrame:
    # Create dataset to store the coordinates
    coordinates = pd.DataFrame(columns=["CPid","i","j"],
                               index=CPfishnet.index)
    coordinates["CPid"] = CPfishnet["newCPid"]
    columns_carreux = coordinates.columns.tolist()
    # Loop into the carreuxEntiers dataset 
    for index,carreu in carreuxEntiers.iterrows():
        # Find the CEid in the CP fisnhet dataframe
        idx_CEs, = np.where(CPfishnet["newCEid"] == carreu["CEid"])
        # Place this values in the coordinates dataframe
        coordinates.iloc[idx_CEs,columns_carreux.index("i")] = carreu["i"]
        coordinates.iloc[idx_CEs,columns_carreux.index("j")] = carreu["j"]
    # Drop nan if it exist
    coordinates = coordinates.dropna(axis="index")
    return coordinates

def get_codes(CPfishnet: gpd.GeoDataFrame)->pd.DataFrame:
    # Get unique CEids
    CEids = np.unique(CPfishnet["newCEid"])
    codes = np.zeros(len(CPfishnet))
    for i in range(len(CEids)):
        # Skip value
        if CEids[i] == 0:
            continue
        idx, = np.where(CPfishnet["newCEid"] == CEids[i])
        codes[idx] = list(range(65,65+len(idx))) 
    # 
    codes = np.delete(codes,0,axis=0)
    return codes

def cumulate_variables(outlet_routes: np.ndarray,
                       pctForet:np.ndarray,
                       pctLacRiviere: list,
                       pctMarais: list)->pd.DataFrame:
    outlet_routes = np.delete(outlet_routes,0,axis=0)
    # 
    # "cumulPctSuperficieLacsAmont",
    #  "cumulPctSuperficieMaraisAmont","cumulPctSuperficieForetAmont"
    # Create the dataset to store the cumulates values
    cumulates = pd.DataFrame(columns=["cumulPctSuperficieLacsAmont",
                                      "cumulPctSuperficieMaraisAmont",
                                      "cumulPctSuperficieForetAmont"],
                             data=np.c_[pctLacRiviere,
                                        pctMarais,
                                        pctForet],
                             index=range(1,len(outlet_routes)+1))
    # Track the upstream cps
    upstreamCPs = []
    for i in range(1,len(outlet_routes)+1):
        if i == 0:
            upstreamCPs.append(0)
            continue
        # Create a copy of the main dataframe
        temp_df = outlet_routes.copy()
        # find the row of the current CP
        idx_row, _ = np.where(temp_df == i)
        # Mask the downstreams values
        temp_df = temp_df[idx_row,:] 
        mask_df = temp_df < i
        # Drop the downstream values
        temp_df = temp_df[~mask_df]
        # Get the unique values
        temp_df = np.unique(temp_df).astype("uint16")
        # Cumulate all the variables
        cumulates.loc[i,"cumulPctSuperficieLacsAmont"] = np.sum(np.array(pctLacRiviere)[temp_df-1])
        cumulates.loc[i,"cumulPctSuperficieForetAmont"] = np.sum(pctForet[temp_df-1])
        cumulates.loc[i,"cumulPctSuperficieMaraisAmont"] = np.sum(np.array(pctMarais)[temp_df-1])
    
    return cumulates

def get_river_geometry(CPfishnet: gpd.GeoDataFrame,
                       rtable: pd.DataFrame)->pd.DataFrame:
    # Create the dataframe to store the data
    river_geometry = pd.DataFrame(columns=["profondeurMin","longueurCoursEauPrincipal",
                                               "largeurCoursEauPrincipal","penteRiviere"],
                                  index = rtable.index)
    # Get the cumulated area in the upstream CPS
    sum_cp_areas = []
    print(rtable.columns)
    for i in range(1,len(rtable)+1):
        # Check if the values in the table are read as string
        if isinstance(rtable.loc[i,"upstreamCPs"],str):
            CP_list = np.array(eval(rtable.loc[i,"upstreamCPs"]),dtype=np.int16)
        else:
            CP_list = np.array(rtable.loc[i,"upstreamCPs"],dtype=np.int16)
        # Drop zero values
        CP_list = np.trim_zeros(CP_list)
        # append the current CP value to sum the area
        CP_list = np.append(CP_list,i)
        # Sum the area and convert m2 to km2
        sum_cp_areas.append(CPfishnet.loc[CP_list,"Area"].sum()*1.0e-6)
    # Compute the values
    # calculated as function of upstream area (units = cm)
    river_geometry["profondeurMin"] = (0.0198*(np.power(sum_cp_areas,0.53)))*100
    # calculated as function of current CP (units = 1/100 km)
    river_geometry["longueurCoursEauPrincipal"] = np.power(CPfishnet["Area"]*1.0e-6,0.5)*10
    # calculated as function of upstream area (units = 1/10m)
    river_geometry["largeurCoursEauPrincipal"] = (0.49*np.power(sum_cp_areas,0.6))*10
    # (units = 1/1000 metres/km)
    river_geometry["penteRiviere"] = 1000 
    
    return river_geometry