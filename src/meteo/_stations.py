from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from osgeo import gdal, ogr
from src.core import projections
from src.core import utils as u


def create_grid_var(ds: xr.Dataset,
                    CEs_df: pd.DataFrame,
                    var_name: str,
                    datenum: list) -> xr.Dataset:
    """_summary_
    Args:
        ds (xr.Dataset): _description_
        CEs_df (pd.DataFrame): _description_
        var_name (str): _description_
        datenum (list): _description_

    Returns:
        xr.Dataset: _description_
    """

    dr = xr.Dataset(
        data_vars={
            var_name: (
                ["pasTemp", "CEid"], np.zeros(
                    [len(datenum), len(CEs_df)]).astype(np.float32)
            )
        },
        coords={
            "CEid": CEs_df["CEid"].values.astype(np.int32),
            "pasTemp": datenum
        }
    )
    # Set indexes as integers
    CEs_df["j"] = CEs_df["j"].values.astype("int")
    CEs_df["i"] = CEs_df["i"].values.astype("int")
    # Filling up the CEQUEAU grid file
    for idx, ce_id in CEs_df.iterrows():
        # Get the row, col indexes
        i = int(ce_id["i"])
        j = int(ce_id["j"])
        dr[var_name][:, idx-1] = ds[var_name].loc[:,j,i].values.astype(np.float32)
    return dr


def interpolation_netCDF(ds: xr.Dataset,
                         CEregrid: str,
                         table: pd.DataFrame,
                         method: str = "linear") -> xr.DataArray:
    """Interpolate

    Args:
        ds (xr.DataArray): _description_
        CEregrid (str): _description_
        table (pd.DataFrame): _description_
        method (str, optional): _description_. Defaults to "linear".

    Returns:
        xr.DataArray: _description_
    """
    # Create objective
    ce_grid = gdal.Open(CEregrid, gdal.GA_ReadOnly)
    ce_array = np.array(ce_grid.GetRasterBand(1).ReadAsArray())
    i_res = ce_array.shape[1]
    j_res = ce_array.shape[0]
    i = np.arange(0, i_res, 1)+10
    j = np.arange(0, j_res, 1)+10
    # Check whether the longitudes need to be corrected
    if ds["lon"].max() > 0:
        ds["lon"] = ds["lon"].values - 360
        # table["lon"] = table["lon"].values - 360
    # Create the ref objects
    # Find the nearest maximun and minimum values
    lon_max_idx = u.find_nearest(ds["lon"].values, table["lon"].max())
    lon_min_idx = u.find_nearest(ds["lon"].values, table["lon"].min())
    lat_max_idx = u.find_nearest(ds["lat"].values, table["lat"].max())
    lat_min_idx = u.find_nearest(ds["lat"].values, table["lat"].min())
    # Slice the dataset
    ds = ds.isel(lat=slice(min(lat_min_idx, lat_max_idx), max(lat_min_idx, lat_max_idx)),
                 lon=slice(min(lon_min_idx, lon_max_idx), max(lon_min_idx, lon_max_idx)))
    # Create reference regridder
    x = np.linspace(min(i), max(i), len(ds["lon"].values))
    y = np.linspace(min(j), max(j), len(ds["lat"].values))
    # y = np.flip(y)
    ds = ds.assign_coords(lat=y)
    ds = ds.assign_coords(lon=x)
    # Create mask
    j = np.flip(j)
    dsi = ds.interp(time=ds["time"], lat=j, lon=i, method=method)
    # Flip the array since gdal reads it upside down
    # mask interpolated dataset
    dsi = dsi.where(ce_array> 0)
    dsi = dsi.rename(lat="j", lon="i")
    dsi = dsi.transpose("time", "j", "i")
    dsi["i"] = i.astype(np.int16)
    dsi["j"] = j.astype(np.int16)
    ds = ds.rename(lat="j", lon="i")
    dsi = _appendCEgrid(dsi, ce_grid)
    dsi = dsi.assign_attrs(
        interpolated=f"Interpolated using the method: {method} from xr.Dataset.interp function")
    return dsi


def get_netCDF_grids(ds: xr.DataArray,
                     CEgrid: gdal.Dataset,
                     watershed: ogr.DataSource) -> np.ndarray:
    """
    This function retrieves the grid point coordinates from an input netcdf
    We check whether this fall or not whitin the watershed boundaries

    Args:
        ds (xr.DataArray): _description_
        CEgrid (gdal.Dataset): _description_
        watershed (ogr.DataSource): _description_

    Returns:
        np.ndarray: _description_
    """
    xtup, ytup, _ = u.GetExtent(CEgrid)
    epsg_dem = projections.get_proj_code(CEgrid)
    # (lon, lat) -> (x, y)
    x, y = projections.utm_to_latlon((np.amin(xtup), np.amax(xtup)),
                                     (np.amin(ytup), np.amax(ytup)),
                                     epsg_dem)
    # TODO: Add another object to deal with the dataset dimension names. For now, the default names are: [time,lat,lon]
    dy = abs(ds["lat"][0].values - ds["lat"][1].values)
    dx = abs(ds["lon"][0].values - ds["lon"][1].values)
    lon_mask = np.arange(x[0], x[1], dx)
    lat_mask = np.arange(y[0], y[1], dy)
    # Get the extent of the vector layer
    shp_layer = watershed.GetLayer()
    xmin, xmax, ymin, ymax = shp_layer.GetExtent()
    # Convert the shp extent to latlon also
    # (lon, lat) -> (x, y)
    x, y = projections.utm_to_latlon((xmin, xmax),
                                     (ymin, ymax),
                                     epsg_dem)
    watershed_extent = (x, y)
    # Check if retrieved points fall in watershed extent
    xypair = u.falls_in_extent(watershed_extent,
                               lon_mask,
                               lat_mask)
    return xypair


def create_station_table(CEregrid: gdal.Dataset,
                         DEM: gdal.Dataset,
                         lon_utm: np.ndarray,
                         lat_utm: np.ndarray,
                         xy_pair: np.ndarray) -> pd.DataFrame:
    """This code generates the station table used to check which grid points
    fall into the watershed boundaries

    Args:
        CEregrid (gdal.Dataset): _description_
        DEM (gdal.Dataset): _description_
        lat_utm (np.ndarray): _description_
        lon_utm (np.ndarray): _description_
        xy_pair (np.ndarray): _description_

    Returns:
        pd.DataFrame: _description_
    """
    ce_grid = CEregrid.ReadAsArray()
    # ce_grid = np.flipud(ce_grid)
    i_res = ce_grid.shape[1]
    j_res = ce_grid.shape[0]
    i = np.linspace(0, i_res, i_res, dtype=int)
    j = np.linspace(0, j_res, j_res, dtype=int)
    # Flip upside down since we are using geographic indexes
    j = np.flipud(j)
    # # Get raster index
    row, col = u.get_index_list(CEregrid,
                                lon_utm,
                                lat_utm)
    # Check if there are values outsite the actual boundaries of the watershed
    # Check first cols
    row = np.array(row, dtype=np.float32)
    col = np.array(col, dtype=np.float32)
    row[row >= j_res] = np.nan
    xy_pair = xy_pair[~np.isnan(row), :]
    col = col[~np.isnan(row)]
    lat_utm = lat_utm[~np.isnan(row)]
    lon_utm = lon_utm[~np.isnan(row)]
    row = row[~np.isnan(row)]

    # then check the rows
    col[col >= i_res] = np.nan
    xy_pair = xy_pair[~np.isnan(col), :]
    row = row[~np.isnan(col)]
    lat_utm = lat_utm[~np.isnan(col)]
    lon_utm = lon_utm[~np.isnan(col)]
    col = col[~np.isnan(col)]

    # Reconvert to integer
    row = np.array(row, dtype=np.int16)
    col = np.array(col, dtype=np.int16)
    # ~np.isnan(xy_pairs).any(axis=1), :
    # Create dataframe and append to the MeteoObject
    stations = pd.DataFrame(data={
        "id": ["NC-grid-"+str(num) for num in range(len(i[col]))],
        "i": i[col]+10,
        "j": j[row]+10,
        "lat": xy_pair[:, 1],
        "lon": xy_pair[:, 0],
        "lat_utm": lat_utm,
        "lon_utm": lon_utm,
        "CEid": ce_grid[row, col],
        "altitude": u.get_altitude_point(DEM, lat_utm, lon_utm)
    }
    )
    return stations


def _appendCEgrid(ds: xr.Dataset,
                  CEregrid: gdal.Dataset) -> xr.DataArray:
    """_summary_

    Args:
        ds (xr.Dataset): _description_
        CEregrid (gdal.Dataset): _description_

    Returns:
        xr.DataArray: _description_
    """
    grid = CEregrid.ReadAsArray().astype(np.float16)
    grid[grid == 0] = np.nan
    dr = xr.Dataset({
        "CE": (
            ("j", "i"),
            grid
        ),
    },
        attrs=dict(units="-",
                   long_name="Whole squares",
                   name="CE")
    )
    return ds.assign(CE=dr.CE)
