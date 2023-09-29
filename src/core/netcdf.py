import xarray as xr
import numpy as np
import pandas as pd
import os


def clip_netcdf(nc_files_path: str,
                bounds: list,
                output_path: str):
    """_summary_

    Args:
        nc_files_path (str): _description_
        bounds (list): _description_
        output_path (str): _description_
    """
    nc_files_list = os.listdir(nc_files_path)
    nc_files_list = [i for i in nc_files_list if i.endswith(".nc")]
    nc_files_list_path = [os.path.join(nc_files_path, i)
                          for i in nc_files_list if i.endswith(".nc")]
    for i, file in enumerate(nc_files_list_path):
        nc_data = xr.open_dataset(file, engine="netcdf4")
        coords_names = list(nc_data.coords)
        # Rename the latitudes-longitudes
        lat = "lat"
        lon = "lon"
        if "longitude" in coords_names:
            nc_data = nc_data.rename({"longitude": "lon",
                                      "latitude": "lat"})
        if nc_data[lon].min() > 0:
            min_lon = bounds[0]+360
            max_lon = bounds[1]+360
        else:
            min_lon = bounds[0]
            max_lon = bounds[1]

        if nc_data[lat][0] > nc_data[lat][-1]:
            min_lat = bounds[2]
            max_lat = bounds[3]
        else:
            min_lat = bounds[3]
            max_lat = bounds[2]
        sliced_nc = nc_data.sel(**{lat: slice(max_lat, min_lat),
                                   lon: slice(min_lon, max_lon)})
        dimensions = list(sliced_nc.dims)
        nc_var_names = list(sliced_nc.variables)
        for dim in dimensions:
            if dim in nc_var_names:
                nc_var_names.remove(dim)
        if "time_bnds" in nc_var_names:
            nc_var_names.remove("time_bnds")
        if nc_files_list[i] == "tmax.nc":
            sliced_nc = sliced_nc.rename({"tmin": "tmax"})
        path = os.path.join(output_path, nc_files_list[i])
        if os.path.isfile(path):
            os.remove(path)
        sliced_nc.to_netcdf(path)


def intermidiate_interpolation(ds: xr.Dataset, scale: int):
    lon = ds["lon"].values
    lat = ds["lat"].values
    dx_original = lon[1] - lon[0]
    dy_original = lat[1] - lat[0]
    # Get the objectvie dx dy
    dx = dx_original/scale
    dy = dy_original/scale
    # Create new vector for lat lon
    lon_objective = np.arange(lon[0], lon[-1], dx)
    lat_objective = np.arange(lat[0], lat[-1], dy)
    dsi = ds.interp(time=ds["time"], lat=lat_objective,
                    lon=lon_objective, method="nearest")
    return dsi
