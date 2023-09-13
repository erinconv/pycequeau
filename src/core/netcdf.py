import xarray as xr
import pandas as pd
import os


def clip_netcdf(nc_files_path: str, bounds: list, output_path: str):
    nc_files_list = os.listdir(nc_files_path)
    nc_files_list = [i for i in nc_files_list if i.endswith(".nc")]
    nc_files_list_path = [os.path.join(nc_files_path, i)
                          for i in nc_files_list if i.endswith(".nc")]
    for i, file in enumerate(nc_files_list_path):
        nc_data = xr.open_dataset(file, engine="netcdf4")
        if nc_data["lon"].min() > 0:
            min_lon = bounds[0]+360
            max_lon = bounds[1]+360
        else:
            min_lon = bounds[0]
            max_lon = bounds[1]

        if nc_data["lat"][0] > nc_data["lat"][-1]:
            min_lat = bounds[2]
            max_lat = bounds[3]
        else:
            min_lat = bounds[3]
            max_lat = bounds[2]
        sliced_nc = nc_data.sel(**{"lat": slice(max_lat, min_lat),
                                   "lon": slice(min_lon, max_lon)})
        sliced_nc.to_netcdf(os.path.join(output_path, nc_files_list[i]))
