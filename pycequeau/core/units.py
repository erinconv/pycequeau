from pyproj import Proj
import xarray as xr
import numpy as np
import pandas as pd


def units_CORDEX(ds: xr.Dataset):
    var_name = list(ds.keys())[0]
    if var_name == "pr":
        # from m d-1 to mm d-1
        ds[var_name].attrs = dict(units="mm d-1",
                                  long_name="Total precipitation")
        # ds[var_name].values = 86400*ds[var_name].values
        ds = ds.rename({var_name: "pTot"})
    if var_name == "rsds":
        # from W m-2 to MJ m2 d-1
        ds[var_name].attrs = dict(units="MJ m-2 d-1",
                                  long_name="Surface solar radiation")
        # ds[var_name].values = 0.0864*ds[var_name].values
        ds = ds.rename({var_name: "rayonnement"})
    if var_name == "sfcWind":
        # from m s-1 to km h-1
        ds[var_name].attrs = dict(units="km h-1",
                                  long_name="Wind speed")
        # ds[var_name].values = 3.6*ds[var_name].values
        ds = ds.rename({var_name: "vitesseVent"})
    if var_name == "clt":
        ds[var_name].attrs = dict(units="0-1",
                                  long_name="Cloud cover")
        ds[var_name].values = ds[var_name].values/100
        ds = ds.rename({var_name: "nebulosite"})
        
    if var_name == "vp":
        # from kPa to mmHg 7.50062
        # ds[var_name] = ds[var_name].transpose("time","lat","lon")
        ds[var_name].attrs = dict(units="mmHg",
                                  long_name="Vapor pressure")
        ds = ds.rename({var_name: "pression"})
        # ds[var_name].attrs = dict(units="mmHg",
        #                            long_name="Vapor pressure")
        # ds[var_name].values = 7.50062e-3*ds[var_name].values
    if var_name == "tasmin" or var_name == "tasmax":
        # from K to degC
        if var_name == "tasmax":
            long_name = "Maximum daily air temperature"
            ds = ds.rename({var_name: "tMax"})
            var_name = "tMax"
        else:
            long_name = "Minimum daily air temperature"
            ds = ds.rename({var_name: "tMin"})
            var_name = "tMin"
        ds[var_name].attrs = dict(units="C",
                                  long_name=long_name)
        # ds[var_name].values = ds[var_name].values - 273.15
    return ds


def units_ERA(ds: xr.Dataset):
    # List the variable names
    var_name = list(ds.keys())[0]
    if var_name == "tp":
        # from m d-1 to mm d-1
        ds[var_name].attrs = dict(units="mm d-1",
                                  long_name="Total precipitation")
        ds[var_name].values = 1e3*ds[var_name].values
        ds = ds.rename({var_name: "pTot"})
    if var_name == "ssr":
        # from W m-2 to MJ m2 d-1
        ds[var_name].attrs = dict(units="MJ m-2 d-1",
                                  long_name="Surface solar radiation")
        ds[var_name].values = 0.0864*ds[var_name].values
        ds = ds.rename({var_name: "rayonnement"})
    if var_name == "wind":
        # from m s-1 to km h-1
        ds[var_name].attrs = dict(units="km h-1",
                                  long_name="Wind speed")
        ds[var_name].values = 3.6*ds[var_name].values
        ds = ds.rename({var_name: "vitesseVent"})
    if var_name == "tcc":
        ds[var_name].attrs = dict(units="0-1",
                                  long_name="Cloud cover")
        ds = ds.rename({var_name: "nebulosite"})
    if var_name == "vp":
        # from kPa to mmHg 7.50062
        ds[var_name].attrs = dict(units="mmHg",
                                  long_name="Vapor pressure")
        ds[var_name].values = 7.50062e-3*ds[var_name].values
        ds = ds.rename({var_name: "pression"})
    if var_name == "tmax" or var_name == "tmin":
        # from K to degC
        if var_name == "tmax":
            long_name = "Maximum daily air temperature"
            ds = ds.rename({var_name: "tMax"})
            var_name = "tMax"
        else:
            long_name = "Minimum daily air temperature"
            ds = ds.rename({var_name: "tMin"})
            var_name = "tMin"
        ds[var_name].attrs = dict(units="C",
                                  long_name=long_name)
        ds[var_name].values = ds[var_name].values - 273.15
    return ds
