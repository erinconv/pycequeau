from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr
from pycequeau.core import projections
from pycequeau.test import tutorial
from pycequeau.physiographic.base import Basin
from pycequeau.core import utils as u
from pycequeau.meteo.meteo_netcdf import StationNetCDF


def main():
    tif_path = "/mnt/c/Users/erinc/Documents/01-PhD/00-Geographic/Melezes/CEgrid.tif"
    watershet_path = "/mnt/c/Users/erinc/Documents/01-PhD/00-Geographic/Melezes/Basin2.shp"
    fishnet_path = "/mnt/c/Users/erinc/Documents/01-PhD/00-Geographic/Melezes/Fishnet_Test.shp"

    # Path where nc files are stored - Historicals
    hist_path = "/mnt/c/Users/erinc/Documents/01-PhD/10-CMIP/04-BasinData/ERA"
    cordex_path = "/mnt/c/Users/erinc/Documents/01-PhD/10-CMIP/04-BasinData/CORDEX/BiasCorrected"

    #

    # Create basin Object
    DEM = gdal.Open(tif_path, gdal.GA_ReadOnly)
    FAC = gdal.Open(tif_path, gdal.GA_ReadOnly)
    CEgrid = gdal.Open(tif_path, gdal.GA_ReadOnly)
    watershed = ogr.Open(watershet_path, gdal.GA_ReadOnly)
    fishnet = ogr.Open(fishnet_path, gdal.GA_Update)
    grid_size = 8000
    # u.rasterize_fishnet(CEgrid, fishnet,8000)
    Melezes = Basin("Melezes", DEM, FAC, CEgrid, fishnet, watershed, grid_size)
    ds = tutorial.open_example("precipitation")
    # MeteoStations = StationNetCDF
    MeteoStations = StationNetCDF.charge_CORDEX_Meteo(cordex_path)
    
    # MeteoStations = StationNetCDF.charge_ERA_Meteo(hist_path)
    MeteoStations.charge_Basin(Melezes)
    dsi= MeteoStations.interpolation("nearest")
    grid = MeteoStations.cequeau_grid(dsi)
    
    # df = pd.DataFrame({"pasTemp":},
    #                       columns=CEs)
    # vitesseVent, 
# tMax, 
# tMin, 
# pTot,
# rayonnement
# nebulosite,
# pression
    
    # dsi = dsi.isel(time=500)
    # ds = ds.isel(time=500)
    # fig, axes = plt.subplots(ncols=2,
    #                          nrows=2,
    #                          figsize=(12, 10))
    # dsi["tmax"].plot(ax=axes[0,0], cmap='Spectral')
    # # axes[0].set_title("Raw data")
    # dsi["tmin"].plot(ax=axes[0,1], cmap='Spectral')
    # dsi["ssr"].plot(ax=axes[1,0], cmap='Spectral')
    # # axes[0].set_title("Raw data")
    # dsi["CE"].plot(ax=axes[1,1], cmap='Spectral')
    # # axes[1].set_title("Interpolated data")
    
    # plt.savefig("interpolated_nearest.png")
    # plt.show()
    # plt.show()
    
    
    # print(ds.i)
    # print(dsi.i)
    # print(dsi)
    # Initialize station object
    # MeteoStations = StationNetCDF(ds,Melezes,CEgrid)
    # first_var = MeteoStations.interpolation("nearest")

    # Meteo = MeteoStation("Melezes",
    #                      grid_size,
    #                      DEM=DEM,
    #                      CEgrid=CEgrid,
    #                      watershed=watershed,
    #                      fishnet=fishnet,
    #                      ds=ds)

    # MeteoStations.stattion_table()

    # print(tif_object,fishnet,grid_size)
    # MeteoStations.stations(CEgrid,fishnet,grid_size)
    # u.loop_zonal_stats(CEgrid,fishnet)
    # print(MeteoStations.altitudes)
    # print(MeteoStations.lon_mask[MeteoStations.lon_idx])
    # print(MeteoStations.lat_mask[MeteoStations.lat_idy])
    # latlon = np.c_[MeteoStations.lat_mask[MeteoStations.lat_idy],
    #                MeteoStations.lon_mask[MeteoStations.lon_idx]]
    # with open("grid_points.csv", 'w') as f:
    #     f.write("x,y \n")
    #     np.savetxt("grid_points.csv",MeteoStations.stations_table(),delimiter=",")
    # print(MeteoStations.watershed_mask)
    #
    # projections.get_proj_code(tif_object)

    # print(type(shp))
    # print("*************************")
    # tif_path = "/mnt/c/Users/erinc/Documents/01-PhD/00-Geographic/Melezes/CEgrid.tif"
    # tif_object = gdal.Open(tif_path)
    # projections.get_proj_json(tif_object)
    # # print(masked_sample)
    # # dsi.values =

    #
    # plt.show()
    # print(dsi)
    # ds_new =

    # print(ds["lon"].values -360)
    # u.find_nearest()

    # ds["lon"].values = ds["lon"].values +correct

    # print(ds["lon"].values)
    # ds = ds.sel(lat=slice(table["lat"].max(),table["lat"].min()),
    #             lon=slice(table["lon"].max()-correct,table["lon"].min()-correct))
    # print(table["lon"]-correct)
    # lat_ref = ds["lat"].values
    # lon_ref = ds["lon"].values
    # print(lon_ref)

    # print(ig)


if __name__ == "__main__":
    main()
