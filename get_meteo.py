from pycequeau.physiographic.base import Basin
from pycequeau.meteo.meteo_netcdf import StationNetCDF
import os


def main():
    # Create the project folder structure.
    # 1- Select folder where the meteo data is stored
    project_folder = "/home/erinconv/01-PhD/River"
    # 2- Select the bassinVersant json file for this basin
    bassinVersant_file = "/home/erinconv/01-PhD/River/results/bassinVersant.json"
    # 3- File list with all the raster and shapes
    files_list = ["Dem1.tif",
                  "FAC_grass.tif",
                  "LC_grass2.tif",
                  "Watershed_GRASS.shp",
                  "sub_basins_grass_clip.shp",
                  "WaterBodies.shp",
                  "WaterBodies.shp"]
    # 4- Create basin Object:
    basin = Basin(project_folder,
                  "Melezes",
                  files_list,
                  bassinVersant_file)
    # 5- Create the Meteo data object
    MeteoStations = StationNetCDF.charge_ERA_Meteo(basin,
                                                   os.path.join(project_folder,"meteo","ERA"))
    # 6- Do the interpolation
    dsi = MeteoStations.interpolation("nearest")
    # 7- Construct the meteo structure for CEQUEAU
    grid = MeteoStations.cequeau_grid(dsi,basin)
    # 8- Save the netcdf with the meteo data
    grid.to_netcdf(os.path.join(project_folder,"meteo","meteo_cequeau.nc"))
    pass
if __name__ == "__main__":
    main()
