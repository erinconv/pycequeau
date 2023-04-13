from pycequeau.physiographic.base import Basin
import os


def main():
    # Create the project folder structure.
    # 1- Select folder
    project_folder = "/home/erinconv/01-PhD/River"
    files_list = ["Dem1.tif",
                  "FAC_grass.tif",
                  "LC_grass2.tif",
                  "Watershed_GRASS.shp",
                  "sub_basins_grass_clip.shp",
                  "WaterBodies.shp",
                  "WaterBodies.shp"]
    # 2- Create basin Object:
    basin = Basin(project_folder,
                  "Melezes",
                  files_list)
    # 2.1 - Select Fisnet dimensions
    basin.set_dimenssions(7500, 7500)
    # 3 Create CE and CP fishnet
    basin.create_CEfishnet()
    basin.create_CPfishnet()
    # 4 - Remove the small CPs in the basin
    basin.polish_CPfishnet()
    # 5 - Do the routing process
    basin.CP_routing()
    # 6 - Create the carreaux entiers structure
    basin.carreauxEntiers_struct()
    # 7 - Create carreux partiels structure
    basin.carreauxPartiels_struct()
    # 8 - Create the bassinVersant structure
    basin.create_bassinVersant_structure()

if __name__ == "__main__":
    main()
