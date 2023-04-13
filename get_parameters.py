from pycequeau.physiographic.base import Basin
from pycequeau.simulations.parameters import Parameters
from pycequeau.simulations._param_examples import send_values_test
import os
import json


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
    # 5- Create parameter object
    params = Parameters(basin)
    # 6- retrieve the example parameter values to populate the structure
    flow_parameters, evapo_parameters, initial_conditions, snow_parameters, simulation_options, transferts, temperature_params = send_values_test()
    # 7- set the different parameter values
    params.set_sol(flow_parameters)
    params.set_solinitial(initial_conditions)
    params.set_transfert(transferts)
    params.set_option(simulation_options)
    params.set_fonte(snow_parameters,1)
    params.set_evapo(evapo_parameters,1)
    params.set_qualite(temperature_params)
    # 8- Create the parameter structure
    params.create_parameter_structure()



if __name__ == "__main__":
    main()
