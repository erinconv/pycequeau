from __future__ import annotations

from osgeo import gdal, ogr
import os
import json
from math import ceil, floor
import pandas as pd
import numpy as np
from src.physiographic import carreauxEntiers as CEs
from src.physiographic import carreauxPartiels as CPs
from src.physiographic import CPfishnet as CPfs
from src.core import utils as u
from src.core import projections as ceqproj
import geopandas as gpd
import sys


class Basin:
    def __init__(self,
                 project_folder: str,
                 basin_name: str,
                 file_list: list,
                 *args) -> None:
        # Create project structure
        self._project_path = project_folder
        self.name = basin_name
        self._project_structure(project_folder, file_list)
        # Create here the fishnet
        self._CEfishnet = os.path.join(
            self._project_path, "geographic", "CE_fishnet.shp")
        self._CPfishnet = os.path.join(
            self._project_path, "geographic", "CP_fishnet.shp")

        # Check if the bassin versant object is an input file
        if len(args) == 1:
            # Check if the provided file is an an actual json file
            try:
                # Open the json file file
                f = open(args[0], "r")
                self.bassinVersant = json.loads(f.read())
            except ValueError:
                raise ValueError("Provided file is not a json file")
        # else:
        #     sys.exit("The number of args must be 1")

    def _project_structure(self,
                           project_folder: str,
                           file_list: list):
        dirs = os.listdir(project_folder)
        if not "geographic" in dirs:
            os.mkdir(os.path.join(project_folder, "geographic"))
        if not "meteo" in dirs:
            os.mkdir(os.path.join(project_folder, "meteo"))
        if not "results" in dirs:
            os.mkdir(os.path.join(project_folder, "results"))
        # Set the file paths in the basin object
        self._set_files_paths(file_list)

    def _set_files_paths(self,
                         file_list: list):
        self._DEM = os.path.join(
            self._project_path, "geographic", file_list[0])
        self._FAC = os.path.join(
            self._project_path, "geographic", file_list[1])
        self._LC = os.path.join(
            self._project_path, "geographic", file_list[2])
        self._Basin = os.path.join(
            self._project_path, "geographic", file_list[3])
        self._SubBasins = os.path.join(
            self._project_path, "geographic", file_list[4])
        self._Waterbodies = os.path.join(
            self._project_path, "geographic", file_list[5])
        self._Wetlands = os.path.join(
            self._project_path, "geographic", file_list[6])

    def set_dimenssions(self, dx: float, dy: float):
        """_summary_

        Args:
            dx (float): _description_
            dy (float): _description_
        """
        ref_dataset = gdal.Open(self._DEM, gdal.GA_ReadOnly)
        transform = ref_dataset.GetGeoTransform()
        pixelWidth = transform[1]
        pixelHeight = -transform[5]
        # Get the number of rows and cols
        self.n_cols = abs(floor(dx/pixelWidth))
        self.n_rows = abs(floor(dy/pixelHeight))
        # Fix the position of the fisnet to match each pixel position
        self._dx = self.n_cols*pixelWidth
        self._dy = self.n_rows*pixelHeight

    def get_dimenssions(self):
        return [self._dx, self._dy]

    def get_EPSG(self):
        return self._epsg

    def _set_EPSG(self):
        self._epsg = ceqproj.get_proj_code(self._DEM)

    @classmethod
    def _create_CEfishnet(cls,
                          basin: str,
                          dx: float,
                          dy: float,
                          fishnet: str,
                          xoffset=0.0,
                          yoffset=0.0) -> None:
        """_summary_

        Args:
            basin (str): _description_
            dx (float): _description_
            dy (float): _description_
            fishnet (str): _description_
        """
        # Open the file from path
        watershed = ogr.Open(basin, gdal.GA_ReadOnly)
        lyr = watershed.GetLayer()
        xmin, xmax, ymin, ymax = lyr.GetExtent()

        # get rows
        rows = ceil((ymax-ymin)/dy)
        # get columns
        cols = ceil((xmax-xmin)/dx)

        # start grid cell envelope
        ringXleftOrigin = xmin - xoffset
        ringXrightOrigin = xmin + dx - xoffset
        ringYtopOrigin = ymax - yoffset
        ringYbottomOrigin = ymax - dy - yoffset

        # create output file
        outDriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(fishnet):
            os.remove(fishnet)
        outDataSource = outDriver.CreateDataSource(fishnet)
        outLayer = outDataSource.CreateLayer(
            fishnet, geom_type=ogr.wkbPolygon)
        featureDefn = outLayer.GetLayerDefn()
        idField = ogr.FieldDefn("CEid", ogr.OFTInteger)
        outLayer.CreateField(idField)
        # create grid cells
        countcols = 0
        id_count = 1
        while countcols < cols:
            countcols += 1

            # reset envelope for rows
            ringYtop = ringYtopOrigin
            ringYbottom = ringYbottomOrigin
            countrows = 0

            while countrows < rows:
                countrows += 1
                ring = ogr.Geometry(ogr.wkbLinearRing)
                ring.AddPoint(ringXleftOrigin, ringYtop)
                ring.AddPoint(ringXrightOrigin, ringYtop)
                ring.AddPoint(ringXrightOrigin, ringYbottom)
                ring.AddPoint(ringXleftOrigin, ringYbottom)
                ring.AddPoint(ringXleftOrigin, ringYtop)
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly = poly.Buffer(0)
                poly.AddGeometry(ring)

                # add new geom to layer
                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(poly)
                outFeature.SetField("CEid", int(id_count))
                outLayer.CreateFeature(outFeature)
                outFeature = None

                # new envelope for next poly
                ringYtop = ringYtop - dy
                ringYbottom = ringYbottom - dy

                id_count += 1

            # new envelope for next poly
            ringXleftOrigin = ringXleftOrigin + dx
            ringXrightOrigin = ringXrightOrigin + dx
        # Create prj file
        proj = lyr.GetSpatialRef()
        prj_path = fishnet.split(os.sep)
        prj_file_name = fishnet.split(os.sep)[-1].replace(".shp", ".prj")
        prj_path[-1] = prj_file_name
        prj_path = os.sep.join(prj_path)
        file = open(prj_path, 'w')
        file.write(proj.ExportToWkt())
        file.close()
        # Clean fishnet
        cls.clean_grids(basin, outDataSource)

    @classmethod
    def clean_grids(cls,
                    watershed: str,
                    CEfishnet: ogr.DataSource):
        """_summary_

        This function deletes the non used CE grid in the created fishnet
        https://gdal.org/doxygen/classOGRGeometry.html
        Args:
            watershed (str): _description_
            CEfishnet (ogr.DataSource): _description_
        """
        shp_watershed = ogr.Open(watershed, gdal.GA_ReadOnly)
        lyr1 = shp_watershed.GetLayer()
        lyr2 = CEfishnet.GetLayer()
        featList2 = range(lyr2.GetFeatureCount())
        feat1 = lyr1.GetFeature(0)
        geom1 = feat1.GetGeometryRef()
        for j in featList2:
            feat2 = lyr2.GetFeature(j)
            geom2 = feat2.GetGeometryRef()
            if not geom2.Within(geom1) and not geom2.Overlaps(geom1):
                lyr2.DeleteFeature(feat2.GetFID())

    @classmethod
    def clip_shps(cls,
                  sub_basins: str,
                  watershed: str,
                  newSubBasins: str):
        gdf1 = gpd.read_file(sub_basins)
        gdf2 = gpd.read_file(watershed)
        new_gdf = gpd.clip(gdf1, gdf2)
        # find and drop allpoints geometry if they exist
        new_gdf = new_gdf.where(new_gdf.geometry.geom_type != "Point")
        new_gdf = new_gdf.dropna()
        # Save the resulting GeoDataFrame as a shapefile
        new_gdf.to_file(newSubBasins)

    @classmethod
    def join_shps(cls,
                  CEfishnet: str,
                  SubBasins: str,
                  CPfishnet: str):
        # Read in the shapefiles
        gdf1 = gpd.read_file(CEfishnet)
        gdf2 = gpd.read_file(SubBasins)
        # Set the CAT Id in the subbasin dataset
        gdf2["CATid"] = range(1, len(gdf2)+1)
        # gdf1 = gdf1.explode()
        # gdf2 = gdf2.explode()
        # Union the shapefiles
        gdf_union = gpd.overlay(gdf1, gdf2, how='union')
        # This is for rasterize it later
        gdf_union["CPid"] = range(1, len(gdf_union)+1)
        if os.path.exists(CPfishnet):
            os.remove(CPfishnet)
        # Save the resulting GeoDataFrame as a shapefile
        gdf_union.to_file(CPfishnet)

    def create_CPfishnet(self):
        # Check if the file already exist
        if os.path.exists(self._CPfishnet):
            os.remove(self._CPfishnet)
        # Clip the sub basins with the watershed file
        self._newSubBasins = os.path.join(
            self._project_path, "geographic", "sub_basins.shp")
        self.clip_shps(self._SubBasins,
                       self._Basin,
                       self._newSubBasins)
        self._SubBasins = self._newSubBasins
        self.join_shps(self._CEfishnet,
                       self._SubBasins,
                       self._CPfishnet)

    def create_CEfishnet(self, xoffset=0.0, yoffset=0.0):
        # Check if the file already exist
        if os.path.exists(self._CEfishnet):
            os.remove(self._CEfishnet)

        self._create_CEfishnet(self._Basin,
                               self._dx,
                               self._dy,
                               self._CEfishnet,
                               xoffset,
                               yoffset)

    def polish_CPfishnet(self, area_th=0.05):
        """_summary_

        Args:
            area_th (float, optional): _description_. Defaults to 0.01.
            flow_th (float, optional): _description_. Defaults to 5e3.
        """
        out_name = os.path.join(
            self._project_path, "geographic", "CP_fishnet_test.shp")
        # Open fishnets as geodataframes
        CEfishnet = gpd.read_file(self._CEfishnet)
        CPfishnet = gpd.read_file(self._CPfishnet)
        # out_name = os.path.join(project_folder, "geographic", "CP_smallCP.shp")
        CPfishnet = CPfs.identify_small_CPs(CEfishnet, CPfishnet, area_th)
        CPfishnet, CEfishnet = CPfs.remove_border_CPs(
            CEfishnet, CPfishnet, self._FAC)
        CPfishnet = CPfs.remove_smallCP(CEfishnet, CPfishnet)
        CPfishnet = CPfs.dissolve_pixels(CEfishnet, CPfishnet, area_th)
        CPfishnet = CPfs.force_4CP(CEfishnet, CPfishnet, area_th)
        # Save the files with all the CP dissolved
        CPfishnet.to_file(self._CPfishnet)

    def CP_routing(self, flow_th=1000):
        # Create the CP grid from the merged shp file
        # This has the same dimensions as the FAC file
        # in order to compare them both to do the routing process
        # The CP grid is already processed and well polished
        CP_array = u.rasterize_shp(self._CPfishnet, self._FAC, "CPid")
        CE_array = u.rasterize_shp(self._CEfishnet, self._FAC, "CEid")

        # Open fishnets as geodataframes
        CEfishnet = gpd.read_file(self._CEfishnet)
        CPfishnet = gpd.read_file(self._CPfishnet)
        # Get the routing table. This table renames the CPs from the outlet
        # up to the last CP. Here the upstream CP are also identify for each
        # individual CP
        self.rtable, CPfishnet = CPfs.routing_table(CPfishnet,
                                                    CEfishnet,
                                                    self._FAC,
                                                    CP_array,
                                                    CE_array,
                                                    self.n_cols,
                                                    self.n_rows)
        # Obtain the downstream CP based on the previous process
        self.rtable = CPfs.get_downstream_CP(self.rtable)
        # self.rtable["newCPid"] = pd.to_numeric(self.rtable["newCPid"])
        # self.rtable["downstreamCPs"] = pd.to_numeric(self.rtable["downstreamCPs"])
        # Compute the route CP by CP
        self.outlet_routes = CPfs.outlet_routes(self.rtable)

        # Renumbering the fishnets
        CPfishnet, CEfishnet = CPfs.renumber_fishnets(
            CPfishnet, CEfishnet, self.rtable)
        # Compute cumulative percentage of surface area
        CPfishnet, upstreamCPs = CPfs.cumulative_areas(
            CPfishnet, CEfishnet, self.outlet_routes)
        # Compute the mean altitudes
        CPfishnet, CEfishnet = CPfs.mean_altitudes(
            CEfishnet, CPfishnet, self._DEM)
        # Add the table to the structure
        # self.rtable = rtable
        # self.outlet_routes = outlet_routes
        # Export the tables as csv into the geographical information
        # self.rtable.to_csv(os.path.join(self._project_path, "geographic", "rtable.csv"),index=False)

        np.savetxt("outlet_routes.csv",
                   self.outlet_routes, delimiter=",", fmt="%1i")
        # outlet_routes(,index=False)

        # Add rtable to the
        # out_name = os.path.join(
        #     self._project_path, "geographic", "CP_fishnet_test.shp")
        # out_name2 = os.path.join(
        #     self._project_path, "geographic", "CE_fishnet_test.shp")
        # if os.path.exists(out_name):
        #     os.remove(out_name)
        #     os.remove(out_name2)

        # Send the fishnets to the basin object to have it available
        self.CPfishnet = CPfishnet
        self.CEfishnet = CEfishnet
        # Save the files
        CPfishnet.to_file(self._CPfishnet)
        CEfishnet.to_file(self._CEfishnet)

    @classmethod
    def get_water_cover(cls,
                        shp_name: str,
                        shp_fishnet: gpd.GeoDataFrame,
                        ref_raster: str, att: str) -> list:
        # Open the waterbody as geopandas
        waterBodies = gpd.read_file(shp_name)
        # Clip the water shp with the fishnet
        clip_shp_water = gpd.clip(waterBodies, shp_fishnet)
        # Store the file in the folder
        # Create a new field to use it as reference raster attribute
        clip_shp_water["mask"] = 1
        # Save to a temporary file
        temp_shp_name = shp_name.replace(".shp", "2.shp")
        clip_shp_water.to_file(temp_shp_name)
        # Save shp as rasterizd tiff
        temp_tif_name = shp_name.replace(".shp", ".tif")
        u.rasterize_shp_as_byte(
            temp_shp_name, ref_raster, "mask", temp_tif_name)
        # Get the x,y vaues of each shp feature
        bounds = shp_fishnet["geometry"].bounds
        shp_fishnet = pd.concat([shp_fishnet, bounds], axis=1)
        ref_dataset = gdal.Open(ref_raster)
        # * This function needs to have the index from 0 to len(df). So, here I fix this isssue.
        # !!Do not change the function since it is also used by other processes in the previous procedures.
        shp_fishnet.index = range(len(shp_fishnet))
        shp_fishnet = CPfs.convert_coords_to_index(shp_fishnet, ref_dataset)
        df1 = pd.DataFrame(shp_fishnet.drop(columns='geometry'))
        # List for storing water bodies percentage
        water = []
        for index, _ in shp_fishnet.iterrows():
            # Get the feature
            feature = shp_fishnet.loc[index]
            raster_feature = u.rasterize_feature(feature, temp_tif_name, att)
            # Count the number of
            water.append(float(np.count_nonzero(
                raster_feature)/raster_feature.size*100.0))

        # Delete the temporary files
        os.remove(shp_name.replace(".shp", "2.shp"))
        os.remove(shp_name.replace(".shp", ".tif"))
        return water

    @classmethod
    def get_land_cover(cls,
                       gdf: gpd.GeoDataFrame,
                       LC: str, att: str) -> tuple:
        # Get the x,y vaues of each shp feature
        bounds = gdf["geometry"].bounds
        gdf = pd.concat([gdf, bounds], axis=1)
        # Open the reference dataset
        LC_dataset = gdal.Open(LC)
        # * This function needs to have the index from 0 to len(df). So, here I fix this isssue.
        # !!Do not change the function since it is also used by other processes in the previous procedures.
        gdf.index = range(len(gdf))
        gdf = CPfs.convert_coords_to_index(gdf, LC_dataset)
        # Loop into each CE
        pctForet = []
        pctSolNu = []
        for index, _ in gdf.iterrows():
            # Get the feature
            feature = gdf.loc[index]
            raster_feature = u.rasterize_feature(feature, LC, att)
            # Count the number of
            pctForet.append(float(np.count_nonzero((raster_feature >= 1) & (
                raster_feature <= 6))/raster_feature.size*100.0))
            pctSolNu.append(float(np.count_nonzero((raster_feature >= 7) & (
                raster_feature <= 13))/raster_feature.size*100.0))
        return pctForet, pctSolNu

    def carreauxEntiers_struct(self):
        # Create the CE grid with the shp dimenssions,
        # not with the reference raster dimensions to save memory and
        # make the processes faster
        CE_Array = self.create_CEgrid()
        # Find i,j coordinates
        coordinates = CEs.find_grid_coordinates(CE_Array)
        # Open the CE and fishnet shp
        # * Temporary. The real instruction is in the previous routine
        # CEfishnet = gpd.read_file(self._CEfishnet)
        # CPfishnet = gpd.read_file(self._CPfishnet)
        # self.CPfishnet = CPfishnet
        # self.CEfishnet = CEfishnet
        # Get the landcover dataset
        pctForet, pctSolNu = self.get_land_cover(self.CEfishnet,
                                                 self._LC,
                                                 "newCEid")
        # Check if there are different files for water bodies
        if self._Waterbodies == self._Wetlands:
            # Get the lakes
            pctLacRiviere = self.get_water_cover(self._Waterbodies,
                                                 self.CEfishnet,
                                                 self._DEM, "newCEid")
            # Get the same value to the marais
            pctMarais = pctLacRiviere.copy()
        else:
            pctLacRiviere = self.get_water_cover(self._Waterbodies,
                                                 self.CEfishnet,
                                                 self._DEM, "newCEid")
            # Get the marshes
            pctMarais = self.get_water_cover(self._Wetlands,
                                             self.CEfishnet,
                                             self._DEM, "newCEid")
        # Scale the percentages to make sure that they all sum up 100%
        CE_shp = np.array(pctLacRiviere) + np.array(pctMarais)
        CE_tif = np.array(pctForet) + np.array(pctSolNu)
        pctSolNu = (np.array(pctSolNu)/CE_tif)*(100-CE_shp)
        pctForet = (np.array(pctForet)/CE_tif)*(100-CE_shp)
        # Fix the nan values
        pctSolNu[np.isnan(pctSolNu)] = 0
        pctForet[np.isnan(pctForet)] = 0
        # Drop the non data CEs
        self.CEfishnet = self.CEfishnet[self.CEfishnet["newCEid"] != 0]
        # Append the i,j to the CE_fishnet
        self.CEfishnet = self.CEfishnet.reindex(
            columns=self.CEfishnet.columns.tolist() + ['i', 'j'])
        # Place the values into the dataset
        self.CEfishnet["i"] = coordinates["i"].values
        self.CEfishnet["j"] = coordinates["j"].values
        # Create the carreauxEntier dataset
        self.carreauxEntiers = pd.DataFrame(columns=["CEid", "i", "j", "pctLacRiviere",
                                                     "pctForet", "pctMarais", "pctSolNu",
                                                     "altitude"],
                                            index=coordinates.index,
                                            data=np.c_[coordinates["CEid"].values,
                                                       coordinates["i"].values,
                                                       coordinates["j"].values,
                                                       np.array(pctLacRiviere).astype(
                                                "float"),
            pctForet.astype("float"),
            np.array(pctMarais).astype(
                                                "float"),
            pctSolNu.astype("float"),
            self.CEfishnet["altitude"].values])
        # Change the values that must be integer types
        self.carreauxEntiers["CEid"] = self.carreauxEntiers["CEid"].astype(
            "uint16")
        self.carreauxEntiers["i"] = self.carreauxEntiers["i"].astype("uint8")
        self.carreauxEntiers["j"] = self.carreauxEntiers["j"].astype("uint8")
        self.CEfishnet.to_file(self._CEfishnet)
        self.carreauxEntiers.to_csv(os.path.join(
            self._project_path, "results", "carreauxEntiers.csv"))

    def carreauxPartiels_struct(self):
        # Start by sorting the values with in the dataframe
        self.CPfishnet = self.CPfishnet.sort_values(by="newCPid")
        # Open the CE and fishnet shp
        # * Temporary. The real instruction is in the previous routine
        # self.CEfishnet  = gpd.read_file(self._CEfishnet)
        # self.CPfishnet = gpd.read_file(self._CPfishnet)
        # Open the route table and outlet
        # # *This is temporaty too
        # self.rtable = pd.read_csv(os.path.join(self._project_path, "geographic", "rtable.csv"))
        # self.outlet_routes = np.genfromtxt(os.path.join(self._project_path, "geographic", "outlet_routes.csv"),delimiter=",")
        # Drop the first row
        # self.rtable = self.rtable.drop(index=0)
        # Get the landcover dataset
        # Find find the coordinates for each carreux partiel
        coordinates = CPs.get_CP_coordinates(self.carreauxEntiers,
                                             self.CPfishnet)
        # short script that gives the new CP code (ie 65,66,67,68) from each CE
        codes = CPs.get_codes(self.CPfishnet)
        # Get the landcover dataset
        pctForet, pctSolNu = self.get_land_cover(self.CPfishnet,
                                                 self._LC,
                                                 "newCPid")
        if self._Waterbodies == self._Wetlands:
            # Get the lakes
            pctLacRiviere = self.get_water_cover(self._Waterbodies,
                                                 self.CPfishnet,
                                                 self._DEM, "newCPid")
            pctMarais = pctLacRiviere.copy()
        else:
            # Get the lakes
            pctLacRiviere = self.get_water_cover(self._Waterbodies,
                                                 self.CPfishnet,
                                                 self._DEM, "newCPid")
            # Get the marshes
            pctMarais = self.get_water_cover(self._Wetlands,
                                             self.CPfishnet,
                                             self._DEM, "newCPid")
        # Scale the percentages to make sure that they all sum up 100%
        CP_shp = np.array(pctLacRiviere) + np.array(pctMarais)
        CP_tif = np.array(pctForet) + np.array(pctSolNu)
        pctSolNu = (np.array(pctSolNu)/CP_tif)*(100-CP_shp)
        pctForet = (np.array(pctForet)/CP_tif)*(100-CP_shp)
        # Fix the nan values (if applicable)
        pctSolNu[np.isnan(pctSolNu)] = 0
        pctForet[np.isnan(pctForet)] = 0
        # Cumulate the variables
        cumulates = CPs.cumulate_variables(
            self.outlet_routes, pctForet, pctLacRiviere, pctMarais)
        # Drop the non data CPs
        # self.CPfishnet = self.CPfishnet[self.CPfishnet["newCPid"] != 0]
        # Get river geometry
        geometry = CPs.get_river_geometry(self.CPfishnet, self.rtable)
        self.CPfishnet = self.CPfishnet.reindex(
            columns=self.CPfishnet.columns.tolist() + ['i', 'j'])
        # Place the values into the dataset
        self.CPfishnet["i"] = coordinates["i"].values
        self.CPfishnet["j"] = coordinates["j"].values
        # Create the carreauxPartiel dataset
        self.carreauxPartiels = pd.DataFrame(columns=["CPid", "i", "j", "code",
                                                      "pctSurface", "idCPAval", "idCPsAmont",
                                                      "idCE", "pctEau", "pctForet", "pctMarais",
                                                      "pctSolNu", "altitudeMoy", "profondeurMin",
                                                      "longueurCoursEauPrincipal", "largeurCoursEauPrincipal",
                                                      "penteRiviere", "cumulPctSuperficieCPAmont", "cumulPctSuperficieLacsAmont",
                                                      "cumulPctSuperficieMaraisAmont", "cumulPctSuperficieForetAmont"],
                                             index=coordinates.index,
                                             data=np.c_[coordinates["CPid"].values,
                                                        coordinates["i"].values,
                                                        coordinates["j"].values,
                                                        np.array(codes),
                                                        self.CPfishnet["pctSurface"].values,
                                                        self.rtable["downstreamCPs"].values,
                                                        self.rtable["upstreamCPs"].values,
                                                        self.CPfishnet["newCEid"].values,
                                                        np.array(
                                                            pctLacRiviere),
                                                        pctForet,
                                                        np.array(pctMarais),
                                                        pctSolNu,
                                                        self.CPfishnet["altitude"].values,
                                                        geometry["profondeurMin"].values,
                                                        geometry["longueurCoursEauPrincipal"].values,
                                                        geometry["largeurCoursEauPrincipal"].values,
                                                        geometry["penteRiviere"].values,
                                                        self.CPfishnet["cumulPctSurf"].values,
                                                        cumulates["cumulPctSuperficieLacsAmont"].values,
                                                        cumulates["cumulPctSuperficieMaraisAmont"].values,
                                                        cumulates["cumulPctSuperficieForetAmont"].values,
                                                        ])
        # Change the values that must be integer types
        self.carreauxPartiels["idCE"] = np.array(
            self.carreauxPartiels["idCE"], dtype=np.int16)
        self.carreauxPartiels["CPid"] = np.array(
            self.carreauxPartiels["CPid"], dtype=np.int16)
        self.carreauxPartiels["idCPAval"] = np.array(
            self.carreauxPartiels["idCPAval"], dtype=np.int16)
        self.carreauxPartiels["i"] = np.array(
            self.carreauxPartiels["i"], dtype=np.int8)
        self.carreauxPartiels["j"] = np.array(
            self.carreauxPartiels["j"], dtype=np.int8)
        self.carreauxPartiels["code"] = np.array(
            self.carreauxPartiels["code"], dtype=np.int8)
        self.carreauxPartiels["penteRiviere"] = np.array(
            self.carreauxPartiels["penteRiviere"], dtype=np.float32)

        self.carreauxPartiels.to_csv(os.path.join(
            self._project_path, "results", "carreauxPartiels.csv"))
        self.CPfishnet.to_file(self._CPfishnet)

    def create_bassinVersant_structure(self):
        # This structure will be stored as json format. This json format
        # will be easily translatet into .mat file for being read by Matlab
        # and also, will serve as one of the main input files in the OpenCEQUEAU
        # system.
        # *Temporary lines to read the csv files.
        # self.carreauxEntiers = pd.read_csv("carreauxEntiers.csv",index_col=0)
        # self.carreauxPartiels = pd.read_csv("carreauxPartiels.csv",index_col=0)
        # Get the column names from the carreux partiels and entiers structures
        columns_CE = self.carreauxEntiers.columns.tolist()
        columns_CP = self.carreauxPartiels.columns.tolist()
        # Create the dictionary structure
        self.bassinVersant = {
            "nbCpCheminLong": [],
            "superficieCE": [],
            "barrage": {},
            "nomBassinVersant": '',
            "carreauxEntiers": {},
            "carreauxPartiels": {}
        }
        # Send the carreuxEntiers values
        for CE_name in columns_CE:
            self.bassinVersant["carreauxEntiers"].update(
                {CE_name: self.carreauxEntiers[CE_name].values.tolist()})
        # Send the carreuxPartiels values
        for CP_name in columns_CP:
            self.bassinVersant["carreauxPartiels"].update(
                {CP_name: self.carreauxPartiels[CP_name].values.tolist()})
        self.bassinVersant["superficieCE"] = self._dx*self._dy*1.0e-6
        self.bassinVersant["nomBassinVersant"] = self.name
        self.bassinVersant["nbCpCheminLong"] = self.outlet_routes.shape[1]
        # Save the files in the results folder
        with open(os.path.join(self._project_path, "results", "bassinVersant.json"), "w") as outfile:
            json.dump(self.bassinVersant, outfile, indent=4)
        pass

    def create_CEgrid(self):
        # Default no data value of the CAT raster
        ref_raster = gdal.Open(self._DEM, gdal.GA_ReadOnly)
        CE_shp = ogr.Open(self._CEfishnet, gdal.GA_ReadOnly)
        lyr = CE_shp.GetLayer()
        # new_field = ogr.FieldDefn('CEid', ogr.OFTInteger)
        # lyr.CreateField(new_field)
        # for i, aa in enumerate(lyr):
        #     feature = lyr.GetFeature(i)
        #     feature.SetField("id2", i+1)
        #     lyr.SetFeature(feature)
        # Get the raster extent
        # xmin, xmax, ymin, ymax, xpixel, ypixel = u.GetExtent(self._DEM)
        # width = self._DEM.RasterXSize
        # height = self._DEM.RasterYSize
        # xmax = xmin + width*xpixel
        # ymin = ymax + height*ypixel
        proj = lyr.GetSpatialRef()
        xmin, xmax, ymin, ymax = lyr.GetExtent()
        # Get the resolution to the new raster
        x_res = ceil((xmax-xmin)/self._dx)
        y_res = ceil((ymax-ymin)/self._dy)
        # Get dimenssions for the CEgrid rasters

        # CEgrid path
        path = os.path.join(self._project_path, "geographic", "CEgrid.tif")
        self._CEgrid = gdal.GetDriverByName('GTiff').Create(
            path, abs(x_res), abs(y_res), 1, gdal.GDT_Int32)
        self._CEgrid.SetProjection(proj.ExportToWkt())
        # self._CEgrid.SetProjection(ref_raster.GetProjection())
        geo_transform = (xmin, self._dx, 0.0, ymax, 0.0, -self._dy)
        self._CEgrid.SetGeoTransform(geo_transform)
        # band.FlushCache()
        band = self._CEgrid.GetRasterBand(1)
        band.SetNoDataValue(0)
        gdal.RasterizeLayer(self._CEgrid, [1], lyr, options=[
                            "ATTRIBUTE=newCEid"])
        band = self._CEgrid.GetRasterBand(1)
        data_array = band.ReadAsArray()
        return data_array

    def create_CPgrid(self):
        # Get the CE and CP hps
        CE_shp = ogr.Open(self._CEfishnet, gdal.GA_ReadOnly)
        CP_shp = ogr.Open(self._CPfishnet, gdal.GA_ReadOnly)
        lyr_CE = CE_shp.GetLayer()
        lyr_CP = CP_shp.GetLayer()
        # Get the CP fishnet projection
        proj = lyr_CP.GetSpatialRef()
        # Here we work  with the CE extention to make sure that both shps have the same extent
        xmin, xmax, ymin, ymax = lyr_CE.GetExtent()
        # Get the resolution to the new raster
        scale = 30
        x_res = ceil((xmax-xmin)/self._dx)*scale
        y_res = ceil((ymax-ymin)/self._dy)*scale
        # Make the Union between the two shps
        path = os.path.join(self._project_path, "geographic", "CPgrid.tif")
        self._CPgrid = gdal.GetDriverByName('GTiff').Create(
            path, x_res, y_res, 1, gdal.GDT_Int32)
        # Seth projection info to the CP tif
        self._CPgrid.SetProjection(proj.ExportToWkt())
        self._CPgrid.SetGeoTransform(
            (xmin, self._dx/scale, 0, ymin, 0, self._dy/scale))
        band = self._CPgrid.GetRasterBand(1)
        band.SetNoDataValue(0)
        gdal.RasterizeLayer(self._CPgrid, [1], lyr_CP, options=[
                            "ATTRIBUTE=newCPid"])
        return band.ReadAsArray()
