from __future__ import annotations
import os
import json
from math import ceil, floor
from osgeo import gdal, ogr
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterstats as rs
import matplotlib.pyplot as plt
import rasterio
from src.physiographic import carreauxEntiers as CEs
from src.physiographic import carreauxPartiels as CPs
from src.physiographic import CPfishnet as CPfs
from src.core import utils as u
from src.core import projections as ceqproj
# import sys
# import rasterio
# from rasterio.features import shapes
# from shapely.geometry import shape


class Basin:
    """_summary_
    """

    def __init__(self,
                 project_folder: str,
                 basin_name: str,
                 file_list: list,
                 *args) -> None:
        """_summary_

        Args:
            project_folder (str): _description_
            basin_name (str): _description_
            file_list (list): _description_

        Raises:
            ValueError: _description_
        """
        # Create project structure
        self._project_path = project_folder
        self.name = basin_name
        self.n_cols = None
        self.n_rows = None
        self._dx = None
        self._dy = None
        self._dx_o = None
        self._dy_o = None
        self.pixel_width = None
        self.pixel_height = None
        self.rtable = None
        self.outlet_routes = None
        self._Basin = None
        self._newSubBasins = None
        self._SubBasins = None
        # Create project structure
        self._project_structure(file_list)
        # Create here the fishnet
        self._CEfishnet = os.path.join(
            self._project_path, "geographic", "CE_fishnet.shp")
        self._CPfishnet = os.path.join(
            self._project_path, "geographic", "CP_fishnet.shp")
        self._CEgrid = os.path.join(
            self._project_path, "geographic", "CEgrid.tif")

        # Check if the bassin versant object is an input file
        if len(args) == 1:
            # Check if the provided file is an an actual json file
            try:
                # Open the json file file
                f = open(args[0], "r", encoding='utf-8')
                self.bassinVersant = json.loads(f.read())
            # except ZeroDivisionError as e
            except ValueError as exc:
                raise ValueError(
                    "Provided file for bassinVerssant structure is not a json file") from exc
        # else:
        #     # Polygonize and process the raster watershed and subbasins file
        #     self.rasterize_maps()

    @property
    def Basin(self):
        return self._Basin

    @property
    def DEM(self):
        return self._DEM

    @property
    def get_CEgrid(self):
        return self._CEgrid

    @property
    def cp_fishnet_name(self):
        return self._CPfishnet

    @property
    def ce_fishnet_name(self):
        return self._CEfishnet

    @property
    def project_path(self):
        return self._project_path

    @property
    def CPfishnet(self):
        return self._CPfihsnet_gdf

    @CPfishnet.setter
    def CPfishnet(self, gdf: gpd.GeoDataFrame):
        self._CPfihsnet_gdf = gdf

    @property
    def CEfishnet(self):
        return self._CEfihsnet_gdf

    @CEfishnet.setter
    def CEfishnet(self, gdf: gpd.GeoDataFrame):
        self._CEfihsnet_gdf = gdf

    @property
    def get_dimenssions(self):
        return [self._dx, self._dy]

    @property
    def get_EPSG(self):
        return self._epsg

    @get_EPSG.setter
    def get_EPSG(self):
        self._epsg = ceqproj.get_proj_code(self._DEM)

    def _project_structure(self, file_list: str):
        """_summary_

        Args:
            file_list (str): _description_
        """
        dirs = os.listdir(self._project_path)
        if not "geographic" in dirs:
            os.mkdir(os.path.join(self._project_path, "geographic"))
        if not "meteo" in dirs:
            os.mkdir(os.path.join(self._project_path, "meteo"))
        if not "results" in dirs:
            os.mkdir(os.path.join(self._project_path, "results"))
        # Set the file paths in the basin object
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
        # self._Waterbodies = os.path.join(
        #     self._project_path, "geographic", file_list[5])
        # self._Wetlands = os.path.join(
        #     self._project_path, "geographic", file_list[6])

    def rasterize_maps(self):
        """_summary_
        """
        # Polygonize both raster
        u.polygonize_raster(self._Basin)
        u.polygonize_raster(self._SubBasins)
        self._Basin = self._Basin.replace(".tif", ".shp")
        self._SubBasins = self._SubBasins.replace(".tif", ".shp")
        # Open shps as geopandas datasets
        watershed = gpd.read_file(self._Basin)
        sub_basins = gpd.read_file(self._SubBasins)
        # Fix geometry
        watershed = u.fix_geometry(watershed)
        sub_basins = u.fix_geometry(sub_basins)
        # Drop non desired values
        watershed = watershed.where(watershed["raster_val"] != 255)
        sub_basins = sub_basins.where(sub_basins["raster_val"] != 255)
        watershed["area"] = watershed.area
        # Select the maximum value
        watershed = watershed.where(
            watershed["area"] == watershed["area"].max())
        watershed = watershed.dropna()
        sub_basins = sub_basins.dropna()
        watershed.to_file(self._Basin)
        sub_basins.to_file(self._SubBasins)

    def set_dimensions(self, dx: float, dy: float):
        """_summary_

        Args:
            dx (float): _description_
            dy (float): _description_
        """
        self._dx_o = dx
        self._dy_o = dy
        ref_dataset = gdal.Open(self._DEM, gdal.GA_ReadOnly)
        transform = ref_dataset.GetGeoTransform()
        self.pixel_width = transform[1]
        self.pixel_height = -transform[5]
        # Get the number of rows and cols
        self.n_cols = abs(floor(dx/self.pixel_width))
        self.n_rows = abs(floor(dy/self.pixel_height))
        # Fix the position of the fisnet to match each pixel position
        self._dx = self.n_cols*self.pixel_width
        self._dy = self.n_rows*self.pixel_height

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
        file = open(prj_path, 'w', encoding='utf-8')
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
        """_summary_

        Args:
            sub_basins (str): _description_
            watershed (str): _description_
            newSubBasins (str): _description_
        """
        gdf1 = gpd.read_file(sub_basins)
        gdf2 = gpd.read_file(watershed)
        new_gdf = gpd.clip(gdf1, gdf2)
        # find and drop all non polygon geometries if they exist
        new_gdf = new_gdf.where(new_gdf.geometry.geom_type != "Point")
        new_gdf = new_gdf.where(
            new_gdf.geometry.geom_type != "MultiLineString")
        new_gdf = new_gdf.dropna()
        # Save the resulting GeoDataFrame as a shapefile
        new_gdf.to_file(newSubBasins)

    @classmethod
    def join_shps(cls,
                  CEfishnet: str,
                  SubBasins: str,
                  CPfishnet: str):
        """_summary_

        Args:
            CEfishnet (str): _description_
            SubBasins (str): _description_
            CPfishnet (str): _description_
        """
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
        """_summary_
        """
        # Check if the file already exist
        # if os.path.exists(self._CPfishnet):
        #     os.remove(self._CPfishnet)
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
        """_summary_

        Args:
            xoffset (float, optional): _description_. Defaults to 0.0.
            yoffset (float, optional): _description_. Defaults to 0.0.
        """
        # Rasterize maps in this step
        self.rasterize_maps()
        # Check if the file already exist
        # if os.path.exists(self._CEfishnet):
        #     os.remove(self._CEfishnet)

        self._create_CEfishnet(self._Basin,
                               self._dx,
                               self._dy,
                               self._CEfishnet,
                               xoffset,
                               yoffset)

    def polish_CPfishnet(self, area_th=0.05):
        """_summary_

        Args:
            area_th (float, optional): _description_. Defaults to 0.05.
            flow_th (float, optional): _description_. Defaults to 5e3.
        """
        # out_name = os.path.join(
        #     self._project_path, "geographic", "CP_fishnet_test.shp")
        # Open fishnets as geodataframes
        CEfishnet = gpd.read_file(self._CEfishnet)
        CPfishnet = gpd.read_file(self._CPfishnet)
        # out_name = os.path.join(project_folder, "geographic", "CP_smallCP.shp")
        CPfishnet = CPfs.identify_small_CPs(CEfishnet, CPfishnet, area_th)
        CPfishnet, CEfishnet = CPfs.remove_border_CPs(
            CEfishnet, CPfishnet, self._FAC)
        CPfishnet = CPfs.remove_smallCP(CEfishnet, CPfishnet, area_th)
        # TODO: This is a temporary solution to drop all the pixel-size features in the CP fishnet. Needs to be corrected in the future
        pixel_size_area = self.pixel_height*self.pixel_width
        # Find the areas this size
        CPfishnet = CPfs.force_4CP(CEfishnet, CPfishnet, area_th)
        CPfishnet = CPfs.dissolve_pixels(CEfishnet, CPfishnet, pixel_size_area)
        # CPfishnet = CPfs.dissolve_pixels(CEfishnet, CPfishnet, pixel_size_area)
        # df = pd.DataFrame(CPfishnet.drop(columns='geometry'))
        # Save the files with all the CP dissolved
        CPfishnet.to_file(self._CPfishnet)
        return None

    def CP_routing(self):
        """
        Create the CP grid from the merged shp file
        This has the same dimensions as the FAC file
        in order to compare them both to do the routing process
        The CP grid is already processed and well polished
        """
        CP_array = u.rasterize_shp(self._CPfishnet, self._FAC, "CPid")
        CE_array = u.rasterize_shp(self._CEfishnet, self._FAC, "CEid")
        self.CEfishnet = gpd.read_file(self._CEfishnet)
        self.CPfishnet = gpd.read_file(self._CPfishnet)
        # Open fishnets as geodataframes
        # Get the routing table. This table renames the CPs from the outlet
        # up to the last CP. Here the upstream CP are also identify for each
        # individual CP
        self.routing = CPfs.routing_table(self.CPfishnet,
                                          self.CEfishnet,
                                          self._FAC,
                                          CP_array,
                                          CE_array,
                                          self.n_cols,
                                          self.n_rows)
        # Get the rtable
        self.rtable, self.CPfishnet = CPfs.get_rtable(self.CPfishnet,
                                                      self.routing)
        # Renumbering the fishnets
        self.CPfishnet, self.CEfishnet = CPfs.renumber_fishnets(self.CPfishnet,
                                                                self.CEfishnet,
                                                                self.rtable,
                                                                self.routing)
        # Obtain the downstream CP based on the previous process
        self.rtable = CPfs.get_downstream_CP(self.rtable)
        # self.rtable["newCPid"] = pd.to_numeric(self.rtable["newCPid"])
        # self.rtable["downstreamCPs"] = pd.to_numeric(self.rtable["downstreamCPs"])
        # Compute the route CP by CP
        # CPfishnet.to_file(self._CPfishnet)
        self.outlet_routes = CPfs.outlet_routes(self.rtable)

        # Compute cumulative percentage of surface area
        self.CPfishnet, upstream_cps = CPfs.cumulative_areas(
            self.CPfishnet, self.CEfishnet, self.outlet_routes)
        # Compute the mean altitudes
        self.CPfishnet, self.CEfishnet = CPfs.mean_altitudes(
            self.CEfishnet, self.CPfishnet, self._DEM)
        # Add the table to the structure
        # self.rtable = rtable
        # self.outlet_routes = outlet_routes
        # Export the tables as csv into the geographical information
        # Save the outlet routes files to have it on hand
        outlet_file_name = os.path.join(
            self._project_path, "results", "outlet_routes.csv")
        np.savetxt(outlet_file_name,
                   self.outlet_routes,
                   delimiter=",", fmt="%1i")
        # Save the upstream CPs file
        upstream_cps_name = os.path.join(
            self._project_path, "results", "upstream_cps.csv")
        np.savetxt(upstream_cps_name,
                   upstream_cps,
                   delimiter=",", fmt="%1i")
        # Save the rtable file
        rtable_name = os.path.join(
            self._project_path, "results", "rtable.csv")
        self.rtable.to_csv(rtable_name)

        # Save the routing file
        routing_name = os.path.join(
            self._project_path, "results", "routing.csv")
        self.routing.to_csv(routing_name)
        # Send the fishnets to the basin object to have it available

        # Save the files
        self.CPfishnet['newCPid'] = self.CPfishnet['newCPid'].astype('int')
        self.CEfishnet['newCEid'] = self.CEfishnet['newCEid'].astype('int')
        self.CPfishnet.to_file(self._CPfishnet)
        self.CEfishnet.to_file(self._CEfishnet)
        # Exprot the files as tif
        self._save_tif_ce_fishnet()

    @classmethod
    def get_water_cover(cls,
                        shp_name: str,
                        gdf: gpd.GeoDataFrame,
                        ref_raster: str,
                        ncols: int,
                        nrows: int) -> list:
        """_summary_

        Args:
            shp_name (str): Water bodies shape file name
            shp_fishnet (gpd.GeoDataFrame): fishnet shape file
            ref_raster (str): reference raster to get the size of the pixel
            ncols (int): _description_
            nrows (int): _description_

        Returns:
            list: _description_
        """
        # Open the waterbody as geopandas
        waterBodies = gpd.read_file(shp_name)
        waterBodies = u.fix_geometry(waterBodies)
        # Clip the water shp with the fishnet
        clip_shp_water = gpd.clip(waterBodies, gdf)
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
        # Compute the zonal sats with the rasterized file
        zonal_stats = rs.zonal_stats(gdf.geometry,
                                     temp_tif_name,
                                     categorical=True)
        # Compute the percentages
        water_df = pd.DataFrame(zonal_stats)
        water_df = water_df.fillna(0)
        water_df = water_df.iloc[:, 0]
        size_ce = ncols*nrows
        water = water_df.values/size_ce*100
        # Delete temporary files
        os.remove(shp_name.replace(".shp", "2.shp"))
        os.remove(shp_name.replace(".shp", ".tif"))
        if np.amax(water) > 100:
            assert False, "Something whent wrong"

        return water

    @classmethod
    def get_land_cover(cls,
                       gdf: gpd.GeoDataFrame,
                       LC_name: str,
                       attr: str,) -> tuple:
        """_summary_

        Args:
            gdf (gpd.GeoDataFrame): _description_
            LC_name (str): _description_
            ncols (int): _description_
            nrows (int): _description_

        Returns:
            tuple: _description_
        """
        #
        gdf = gdf.sort_values(by=attr)
        zonal_stats = rs.zonal_stats(gdf.geometry,
                                     LC_name,
                                     categorical=True)
        # Get the pixel size of the land cover raster
        raster = gdal.Open(LC_name)
        gt = raster.GetGeoTransform()
        pixelSizeX = gt[1]
        pixelSizeY = gt[5]
        # Get the the coefficient to convert from number of pixels to area in m2
        pixel_area_coeff = abs(pixelSizeX*pixelSizeY)
        # Find the values of each land category
        df_stats = pd.DataFrame(zonal_stats)
        col_vals = df_stats.columns.values
        df_stats = df_stats.loc[:, col_vals]
        df_stats = df_stats.fillna(0)
        # Get the indexes of the forest and bare soil
        idx_forest = (col_vals >= 1) & (col_vals <= 6)
        idx_sol_nu = ((col_vals >= 7) & (col_vals <= 13)) | ((
            col_vals >= 15) & (col_vals <= 17))
        idx_wetland = col_vals == 14
        # Label 0 is water that is not in the land (mostly sea water and stuaries)
        idx_water = (col_vals == 0) | (col_vals == 18)
        # Get the impermeable surface
        idx_tri_s = col_vals == 16
        forest_df = df_stats.loc[:, idx_forest]
        sol_nu_df = df_stats.loc[:, idx_sol_nu]
        wetlands_df = df_stats.loc[:, idx_wetland]
        water_df = df_stats.loc[:, idx_water]
        tri_df = df_stats.loc[:, idx_tri_s]
        # Compute the percentage of forest and bare soil
        size_carreaux = gdf.area.values
        pctForet = pixel_area_coeff * \
            forest_df.sum(axis=1).values/size_carreaux*100
        pctSolNu = pixel_area_coeff * \
            sol_nu_df.sum(axis=1).values/size_carreaux*100
        pctWater = pixel_area_coeff * \
            water_df.sum(axis=1).values/size_carreaux*100
        pctWetlands = pixel_area_coeff * \
            wetlands_df.sum(axis=1).values/size_carreaux*100
        pctImpermeable = pixel_area_coeff * \
            tri_df.sum(axis=1).values/size_carreaux*100
        return pctForet, pctSolNu, pctWater, pctWetlands, pctImpermeable

    def carreauxEntiers_struct(self) -> None:
        """_summary_
        """
        # Create the CE grid with the shp dimenssions,
        # not with the reference raster dimensions to save memory and
        # make the processes faster
        ds = gdal.Open(self._CEgrid, gdal.GA_ReadOnly)
        ce_array = np.array(ds.GetRasterBand(1).ReadAsArray())
        # Find i,j coordinates
        coordinates = CEs.find_grid_coordinates(ce_array)
        # Get the landcover dataset
        self.CEfishnet.index = self.CEfishnet["newCEid"].values
        pctForet, pctSolNu, pctWater, pctWetlands, pctImpermeable = self.get_land_cover(self.CEfishnet,
                                                                                        self._LC,
                                                                                        "newCEid")
        # Scale the percentages to make sure that they all sum up 100%
        total = pctForet + pctSolNu + pctWater + pctWetlands  # This is 100%
        pctForet = np.divide(pctForet, total)*100
        pctSolNu = np.divide(pctSolNu, total)*100
        pctWater = np.divide(pctWater, total)*100
        pctWetlands = np.divide(pctWetlands, total)*100
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
                                                     "altitude", "pctImpermeable"],
                                            index=coordinates.index,
                                            data=np.c_[coordinates["CEid"].values,
                                                       coordinates["i"].values,
                                                       coordinates["j"].values,
                                                       pctWater,
                                                       pctForet,
                                                       pctWetlands,
                                                       pctSolNu,
                                                       self.CEfishnet["altitude"].values,
                                                       pctImpermeable]
                                            )
        # Change the values that must be integer types
        self.carreauxEntiers["CEid"] = self.carreauxEntiers["CEid"].astype(
            "uint16")
        self.carreauxEntiers["i"] = self.carreauxEntiers["i"].astype("uint8")
        self.carreauxEntiers["j"] = self.carreauxEntiers["j"].astype("uint8")
        self.CEfishnet.to_file(self._CEfishnet)
        self.carreauxEntiers.to_csv(os.path.join(
            self._project_path, "results", "carreauxEntiers.csv"))

    def carreauxPartiels_struct(self) -> None:
        """_summary_
        """
        # Start by sorting the values with in the dataframe
        self.CPfishnet = self.CPfishnet.sort_values(by="newCPid")
        self.CPfishnet.index = self.CPfishnet["newCPid"].values
        # Find find the coordinates for each carreux partiel
        coordinates = CPs.get_CP_coordinates(self.carreauxEntiers,
                                             self.CPfishnet)
        # short script that gives the new CP code (ie 65,66,67,68) from each CE
        codes = CPs.get_codes(self.CPfishnet)
        # Get the landcover dataset
        pctForet, pctSolNu, pctWater, pctWetlands, pctImpermeable = self.get_land_cover(self.CPfishnet,
                                                                                        self._LC,
                                                                                        "newCPid")
        # Scale the percentages to make sure that they all sum up 100%
        total = pctForet + pctSolNu + pctWater + pctWetlands  # This is 100%
        remain = 100 - total
        remain[remain < 0] = 0
        pctForet = np.divide(pctForet, total)*100
        pctSolNu = np.divide(pctSolNu, total)*100
        pctWater = np.divide(pctWater, total)*100
        pctWetlands = np.divide(pctWetlands, total)*100
        # Scale the percentages to make sure that they all sum up 100%
        # Cumulate the variables
        cumulates = CPs.cumulate_variables(
            self.outlet_routes, pctForet, pctWater, pctWetlands)
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
                                                      "cumulPctSuperficieMaraisAmont", "cumulPctSuperficieForetAmont",
                                                      "cumulArea", "pctImpermeable"],
                                             index=coordinates.index,
                                             data=np.c_[coordinates["CPid"].values,
                                                        coordinates["i"].values,
                                                        coordinates["j"].values,
                                                        np.array(codes),
                                                        self.CPfishnet["pctSurface"].values,
                                                        self.rtable["downstreamCPs"].values,
                                                        self.rtable["upstreamCPs"].values,
                                                        self.CPfishnet["newCEid"].values,
                                                        pctWater,
                                                        pctForet,
                                                        pctWetlands,
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
                                                        self.CPfishnet["cumulArea"].values,
                                                        pctImpermeable
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
        """_summary_
        """
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

    def plot_routing(self, area_th=0.01):
        """_summary_

        Args:
            area_th (float, optional): _description_. Defaults to 0.01.
        """
        # Centroinds
        CP_fishnet = gpd.read_file(self._CPfishnet)
        cp_struct_name = os.path.join(
            self._project_path, "results", "carreauxPartiels.csv")
        carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
        CP_fishnet["x_c"] = CP_fishnet.centroid.x
        CP_fishnet["y_c"] = CP_fishnet.centroid.y
        # Define the pair of coordinates
        p1 = np.zeros((len(CP_fishnet), 2))
        p2 = p1.copy()
        u_c = np.zeros(len(CP_fishnet))
        v_c = u_c.copy()

        fig, ax = plt.subplots(nrows=1, ncols=1, layout='constrained')
        CP_fishnet.plot(ax=ax, column="cumulArea",)
        max_area = CP_fishnet["cumulArea"].max()
        idx_small_area = CP_fishnet["cumulArea"] > area_th*max_area
        CP_fishnet.index = CP_fishnet["newCPid"].values
        carreaux_partiels.index = CP_fishnet["newCPid"].values
        for i, _ in CP_fishnet.iterrows():
            cp_aval = carreaux_partiels.loc[i, "idCPAval"]
            if cp_aval == 0:
                p1[i-1, :] = [CP_fishnet.loc[i, "x_c"],
                              CP_fishnet.loc[i, "y_c"]]
                p2[i-1, :] = [CP_fishnet.loc[i, "x_c"],
                              CP_fishnet.loc[i, "y_c"]]
            else:
                p1[i-1, :] = [CP_fishnet.loc[i, "x_c"],
                              CP_fishnet.loc[i, "y_c"]]
                p2[i-1, :] = [CP_fishnet.loc[cp_aval, "x_c"],
                              CP_fishnet.loc[cp_aval, "y_c"]]

        # Computing vector components
        u_c = p2[:, 0] - p1[:, 0]
        v_c = p2[:, 1] - p1[:, 1]

        ax.quiver(p1[idx_small_area, 0], p1[idx_small_area, 1],
                  u_c[idx_small_area], v_c[idx_small_area],
                  scale_units='xy',
                  angles='xy',
                  scale=1)
        plt.show()

    def _save_tif_ce_fishnet(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        # Get the fishnet vector file
        CE_shp = ogr.Open(self._CEfishnet, gdal.GA_ReadOnly)
        lyr = CE_shp.GetLayer()
        # Open shp a geopandas
        vector = gpd.read_file(self.ce_fishnet_name)
        geom_value = ((geom, value)
                      for geom, value in zip(vector.geometry, vector['newCEid']))
        xmin, xmax, ymin, ymax = lyr.GetExtent()
        # Get the resolution to the new raster
        pixel_size = self._dx
        x_res = floor((xmax-xmin)/pixel_size)
        y_res = ceil((ymax-ymin)/pixel_size)
        # Get the rasterio transformer and rasterize
        transform = rasterio.Affine(pixel_size, 0.0, xmin, 0.0,
                                    (-1)*pixel_size, ymax)
        rasterized = rasterio.features.rasterize(geom_value,
                                                 out_shape=(y_res, x_res),
                                                 fill=0,
                                                 out=None,
                                                 transform=transform,
                                                 all_touched=False,
                                                 default_value=1,
                                                 dtype=None)
        # Save file
        self._CEgrid = os.path.join(
            self._project_path, "geographic", "CEgrid.tif")
        with rasterio.open(self._CEgrid, "w",
                           driver="GTiff",
                           crs=vector.crs,
                           transform=transform,
                           dtype=rasterio.uint32,
                           count=1,
                           width=x_res,
                           nodata=0,
                           height=y_res) as dst:
            dst.write(rasterized, indexes=1)

    def create_CPgrid(self):
        """_summary_

        Returns:
            _type_: _description_
        """
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
