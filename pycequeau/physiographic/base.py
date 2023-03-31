from __future__ import annotations

from osgeo import gdal, ogr
import os
from math import ceil
import pandas as pd
from pycequeau.physiographic import CPfishnet as CPS
from pycequeau.core import utils as u
from pycequeau.core import projections as proj
import geopandas as gpd
import sys


class Basin:
    def __init__(self,
                 project_folder: str,
                 basin_name: str,
                 file_list: list) -> None:
        # Create project structure
        self._project_path = project_folder
        self.name = basin_name
        self._project_structure(project_folder, file_list)
        # Create here the fishnet
        self._CEfishnet = os.path.join(
            self._project_path, "geographic", "CE_fishnet.shp")
        self._CPfishnet = os.path.join(
            self._project_path, "geographic", "CP_fishnet.shp")

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
        self._DIR = os.path.join(
            self._project_path, "geographic", file_list[2])
        self._Basin = os.path.join(
            self._project_path, "geographic", file_list[3])
        self._SubBasins = os.path.join(
            self._project_path, "geographic", file_list[4])

    def set_dimenssions(self, dx: float, dy: float):
        self._dx = dx
        self._dy = dy

    def get_dimenssions(self):
        return [self._dx, self._dy]

    def get_EPSG(self):
        return self._epsg

    def _set_EPSG(self):
        self._epsg = proj.get_proj_code(self._DEM)

    @classmethod
    def _create_CEfishnet(cls,
                          basin: str,
                          dx: float,
                          dy: float,
                          fishnet: str) -> None:
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
        ringXleftOrigin = xmin
        ringXrightOrigin = xmin + dx
        ringYtopOrigin = ymax
        ringYbottomOrigin = ymax-dy

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
            pass
        # Execute the merge
        else:
            pass
        self.join_shps(self._CEfishnet,
                    self._SubBasins,
                    self._CPfishnet)

    def create_CEfishnet(self):
        
        # Check if the file already exist
        if os.path.exists(self._CEfishnet):
            pass
        else:
            pass
        self._create_CEfishnet(self._Basin,
                            self._dx,
                            self._dy,
                            self._CEfishnet)

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
        CPfishnet = CPS.identify_small_CPs(CEfishnet, CPfishnet, area_th)
        CPfishnet, CEfishnet = CPS.remove_border_CPs(CEfishnet, CPfishnet, self._FAC)
        CPfishnet = CPS.remove_smallCP(CEfishnet,CPfishnet)
        CPfishnet = CPS.dissolve_pixels(CEfishnet,CPfishnet,area_th)
        CPfishnet = CPS.force_4CP(CEfishnet,CPfishnet,area_th)
        CPfishnet.to_file(self._CPfishnet)
        # if os.path.exists(out_name):
        #     os.remove(out_name)
        # CPfishnet.to_file(out_name)

    def CP_routing(self, flow_th=1000):
        # Create the CP grid from the merged shp file
        # This has the same dimensions as the FAC file 
        # in order to compare them both to do the routing process
        # The CP grid is already processed and well polished
        CP_array = u.rasterize_shp(self._CPfishnet,self._FAC,"CPid")
        # Open fishnets as geodataframes
        CEfishnet = gpd.read_file(self._CEfishnet)
        CPfishnet = gpd.read_file(self._CPfishnet)
        # Get the routing table. This table renames the CPs from the outlet
        # up to the last CP. Here the upstream CP are also identify for each
        # individual CP
        rtable = CPS.routing_table(CEfishnet,
                                     CPfishnet,
                                     self._FAC,
                                     self._DIR,
                                     CP_array,
                                     flow_th)
        # Obtain the downstream CP based on the previous process
        rtable = CPS.get_downstream_CP(rtable)
        rtable["newCPid"] = pd.to_numeric(rtable["newCPid"])
        rtable["downstreamCPs"] = pd.to_numeric(rtable["downstreamCPs"])
        # Compute the route CP by CP
        outlet_routes = CPS.outlet_routes(rtable)
        # Renumbering the fishnets
        CPfishnet,CEfishnet = CPS.renumber_fishnets(CPfishnet,CEfishnet,rtable)
        # Compute cumulative percentage of surface area
        CPfishnet, upstreamCPs = CPS.cumulative_areas(CPfishnet,CEfishnet,outlet_routes)
        
        
        # Add rtable to the 
        # out_name = os.path.join(
        #     self._project_path, "geographic", "CP_fishnet_test.shp")
        # out_name2 = os.path.join(
        #     self._project_path, "geographic", "CE_fishnet_test.shp")
        # if os.path.exists(out_name):
        #     os.remove(out_name)
        #     os.remove(out_name2)
        
        # Save the shp files in the hard drive
        CPfishnet.to_file(self._CPfishnet)
        CEfishnet.to_file(self._CEfishnet)

    def find_grid_coordinates(self):
        # Create the CE grid with the shp dimenssions,
        # not with the reference raster dimensions to save memory and
        # make the processes faster
        CE_Array = self.create_CEgrid()
        
        pass

    def create_CEgrid(self):
        # Default no data value of the CAT raster
        ref_raster = gdal.Open(self._DEM,gdal.GA_ReadOnly)
        nodata = ref_raster.GetRasterBand(1).GetNoDataValue()
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
        x_res = int((xmax-xmin)/self._dx)
        y_res = int((ymax-ymin)/self._dy)
        # Get dimenssions for the CEgrid rasters

        # CEgrid path
        self._CEgrid = gdal.GetDriverByName('MEM').Create(
            '', abs(x_res), abs(y_res), 1, gdal.GDT_Int16)
        self._CEgrid.SetProjection(proj.ExportToWkt())
        self._CEgrid.SetGeoTransform((xmin, self._dx, 0, ymin, 0, self._dy))
        band = self._CEgrid.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        gdal.RasterizeLayer(self._CEgrid, [1], lyr, options=["ATTRIBUTE=newCEid"])
        data_array = band.ReadAsArray()
        return data_array

    def create_CPgrid(self):
        # Default no data value of the CAT raster
        nodata = self._CAT.GetRasterBand(1).GetNoDataValue()
        lyr = self._CPfishnet.GetLayer()
        proj = lyr.GetSpatialRef()
        xmin, xmax, ymin, ymax = lyr.GetExtent()
        # Get the resolution to the new raster
        scale = 25
        x_res = int((xmax-xmin)/self._dx)*scale
        y_res = int((ymax-ymin)/self._dy)*scale
        # Make the Union between the two shps
        path = os.path.join(self._project_path, "geographic", "CPgrid.tif")
        self._CPgrid = gdal.GetDriverByName('GTiff').Create(
            path, abs(x_res), abs(y_res), 1, gdal.GDT_Int32)
        # print(self._CPgrid)
        # print(path)
        print(abs(x_res), abs(y_res))
        self._CPgrid.SetProjection(proj.ExportToWkt())
        self._CPgrid.SetGeoTransform(
            (xmin, self._dx/scale, 0, ymin, 0, self._dy/scale))
        band = self._CPgrid.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        gdal.RasterizeLayer(self._CPgrid, [1], lyr, options=["ATTRIBUTE=CPid"])
