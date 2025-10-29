from __future__ import annotations
import sys
from math import floor
import itertools
import geopandas as gpd
import pandas as pd
from osgeo import gdal, ogr, osr
import numpy as np
from shapely.geometry import Polygon
from shapely.validation import make_valid, explain_validity
from . import projections
# from shapely.validation import make_valid, explain_validity
# import xarray as xr
# MultiPolygon
# import matplotlib.pyplot as plt
# from src.meteo.base import MeteoStation
# from src.core import projections
# from src.physiographic.base import Basin
# import rasterio
# from rasterio.features import shapes


# def polygonize_raster(raster_name: str):
#     # https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons
#     # Polygonize both rasters and change the path name to the shp
#     mask = None
#     with rasterio.Env():
#         with rasterio.open(raster_name) as src:
#             image = src.read(1)  # first band
#             results = (
#                 {'properties': {'raster_val': v}, 'geometry': s}
#                 for i, (s, v) in enumerate(shapes(image, mask=mask, transform=src.transform)))
#     geoms = list(results)
#     # Convert the list into a gpd
#     gdf = gpd.GeoDataFrame.from_features(geoms)
#     # Open the tif file to get the epsg
#     raster_ref = gdal.Open(raster_name, gdal.GA_ReadOnly)
#     epsg = projections.get_proj_code(raster_ref)
#     gdf = gdf.set_crs(epsg)
#     gdf['validity'] = gdf.apply(
#         lambda row: explain_validity(row.geometry), axis=1)
#     gdf.geometry = gdf.apply(
#         lambda row: make_valid(row.geometry), axis=1)
# print(gdf)
# gdf.geometry.make_valid()
# print(make_valid(gdf["geometry"]))
# gdf.geometry = make_valid(gdf.geometry)
# return gdf
# gdf.to_file(raster_name.replace("tif","shp"))
# print(gpd_d)
# self._Basin = self._Basin.replace(".tif", ".shp")
# self._SubBasins = self._SubBasins.replace(".tif", ".shp")
# gpd_d.to_file(self._Basin)
def polygonize_raster(raster_name: str):
    """_summary_

    Args:
        raster_name (str): _description_
    """
    # Open the file
    src_ds = gdal.Open(raster_name)
    # Define shp file attributes
    srcband = src_ds.GetRasterBand(1)
    dst_layername = 'polygonized'
    drv = ogr.GetDriverByName("ESRI Shapefile")
    # Set the export name
    dst_ds = drv.CreateDataSource(raster_name.replace(".tif", ".shp"))
    # Define the spatial reference based on the TIF attributes
    sp_ref = osr.SpatialReference()
    epsg = projections.get_proj_code(src_ds)
    sp_ref.SetFromUserInput(epsg)
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)
    # Rasterize the object
    fld = ogr.FieldDefn("raster_val", ogr.OFTInteger)
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("raster_val")
    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)


def fix_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """_summary_
    https://gis.stackexchange.com/questions/430384/using-shapely-methods-explain-validity-and-make-valid-on-shapefile
    Args:
        gdf (gpd.GeoDataFrame): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    gdf['validity'] = gdf.apply(
        lambda row: explain_validity(row.geometry), axis=1)
    gdf.geometry = gdf.apply(lambda row: make_valid(
        row.geometry) if not row.geometry.is_valid == 'Valid Geometry' else row.geometry, axis=1)
    return gdf
# def drop_duplicated_geometries(geoseries: gpd.GeoSeries):
#     """
#     taken from:
#     https://ml-gis-service.com/index.php/2021/09/24/toolbox-drop-duplicated-geometries-from-geodataframe/
#     Function drops duplicated geometries from a geoseries. It works as follow:
#         1. Take record from the dataset. Check it's index against list of indexes-to-skip. If it's not there then move to the next step.
#         2. Store record's index in the list of processed indexes (to re-create geoseries without duplicates) and in the list of indexes-to-skip.
#         3. Compare this record to all other records. If any of them is a duplicate then store its index in the indexes-to-skip.
#         4. If all records are checked then re-create dataframe without duplicates based on the list of processed indexes.
#     INPUT:

#     :param geoseries: (gpd.GeoSeries)

#     OUTPUT:

#     :returns: (gpd.GeoDataFrame)
#     """

#     indexes_to_skip = []
#     processed_indexes = []

#     for index, geom in geoseries.items():
#         if index not in indexes_to_skip:
#             processed_indexes.append(index)
#             indexes_to_skip.append(index)
#             for other_index, other_geom in geoseries.items():
#                 if other_index in indexes_to_skip:
#                     pass
#                 else:
#                     if geom.equals(other_geom):
#                         indexes_to_skip.append(other_index)
#                     else:
#                         pass
#     output_gs = geoseries[processed_indexes].copy()
#     return indexes_to_skip, processed_indexes


def rasterize_shp(grid_shp: str,
                  ref_name: str,
                  field: str) -> np.ndarray:
    """_summary_

    Args:
        grid_shp (str): _description_
        ref_name (str): _description_
        field (str): _description_

    Returns:
        np.ndarray: _description_
    """
    # Get raster georeference info
    raster = gdal.Open(ref_name, gdal.GA_ReadOnly)
    # transform = raster.GetGeoTransform()
    # xOrigin = transform[0]
    # yOrigin = transform[3]
    # pixelWidth = transform[1]
    # pixelHeight = -transform[5]
    x_res = raster.RasterXSize
    y_res = raster.RasterYSize
    shp = ogr.Open(grid_shp, gdal.GA_ReadOnly)
    lyr = shp.GetLayer()
    # Get shp info
    # proj = lyr.GetSpatialRef()
    # xmin, xmax, ymin, ymax = lyr.GetExtent()
    # Get the resolution to the new raster
    # Create tiff format file
    # grid_raster = gdal.GetDriverByName('GTiff').Create(
    #     grid_shp.replace(".shp", ".tif"),
    #     x_res, y_res, 1, gdal.GDT_Int32)
    # Create the raster in the memory
    grid_raster = gdal.GetDriverByName('MEM').Create(
        '', x_res, y_res, 1, gdal.GDT_Int32)
    grid_raster.SetGeoTransform(raster.GetGeoTransform())
    grid_raster.SetProjection(raster.GetProjection())
    band = grid_raster.GetRasterBand(1)
    band.SetNoDataValue(0)
    gdal.RasterizeLayer(grid_raster, [1], lyr, options=[
                        "ATTRIBUTE="+field])
    band = grid_raster.GetRasterBand(1)
    data_array = band.ReadAsArray()
    return data_array


def rasterize_shp_as_byte(grid_shp: str,
                          ref_name: str, field: str,
                          name: str) -> None:
    """_summary_

    Args:
        grid_shp (str): _description_
        ref_name (str): _description_
        field (str): _description_
        name (str): _description_

    Returns:
        _type_: _description_
    """
    # Get raster georeference info
    raster = gdal.Open(ref_name, gdal.GA_ReadOnly)
    # transform = raster.GetGeoTransform()
    # xOrigin = transform[0]
    # yOrigin = transform[3]
    # pixelWidth = transform[1]
    # pixelHeight = -transform[5]
    x_res = raster.RasterXSize
    y_res = raster.RasterYSize
    # Get shp info
    shp = ogr.Open(grid_shp, gdal.GA_ReadOnly)
    lyr = shp.GetLayer()
    # proj = lyr.GetSpatialRef()
    # xmin, xmax, ymin, ymax = lyr.GetExtent()
    # Get the resolution to the new raster
    # Create tiff format file
    # grid_raster = gdal.GetDriverByName('GTiff').Create(
    #     grid_shp.replace(".shp",".tif"),
    #     x_res, y_res, 1, gdal.GDT_Int32)
    # Create the raster in the memory
    grid_raster = gdal.GetDriverByName('GTiff').Create(
        name, x_res, y_res, 1, gdal.GDT_Byte)
    grid_raster.SetGeoTransform(raster.GetGeoTransform())
    grid_raster.SetProjection(raster.GetProjection())
    band = grid_raster.GetRasterBand(1)
    band.SetNoDataValue(0)
    gdal.RasterizeLayer(grid_raster, [1], lyr, options=[
                        "ATTRIBUTE="+field])
    band = grid_raster.GetRasterBand(1)
    data_array = band.ReadAsArray()
    return data_array


def get_altitude_point(DEM: gdal.Dataset,
                       lat_utm: np.array,
                       lon_utm: np.array):
    """_summary_

    Args:
        DEM (gdal.Dataset): _description_
        lat_utm (np.array): _description_
        lon_utm (np.array): _description_

    Returns:
        _type_: _description_
    """
    xtup, ytup, ptup = GetExtent(DEM)
    # (xmin, xmax), (ymin, ymax), (xpixel, ypixel)
    # xmin, xmax, ymin, ymax, xpixel, ypixel = GetExtent(DEM)
    DEM_array = DEM.ReadAsArray()
    altitudes = []
    for y, x in zip(lat_utm, lon_utm):
        col = int((x-np.amin(xtup))/ptup[0])
        row = int((np.amax(ytup)-y)/(-ptup[1]))
        if col < 0 or row < 0:
            altitudes.append(-9999)
        elif row >= DEM_array.shape[0] or col > DEM_array.shape[1]:
            altitudes.append(-9999)
        else:
            altitudes.append(DEM_array[row, col])
    return altitudes


def GetExtent(raster: gdal.Dataset):
    """ Return list of corner coordinates from a gdal Dataset """
    xmin, xpixel, _, ymax, _, ypixel = raster.GetGeoTransform()
    width, height = raster.RasterXSize, raster.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return (xmin, xmax), (ymin, ymax), (xpixel, ypixel)
    # return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)


def regrid_CE(raster: gdal.Dataset,
              shp: ogr.DataSource,
              grid_size: int) -> gdal.Dataset:
    """_summary_

    Args:
        raster (gdal.Dataset): _description_
        shp (ogr.DataSource): _description_
        grid_size (int): _description_

    Returns:
        gdal.Dataset: _description_
    """
    xmin, xpixel, _, ymax, _, ypixel = raster.GetGeoTransform()
    width, height = raster.RasterXSize, raster.RasterYSize
    # transform = raster.GetGeoTransform()
    # xOrigin = transform[0]
    # yOrigin = transform[3]
    # pixelWidth = transform[1]
    # pixelHeight = transform[5]
    # xmin = transform[0]
    # ymin = transform[3]
    xmax = xmin + width*xpixel
    ymin = ymax + height*ypixel
    lyr = shp.GetLayer()
    proj = lyr.GetSpatialRef()
    xmin, xmax, ymin, ymax = lyr.GetExtent()
    # Get the resolution to the new raster
    xres = int((xmax-xmin)/grid_size)
    yres = int((ymax-ymin)/grid_size)
    # yres = raster.RasterYSize
    # print(xres,yres)
    # GTiff
    mem_ds = gdal.GetDriverByName("GTiff").Create(
        'test.tif', abs(xres), abs(yres), 1, gdal.GDT_Int32)
    mem_ds.SetProjection(proj.ExportToWkt())
    mem_ds.SetGeoTransform((xmin, grid_size, 0, ymin, 0, grid_size))
    band = mem_ds.GetRasterBand(1)
    NoData_value = 0
    band.SetNoDataValue(NoData_value)
    gdal.RasterizeLayer(mem_ds, [1], lyr, options=["ATTRIBUTE=newid"])
    return mem_ds


def rasterize_feature(gdf: gpd.GeoDataFrame,
                      raster_name: str,
                      att: str) -> np.ndarray:
    """_summary_

    Args:
        gdf (gpd.GeoDataFrame): _description_
        raster_name (str): _description_
        att (str): _description_

    Returns:
        np.ndarray: _description_
    """
    # Get raster georeference info
    raster = gdal.Open(raster_name, gdal.GA_ReadOnly)
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    # yOrigin1 = yOrigin + pixelHeight[5] * raster.RasterYSize
    # Get no data value
    # no_data = raster.GetRasterBand(1).GetNoDataValue()
    # xmin, ymin, xmax, ymax = gdf.bounds

    # Specify offset and rows and columns to read
    xoff = floor((gdf['minx'] - xOrigin)/pixelWidth)
    yoff = floor((yOrigin - gdf['maxy'])/pixelHeight)
    xcount = int((gdf['maxx'] - gdf['minx'])/pixelWidth)
    ycount = int((gdf['maxy'] - gdf['miny'])/pixelHeight)

    # Get the projection
    proj = osr.SpatialReference(wkt=raster.GetProjection())
    # create the spatial reference system, WGS84
    srs = osr.SpatialReference()
    epsg = int(proj.GetAttrValue('AUTHORITY', 1))
    srs.ImportFromEPSG(epsg)

    # Create an in-memory vector dataset from the GeoDataFrame's geometry column
    outDriver = ogr.GetDriverByName("MEMORY")
    tempSHP = outDriver.CreateDataSource("")
    if gdf.geometry.geom_type == "Polygon":
        temp_layer = tempSHP.CreateLayer('feat', srs, ogr.wkbPolygon)
    elif gdf.geometry.geom_type == "MultiPolygon":
        temp_layer = tempSHP.CreateLayer('feat', srs, ogr.wkbMultiPolygon)
    else:
        sys.exit("Geometry type not valid. Fix it")

    # Add an ID field
    idField = ogr.FieldDefn(att, ogr.OFTInteger)
    geom = ogr.CreateGeometryFromWkb(gdf['geometry'].wkb)
    temp_layer.CreateField(idField)

    # Create the feature and set values
    featureDefn = temp_layer.GetLayerDefn()
    feat = ogr.Feature(featureDefn)
    feat.SetGeometry(geom)
    feat.SetField(att, 1)
    temp_layer.CreateFeature(feat)

    target_ds = gdal.GetDriverByName('MEM').Create(
        "", xcount, ycount, gdal.GDT_Byte)
    target_ds.SetGeoTransform((gdf['minx'], pixelWidth, 0,
                               gdf['miny'], 0, pixelHeight,))

    # Rasterize the in-memory vector dataset onto the mask array
    gdal.RasterizeLayer(target_ds, [1],
                        temp_layer, options=["ATTRIBUTE="+att])

    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray()

    banddataraster = raster.GetRasterBand(1)
    # no_data = raster.GetRasterBand(1).GetNoDataValue()
    # Mask array based on the nondata value
    dataraster = banddataraster.ReadAsArray(
        xoff, yoff, xcount, ycount)
    # masked_dataraster = np.ma.masked_where(dataraster==no_data,dataraster)
    # Change
    # np.set_printoptions(threshold=sys.maxsize)
    # print(datamask)
    # print(dataraster)
    # print(datamask.shape)
    return datamask*dataraster


# def rasterize_shp2(gdf: gpd.GeoDataFrame,
#                    raster_name: str) -> np.ndarray:
#     # Get affine transformation matrix from the raster dataset
#     # Get raster georeference info
#     raster = gdal.Open(raster_name, gdal.GA_ReadOnly)
#     transform = raster.GetGeoTransform()
#     xOrigin = transform[0]
#     yOrigin = transform[3]
#     pixelWidth = transform[1]
#     pixelHeight = -transform[5]
#     # Get no data value
#     # no_data = raster.GetRasterBand(1).GetNoDataValue()
#     xmin, ymin, xmax, ymax = gdf.geometry.bounds

#     # Specify offset and rows and columns to read
#     xoff = ceil((xmin - xOrigin)/pixelWidth)
#     yoff = int((yOrigin - ymax)/pixelHeight)
#     xcount = int((xmax - xmin)/pixelWidth)
#     ycount = int((ymax - ymin)/pixelHeight)

#     proj = osr.SpatialReference(wkt=raster.GetProjection())
#     # create the spatial reference system, WGS84
#     srs = osr.SpatialReference()
#     epsg = int(proj.GetAttrValue('AUTHORITY', 1))
#     srs.ImportFromEPSG(epsg)

#     # Create an in-memory vector dataset from the GeoDataFrame's geometry column
#     driver = ogr.GetDriverByName('Memory')

#     temp_ds = driver.CreateDataSource('')
#     temp_layer = temp_ds.CreateLayer('temp', srs, ogr.wkbMultiPolygon)
#     feature = ogr.Feature(temp_layer.GetLayerDefn())
#     idField = ogr.FieldDefn("temp", ogr.OFTInteger)
#     temp_layer.CreateField(idField)
#     feature.SetField("temp", 1)
#     # print(feature.GetField("temp"))
#     # Make a geometry, from Shapely object
#     geom = ogr.CreateGeometryFromWkb(gdf['geometry'].wkb)
#     feature.SetGeometry(geom)
#     temp_layer.CreateFeature(feature)

#     # Create an in-memory temporal raster
#     target_ds = gdal.GetDriverByName('MEM').Create(
#         '', xcount, ycount, ogr.OFTInteger)
#     target_ds.SetGeoTransform((
#         xmin, pixelWidth, 0,
#         ymax, 0, pixelHeight,
#     ))

#     bandmask = target_ds.GetRasterBand(1)
#     bandmask.SetNoDataValue(0)
#     bandmask.Fill(0)
#     # Create for target raster the same projection as for the value raster
#     raster_srs = osr.SpatialReference()
#     raster_srs.ImportFromWkt(raster.GetProjectionRef())
#     target_ds.SetProjection(raster_srs.ExportToWkt())

#     # Rasterize the in-memory vector dataset onto the mask array
#     gdal.RasterizeLayer(target_ds, [1], temp_layer, options=["ATTRIBUTE=temp"])

#     bandmask = target_ds.GetRasterBand(1)
#     datamask = bandmask.ReadAsArray()

#     banddataraster = raster.GetRasterBand(1)

#     dataraster = banddataraster.ReadAsArray(
#         xoff, yoff, xcount, ycount)

#     pass


# def zonal_statistics(feat,
#                      raster: gdal.Dataset,
#                      shp: ogr.DataSource,
#                      field: str):
#     # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#calculate-zonal-statistics
#     # Get raster georeference info
#     transform = raster.GetGeoTransform()
#     xOrigin = transform[0]
#     yOrigin = transform[3]
#     pixelWidth = transform[1]
#     pixelHeight = transform[5]
#     # Get no data value
#     no_data = raster.GetRasterBand(1).GetNoDataValue()
#     burn_value = int(feat.GetField(field))
#     lyr = shp.GetLayer()
#     # Get extent of feat
#     geom = feat.GetGeometryRef()
#     if (geom.GetGeometryName() == 'MULTIPOLYGON'):
#         count = 0
#         pointsX = []
#         pointsY = []
#         for polygon in geom:
#             geomInner = geom.GetGeometryRef(count)
#             ring = geomInner.GetGeometryRef(0)
#             numpoints = ring.GetPointCount()
#             for p in range(numpoints):
#                 lon, lat, z = ring.GetPoint(p)
#                 pointsX.append(lon)
#                 pointsY.append(lat)
#             count += 1
#     elif (geom.GetGeometryName() == 'POLYGON'):
#         ring = geom.GetGeometryRef(0)
#         numpoints = ring.GetPointCount()
#         pointsX = []
#         pointsY = []
#         for p in range(numpoints):
#             lon, lat, z = ring.GetPoint(p)
#             pointsX.append(lon)
#             pointsY.append(lat)
#     else:
#         sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")

#     xmin = min(pointsX)
#     xmax = max(pointsX)
#     ymin = min(pointsY)
#     ymax = max(pointsY)

#     # Specify offset and rows and columns to read
#     xoff = ceil((xmin - xOrigin)/pixelWidth)
#     yoff = int((yOrigin - ymax)/pixelWidth)
#     xcount = int((xmax - xmin)/pixelWidth)
#     ycount = int((ymax - ymin)/pixelWidth)

#     # print(lyr.GetExtent())
#     # print(xmin,xmax,ymin,ymax)
#     # print(xoff,yoff,ycount,xcount)
#     # Create memory target raster
#     # GTiff
#     target_ds = gdal.GetDriverByName('MEM').Create(
#         '', xcount, ycount, gdal.GDT_Int16)
#     target_ds.SetGeoTransform((
#         xmin, pixelWidth, 0,
#         ymax, 0, pixelHeight,
#     ))
#     bandmask = target_ds.GetRasterBand(1)
#     bandmask.SetNoDataValue(0)
#     bandmask.Fill(0)
#     # Create for target raster the same projection as for the value raster
#     raster_srs = osr.SpatialReference()
#     raster_srs.ImportFromWkt(raster.GetProjectionRef())
#     target_ds.SetProjection(raster_srs.ExportToWkt())
#     # bandmask = target_ds.GetRasterBand(1)
#     # Rasterize zone polygon to raster
#     # query = field + " = '" + str(burn_value) +"'"
#     query = "INTERSECTS(geometry, POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s)))" % \
#         (xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin)
#     # Rasterize the shapefile with a value of 1 for all pixels inside the selected polygons
#     gdal.RasterizeLayer(target_ds, [1], lyr, options=[
#                         'WHERE="%s"' % query], burn_values=[1])
#     # print(shp)
#     # gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])
#     # Read raster as arrays
#     banddataraster = raster.GetRasterBand(1)

#     dataraster = banddataraster.ReadAsArray(
#         xoff, yoff, xcount, ycount)
#     # bandmask = target_ds.GetRasterBand(1)
#     datamask = bandmask.ReadAsArray()
#     print(datamask)
#     zoneraster = np.ma.masked_array(
#         dataraster,  np.logical_not(datamask))
#     # print(dataraster)
#     # row, col = centroid_index(feat, target_ds)
#     # if zoneraster[abs(row), abs(col)] == 0:
#     #     return 0
#     # else:
#     #     # Calculate statistics of zonal raster
#     #     return zoneraster[abs(row), abs(col)]
#     # np.mean(zoneraster), np.median(zoneraster), np.std(zoneraster), np.var(zoneraster)
#     return zoneraster


# def centroid_index(feat, raster: gdal.Dataset):
#     transform = raster.GetGeoTransform()
#     xOrigin = transform[0]
#     yOrigin = transform[3]
#     pixelWidth = transform[1]
#     pixelHeight = transform[5]
#     # Get the centroid
#     geom_poly = feat.GetGeometryRef()
#     centroid = geom_poly.Centroid()
#     point = centroid.GetPoint(0)
#     col = int((point[0] - xOrigin) / pixelWidth)
#     row = int((yOrigin - point[1]) / pixelHeight)
#     return row, col


def get_index_list(raster: gdal.Dataset,
                   x: np.ndarray,
                   y: np.ndarray) -> tuple:
    """_summary_

    Args:
        raster (gdal.Dataset): _description_
        x (np.ndarray): _description_
        y (np.ndarray): _description_

    Returns:
        tuple: _description_
    """
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    col = np.abs(((x - xOrigin) / pixelWidth)).astype(int)
    row = np.abs(((yOrigin - y) / pixelHeight)).astype(int)
    # xmin, xpixel, _, ymax, _, ypixel = raster.GetGeoTransform()
    # width, height = raster.RasterXSize, raster.RasterYSize
    # xmax = xmin + width*xpixel
    # ymin = ymax + height*ypixel
    # col = np.abs(((x - xmin) / width)).astype(int)
    # row = np.abs(((ymax - y) / height)).astype(int)
    return row, col


# def loop_zonal_stats(raster: gdal.Dataset,
#                      shp: ogr.DataSource):
#     lyr = shp.GetLayer()
#     featList = range(lyr.GetFeatureCount())
#     statDict = {}
#     # idField = ogr.FieldDefn("newid", ogr.OFTInteger)
#     # lyr.CreateField(idField)
#     for i, FID in enumerate(featList):
#         feat = lyr.GetFeature(FID)

#         # Get extent of feat
#         CEval = zonal_statistics(feat, raster, shp)
#         feat.SetField("newid", int(CEval))
#         lyr.SetFeature(feat)
#     return statDict


def falls_in_extent(extent: tuple,
                    x: list,
                    y: list) -> np.ndarray:
    """_summary_

    Args:
        extent (tuple): _description_
        x (list): _description_
        y (list): _description_

    Returns:
        np.ndarray: _description_
    """
    # Extract the extent from the given values
    x_ext, y_ext = extent
    # (lon, lat) -> (x, y)
    xy_pairs = np.c_[list(itertools.product(x, y))]
    idx, = np.where((xy_pairs[:, 0] <= np.amin(x_ext)) |
                    (xy_pairs[:, 0] >= np.amax(x_ext)) |
                    (xy_pairs[:, 1] <= np.amin(y_ext)) |
                    (xy_pairs[:, 1] >= np.amax(y_ext)))
    xy_pairs[idx] = np.nan
    return xy_pairs[~np.isnan(xy_pairs).any(axis=1), :]


# def remap_CEgrid(CEgrid: gdal.Dataset,
#                  fishnet: ogr.DataSource):
#     array = CEgrid.ReadAsArray()
#     dx, dy = np.where(array == 699)
#     pass


def find_nearest(array: np.ndarray, value: float) -> np.ndarray:
    """_summary_
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    Args:
        array (np.ndarray): _description_
        value (float): _description_

    Returns:
        np.ndarray: _description_
    """
    # array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()
    return (np.abs(array - value)).argmin()


def convert_multi_to_poly(geometry: ogr.Geometry) -> ogr.Geometry:
    """_summary_

    Args:
        geometry (ogr.Geometry): _description_

    Returns:
        ogr.Geometry: _description_
    """
    # Convert the Multipolygon to a Polygon
    # if geometry.GetGeometryCount() == 1:
    #     polygon = geometry.GetGeometryRef(0)
    # else:
    polygon = ogr.Geometry(ogr.wkbPolygon)
    for i in range(geometry.GetGeometryCount()):
        poly = geometry.GetGeometryRef(i)
        for ring in range(poly.GetGeometryCount()):
            polygon.AddGeometry(poly.GetGeometryRef(ring))

    return polygon


def ogr_to_gpd(shp: ogr.DataSource) -> gpd.GeoDataFrame:
    """_summary_
    This function converst a input OGR object into a GeoDataFrame
    The results is taken to be used in the merge function.
    Or any other funtion that needs geopandas
    Args:
        shp (ogr.DataSource): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    layer = shp.GetLayer()
    srs = layer.GetSpatialRef()
    col_names = [field.name for field in layer.schema]
    col_names.append("geometry")
    # Convert the layer to a pandas DataFrame
    df = pd.DataFrame(columns=col_names,
                      index=range(len(layer)))
    for i, feature in enumerate(layer):
        geom2 = feature.GetGeometryRef()
        geom2 = geom2.MakeValid()
        if geom2.GetGeometryName() == 'MULTIPOLYGON':

            polygon = convert_multi_to_poly(geom2)
            print(polygon)
            #
            # num_rings = polygon.GetGeometryCount()
            exterior_ring = polygon.GetGeometryRef(0)
        else:
            # num_rings = geom2.GetGeometryCount()
            exterior_ring = geom2.GetGeometryRef(0)

        num_points = exterior_ring.GetPointCount()
        points = []
        # Getting vextexes of each polygon
        for j in range(num_points):
            point = exterior_ring.GetPoint(j)
            points.append(point)

        for field in layer.schema:
            df[field.name][i] = feature.GetField(field.name)
        # print(geom2.GetGeometryName())
        # if geom2.GetGeometryName()=="MULTIPOLYGON":
        #     df["geometry"][i] = MultiPolygon(points)
        # else:
        df["geometry"][i] = Polygon(points)

    EPSG_code = srs.GetAttrValue('AUTHORITY', 1)
    EPSG = f"EPSG:{EPSG_code}"
    gdf = gpd.GeoDataFrame(df, crs=EPSG)
    print(gdf)
    return gdf


def saveGTIFF(ref_TIF: str,
              data_array: np.ndarray,
              output_name: str):
    ds = gdal.Open(ref_TIF, gdal.GA_ReadOnly)
    [rows, cols] = data_array.shape
    # Save the file depending on the file type name
    type_name = data_array.dtype.base.name
    if type_name in ["uint8"]:
        gdal_type = gdal.GDT_Byte
    if type_name in ["int16"]:
        gdal_type = gdal.GDT_UInt16
    if type_name in ["int32", "int64"]:
        gdal_type = gdal.GDT_UInt32
    elif type_name in ["float32", "float64"]:
        gdal_type = gdal.GDT_Float32
    # type_name = data_type.base.name
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(output_name, cols, rows, 1, gdal_type)
    # sets same geotransform as input
    outdata.SetGeoTransform(ds.GetGeoTransform())
    outdata.SetProjection(ds.GetProjection())  # sets same projection as input
    outdata.GetRasterBand(1).WriteArray(data_array)
    # if you want these values transparent
    outdata.GetRasterBand(1).SetNoDataValue(-9999)
    outdata.FlushCache()  # saves to disk!!


def convert_slope(slope_file_path: str) -> np.ndarray:
    ds = gdal.Open(slope_file_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    nan_val = band.GetNoDataValue()
    array_d = band.ReadAsArray()
    # convert the slope from degrees to m/m -> atan(slope)
    array_d[array_d == nan_val] = np.nan
    array_d = np.deg2rad(array_d)
    array_d = np.arctan(array_d)
    return array_d


def reclassify_landcover(LC_file_path: str, classes_idx: list, output_tif_name: str):
    ds = gdal.Open(LC_file_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    nan_val = band.GetNoDataValue()
    if nan_val == None:
        nan_val = 0
    array_d = band.ReadAsArray()
    array_d[array_d == 0] = nan_val
    reclass_array = np.copy(array_d)
    # This is class name for cequeau. Just for reference
    # classes_name = ["forest","bare","wetlands","water","urban"]
    # Loop into each of the classes from the land cover
    for land_type, ceq_class in zip(classes_idx, [1, 7, 14, 18, 16]):
        # idx = np.stack([land_type])
        for idx in land_type:
            reclass_array[array_d == idx] = ceq_class
        # np.delete(a, idx[:, 1:])
    saveGTIFF(LC_file_path, reclass_array, output_tif_name)
    return 0

# Function to retrieve the outlet point of the basin
def get_outlet_point(FAC_path: str):
    ds = gdal.Open(FAC_path, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    array_d = band.ReadAsArray()
    i, j = np.where(array_d == np.amax(array_d))
    xoff, a, b, yoff, d, e = ds.GetGeoTransform()
    xp = a * j[0] + b * i[0] + xoff
    yp = d * j[0] + e * i[0] + yoff
    return xp, yp
