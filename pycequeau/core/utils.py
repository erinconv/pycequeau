from __future__ import annotations
import sys
import numpy as np
import xarray as xr
from osgeo import gdal, ogr, osr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid
from math import ceil
# from pycequeau.meteo.base import MeteoStation
# from pycequeau.core import projections
# from pycequeau.physiographic.base import Basin
import itertools
import matplotlib.pyplot as plt

def rasterize_shp(grid_shp: str,
                   ref_name: str,
                   field: str)-> np.ndarray:
    # Get raster georeference info
    raster = gdal.Open(ref_name, gdal.GA_ReadOnly)
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    x_res = raster.RasterXSize
    y_res = raster.RasterYSize
    # Get shp info
    shp = ogr.Open(grid_shp,gdal.GA_ReadOnly)
    lyr = shp.GetLayer()
    proj = lyr.GetSpatialRef()
    xmin, xmax, ymin, ymax = lyr.GetExtent()
    # Get the resolution to the new raster
    # Create tiff format file
    # grid_raster = gdal.GetDriverByName('GTiff').Create(
    #     grid_shp.replace(".shp",".tif"), 
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


def get_i_j_CEQUEAU_grid(CEgrid: gdal.Dataset):
    pass


def get_altitude_point(DEM: gdal.Dataset,
                       lat_utm: np.array,
                       lon_utm: np.array):
    xmin, xmax, ymin, ymax, xpixel, ypixel = GetExtent(DEM)
    DEM_array = DEM.ReadAsArray()
    altitudes = []
    for x, y in zip(lat_utm, lon_utm):
        col = int((x-xmin)/xpixel)
        row = int((ymax-y)/-ypixel)
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

    return (xmin, xmax, ymin, ymax, xpixel, ypixel)
    # return (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)


def regrid_CE(raster: gdal.Dataset,
              shp: ogr.DataSource,
              grid_size: int) -> gdal.Dataset:
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
                      raster_name: str) -> np.ndarray:
    # Get raster georeference info
    raster = gdal.Open(raster_name,gdal.GA_ReadOnly)
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    # yOrigin1 = yOrigin + pixelHeight[5] * raster.RasterYSize
    # Get no data value
    no_data = raster.GetRasterBand(1).GetNoDataValue()
    xmin, ymin, xmax, ymax = gdf.geometry.bounds

    # Specify offset and rows and columns to read
    xoff = ceil((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelHeight)
    xcount = int((xmax - xmin)/pixelWidth)
    ycount = int((ymax - ymin)/pixelHeight)

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
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    geom = ogr.CreateGeometryFromWkb(gdf['geometry'].wkb)
    temp_layer.CreateField(idField)

    # Create the feature and set values
    featureDefn = temp_layer.GetLayerDefn()
    feat = ogr.Feature(featureDefn)
    feat.SetGeometry(geom)
    feat.SetField("id", 1)
    temp_layer.CreateFeature(feat)

    target_ds = gdal.GetDriverByName('MEM').Create(
        "", xcount, ycount, gdal.GDT_Byte)
    target_ds.SetGeoTransform((xmin, pixelWidth, 0,
                               ymin, 0, pixelHeight,))

    # Rasterize the in-memory vector dataset onto the mask array
    gdal.RasterizeLayer(target_ds, [1], 
                        temp_layer, options=["ATTRIBUTE=id"])

    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray()
    
    banddataraster = raster.GetRasterBand(1)

    dataraster = banddataraster.ReadAsArray(
        xoff, yoff, xcount, ycount)
    print(datamask)
    print(dataraster)
    print(datamask.shape)
    pass


def rasterize_shp2(gdf: gpd.GeoDataFrame,
                  raster_name: str) -> np.ndarray:
    # Get affine transformation matrix from the raster dataset
    # Get raster georeference info
    raster = gdal.Open(raster_name,gdal.GA_ReadOnly)
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    # Get no data value
    no_data = raster.GetRasterBand(1).GetNoDataValue()
    xmin, ymin, xmax, ymax = gdf.geometry.bounds

    # Specify offset and rows and columns to read
    xoff = ceil((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelHeight)
    xcount = int((xmax - xmin)/pixelWidth)
    ycount = int((ymax - ymin)/pixelHeight)

    proj = osr.SpatialReference(wkt=raster.GetProjection())
    # create the spatial reference system, WGS84
    srs = osr.SpatialReference()
    epsg = int(proj.GetAttrValue('AUTHORITY', 1))
    srs.ImportFromEPSG(epsg)

    # Create an in-memory vector dataset from the GeoDataFrame's geometry column
    driver = ogr.GetDriverByName('Memory')

    temp_ds = driver.CreateDataSource('')
    temp_layer = temp_ds.CreateLayer('temp', srs, ogr.wkbMultiPolygon)
    feature = ogr.Feature(temp_layer.GetLayerDefn())
    idField = ogr.FieldDefn("temp", ogr.OFTInteger)
    temp_layer.CreateField(idField)
    feature.SetField("temp", 1)
    # print(feature.GetField("temp"))
    # Make a geometry, from Shapely object
    geom = ogr.CreateGeometryFromWkb(gdf['geometry'].wkb)
    feature.SetGeometry(geom)
    temp_layer.CreateFeature(feature)

    # Create an in-memory temporal raster
    target_ds = gdal.GetDriverByName('MEM').Create(
        '', xcount, ycount, ogr.OFTInteger)
    target_ds.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
    ))

    bandmask = target_ds.GetRasterBand(1)
    bandmask.SetNoDataValue(0)
    bandmask.Fill(0)
    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())

    # Rasterize the in-memory vector dataset onto the mask array
    gdal.RasterizeLayer(target_ds, [1], temp_layer, options=["ATTRIBUTE=temp"])

    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray()

    banddataraster = raster.GetRasterBand(1)

    dataraster = banddataraster.ReadAsArray(
        xoff, yoff, xcount, ycount)

    pass


def zonal_statistics(feat,
                     raster: gdal.Dataset,
                     shp: ogr.DataSource,
                     field: str):
    # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#calculate-zonal-statistics
    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    # Get no data value
    no_data = raster.GetRasterBand(1).GetNoDataValue()
    burn_value = int(feat.GetField(field))
    # TODO: Check for projection matching
    lyr = shp.GetLayer()
    # Get extent of feat
    geom = feat.GetGeometryRef()
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []
        pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []
        pointsY = []
        for p in range(numpoints):
            lon, lat, z = ring.GetPoint(p)
            pointsX.append(lon)
            pointsY.append(lat)
    else:
        sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")

    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    # Specify offset and rows and columns to read
    xoff = ceil((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)
    ycount = int((ymax - ymin)/pixelWidth)

    # print(lyr.GetExtent())
    # print(xmin,xmax,ymin,ymax)
    # print(xoff,yoff,ycount,xcount)
    # Create memory target raster
    # GTiff
    target_ds = gdal.GetDriverByName('MEM').Create(
        '', xcount, ycount, gdal.GDT_Int16)
    target_ds.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
    ))
    bandmask = target_ds.GetRasterBand(1)
    bandmask.SetNoDataValue(0)
    bandmask.Fill(0)
    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())
    # bandmask = target_ds.GetRasterBand(1)
    # Rasterize zone polygon to raster
    # query = field + " = '" + str(burn_value) +"'"
    query = "INTERSECTS(geometry, POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s)))" % \
        (xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin)
    # Rasterize the shapefile with a value of 1 for all pixels inside the selected polygons
    gdal.RasterizeLayer(target_ds, [1], lyr, options=[
                        'WHERE="%s"' % query], burn_values=[1])
    # print(shp)
    # gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])
    # Read raster as arrays
    banddataraster = raster.GetRasterBand(1)

    dataraster = banddataraster.ReadAsArray(
        xoff, yoff, xcount, ycount)
    # bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray()
    print(datamask)
    zoneraster = np.ma.masked_array(
        dataraster,  np.logical_not(datamask))
    # print(dataraster)
    # row, col = centroid_index(feat, target_ds)
    # if zoneraster[abs(row), abs(col)] == 0:
    #     return 0
    # else:
    #     # Calculate statistics of zonal raster
    #     return zoneraster[abs(row), abs(col)]
    # np.mean(zoneraster), np.median(zoneraster), np.std(zoneraster), np.var(zoneraster)
    return zoneraster


def centroid_index(feat, raster: gdal.Dataset):
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    # Get the centroid
    geom_poly = feat.GetGeometryRef()
    centroid = geom_poly.Centroid()
    point = centroid.GetPoint(0)
    col = int((point[0] - xOrigin) / pixelWidth)
    row = int((yOrigin - point[1]) / pixelHeight)
    return row, col


def get_index_list(raster: gdal.Dataset,
                   x: np.array,
                   y: np.array):
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


def loop_zonal_stats(raster: gdal.Dataset,
                     shp: ogr.DataSource):
    lyr = shp.GetLayer()
    featList = range(lyr.GetFeatureCount())
    statDict = {}
    # TODO: Check if feature already exist
    # idField = ogr.FieldDefn("newid", ogr.OFTInteger)
    # lyr.CreateField(idField)
    for i, FID in enumerate(featList):
        feat = lyr.GetFeature(FID)

        # Get extent of feat
        CEval = zonal_statistics(feat, raster, shp)
        feat.SetField("newid", int(CEval))
        lyr.SetFeature(feat)
    return statDict


def falls_in_extent(fishnet: ogr.DataSource,
                    x: list,
                    y: list):
    # Fishnet must alway be in lat lon.
    # TODO: Generalize method to deal with coord conversion
    shp_layer = fishnet.GetLayer()
    xmin, xmax, ymin, ymax = shp_layer.GetExtent()
    xy_pairs = np.c_[list(itertools.product(x, y))]
    idx, = np.where((xy_pairs[:, 0] <= xmin) |
                    (xy_pairs[:, 0] >= xmax) |
                    (xy_pairs[:, 1] <= ymin) |
                    (xy_pairs[:, 1] >= ymax))
    xy_pairs[idx] = np.nan
    return xy_pairs[~np.isnan(xy_pairs).any(axis=1), :]


def remap_CEgrid(CEgrid: gdal.Dataset,
                 fishnet: ogr.DataSource):
    array = CEgrid.ReadAsArray()
    dx, dy = np.where(array == 699)
    pass


def find_nearest(array: np.ndarray, value: float):
    # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    # array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()
    return (np.abs(array - value)).argmin()


def convert_multi_to_poly(geometry: ogr.Geometry) -> ogr.Geometry:
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
