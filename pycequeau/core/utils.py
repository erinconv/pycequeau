from __future__ import annotations

import numpy as np
import xarray as xr
from osgeo import gdal, ogr, osr
from pycequeau.meteo.base import MeteoStation
from pycequeau.core import projections
import sys
from pycequeau.physiographic.base import Basin
import itertools


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


def zonal_statistics(feat,
                     raster: gdal.Dataset,
                     shp: ogr.DataSource):
    # https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#calculate-zonal-statistics
    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
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
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1

    # Create memory target raster
    # GTiff
    target_ds = gdal.GetDriverByName('MEM').Create(
        '', xcount, ycount, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
    ))

    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())

    # Rasterize zone polygon to raster
    gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])

    # Read raster as arrays
    banddataraster = raster.GetRasterBand(1)
    dataraster = banddataraster.ReadAsArray(
        xoff, yoff, xcount, ycount).astype(np.float)

    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(np.float)

    # Mask zone of raster
    zoneraster = np.ma.masked_array(
        dataraster,  np.logical_not(datamask))
    row, col = centroid_index(feat, target_ds)
    if zoneraster[abs(row), abs(col)] == 0:
        return 0
    else:
        # Calculate statistics of zonal raster
        return zoneraster[abs(row), abs(col)]
    # return np.average(zoneraster), np.mean(zoneraster), np.median(zoneraster), np.std(zoneraster), np.var(zoneraster)


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

