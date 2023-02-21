from osgeo import ogr
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pycequeau.core import projections
# Extract first layer of features from shapefile using OGR
ds = ogr.Open('/mnt/c/Users/erinc/Documents/01-PhD/00-Geographic/Melezes/Basin2.shp')
nlay = ds.GetLayerCount()
lyr = ds.GetLayer(0)

# Get extent and calculate buffer size
ext = lyr.GetExtent()
xmin,xmax,ymin,ymax = lyr.GetExtent()
print(xmax,xmin,ymax,ymin)
ymin, xmin = projections.latlon_to_utm(xmin, ymin,"EPSG:32198")
ymax, xmax = projections.latlon_to_utm(xmax, ymax,"EPSG:32198")
# all_y, all_x = projections.latlon_to_utm(all_x, all_y,"EPSG:32198")
print(xmax,xmin,ymax,ymin)
# xoff = (xmin-xmax)
# yoff = (ymin-ymax)

# Prepare figure
fig = plt.figure()
ax = fig.add_subplot(111)
# ax.set_xlim(xmin-xoff,xmax+xoff)
# ax.set_ylim(ymin-yoff,ymax+yoff)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

paths = []
lyr.ResetReading()

# Read all features in layer and store as paths
for feat in lyr:
    geom = feat.geometry()
    codes = []
    all_x = []
    all_y = []
    for i in range(geom.GetGeometryCount()):
        # Read ring geometry and create path
        r = geom.GetGeometryRef(i)
        x = [r.GetX(j) for j in range(r.GetPointCount())]
        y = [r.GetY(j) for j in range(r.GetPointCount())]
        y,x = projections.latlon_to_utm(x, y,"EPSG:32198")
        # skip boundary between individual rings
        # 
        codes += [mpath.Path.MOVETO] + \
                     (len(x)-1)*[mpath.Path.LINETO]
        all_x += x
        all_y += y
    # all_y, all_x = projections.latlon_to_utm(all_x, all_y,"EPSG:32198")
    path = mpath.Path(np.column_stack((all_x,all_y)), codes)
    paths.append(path)

# x,y = projections.latlon_to_utm(x,y,"EPSG:32198")
# print(paths)
# Add paths as patches to axes
for path in paths:
    patch = mpatches.PathPatch(path, \
            facecolor='none', edgecolor='black')
    ax.add_patch(patch)

ax.set_aspect(1.0)
plt.show()