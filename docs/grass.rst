Getting the geographical data with grass
========================================

The pycequeau library requires several ``.shp`` and ``.tif`` files to extract the physiographical. The toolbox is designed to accept data formats from GIS-based software such as `QGIS`_, `ArcGIS`_, `SAGA-GIS`_, or `GRASS-GIS`_.

.. _QGIS: https://www.qgis.org/en/site/forusers/download.html
.. _ArcGIS: https://pro.arcgis.com/en/pro-app/latest/get-started/download-arcgis-pro.htm
.. _GRASS-GIS: https://grass.osgeo.org/download/
.. _SAGA-GIS: https://sourceforge.net/projects/saga-gis/files/


Download the DEM
----------------

We will use QGIS and GRASS-GIS as the base software to obtain all the files in this case. We will start by downloading the DEM using the SRTM-downloader plugin available through QGIS. Then, in the plugins manager, click on ``Manage and Install plugins`` and write ``SRTM Downloader``. Once you have successfully installed the plugin, we can download the DEM for our basin.

Downloading the Mélèzes River data
----------------------------------
To download the DEM we will need to place our QGIS canvas in the extent that covers the whole basin area. To find out where exactly your basin is located, you will require to have the outlet point coordinates and a basemap to help you to find the boundaries (in our case, we will use OpenMaps as basemap) as shown in the image:

.. image:: figures/GRASS-tutorial/01-Extent.png
  :width: 400
  :alt: Getting the extent of the basin
  :align: center

Now, open the SRTM-downloader plugin, and when you click on the ``Set canvas extent`` option it will automatically fill up the extent area to download the DEM files. When you click on download, you will require to register and input your credentilas from the `EARTHDATA`_ portal

.. _EARTHDATA: https://urs.earthdata.nasa.gov//users/new

.. image:: figures/GRASS-tutorial/02-Download.png
  :width: 300
  :alt: Download the DEM
  :align: center

In this step, you will have several raster files that you need to merge. To do so, we will create a ``virtual raster`` to mosaic all the images by clicking on the ``Raster -> Miscellaneous -> Build virtual Raster`` . Select all the images and a directory to store the mosaic and run.

.. image:: figures/GRASS-tutorial/03-Mosaic.png
  :width: 300
  :alt: Mosaic image
  :align: center

Now we need to discern the extent of the actual basin. To do so, we must trace the approximate watershed extent using a shapefile. To create a new shapefile, go to ``New shape file layer`` and fill up the options as follows:

.. image:: figures/GRASS-tutorial/04-Extent-shp.png
  :width: 300
  :alt: Mosaic image
  :align: center

Now, let's crop this mosaic to remove any undesired areas. To do this, right click on the mosaic layer in the layer panel in QGIS, select ``Export`` and then ``Save as...``. In this step, two important things need to be defined.

#. **Define the reference system**: This is fundamental because the pycequeau toolbox works only in projected coordinates. For this case, the ``EPSG:32198`` is required. The selection of the Projected system is case-dependent, and you need to find out what is your EPSG code based on the location of your basin.
#. **The extent to clip the basin**: You can calculate the extent to clip the mosaic by using the previously created polygon. Click in ``Calculate from layer`` and select the shpe layer to have the same extent.

.. image:: figures/GRASS-tutorial/05-Export.png
  :width: 300
  :alt: Mosaic image
  :align: center

Now we are all set to delineate the watershed using this DEM.


Delineating the watershed
-------------------------

Once obtained the DEM file, let's open GRASS-GIS and create a new location based on the geographical metadata of our DEM file. In this case we will, I set up the project as follows:

.. image:: figures/GRASS-tutorial/06-OpenGRASS.png
  :width: 400
  :alt: Grass start
  :align: center

Once create the region, in the right panel clic on new to create a new mapset and define your own basin name. Then double click to start working

.. image:: figures/GRASS-tutorial/07-NewMapset.png
  :width: 400
  :alt: new mapset
  :align: center

We need to start by uploading our DEM file. To do that, go to the console option and run this command line, replacing ``path/to/your/project/folder/DEM.tif`` with your own folder

.. code-block:: bash
  
  r.in.gdal input=path/to/your/project/folder/DEM.tif output=DEM

Once imported, set the region by runnig the next command line:

.. code-block:: bash  
  
  g.region -d raster=DEM

.. note::

    OPTIONAL: If we have an existent stream nerwork, we can burn this vector file to make our process even more exact by running the next command.

    .. code-block:: bash 

      v.in.ogr input=path/to/your/project/folder/streams.shp output=streams
    
    .. code-block:: bash 

      r.carve -n --overwrite raster=DEM vector=streams output=BurnedDEM

    Doing this step is highly recommended in basins where the surface is mostly flat and the flow direction must be corrected. If you do not have any reference stream network, we recommend download it form the following global dataset: https://www.hydrosheds.org/



Now, you can see the DEM raster from the display window:

.. image:: figures/GRASS-tutorial/08-DEM-GRASS.png
  :width: 400
  :alt: new mapset
  :align: center

Then, let's compute the the DIR file

.. code-block:: bash
  
  r.watershed --overwrite elevation=DEM drainage=DIR


To process the corrected files, we need to install one grass extention. You can do so by running the next command into the console:

.. code-block:: bash
  
  g.extension extension=r.accumulate 

Once finished, let's obtain the flow accumulation, the subbasins and the stream network by running this line

.. code-block:: bash
  
  r.accumulate --verbose --overwrite direction=DIR format=auto accumulation=FAC subwatershed=Watershed stream=Streams threshold=THERESHOLD coordinates=X,Y

where ``THRESHOLD`` is the minimum flow accumulation to be considered as river stream, ``X,Y`` are the outlet coordinates.

Correcting the watershed delineation
------------------------------------

It is possible to obtain a wrong delineation of the watershed because of give ``X,Y`` do not fall into an actual ``Flow Accumulation`` pixel. To correct this, use the display window of GRASS to find the coordinates where this X,Y pair that fall in the FAC map:

.. image:: figures/GRASS-tutorial/09-outlet_point.png
  :width: 400
  :alt: new mapset
  :align: center

Now, let's run the following command line to obtain the watershed delineation correctly:

.. code-block:: bash
  
  r.water.outlet --overwrite input=DIR output=Watershed coordinates=Xcorrected,Ycorrected

Retrieve the subbasin raster file.
----------------------------------

To retrieve the subbasin raster, run the following command line

.. code-block:: bash
  
  r.watershed --overwrite elevation=DEM threshold=THEREDHOLD drainage=DIR basin=CAT

.. warning::

  If the previous instruction does not provide the correct subbasin strcuture, the following command lines can solve the issue:
  
  .. code-block:: bash

    g.extension extension=r.stream.basins

  .. code-block:: bash

    r.stream.extract --overwrite elevation=DEM threshold=THEREDHOLD stream_raster=streams_r stream_vector=streams_v direction=DIR

  .. code-block:: bash

    r.stream.basins --overwrite direction=DIR stream_rast=streams_r basins=CAT



and this is the result

.. image:: figures/GRASS-tutorial/09-Subbasins.png
  :width: 400
  :alt: new mapset
  :align: center

Now, let's mask all the results using the obtained watershed delineation as follows:

.. code-block:: bash

  r.mask raster=Watershed

We can translate the CAT file from `tif` format into `shp` as follows:

.. code-block:: bash

  r.to.vect -s input=CAT output=CAT type=area


Now, export the raster as standard TIF formats.


.. code-block:: bash

  r.out.gdal input=FilledDEM output=path/to/your/project/folder/DEM_Filled.tif
  r.out.gdal input=DIR output=path/to/your/project/folder/DIR.tif
  r.out.gdal input=FAC output=path/to/your/project/folder/FAC.tif
  r.out.gdal input=CAT output=path/to/your/project/folder/CAT.tif
  r.out.gdal input=Watershed output=path/to/your/project/folder/Watershed.tif

You can now open those files in your favourite GIS-based software and continue the following steps.