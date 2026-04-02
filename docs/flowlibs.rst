Getting the Geographical Data Using FlowLibs
============================================

As part of the development of Pycequeau, an additional repository provides helper libraries for computationally
intensive hydrological terrain processing. These utilities are distributed in the following GitHub repository:

https://github.com/erinconv/FlowLibs


What is FlowLibs?
-----------------

FlowLibs is a collection of C++ tools for hydrological terrain analysis, distributed primarily through Docker and
wrapped by small Python scripts included in the repository. The current workflow supports:

- Complete breaching of a DEM
- D8 flow-direction generation in ESRI/ArcGIS format
- D8 flow-accumulation computation
- Watershed delineation from pour points
- Subcatchment delineation
- Optional masking of outputs by watershed extent
- Export of the river network to a GeoPackage

FlowLibs can be used as a standalone preprocessing tool, but its main purpose is to support the Pycequeau workflow.


How to use FlowLibs
-------------------

FlowLibs is easiest to use through Docker. Install Docker first:

https://www.docker.com/

Then either build the image locally from the repository root:

.. code-block:: bash

    docker build -t erinconv/flowlibs:latest -f docker/Dockerfile .

or pull the latest published image:

.. code-block:: bash

    docker pull erinconv/flowlibs:latest

The repository provides Python wrappers in the ``python/`` directory. In the current documented workflow, the steps
are run separately:

1. carve the DEM
2. delineate watersheds and subcatchments from pour points
3. optionally mask outputs by watershed extent
4. export the river network to a GeoPackage

By default, the carving wrapper writes outputs in a folder named ``dephier_outputs`` located next to the input DEM.
The later wrappers use that folder as their working directory.

The Python wrapper scripts are available in the FlowLibs GitHub repository:
`python folder <https://github.com/erinconv/FlowLibs/tree/main/python>`_


Running the documented workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Carve the DEM:

.. code-block:: bash

    python python/run_carving.py path/to/your/DEM.tif

Script:
`run_carving.py <https://github.com/erinconv/FlowLibs/blob/main/python/run_carving.py>`_

This creates the ``dephier_outputs`` folder next to the DEM.

The carving step writes:

- ``run1_DEM_breached.tif``: breached DEM
- ``run1_DIR.tif``: D8 flow directions in ESRI/ArcGIS format
- ``run1_FAC.tif``: D8 flow accumulation

2. Delineate watersheds from one or more pour points:

.. code-block:: bash

    python python/run_watersheds.py path/to/dephier_outputs --pour_points X Y id basin_name

Script:
`run_watersheds.py <https://github.com/erinconv/FlowLibs/blob/main/python/run_watersheds.py>`_

Where ``X`` and ``Y`` correspond to the point coordinate, ``id`` is the label that will be assigned to the delineated watershed, and ``basin_name`` is 
the name of the basin. This command will delineate watersheds based on the provided pour point coordinates.

Useful options include:

- ``--dem_name`` to override the default breached DEM name
- ``--flow_dirs_name`` to override the default flow-direction raster name
- ``--accumulation_name`` to override the default accumulation raster name
- ``--output_prefix`` to change the default prefix ``run1``
- ``--stream_threshold`` to define the stream and subcatchment threshold in cells

The watershed step writes:

- ``run1_RIV.tif``: binary stream raster
- ``run1_RIV_ORD.tif``: Horton-Strahler stream order
- ``run1_WAT.tif``: watershed labels
- ``run1_CAT.tif``: subcatchment labels

3. Optionally mask the generated rasters using the watershed extent:

.. code-block:: bash

    python python/run_mask_outputs_by_watershed.py path/to/dephier_outputs

Script:
`run_mask_outputs_by_watershed.py <https://github.com/erinconv/FlowLibs/blob/main/python/run_mask_outputs_by_watershed.py>`_

Useful options include:

- ``--watershed-name`` to select the watershed raster used as the mask
- ``--suffix`` to control the suffix added to masked rasters

With the default suffix ``_masked``, the masking step writes files such as:

- ``run1_DEM_breached_masked.tif``
- ``run1_DIR_masked.tif``
- ``run1_FAC_masked.tif``
- ``run1_RIV_masked.tif``
- ``run1_RIV_ORD_masked.tif``
- ``run1_CAT_masked.tif``

The watershed raster itself remains ``run1_WAT.tif`` and is reused as the mask reference.

4. Export the stream network to a GeoPackage:

.. code-block:: bash

    python python/run_stream_vectorize.py path/to/dephier_outputs

Script:
`run_stream_vectorize.py <https://github.com/erinconv/FlowLibs/blob/main/python/run_stream_vectorize.py>`_

Useful options include:

- ``--dem-name``
- ``--flow-dirs-name``
- ``--accumulation-name``
- ``--streams-name``
- ``--stream-order-name``
- ``--watersheds-name``
- ``--subcatchments-name``
- ``--output-name``

The vectorization step writes a GeoPackage, by default:

- ``run1_river_network.gpkg``

If you vectorize masked rasters and choose a masked output name, a typical result is:

- ``run1_river_network_masked.gpkg``

The GeoPackage contains two layers:

- ``river_segments``
- ``river_nodes``

These layers store hydrological attributes such as stream order, contributing accumulation, watershed and
subcatchment IDs, segment length, elevation drop, and slope.
