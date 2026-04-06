Conditioning the DEM
====================

One crucial step to improve the accuracy of the flow directions derived from the DEM is to condition the flow lines
to follow a pre-defined path set by a reference hydrological network or by the existence of water body areas.
This process produces the so-called Conditioned DEM, which is a modified version of the original elevation model that
no longer represent reality but is very powerful to accurately derive rivers form the DEM.

Pycequeau provides built in functions to use an existent waterbody mask to perform such conditioning or even though, a reference
stream network that can be burnt in the DEM to pre set the flow paths, thus improving the results obtained by the hydrological simulations.

The hydrological conditioning perform by Pycequeau is very subtle and only makes sure that no flatten water surfaces are present
in the condtioned DEM. This is very important because the derivation of the Flo Direction Model can be very impredictible in flat
areas and and normallly the algorithms deliver very poor results.

Condition the DEM using an existent WBM
---------------------------------------

Although this demonstration will use the already derived raster from Copernicus, this can be applied to any DEM and WBM that are well aligned
with each other. To perform it with the maps we downlaoded in :doc:`the previous step <cop_dem>`, you must use the following code:

.. code-block:: python

    from pycequeau.core import CopernicusDEMProcessor

    dem_path = "path/to/your/dem.tif"
    wbm_path = "path/to/your/wbm.tif"
    cop_dem = CopernicusDEMProcessor.basic_dem_condition(dem_path, wbm_path)

As a result of this, you wull obtain a raster called ``DEM_condidioned.tif`` which is the one you can use as input fro the
:doc:`GRASS tools <grass>` that were previously presented in the documentation.

Handling Depressions (Sinks) in DEM Conditioning
------------------------------------------------

An important step in DEM preprocessing is the treatment of depressions (commonly referred to as *sinks*).

A sink is defined as a grid cell (or group of cells) whose surrounding neighbouring cells all have higher elevation values.
In raster-based flow routing, such cells prevent the proper computation of drainage directions because no downslope path exists.

Sinks may originate from:

- Artefacts introduced during DEM generation (e.g., interpolation errors, radar shadow, vegetation bias, or digitization artefacts),
- Resolution limitations,
- Natural landscape features such as closed basins or lakes.

Because flow direction algorithms (e.g., D8, D∞) require a continuous downslope gradient, untreated sinks can disrupt drainage network extraction
and flow accumulation modelling.

There are different approaches to correct the presence of sinks in a DEM. One of the most common methods is **depression filling**, in which the
elevation of sink cells is raised to match the minimum surrounding elevation required to create a valid downslope flow path.
One of the main drawbacks of this approach is that it distorts the original terrain structure by creating artificial flat areas around the filled sinks.
These flat surfaces may introduce ambiguities when computing the Flow Direction model and can alter local slope gradients.

Another approach to address the existence of sinks is **depression carving (breaching)**.
Instead of raising the elevation within the sink, the carving process lowers selected cells along a path connecting the depression to a valid downslope outlet.
This approach better preserves realistic terrain gradients and avoids the creation of artificial flat areas.
Carving is particularly appropriate when sinks in the DEM result from noise or production errors rather than representing natural landscape features.
Nonetheless, carving may not be the most appropriate method in landscapes where depressions correspond to real hydrological features (e.g., endorheic basins).


Integration with the Pycequeau workflow
---------------------------------------

The Copernicus DEM products used in this workflow are distributed in a geographic coordinate reference system. After applying the basic conditioning step, 
the resulting raster must be reprojected to a projected CRS before it can be used in Pycequeau.
This requirement is important because several downstream preprocessing operations in Pycequeau assume planar coordinates and linear units. 
In particular, grid dimensions, flow-path lengths, polygon areas, and proximity-based calculations are expected to be expressed in projected units rather than angular coordinates. 
Using a DEM in geographic coordinates would therefore lead to inconsistent spatial measurements and unreliable hydrological outputs.

For this reason, the conditioned DEM should be reprojected with GIS software such as QGIS or ArcGIS before continuing with the Pycequeau workflow. 
The selected projected CRS should be appropriate for the study area and should be consistent with the other rasters.
