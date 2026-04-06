Welcome to pycequeau
====================

**pycequeau** is a Python library that prepares the inputs needed to build and run the **CEQUEAU** hydrological and water temperature model. It helps you obtain and process physiographic data, meteorological inputs, and parameters for your watershed.

.. note::

   This project is under active development. Feedback, suggestions, and contributions are welcome.

What you can do
---------------

- **Physiographic data**: Extract or prepare basin geometry and terrain-derived inputs.
- **Meteorological data**: Get and format weather inputs (e.g. from NetCDF).
- **Parameters**: Generate or manage CEQUEAU model parameters structure.

The library is designed to work with common geospatial formats and tools (e.g. QGIS, GRASS GIS, or other GIS software).

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   Installation <readme>

.. toctree::
   :maxdepth: 2
   :caption: Physiographic

   physiographic/cop_dem
   physiographic/condition_dem
   physiographic/grass
   physiographic/flowlibs
   physiographic/get_physio.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Meteorological

   meteo/get_meteo.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Simulation

   simulation/get_params.ipynb

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
