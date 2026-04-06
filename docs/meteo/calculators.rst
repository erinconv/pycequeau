Using the Meteorological Calculators
====================================

The meteorological calculators are small preprocessing utilities that help
derive variables from already-downloaded NetCDF meteorological inputs. They
are useful when the source dataset does not directly provide the exact daily
variable expected by the pycequeau workflow.

At the moment, the calculator interface includes tools for:

- computing wind speed from ``u10`` and ``v10``
- computing vapor pressure from dewpoint temperature


When to use them
----------------

Use the calculators after downloading and organizing your meteorological
NetCDF files, but before loading them into
:class:`pycequeau.meteo.meteo_netcdf.NetCDFMeteo`.

The source variables still need valid unit metadata. For the supported input
unit strings and the canonical units expected by pycequeau, see :doc:`units`.

In practice, the workflow is usually:

1. Download ERA5 or another meteorological dataset.
2. Aggregate the source fields to daily values if needed.
3. Use a calculator when a required pycequeau variable must be derived from
   another source field.
4. Save the derived variable as its own NetCDF file.
5. Load the prepared set of NetCDF files into the pycequeau meteorological
   workflow.


Main entry points
-----------------

The main calculator entry points are:

- :meth:`pycequeau.meteo.calculators.base.MeteoCalculator.create_variable_dataset`
  to return an :class:`xarray.Dataset`
- :meth:`pycequeau.meteo.calculators.base.MeteoCalculator.create_variable_file`
  to write the derived variable directly to a NetCDF file

These entry points can work with:

- a folder containing NetCDF files
- a single NetCDF file
- a list of NetCDF files

You can inspect the currently available derivations with:

.. code-block:: python

   from pycequeau.meteo.calculators import MeteoCalculator

   MeteoCalculator.available_derivations()


Dataset example
---------------

.. code-block:: python

   from pycequeau.meteo.calculators import MeteoCalculator

   ds = MeteoCalculator.create_variable_dataset(
       r"path/to/netcdf_folder",
       variable="wind_speed",
   )

This call asks pycequeau to:

- resolve the calculator associated with ``wind_speed``
- look for the required input variables in the provided folder
- build the derived dataset


Wind speed example
------------------

To compute wind speed from daily mean wind components:

.. code-block:: python

   from pycequeau.meteo.calculators import MeteoCalculator

   wind = MeteoCalculator.create_variable_dataset(
       r"path/to/netcdf_folder",
       variable="wind_speed",
   )

By default, this uses
:class:`pycequeau.meteo.calculators.wind_speed.WindSpeedCalculator` and looks
for either:

- ``u10`` or ``10m_u_component_of_wind``
- ``v10`` or ``10m_v_component_of_wind``

To write the result directly to a NetCDF file:

.. code-block:: python

   from pycequeau.meteo.calculators import MeteoCalculator

   output_path = MeteoCalculator.create_variable_file(
       r"path/to/netcdf_folder",
       variable="wind_speed",
       output_path=r"path/to/output/wind.nc",
   )


Vapor pressure example
----------------------

To compute vapor pressure from dewpoint temperature:

.. code-block:: python

   from pycequeau.meteo.calculators import MeteoCalculator

   vp = MeteoCalculator.create_variable_dataset(
       r"path/to/netcdf_folder",
       variable="vapor_pressure",
       vapor_pressure_method="murray_1967",
   )

This uses
:class:`pycequeau.meteo.calculators.vapor_pressure.VaporPressureCalculator`.
By default, it looks for ``d2m`` or ``dewpoint_temperature``.

Two formulations are currently available through the
``vapor_pressure_method`` argument:

- ``lowe_1977``
- ``murray_1967``


Overriding source variable names
--------------------------------

If your NetCDF files use different variable names, you can provide your own
names explicitly.

For a single-source derivation:

.. code-block:: python

   ds = MeteoCalculator.create_variable_dataset(
       r"path/to/netcdf_folder",
       variable="vapor_pressure",
       source_variable="tdps",
   )

For a multi-source derivation:

.. code-block:: python

   ds = MeteoCalculator.create_variable_dataset(
       r"path/to/netcdf_folder",
       variable="wind_speed",
       source_variable=("u_wind", "v_wind"),
   )
