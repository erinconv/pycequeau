Getting Meteorological Data from ERA5
=====================================

For the current pycequeau meteorological workflow, the recommended starting
point is **ERA5** from the Copernicus Climate Data Store (CDS). ERA5 is a
global atmospheric reanalysis distributed by the Copernicus Climate Change
Service and accessed through the public CDS web portal and API.

If you are starting a new project, ERA5 is usually the simplest source to use
because it is global, well documented, and already distributed in NetCDF-ready
formats that fit naturally with pycequeau workflow.


Why ERA5?
---------

ERA5 is a practical choice for pycequeau because it provides:

- global spatial coverage
- long historical records
- hourly products that can be aggregated to daily values
- easy access through the CDS web interface and API

The pycequeau meteorological workflow is built around **gridded NetCDF
datasets**. In practice, that means ERA5 can be downloaded, preprocessed into
daily fields, and then ingested into
:class:`pycequeau.meteo.meteo_netcdf.NetCDFMeteo`.


Where to get the data
---------------------

ERA5 data can be downloaded from the official Copernicus Climate Data Store:

- Climate Data Store documentation:
  `<https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+documentation>`_
- ERA5 hourly data on single levels:
  `<https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download>`_
- ERA5 hourly data on pressure levels:
  `<https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download>`_

For most pycequeau applications, the **single-level ERA5 product** is the one
to start with. Pressure-level products are usually only needed for more
specialized atmospheric analyses and are not the usual entry point for the
standard pycequeau workflow.


How to access ERA5 with the CDS API
-----------------------------------

The recommended way to automate downloads is the official `cdsapi` client.
This documentation intentionally keeps that part brief, because account setup,
API credentials, licence acceptance, and request syntax are maintained by the
Copernicus team and may evolve over time.

Please follow the official CDS guidance for:

- creating an account
- accepting dataset licences
- installing ``cdsapi``
- configuring your API key
- writing download requests

Official references:

- CDS documentation:
  `<https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+documentation>`_
- CDS migration and quick-guide notice:
  `<https://confluence.ecmwf.int/display/CKB/Please+read%3A+CDS+and+ADS+migrating+to+new+infrastructure%3A+Common+Data+Store+%28CDS%29+Engine>`_

In other words, pycequeau does **not** replace the Copernicus download step.
You should first retrieve the raw ERA5 files from CDS, then prepare those
files for the pycequeau meteorological workflow.


What pycequeau expects
----------------------

The pycequeau meteorological workflow works with **daily gridded NetCDF
variables**. The ingestion layer normalizes supported variables through the
meteorological schema, and the NetCDF workflow then interpolates those fields
to the CE grid.

Some of the main variable names currently recognized by pycequeau include:

- ``tasmax`` or ``tMax`` for daily maximum temperature
- ``tasmin`` or ``tMin`` for daily minimum temperature
- ``tp``, ``pr`` or ``pTot`` for daily precipitation
- ``ssrd``, ``ssr`` or ``rsds`` for shortwave radiation
- ``strd`` or ``msdwlwrf`` for longwave radiation
- ``d2m`` for dewpoint temperature
- ``tcc`` or ``clt`` for cloud cover
- ``sp`` for surface pressure
- ``wind`` or ``sfcWind`` for scalar wind speed

At the moment, the workflow is centered on **daily** inputs, not raw hourly
ERA5 time steps. That means the CDS download is only the first part of the
workflow. The ERA5 fields still need to be aggregated or converted into the
daily variables expected by pycequeau.


Units
-----

The meteorological inputs must include a valid ``units`` attribute in the
NetCDF variables. pycequeau can normalize several common ERA5 and meteorology
unit conventions during preprocessing, but the metadata still needs to be
present and correct.

For the full list of supported unit strings and the canonical units expected
by the workflow, see :doc:`units`.


Recommended ERA5 variables to download first
--------------------------------------------

A practical starting set is:

- 2 m temperature
- total precipitation
- surface solar radiation downwards
- surface thermal radiation downwards
- 2 m dewpoint temperature
- total cloud cover
- surface pressure
- 10 m u-component of wind
- 10 m v-component of wind

This combination gives you enough information to prepare most of the variables
used in the current gridded workflow.


A simple preparation strategy
-----------------------------

In practice, a typical workflow looks like this:

1. Download ERA5 data from CDS for the basin extent and study period.
2. Store the downloaded fields as NetCDF files.
3. Aggregate hourly variables into daily products where needed.
4. Convert derived variables before ingestion when appropriate.
   For example:

   - build daily wind speed from ``u10`` and ``v10`` with
     :class:`pycequeau.meteo.calculators.wind_speed.WindSpeedCalculator`
   - derive vapor pressure from dewpoint temperature, if that variable is
     needed, with
     :class:`pycequeau.meteo.calculators.vapor_pressure.VaporPressureCalculator`
5. Load the prepared files into :class:`pycequeau.meteo.meteo_netcdf.NetCDFMeteo`.
6. Interpolate the meteorological fields to the CE grid and export them to the
   CEQUEAU-ready format.

The preprocessing utilities in :mod:`pycequeau.meteo.calculators` can help
with some of these derivations once the raw NetCDF files have already been
downloaded. A separate page will document that workflow in more detail.


What this page does not cover
-----------------------------

This page does not try to reproduce the official CDS API manual or provide a
full ERA5 request cookbook. Those details are better maintained in the
official Copernicus documentation.

The goal here is simply to document the recommended **data source** and show
how it connects to the pycequeau meteorological workflow.
