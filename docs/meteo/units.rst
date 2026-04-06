Meteorological Units
====================

This page documents the units expected by the pycequeau meteorological
workflow and the source-unit strings that can currently be normalized during
preprocessing.

.. important::

   All meteorological variables loaded from NetCDF must include a ``units``
   attribute. If that metadata is missing, pycequeau cannot normalize the
   input.


Canonical units
---------------

The main internal target units are:

- temperature and dewpoint temperature: ``C``
- precipitation: ``mm d-1``
- shortwave and longwave radiation: ``MJ m-2 d-1``
- cloud cover: ``0-1``
- wind speed: ``km h-1``
- relative humidity: ``%``
- vapor pressure: ``mmHg``
- surface pressure: ``Pa``


Reference Table
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 25 12 28

   * - Variable
     - Supported names
     - Canonical unit
     - Accepted source units
   * - Temperature
     - ``tasmax``, ``tMax``;

       ``tasmin``, ``tMin``;

       ``d2m``
     - ``C``
     - ``C``, ``K``
   * - Precipitation
     - ``tp``, ``pr``, ``pTot``
     - ``mm d-1``
     - ``mm d-1``, ``mm``,

       ``m d-1``, ``m``,

       ``kg m-2 s-1``
   * - Shortwave radiation
     - ``ssrd``, ``ssr``, ``rsds``, ``rayonnement``
     - ``MJ m-2 d-1``
     - ``MJ m-2 d-1``,

       ``MJ m-2``,

       ``J m-2 d-1``,

       ``J m-2``,

       ``W m-2``
   * - Longwave radiation
     - ``strd``, ``msdwlwrf``, ``longwaveRad``
     - ``MJ m-2 d-1``
     - ``MJ m-2 d-1``,

       ``MJ m-2``,

       ``J m-2 d-1``,

       ``J m-2``,

       ``W m-2``
   * - Cloud cover
     - ``tcc``, ``clt``, ``nebulosite``
     - ``0-1``
     - ``0-1``, ``fraction``,

       ``1``, ``%``,

       ``percent``
   * - Wind speed
     - ``wind``, ``sfcWind``, ``vitesseVent``
     - ``km h-1``
     - ``km h-1``, ``m s-1``
   * - Relative humidity
     - ``hurs``, ``humiditeRelative``
     - ``%``
     - ``%``, ``percent``,

       ``0-1``, ``fraction``,

       ``1``
   * - Vapor pressure
     - ``vp``, ``pression``
     - ``mmHg``
     - ``mmHg``, ``Pa``, ``kPa``
   * - Surface pressure
     - ``sp``, ``surfacePressure``
     - ``Pa``
     - ``Pa``, ``kPa``
