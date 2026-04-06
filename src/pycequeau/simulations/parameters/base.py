from __future__ import annotations

import os
from abc import ABC, abstractmethod

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from ...core import projections
from ...physiographic.base import Basin


class SimulationParameterBase(ABC):
    """
    Shared base class for simulation-parameter sections.

    Each parameter family is responsible for exporting its own CEQUEAU-compatible
    dictionary fragment through :meth:`to_dict`.
    """

    def __init__(self, basin_structure: Basin) -> None:
        self.basin_structure = basin_structure

    @abstractmethod
    def to_dict(self) -> dict:
        """Return the CEQUEAU-compatible representation of the section."""


class BasinParameterBase(SimulationParameterBase):
    """
    Shared helper layer for parameter sections that depend on basin data.

    This class groups the basin-derived calculations shared by simulation
    parameter sections:

    - conversion of the watershed centroid latitude into the CEQUEAU ``xla`` code,
    - determination of the maximum-insolation day from the meteorological file,
    - empirical time-of-concentration estimators used to populate the transfer block.
    """

    @staticmethod
    def normalize_numeric_for_mex(obj):
        """Recursively convert numeric values to float64-compatible structures."""
        if isinstance(obj, dict):
            return {
                key: BasinParameterBase.normalize_numeric_for_mex(value)
                for key, value in obj.items()
            }

        if isinstance(obj, (list, tuple)):
            converted = [BasinParameterBase.normalize_numeric_for_mex(value) for value in obj]
            if converted and all(
                isinstance(value, (int, float, np.integer, np.floating, bool, np.bool_))
                for value in converted
            ):
                return np.asarray(converted, dtype=np.float64)
            return converted

        if isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.number):
                return obj.astype(np.float64)
            return obj

        if isinstance(obj, (np.integer, np.floating, int, float, bool, np.bool_)):
            return float(obj)

        return obj

    def day_max_insolation(self, meteo_file_name: str) -> int:
        """
        Compute the CEQUEAU insolation day index from the meteorological forcing.

        The CEQUEAU option structure uses ``jonei`` and ``joeva`` to represent the
        reference day of maximum insolation. This value is estimated from the day
        of year with the highest multi-annual mean ``tMax``.
        """
        meteo_file_name = os.path.join(
            self.basin_structure.project_path,
            "meteo",
            meteo_file_name,
        )
        nc_meteo = xr.open_dataset(meteo_file_name)
        tmax = nc_meteo["tMax"].to_pandas()
        time = pd.to_datetime(nc_meteo["pasTemp"].values - 719529, unit="D")
        tmax.index = time
        multi_annual = tmax.groupby(tmax.index.day_of_year).mean()
        idx_max_temp = multi_annual.idxmax()
        jour = multi_annual.index[idx_max_temp].values
        return int(np.mean(jour - 365 / 4))

    def compute_xla(self) -> int:
        """
        Compute the CEQUEAU ``xla`` code from the basin centroid latitude.

        CEQUEAU expects latitude in the compact integer form ``XXxx`` instead of a
        decimal degree value ``XX.xx``. For example, ``46.73`` becomes ``4673``.
        """
        watershed = gpd.read_file(self.basin_structure.watershed_shapefile)
        centroid = watershed.centroid
        epsg_code = centroid.crs.srs
        x = centroid.x.values.tolist()
        y = centroid.y.values.tolist()
        _, y = projections.utm_to_latlon(x, y, epsg_code)
        y_str = str(round(y[0], 2)).replace(".", "")
        if len(y_str) < 4:
            y_str = y_str + "0"
        return int(y_str)

    def compute_tc(self) -> tuple[float, dict[str, float]]:
        """
        Compute the basin time of concentration and return the contributing methods.

        The exported transfer block stores both the average time of concentration
        and the individual empirical estimates used to derive it.
        """
        cp_struct_name = os.path.join(
            self.basin_structure.project_path,
            "results",
            "carreauxPartiels.csv",
        )
        carreaux_partiels = pd.read_csv(cp_struct_name, index_col=0)
        hmin = carreaux_partiels["altitudeMoy"].min()
        hmax = carreaux_partiels["altitudeMoy"].max()
        hmean = carreaux_partiels["altitudeMoy"].mean()
        basin_area = (
            self.basin_structure.bassinVersant.get("superficieCE")
            * carreaux_partiels["pctSurface"].sum()
            / 100
        )
        main_channel_length = self.basin_structure.bassinVersant.get("longCanalPrincipal")
        if main_channel_length in (None, 0):
            raise ValueError(
                "The basin structure does not contain 'longCanalPrincipal'. "
                "Regenerate the basin structure before computing simulation parameters."
            )
        methods = {
            "Kirpich": self.tc_kirpich(main_channel_length, hmax - hmin),
            "Giandotti": self.tc_giandotti(
                basin_area,
                main_channel_length,
                hmean - hmin,
            ),
            "Department_of_public_work": self.tc_pw(main_channel_length, hmax - hmin),
        }
        avg_tc = pd.DataFrame.from_dict(methods, orient="index")[0].mean(axis=0)
        return avg_tc, methods

    @staticmethod
    def tc_pw(main_channel_length: float, height_differences: float) -> float:
        r"""
        This method uses the Department of Public Works (1995) formulation to
        compute the time of concentration:

        .. math::

            T_{c} = 60 \left ( 11.9 \frac{L^{3}}{H} \right )^{0.385}

        where:

        - :math:`L` is the length of the main channel in :math:`mi`
        - :math:`H` is the maximum elevation difference in :math:`ft`
        - :math:`T_{c}` is the time of concentration in :math:`min`

        The returned value is converted to days.
        """
        tc = 60 * (
            11.9 * (0.000621371 * main_channel_length) ** 3 / (3.28084 * height_differences)
        ) ** 0.385
        return tc / 1440

    @staticmethod
    def tc_giandotti(
        basin_area: float,
        main_channel_length: float,
        height_differences: float,
    ) -> float:
        r"""
        This method uses the Giandotti (1934) formulation to compute the time of
        concentration:

        .. math::

            T_{c} = \frac{4\sqrt{A} + 1.5L}{0.8\sqrt{H}}

        where:

        - :math:`A` is the area of the basin in :math:`\text{km}^{2}`
        - :math:`L` is the length of the main channel in :math:`\text{km}`
        - :math:`H` is the difference between the mean basin elevation and the
          outlet elevation
        - :math:`T_{c}` is the time of concentration in hours

        This formula was calibrated on 12 basins with drainage areas between
        170 and 70 000 :math:`\text{km}^{2}`.

        The returned value is converted to days.
        """
        tc = (4 * np.sqrt(basin_area) + 1.5 * main_channel_length / 1e3) / (
            0.8 * np.sqrt(height_differences)
        )
        return tc / 24

    @staticmethod
    def tc_kirpich(main_channel_length: float, height_differences: float) -> float:
        r"""
        This method uses the Kirpich (1940) formulation to compute the time of
        concentration:

        .. math::

            T_{c} = 0.0078 L^{0.77}S^{-0.385}

        where :math:`L` is the length of the main channel in :math:`ft` and
        :math:`S` is the mean basin slope in :math:`\text{m m}^{-1}` computed as:

        .. math::

            S = \frac{H_{max} - H_{min}}{L}

        :math:`T_{c}` is obtained in minutes and then converted to days.
        """
        s_basin = height_differences / main_channel_length
        tc = 0.0078 * (3.28084 * main_channel_length) ** 0.77 * (s_basin ** (-0.385))
        return tc / 1440

    def to_dict(self) -> dict:
        raise NotImplementedError("Context helpers do not export a standalone section.")
