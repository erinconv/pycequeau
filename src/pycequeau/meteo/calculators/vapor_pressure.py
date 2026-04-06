from __future__ import annotations

import numpy as np
import xarray as xr

from ...core import UnitHandler
from .base import MeteoCalculator

SATURATED_VAPOR_PRESSURE_REFERENCES = {
    "lowe_1977": {
        "label": "Lowe (1977)",
        "url": (
            "https://journals.ametsoc.org/view/journals/apme/16/1/1520-0450_1977_016_0100_aapftc_2_0_co_2.xml"
        ),
    },
    "murray_1967": {
        "label": "Murray (1967)",
        "url": (
            "https://journals.ametsoc.org/view/journals/apme/6/1/1520-0450_1967_006_0203_otcosv_2_0_co_2.xml"
        ),
    },
}


class VaporPressureCalculator(MeteoCalculator):
    """Compute saturated vapor pressure from dewpoint temperature.

    The saturated vapor pressure is evaluated at the dewpoint temperature.
    Under that condition, the saturated vapor pressure is equal to the actual
    vapor pressure.
    """

    variable_name = "vapor_pressure"
    default_output_name = "vp"
    source_variable_groups = (("d2m", "dewpoint_temperature"),)

    @classmethod
    def _build_output_dataset(
        cls,
        source_dataarrays: dict[str, xr.DataArray],
        *,
        output_name: str,
        vapor_pressure_method: str = "lowe_1977",
        **kwargs,
    ) -> xr.Dataset:
        dewpoint_temperature = source_dataarrays[next(iter(source_dataarrays))]
        vapor_pressure = cls.vapor_pressure_from_dewpoint_dataarray(
            dewpoint_temperature,
            method=vapor_pressure_method,
        ).rename(output_name)

        result = xr.Dataset(
            data_vars={output_name: vapor_pressure},
            coords=dewpoint_temperature.coords,
            attrs={},
        )
        result.attrs["derivation_method"] = vapor_pressure_method
        result.attrs["derivation_variable"] = cls.variable_name
        result.attrs["derivation_source_variable"] = (
            dewpoint_temperature.name or "dewpoint_temperature"
        )
        result.attrs["derivation_reference"] = cls._get_reference(vapor_pressure_method)["label"]
        result.attrs["derivation_reference_url"] = cls._get_reference(vapor_pressure_method)["url"]
        return result

    @classmethod
    def vapor_pressure_from_dewpoint_dataarray(
        cls,
        dewpoint_temperature: xr.DataArray,
        method: str = "lowe_1977",
    ) -> xr.DataArray:
        """Compute vapor pressure from a dewpoint-temperature data array.

        This method first converts the dewpoint temperature to degrees Celsius,
        then computes the saturated vapor pressure in :math:`Pa`, and finally
        converts the result to :math:`mmHg`.

        The resulting vapor pressure is:

        .. math::

            e_{a} = e_{s}(T_{dew})

        Where:
            - :math:`e_{a}` is the actual vapor pressure.
            - :math:`e_{s}` is the saturated vapor pressure function.
            - :math:`T_{dew}` is the dewpoint temperature.
        """
        source_unit = str(dewpoint_temperature.attrs.get("units", "")).strip()
        if not source_unit:
            raise ValueError(
                f"Variable '{dewpoint_temperature.name}' is missing the 'units' attribute."
            )

        dewpoint_celsius = UnitHandler.convert_temperature_to_celsius(
            dewpoint_temperature,
            source_unit,
        )
        vapor_pressure_pa = cls.saturated_vapor_pressure(
            dewpoint_celsius,
            method=method,
        )
        vapor_pressure_mmhg = UnitHandler.convert_vapor_pressure_to_mmhg(
            vapor_pressure_pa,
            "Pa",
        )

        reference = cls._get_reference(method)
        vapor_pressure_mmhg = vapor_pressure_mmhg.rename("vapor_pressure")
        vapor_pressure_mmhg.attrs = {
            "units": "mmHg",
            "long_name": "Saturated vapor pressure derived from dewpoint temperature",
            "source_variable": dewpoint_temperature.name or "dewpoint_temperature",
            "source_units": source_unit,
            "derivation_method": method,
            "derivation_reference": reference["label"],
            "derivation_reference_url": reference["url"],
        }
        return vapor_pressure_mmhg

    @classmethod
    def vapor_pressure_from_dewpoint_array(
        cls,
        dewpoint_temperature: np.ndarray,
        method: str = "lowe_1977",
    ) -> np.ndarray:
        """Compute vapor pressure from a dewpoint-temperature array.

        The vapor pressure is obtained by evaluating the saturated vapor
        pressure at the dewpoint temperature:

        .. math::

            e_{a} = e_{s}(T_{dew})

        The returned values are expressed in :math:`mmHg`.
        """
        vapor_pressure_pa = cls.saturated_vapor_pressure(
            dewpoint_temperature,
            method=method,
        )
        return UnitHandler.convert_vapor_pressure_to_mmhg(vapor_pressure_pa, "Pa")

    @classmethod
    def saturated_vapor_pressure(
        cls,
        temperature_celsius: np.ndarray | xr.DataArray,
        method: str = "lowe_1977",
    ) -> np.ndarray | xr.DataArray:
        """Compute saturated vapor pressure from air temperature.

        The formulation is selected with ``method``:

            - ``lowe_1977`` uses the polynomial approximation from Lowe (1977).
            - ``murray_1967`` uses the exponential form used in the Murray
              implementation adopted in the legacy CORDEX workflow.

        The result is returned in :math:`Pa`.
        """
        if method == "lowe_1977":
            return cls._saturated_vapor_pressure_lowe_1977(temperature_celsius)
        if method == "murray_1967":
            return cls._saturated_vapor_pressure_murray_1967(temperature_celsius)
        supported = ", ".join(sorted(SATURATED_VAPOR_PRESSURE_REFERENCES))
        raise ValueError(
            f"Unsupported vapor-pressure method '{method}'. Supported methods are: {supported}."
        )

    @classmethod
    def _saturated_vapor_pressure_lowe_1977(
        cls,
        temperature_celsius: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Compute saturated vapor pressure with the Lowe (1977) formulation.

        The saturated vapor pressure is computed with two sixth-order
        polynomials, one over liquid water and one over ice:

        .. math::

            e_{s,w}(T) = 100 \left(a_{0} + T\left(a_{1} + T\left(a_{2} +
            T\left(a_{3} + T\left(a_{4} + T\left(a_{5} + a_{6}T\right)\right)\right)\right)\right)\right)

        .. math::

            e_{s,i}(T) = 100 \left(a_{0} + T\left(a_{1} + T\left(a_{2} +
            T\left(a_{3} + T\left(a_{4} + T\left(a_{5} + a_{6}T\right)\right)\right)\right)\right)\right)

        Where:
            - :math:`T` is the temperature in :math:`^\circ C`.
            - :math:`e_{s,w}` is the saturated vapor pressure over water in :math:`Pa`.
            - :math:`e_{s,i}` is the saturated vapor pressure over ice in :math:`Pa`.

        The water polynomial is used for :math:`T \\ge 0`, and the ice
        polynomial is used for :math:`T < 0`.
        """
        values = temperature_celsius
        over_water = cls._lowe_over_water(values)
        over_ice = cls._lowe_over_ice(values)
        if isinstance(values, xr.DataArray):
            return xr.where(values >= 0.0, over_water, over_ice)
        return np.where(values >= 0.0, over_water, over_ice)

    @classmethod
    def _saturated_vapor_pressure_murray_1967(
        cls,
        temperature_celsius: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Compute saturated vapor pressure with the Murray (1967) formulation.

        Murray gives the convenient exponential form:

        .. math::

            e_{s,w}(T_{k}) = 100 \times 6.1078 \exp\left(
            17.2693882 \frac{T_{k} - 273.16}{T_{k} - 35.86}\right)

        .. math::

            e_{s,i}(T_{k}) = 100 \times 6.1078 \exp\left(
            21.8745584 \frac{T_{k} - 273.16}{T_{k} - 7.66}\right)

        Where:
            - :math:`T_{k}` is the temperature in :math:`K`.
            - :math:`e_{s,w}` is the saturated vapor pressure over water in :math:`Pa`.
            - :math:`e_{s,i}` is the saturated vapor pressure over ice in :math:`Pa`.
            - over water, :math:`u = 7.5` and :math:`v = 237.3`.
            - over ice, :math:`u = 9.5` and :math:`v = 265.5`.
            - the transformed coefficients satisfy :math:`a = u \ln 10` and
              :math:`b = 273.16 - v`.

        The water branch is used for :math:`T \\ge 0`, and the ice branch is
        used for :math:`T < 0`.
        """
        temperature_kelvin = temperature_celsius + 273.15
        over_water = cls._murray_over_water(temperature_kelvin)
        over_ice = cls._murray_over_ice(temperature_kelvin)
        values = temperature_celsius
        if isinstance(values, xr.DataArray):
            return xr.where(values >= 0.0, over_water, over_ice)
        return np.where(values >= 0.0, over_water, over_ice)

    @staticmethod
    def _murray_over_water(
        temperature_kelvin: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Murray (1967) branch for temperatures over water."""
        a = 17.2693882
        b = 35.86
        vapor_pressure_hpa = 6.1078 * np.exp(
            a * (temperature_kelvin - 273.16) / (temperature_kelvin - b)
        )
        return 100.0 * vapor_pressure_hpa

    @staticmethod
    def _murray_over_ice(
        temperature_kelvin: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Murray (1967) branch for temperatures over ice."""
        a = 21.8745584
        b = 7.66
        vapor_pressure_hpa = 6.1078 * np.exp(
            a * (temperature_kelvin - 273.16) / (temperature_kelvin - b)
        )
        return 100.0 * vapor_pressure_hpa

    @staticmethod
    def _lowe_over_water(
        temperature_celsius: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Lowe (1977) branch for temperatures over water."""
        a0 = 6.107799961
        a1 = 0.4436518521
        a2 = 0.01428945805
        a3 = 0.0002650648471
        a4 = 3.031240936e-6
        a5 = 2.034080948e-8
        a6 = 6.136820929e-11
        vapor_pressure_hpa = a0 + temperature_celsius * (
            a1
            + temperature_celsius * (
                a2
                + temperature_celsius * (
                    a3
                    + temperature_celsius * (a4 + temperature_celsius * (a5 + a6 * temperature_celsius))
                )
            )
        )
        return 100.0 * vapor_pressure_hpa

    @staticmethod
    def _lowe_over_ice(
        temperature_celsius: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Lowe (1977) branch for temperatures over ice."""
        a0 = 6.109177956
        a1 = 0.503469897
        a2 = 0.01886013408
        a3 = 0.0004176223716
        a4 = 5.82472028e-6
        a5 = 4.838803174e-8
        a6 = 1.838826904e-10
        vapor_pressure_hpa = a0 + temperature_celsius * (
            a1
            + temperature_celsius * (
                a2
                + temperature_celsius * (
                    a3
                    + temperature_celsius * (a4 + temperature_celsius * (a5 + a6 * temperature_celsius))
                )
            )
        )
        return 100.0 * vapor_pressure_hpa

    @staticmethod
    def _get_reference(method: str) -> dict[str, str]:
        try:
            return SATURATED_VAPOR_PRESSURE_REFERENCES[method]
        except KeyError as exc:
            supported = ", ".join(sorted(SATURATED_VAPOR_PRESSURE_REFERENCES))
            raise ValueError(
                f"Unsupported vapor-pressure method '{method}'. Supported methods are: {supported}."
            ) from exc
