from __future__ import annotations

import numpy as np
import xarray as xr

from .base import MeteoCalculator


class WindSpeedCalculator(MeteoCalculator):
    """Compose daily wind speed from daily mean u10 and v10 components."""

    variable_name = "wind_speed"
    default_output_name = "wind"
    source_variable_groups = (
        ("u10", "10m_u_component_of_wind"),
        ("v10", "10m_v_component_of_wind"),
    )

    @classmethod
    def _build_output_dataset(
        cls,
        source_dataarrays: dict[str, xr.DataArray],
        *,
        output_name: str,
        **kwargs,
    ) -> xr.Dataset:
        """Build a wind-speed dataset from two daily wind-component inputs."""
        source_names = list(source_dataarrays)
        u10 = source_dataarrays[source_names[0]]
        v10 = source_dataarrays[source_names[1]]
        cls._validate_compatible_inputs(u10, v10)
        wind_speed = cls.wind_speed_from_components_dataarray(
            u10,
            v10,
            output_name=output_name,
        )
        result = xr.Dataset(
            data_vars={output_name: wind_speed},
            coords=u10.coords,
            attrs={},
        )
        result.attrs["derivation_variable"] = cls.variable_name
        result.attrs["derivation_source_variable"] = (
            f"{u10.name or 'u10'}, {v10.name or 'v10'}"
        )
        result.attrs["derivation_method"] = "daily_vector_magnitude_from_daily_mean_components"
        return result

    @classmethod
    def wind_speed_from_components_dataarray(
        cls,
        u10: xr.DataArray,
        v10: xr.DataArray,
        *,
        output_name: str = "wind",
    ) -> xr.DataArray:
        """Compose wind speed from two data-array wind components."""
        speed = cls.wind_speed_from_components_array(u10, v10)
        speed = speed.rename(output_name)
        speed.attrs = dict(u10.attrs)
        speed.attrs["long_name"] = "Daily 10 m wind speed"
        speed.attrs["standard_name"] = "wind_speed"
        speed.attrs["source_variables"] = f"{u10.name or 'u10'}, {v10.name or 'v10'}"
        speed.attrs["source_time_step"] = "daily"
        speed.attrs["aggregation"] = "daily_vector_magnitude_from_daily_mean_components"
        speed.attrs["conversion_method"] = "wind speed composed from daily mean u10 and v10"
        speed.attrs["conversion_equation"] = "wind = sqrt(u10^2 + v10^2)"
        return speed

    @staticmethod
    def wind_speed_from_components_array(
        u10: np.ndarray | xr.DataArray,
        v10: np.ndarray | xr.DataArray,
    ) -> np.ndarray | xr.DataArray:
        """Compose wind speed from array-like wind components."""
        return np.sqrt(u10**2 + v10**2)

    @staticmethod
    def _validate_compatible_inputs(u10: xr.DataArray, v10: xr.DataArray) -> None:
        """Validate that the wind-component inputs are dimensionally compatible."""
        if u10.dims != v10.dims:
            raise ValueError(f"Wind component dimensions do not match: {u10.dims} vs {v10.dims}")

        for dim in u10.dims:
            if dim not in v10.coords:
                raise ValueError(f"Coordinate '{dim}' is missing from the second wind component.")
            if not np.array_equal(u10[dim].values, v10[dim].values):
                raise ValueError(f"Coordinate '{dim}' does not match between wind components.")
