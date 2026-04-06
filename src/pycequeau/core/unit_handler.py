"""Utilities for normalizing meteorological units."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from ..meteo.schema import VariableSpec


class UnitHandler:
    """Convert meteorological variables to their canonical units."""

    @classmethod
    def normalize_unit_text(cls, unit: str) -> str:
        unit = unit.strip()
        unit = unit.replace("**", "^")
        unit = unit.replace("day-1", "d-1")
        unit = unit.replace("/day", " d-1")
        unit = unit.replace("/d", " d-1")
        unit = unit.replace("degc", "c")
        unit = unit.replace("(", "").replace(")", "")
        unit = " ".join(unit.split())
        unit = unit.replace(" - ", "-")
        unit = unit.replace("- ", "-")
        unit = unit.replace(" -", "-")
        return unit.lower()

    @classmethod
    def convert_dataarray_to_canonical_units(
        cls,
        data_array: xr.DataArray,
        spec: VariableSpec,
    ) -> xr.DataArray:
        source_unit = str(data_array.attrs.get("units", "")).strip()
        if not source_unit:
            raise ValueError(f"Variable '{data_array.name}' is missing the 'units' attribute.")

        converted = cls.convert_array_to_canonical_units(
            data_array,
            source_unit,
            spec,
        )

        converted.attrs = dict(data_array.attrs)
        converted.attrs["units"] = spec.canonical_unit
        converted.attrs["source_units"] = source_unit
        return converted

    @classmethod
    def convert_array_to_canonical_units(
        cls,
        values: np.ndarray | xr.DataArray,
        source_unit: str,
        spec: VariableSpec,
    ) -> np.ndarray | xr.DataArray:
        canonical_name = spec.canonical_name
        if canonical_name in {
            "temperature_max",
            "temperature_min",
            "dewpoint_temperature",
        }:
            return cls._convert_temperature_to_celsius(values, source_unit)
        if canonical_name == "precipitation":
            return cls._convert_precipitation_to_mm_per_day(values, source_unit)
        if canonical_name in {"shortwave_radiation", "longwave_radiation"}:
            return cls._convert_radiation_to_mj_per_m2_day(values, source_unit)
        if canonical_name == "cloud_cover":
            return cls._convert_cloud_cover_to_fraction(values, source_unit)
        if canonical_name == "wind_speed":
            return cls._convert_wind_to_km_per_hour(values, source_unit)
        if canonical_name == "relative_humidity":
            return cls._convert_relative_humidity(values, source_unit)
        if canonical_name == "vapor_pressure":
            return cls._convert_vapor_pressure_to_mmhg(values, source_unit)
        if canonical_name == "surface_pressure":
            return cls._convert_surface_pressure_to_pa(values, source_unit)
        return values

    @classmethod
    def convert_temperature_to_celsius(
        cls,
        values: np.ndarray | xr.DataArray,
        source_unit: str,
    ) -> np.ndarray | xr.DataArray:
        return cls._convert_temperature_to_celsius(values, source_unit)

    @classmethod
    def convert_vapor_pressure_to_mmhg(
        cls,
        values: np.ndarray | xr.DataArray,
        source_unit: str,
    ) -> np.ndarray | xr.DataArray:
        return cls._convert_vapor_pressure_to_mmhg(values, source_unit)

    @classmethod
    def _convert_temperature_to_celsius(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"c", "°c"}:
            return values
        if normalized in {"k", "kelvin"}:
            return values - 273.15
        raise ValueError(f"Unsupported temperature unit '{unit}'. Expected C or K.")

    @classmethod
    def _convert_precipitation_to_mm_per_day(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"mm d-1", "mm"}:
            return values
        if normalized in {"m d-1", "m", "meter", "metre", "m of water equivalent"}:
            return values * 1000.0
        if normalized in {"kg m-2 s-1"}:
            return values * 86400.0
        raise ValueError(
            f"Unsupported precipitation unit '{unit}'. Expected mm d-1, m d-1, or kg m-2 s-1."
        )

    @classmethod
    def _convert_radiation_to_mj_per_m2_day(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"mj m-2 d-1", "mj m-2"}:
            return values
        if normalized in {"j m-2 d-1", "j m-2"}:
            return values / 1_000_000.0
        if normalized in {"w m-2", "w m^-2", "w m^(-2)"}:
            return values * 0.0864
        raise ValueError(
            f"Unsupported radiation unit '{unit}'. Expected MJ m-2 d-1, J m-2, or W m-2."
        )

    @classmethod
    def _convert_cloud_cover_to_fraction(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"0-1", "fraction", "1"}:
            return values
        if normalized in {"%", "percent"}:
            return values / 100.0
        raise ValueError(f"Unsupported cloud-cover unit '{unit}'. Expected 0-1 or %.")

    @classmethod
    def _convert_wind_to_km_per_hour(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"km h-1", "km h^-1", "km/h"}:
            return values
        if normalized in {"m s-1", "m s^-1", "m/s"}:
            return values * 3.6
        raise ValueError(f"Unsupported wind unit '{unit}'. Expected km h-1 or m s-1.")

    @classmethod
    def _convert_relative_humidity(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"%", "percent"}:
            return values
        if normalized in {"0-1", "fraction", "1"}:
            return values * 100.0
        raise ValueError(f"Unsupported relative-humidity unit '{unit}'. Expected % or 0-1.")

    @classmethod
    def _convert_vapor_pressure_to_mmhg(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"mmhg"}:
            return values
        if normalized in {"pa"}:
            return values * 0.00750062
        if normalized in {"kpa"}:
            return values * 7.50062
        raise ValueError(f"Unsupported vapor-pressure unit '{unit}'. Expected mmHg, Pa, or kPa.")

    @classmethod
    def _convert_surface_pressure_to_pa(
        cls,
        values: np.ndarray | xr.DataArray,
        unit: str,
    ) -> np.ndarray | xr.DataArray:
        normalized = cls.normalize_unit_text(unit)
        if normalized in {"pa"}:
            return values
        if normalized in {"kpa"}:
            return values * 1000.0
        raise ValueError(f"Unsupported surface-pressure unit '{unit}'. Expected Pa or kPa.")
