from __future__ import annotations

import os
from abc import ABC, abstractmethod

import xarray as xr


class MeteoCalculator(ABC):
    """Base class for explicit meteorological preprocessing calculations."""

    registry: dict[str, type["MeteoCalculator"]] = {}
    variable_name: str | None = None
    default_output_name: str | None = None
    source_variable_groups: tuple[tuple[str, ...], ...] = ()

    def __init_subclass__(cls, **kwargs) -> None:
        """Register concrete calculator subclasses by derived variable name."""
        super().__init_subclass__(**kwargs)
        if cls.variable_name:
            MeteoCalculator.registry[cls.variable_name] = cls

    @classmethod
    def available_derivations(cls) -> tuple[str, ...]:
        """Return the derived meteorological variables supported by the registry."""
        return tuple(cls.registry)

    @classmethod
    def get_calculator_class(cls, variable: str) -> type["MeteoCalculator"]:
        """Resolve the calculator class that is responsible for a derived variable."""
        try:
            return cls.registry[variable]
        except KeyError as exc:
            supported = ", ".join(cls.available_derivations())
            raise ValueError(
                f"Unsupported derived variable '{variable}'. Supported derived variables are: {supported}."
            ) from exc

    @classmethod
    def create_variable_dataset(
        cls,
        inputs: str | list[str] | tuple[str, ...],
        variable: str,
        *,
        source_variable: str | tuple[str, ...] | list[str] | None = None,
        output_name: str | None = None,
        **kwargs,
    ) -> xr.Dataset:
        """Build a derived-variable dataset from one or more NetCDF inputs."""
        calculator_class = cls.get_calculator_class(variable)
        source_dataarrays = calculator_class._load_required_sources(
            inputs,
            calculator_class._resolve_required_source_groups(source_variable),
        )
        target_name = output_name or calculator_class.default_output_name
        if target_name is None:
            raise ValueError(f"Calculator '{variable}' does not define a default output name.")
        return calculator_class._build_output_dataset(
            source_dataarrays,
            output_name=target_name,
            **kwargs,
        )

    @classmethod
    def create_variable_file(
        cls,
        inputs: str | list[str] | tuple[str, ...],
        variable: str,
        output_path: str | None = None,
        *,
        source_variable: str | tuple[str, ...] | list[str] | None = None,
        output_name: str | None = None,
        **kwargs,
    ) -> str:
        """Compute a derived variable and write it to a NetCDF file."""
        ds = cls.create_variable_dataset(
            inputs,
            variable,
            source_variable=source_variable,
            output_name=output_name,
            **kwargs,
        )
        destination = output_path or cls._default_output_path(inputs, ds, output_name)
        ds.to_netcdf(destination)
        return destination

    @classmethod
    def _resolve_required_source_groups(
        cls,
        source_variable: str | tuple[str, ...] | list[str] | None,
    ) -> tuple[tuple[str, ...], ...]:
        """Resolve source-variable overrides against the calculator requirements."""
        if source_variable is None:
            return cls.source_variable_groups

        if isinstance(source_variable, str):
            if len(cls.source_variable_groups) != 1:
                raise ValueError(
                    "A single 'source_variable' override can only be used for "
                    "single-source derivations."
                )
            return ((source_variable,),)

        override_values = tuple(source_variable)
        if len(cls.source_variable_groups) == 1:
            return (override_values,)

        if len(override_values) != len(cls.source_variable_groups):
            raise ValueError(
                "The 'source_variable' override must provide one source name for each "
                "required source input."
            )

        return tuple((name,) for name in override_values)

    @classmethod
    def _load_required_sources(
        cls,
        inputs: str | list[str] | tuple[str, ...],
        required_sources: tuple[tuple[str, ...], ...],
    ) -> dict[str, xr.DataArray]:
        """Load the source data arrays required to compute a derived variable."""
        loaded_sources: dict[str, xr.DataArray] = {}

        for candidate_group in required_sources:
            data_array = None
            resolved_name = None

            if isinstance(inputs, str) and os.path.isdir(inputs):
                file_path, resolved_name = cls._find_variable_file(inputs, candidate_group)
                with xr.open_dataset(file_path, engine="netcdf4") as ds:
                    data_array = ds[resolved_name].load()
            elif isinstance(inputs, str):
                with xr.open_dataset(inputs, engine="netcdf4") as ds:
                    resolved_name = cls._resolve_variable_from_candidates(ds, candidate_group)
                    data_array = ds[resolved_name].load()
            else:
                for file_path in inputs:
                    with xr.open_dataset(file_path, engine="netcdf4") as ds:
                        try:
                            resolved_name = cls._resolve_variable_from_candidates(ds, candidate_group)
                        except ValueError:
                            continue
                        data_array = ds[resolved_name].load()
                        break

            if data_array is None or resolved_name is None:
                requested = ", ".join(candidate_group)
                raise ValueError(
                    f"Could not find a NetCDF input containing one of the source variables: {requested}."
                )

            loaded_sources[resolved_name] = data_array

        return loaded_sources

    @classmethod
    def _default_output_path(
        cls,
        inputs: str | list[str] | tuple[str, ...],
        ds: xr.Dataset,
        output_name: str | None,
    ) -> str:
        """Infer a default NetCDF output path from the source input location."""
        target_name = next(iter(ds.data_vars), output_name or cls.default_output_name or "output")
        if isinstance(inputs, str) and os.path.isdir(inputs):
            return os.path.join(inputs, f"{target_name}.nc")

        if isinstance(inputs, str):
            folder = os.path.dirname(inputs)
            extension = os.path.splitext(inputs)[1] or ".nc"
            return os.path.join(folder, f"{target_name}{extension}")

        first_path = inputs[0]
        folder = os.path.dirname(first_path)
        extension = os.path.splitext(first_path)[1] or ".nc"
        return os.path.join(folder, f"{target_name}{extension}")

    @staticmethod
    def _find_variable_file(
        folder_path: str,
        variable_names: tuple[str, ...],
    ) -> tuple[str, str]:
        """Find the first NetCDF file in a folder containing one of the requested variables."""
        for file_name in sorted(os.listdir(folder_path)):
            if not file_name.endswith(".nc"):
                continue
            file_path = os.path.join(folder_path, file_name)
            with xr.open_dataset(file_path, engine="netcdf4") as ds:
                try:
                    resolved_name = MeteoCalculator._resolve_variable_from_candidates(
                        ds,
                        variable_names,
                    )
                except ValueError:
                    continue
                return file_path, resolved_name

        raise ValueError(
            "Could not find a NetCDF file containing any of the variables "
            f"{', '.join(variable_names)} in '{folder_path}'."
        )

    @staticmethod
    def _resolve_variable_from_candidates(
        ds: xr.Dataset,
        variable_names: tuple[str, ...],
    ) -> str:
        """Resolve the first matching variable name found in a dataset."""
        for variable_name in variable_names:
            if variable_name in ds.data_vars:
                return variable_name
        available = ", ".join(ds.data_vars)
        requested = ", ".join(variable_names)
        raise ValueError(
            f"Dataset does not contain any of the requested variables ({requested}). "
            f"Available variables are: {available}."
        )

    @classmethod
    @abstractmethod
    def _build_output_dataset(
        cls,
        source_dataarrays: dict[str, xr.DataArray],
        *,
        output_name: str,
        **kwargs,
    ) -> xr.Dataset:
        """Build the derived-variable dataset for a concrete calculator."""
        raise NotImplementedError
