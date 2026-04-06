from __future__ import annotations

import xarray as xr

from ..physiographic.base import Basin
from .schema import DEFAULT_METEO_SCHEMA, MeteoSchema


class Meteo:
    """Base class for meteorological datasets used by pycequeau."""

    def __init__(
        self,
        basin_struct: Basin,
        ds: xr.Dataset,
        schema: MeteoSchema | None = None,
    ) -> None:
        """Store the basin context and the meteorological dataset."""
        self.basin_struct = basin_struct
        self.schema = schema or DEFAULT_METEO_SCHEMA
        self.ds = ds

    @property
    def variables(self) -> list[str]:
        """Return the meteorological variable names available in the dataset."""
        return list(self.ds.data_vars)

    @classmethod
    def cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        """Convert an interpolated meteorological dataset to the CEQUEAU layout."""
        return cls._cequeau_grid(ds, basin_struct)

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        """Implement the CEQUEAU export for a concrete meteorological source."""
        raise NotImplementedError
