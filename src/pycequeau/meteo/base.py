from __future__ import annotations

import xarray as xr

from ..physiographic.base import Basin


class Meteo:
    """Shared basin context for meteorological data processors."""

    def __init__(self, basin_struct: Basin) -> None:
        self.basin_struct = basin_struct

    @classmethod
    def cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        """Convert an interpolated meteorological dataset to the CEQUEAU layout."""
        return cls._cequeau_grid(ds, basin_struct)

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        raise NotImplementedError
