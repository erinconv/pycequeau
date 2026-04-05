from __future__ import annotations

import xarray as xr

from ..physiographic.base import Basin
from .base import Meteo


class StationMeteo(Meteo):
    """Placeholder for future surface-station meteorological support."""

    def __init__(self, basin_struct: Basin, ds: xr.Dataset) -> None:
        super().__init__(basin_struct, ds)
        raise NotImplementedError(
            "StationMeteo is not implemented yet. "
            "Surface-station support is planned for a future refactor."
        )

    @classmethod
    def from_dataset(cls, basin_struct: Basin, ds: xr.Dataset) -> "StationMeteo":
        raise NotImplementedError(
            "StationMeteo.from_dataset() is not implemented yet. "
            "Surface-station support is planned for a future refactor."
        )

    @classmethod
    def from_files(cls, basin_struct: Basin, data_path: str) -> "StationMeteo":
        raise NotImplementedError(
            "StationMeteo.from_files() is not implemented yet. "
            "Surface-station support is planned for a future refactor."
        )

    @classmethod
    def _cequeau_grid(cls, ds: xr.Dataset, basin_struct: Basin) -> xr.Dataset:
        raise NotImplementedError(
            "StationMeteo CEQUEAU export is not implemented yet."
        )
