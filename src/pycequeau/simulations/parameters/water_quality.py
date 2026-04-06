from __future__ import annotations

import numpy as np

from .base import BasinParameterBase


class WaterTemperatureParameters(BasinParameterBase):
    """
    Water-temperature-specific parameters under the water-quality family.

    At the moment the water-quality branch only contains the CEQUEAU water
    temperature model, but the container is intentionally scoped so additional
    quality modules can be added later without reshaping the top-level API.
    """

    def __init__(self, basin_structure) -> None:
        super().__init__(basin_structure)
        self.longwave_radiation = self._default_longwave_radiation()
        self.model_parameters: dict | None = None

    @classmethod
    def from_values(
        cls,
        basin_structure,
        values: np.ndarray,
        model: int = 1,
    ) -> "WaterTemperatureParameters":
        """Build a water-temperature parameter section from a raw value vector."""
        obj = cls(basin_structure)
        obj.set_model_parameters(values, model=model)
        return obj

    def _default_longwave_radiation(self) -> dict:
        """Return the default DLI parameter sets used by the temperature model."""
        return {
            "m1": {"u": 0, "v": 1, "a": 0.8171},
            "m2": {"u": 0.2184, "v": 2.4311, "a": 9.3645e-6},
            "m3": {"u": 0.17, "v": 4, "a": 0.2610, "b": -7.77e-2},
            "m4": {"u": 0.17, "v": 2, "a": 0.74, "b": 0.0065},
            "m5": {"u": 0.2187, "v": 1.6689, "a": 1.24, "b": 0.1429},
            "m6": {"u": 0.1296, "v": 4, "a": 1.08, "b": 2016},
            "m7": {"u": 0.1827, "v": 3.4472, "a": 46.5, "b": 1.2, "c": 3, "d": 0.5},
            "m8": {"u": 0.1533, "v": 4, "a": 0.72, "b": 0.0090, "c": 0.078},
        }

    def set_model_parameters(self, values: np.ndarray, model: int = 1) -> None:
        """
        Populate the water-temperature parameter block.

        The exported layout matches the ``qualite.cequeau.temperat`` structure
        expected by the CEQUEAU workflow.
        """
        if model != 1:
            raise ValueError(
                "No water temperature model with this label has ben yet added to the CEQUEAU model"
            )
        temperat = {
            "crayso": values[2],
            "crayin": values[3],
            "cevapo": values[4],
            "cconve": values[5],
            "crigel": values[6],
            "tnap": values[7],
            "bassol": values[8],
            "corsol": values[9],
            "panap": 1,
            "tinit": 0,
        }
        self.model_parameters = {
            "cequeau": {
                "coprom": values[0],
                "colarg": values[1],
                "temperat": temperat,
            }
        }

    def to_dict(self) -> dict:
        return {
            "qualite": self.model_parameters,
            "dli": self.longwave_radiation,
        }


class WaterQualityParameters(BasinParameterBase):
    """
    Top-level water-quality parameter group.

    Water temperature is currently the only implemented quality module, but it is
    kept as a child structure so future quality components can be added cleanly.
    """

    def __init__(self, basin_structure) -> None:
        super().__init__(basin_structure)
        self.temperature = WaterTemperatureParameters(basin_structure)

    @classmethod
    def from_values(
        cls,
        basin_structure,
        temperature_parameters: np.ndarray,
        temperature_model: int = 1,
    ) -> "WaterQualityParameters":
        """Build the water-quality parameter group from raw model values."""
        obj = cls(basin_structure)
        obj.set_temperature(temperature_parameters, model=temperature_model)
        return obj

    @property
    def qualite(self) -> dict | None:
        return self.temperature.model_parameters

    @property
    def dli(self) -> dict:
        return self.temperature.longwave_radiation

    def set_temperature(self, values: np.ndarray, model: int = 1) -> None:
        self.temperature.set_model_parameters(values, model=model)

    def to_dict(self) -> dict:
        return self.temperature.to_dict()
