from __future__ import annotations

import json
import os

from scipy.io import savemat

from ...core import matlab_utils as matu
from ...physiographic.base import Basin
from .base import BasinParameterBase
from .hydrology import HydrologicalParameters
from .water_quality import WaterQualityParameters


class PycequeauParams(BasinParameterBase):
    """Top-level CEQUEAU parameter facade.

    This class assembles pre-built hydrology and water-quality sections into the
    final CEQUEAU parameter dictionary and handles export to JSON/MAT formats.
    """

    def __init__(
        self,
        bassinVersant: Basin,
        hydrology: HydrologicalParameters | None = None,
        water_quality: WaterQualityParameters | None = None,
        *,
        ctp: float = 0.0,
        lac: float = 0.0,
        surface: float = 0.0,
    ) -> None:
        super().__init__(bassinVersant)
        self.ctp = ctp
        self.lac = lac
        self.surface = surface
        self.hydrology = hydrology or HydrologicalParameters(bassinVersant)
        self.water_quality = water_quality or WaterQualityParameters(bassinVersant)
        self.parametres = None

    def create_parameter_structure(self):
        """
        Assemble and return the final CEQUEAU parameter dictionary.
        """
        sections = {}
        sections.update(self.hydrology.to_dict())
        sections.update(self.water_quality.to_dict())
        sections.update(
            {
                "ctp": self.ctp,
                "lac": self.lac,
                "surface": self.surface,
            }
        )
        self.parametres = self.normalize_numeric_for_mex(sections)
        return self.parametres

    def export_parameter_structure_json(self, file_name: str = "parameters.json"):
        """Export the current parameter dictionary to JSON."""
        if self.parametres is None:
            raise ValueError(
                "Parameter structure not built yet. Call create_parameter_structure() first."
            )
        out_json = os.path.join(self.basin_structure.project_path, "results", file_name)
        with open(out_json, "w", encoding="utf-8") as outfile:
            json.dump(self.parametres, outfile, indent=4, default=tuple)

    def export_parameter_structure_mat(self, file_name: str = "parameters.mat"):
        """Export current parameter dictionary to MATLAB/Octave .mat format."""
        if self.parametres is None:
            raise ValueError(
                "Parameter structure not built yet. Call create_parameter_structure() first."
            )
        out_mat = os.path.join(self.basin_structure.project_path, "results", file_name)
        savemat(
            out_mat,
            {"parametres": matu.to_mat_compatible(self.parametres)},
            long_field_names=True,
            do_compression=True,
        )
