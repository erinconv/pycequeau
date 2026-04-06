from __future__ import annotations

import os
from enum import IntEnum

import numpy as np
import pandas as pd

from .base import BasinParameterBase


class SnowmeltModel(IntEnum):
    """Snowmelt model ids aligned with the CEQUEAU source code."""

    CEQUEAU = 1
    CEMA_NEIGE = 2
    UEB = 3


SNOWMELT_MODEL_NAMES = {
    SnowmeltModel.CEQUEAU: "CEQUEAU",
    SnowmeltModel.CEMA_NEIGE: "CemaNeige",
    SnowmeltModel.UEB: "UEB",
}

SNOWMELT_EXPORT_KEYS = {
    SnowmeltModel.CEQUEAU: "cequeau",
    SnowmeltModel.CEMA_NEIGE: "cemaNeige",
    SnowmeltModel.UEB: "UEB",
}


class EvapotranspirationModel(IntEnum):
    """Evapotranspiration model ids aligned with ``Simulation.h``."""

    CEQUEAU = 1
    KPENMAN = 2
    PRIESTLEYTAYLOR = 3
    MCGUINNESS = 4
    PENMONT = 5
    MORTON = 6


EVAPOTRANSPIRATION_MODEL_NAMES = {
    EvapotranspirationModel.CEQUEAU: "Thornthwaite",
    EvapotranspirationModel.KPENMAN: "Kimberley-Penman",
    EvapotranspirationModel.PRIESTLEYTAYLOR: "Priestley-Taylor",
    EvapotranspirationModel.MCGUINNESS: "McGuinness",
    EvapotranspirationModel.PENMONT: "Penman-Monteith",
    EvapotranspirationModel.MORTON: "Morton",
}

EVAPOTRANSPIRATION_EXPORT_KEYS = {
    EvapotranspirationModel.CEQUEAU: "cequeau",
    EvapotranspirationModel.KPENMAN: "KPenman",
    EvapotranspirationModel.PRIESTLEYTAYLOR: "PriestleyTaylor",
    EvapotranspirationModel.MCGUINNESS: "McGuinness",
    EvapotranspirationModel.PENMONT: "PenmanMonteith",
    EvapotranspirationModel.MORTON: "Morton",
}


class SnowmeltParameters(BasinParameterBase):
    """
    Snowmelt parameter structure for the supported hydrological models.

    The CEQUEAU parameter file stores several snowmelt model blocks in the
    ``fonte`` section. This class keeps the snow-specific defaults and
    model-dependent parameters in one place.
    """

    def __init__(self, basin_structure) -> None:
        super().__init__(basin_structure)
        self.data: dict | None = None

    def set_parameters(self, values: np.ndarray) -> None:
        """
        Populate the snowmelt parameter dictionaries.

        Parameters
        ----------
        values
            Ordered parameter vector used to populate the snowmelt structures.
            The first seven values correspond to the degree-day CEQUEAU model:

            ``strne_s, tfc_s, tfd_s, tsc_s, tsd_s, ttd, tts_s``.
        The exported structure contains the CEQUEAU, UEB, and CemaNeige blocks
        regardless of the active model selected in the simulation options.
        """
        values = np.asarray(values, dtype=float)
        builder_map = {
            SnowmeltModel.CEQUEAU: self._build_cequeau_parameters,
            SnowmeltModel.CEMA_NEIGE: self._build_cema_neige_parameters,
            SnowmeltModel.UEB: self._build_ueb_parameters,
        }
        self.data = {
            SNOWMELT_EXPORT_KEYS[snowmelt_model]: builder(values)
            for snowmelt_model, builder in builder_map.items()
        }

    @staticmethod
    def _build_cequeau_parameters(values: np.ndarray) -> dict:
        return {
            "strne_s": values[0],
            "tfc_s": values[1],
            "tfd_s": values[2],
            "tsc_s": values[3],
            "tsd_s": values[4],
            "ttd": values[5],
            "tts_s": values[6],
        }

    @staticmethod
    def _build_ueb_parameters(values: np.ndarray) -> dict:
        return {
            "strne_s": values[0],
            "K_s": 0.15,
            "z0": 0.003,
            "aep": 0.2,
            "K_sat": 350,
            "rho_s": 450,
            "melt_frac": 0.99,
            "melt_thr": 0,
            "hours": 16,
            "z": 2.0,
            "avo": 0.8,
            "airo": 0.6,
            "Lc": 0.05,
            "fstab": 1,
            "D": 0.001,
            "de": 0.4,
            "snow_temp_method": 1,
            "w": 0,
            "ub": -2000,
            "E": 0,
            "tausn": 0,
            "tsurf": 0,
            "tave": 0,
            "Mr": 0,
            "albedo": 0.25,
        }

    @staticmethod
    def _build_cema_neige_parameters(values: np.ndarray) -> dict:
        return {
            "strne": values[0],
            "Kf": 15,
            "Tf": 0,
            "CTg": 0.85,
            "theta": 0.8,
            "Gseuil": 250 * 0.9,
            "Vmin": 0.5,
            "Zmed": 300,
            "eTg": 0.1,
            "G": 0,
        }

    def apply_option_defaults(self, option: dict | None, sol_initial: dict | None) -> None:
        """Inject option- and state-dependent values into the CEQUEAU snow block."""
        if self.data is None or option is None or sol_initial is None:
            return
        cequeau_key = SNOWMELT_EXPORT_KEYS[SnowmeltModel.CEQUEAU]
        self.data[cequeau_key]["jonei"] = option["jonei"]
        self.data[cequeau_key]["tmur"] = sol_initial["tmur"]
        self.data[cequeau_key]["tstock"] = sol_initial["tstock"]

    def to_dict(self) -> dict:
        return self.data or {}


class EvapotranspirationParameters(BasinParameterBase):
    """
    Evapotranspiration parameter structure for supported models.

    The selected model controls which keys are exported in the ``evapo`` block.
    """

    def __init__(self, basin_structure) -> None:
        super().__init__(basin_structure)
        self.data: dict | None = None

    def set_parameters(self, values: np.ndarray) -> None:
        """
        Populate the evapotranspiration section for all supported CEQUEAU models.

        The exported structure contains one parameter block per supported model.
        The active evapotranspiration model is controlled separately through the
        simulation option block.
        """
        values = np.asarray(values, dtype=float)
        builder_map = {
            EvapotranspirationModel.CEQUEAU: self._build_cequeau_parameters,
            EvapotranspirationModel.KPENMAN: self._build_kpenman_parameters,
            EvapotranspirationModel.PRIESTLEYTAYLOR: self._build_priestley_taylor_parameters,
            EvapotranspirationModel.MCGUINNESS: self._build_mcguinness_parameters,
            EvapotranspirationModel.PENMONT: self._build_penman_monteith_parameters,
            EvapotranspirationModel.MORTON: self._build_morton_parameters,
        }
        self.data = {
            EVAPOTRANSPIRATION_EXPORT_KEYS[evapo_model]: builder(values)
            for evapo_model, builder in builder_map.items()
        }

    @staticmethod
    def _build_cequeau_parameters(values: np.ndarray) -> dict:
        return {
            "evnap": values[0],
            "xaa": values[1],
            "xit": values[2],
        }

    @staticmethod
    def _build_kpenman_parameters(values: np.ndarray) -> dict:
        return {
            "evnap": values[0],
        }

    @staticmethod
    def _build_priestley_taylor_parameters(values: np.ndarray) -> dict:
        return {
            "alpha": values[0],
            "evnap": values[1],
        }

    @staticmethod
    def _build_mcguinness_parameters(values: np.ndarray) -> dict:
        return {
            "evnap": values[0],
        }

    @staticmethod
    def _build_penman_monteith_parameters(values: np.ndarray) -> dict:
        return {
            "evnap": values[0],
        }

    @staticmethod
    def _build_morton_parameters(values: np.ndarray) -> dict:
        return {
            "alpha": values[0],
            "evnap": values[1],
        }

    def apply_option_defaults(self, option: dict | None) -> None:
        """Propagate option-dependent values into the evapotranspiration blocks."""
        if self.data is None or option is None:
            return
        cequeau_key = EVAPOTRANSPIRATION_EXPORT_KEYS[EvapotranspirationModel.CEQUEAU]
        self.data[cequeau_key]["joeva"] = option["joeva"]

    def to_dict(self) -> dict:
        return self.data or {}


class HydrologicalParameters(BasinParameterBase):
    """
    Container for the hydrological part of the CEQUEAU parameter structure.

    This section groups the option, soil, initial-state, transfer, snowmelt, and
    evapotranspiration blocks that are required by the hydrological model.
    """

    def __init__(self, basin_structure) -> None:
        super().__init__(basin_structure)
        self.option: dict | None = None
        self.sol: dict | None = None
        self.solInitial: dict | None = None
        self.transfert: dict | None = None
        self.snowmelt = SnowmeltParameters(basin_structure)
        self.evapotranspiration = EvapotranspirationParameters(basin_structure)
        self.time_of_concentrations: dict[str, float] | None = None

    @classmethod
    def from_values(
        cls,
        basin_structure,
        flow_parameters: np.ndarray,
        initial_conditions: np.ndarray,
        transfer_parameters: np.ndarray,
        simulation_options: np.ndarray,
        snow_parameters: np.ndarray,
        evapotranspiration_parameters: np.ndarray,
        meteo_file_name: str,
    ) -> "HydrologicalParameters":
        """
        Build the hydrological parameter group from the raw CEQUEAU input vectors.

        This constructor hides the ordering dependencies of the setter-based API
        by computing the insolation day and propagating dependent values
        internally before returning the fully initialized section.
        """
        obj = cls(basin_structure)
        jonei = obj.day_max_insolation(meteo_file_name)
        obj.set_soil(flow_parameters)
        obj.set_initial_soil_conditions(initial_conditions)
        obj.set_transfer(transfer_parameters)
        obj.set_option(simulation_options, jonei)
        obj.set_snowmelt(snow_parameters)
        obj.set_evapotranspiration(evapotranspiration_parameters)
        return obj

    def set_option(self, values: np.ndarray, jonei: int | None) -> None:
        """
        Populate the simulation option block.

        The ``jonei`` and ``joeva`` entries are computed from the basin
        meteorological forcing through :meth:`day_max_insolation`.
        """
        if jonei is None:
            raise ValueError(
                "The maximum insolation day has not been computed yet. "
                "Call day_max_insolation(...) before set_option(...)."
            )
        self.option = {
            "ipassim": float(24),
            "moduleFonte": float(values[0]),
            "moduleEvapo": float(values[1]),
            "moduleOmbrage": float(0),
            "moduleDLI": float(1),
            "calculQualite": float(values[2]),
            "jonei": float(jonei),
            "joeva": float(jonei),
        }
        self.snowmelt.apply_option_defaults(self.option, self.solInitial)
        self.evapotranspiration.apply_option_defaults(self.option)

    def set_soil(self, values: np.ndarray) -> None:
        """
        Populate the ``sol`` block.

        The exported ``xla`` value is derived from the watershed centroid because
        CEQUEAU stores latitude as an integer code rather than decimal degrees.
        """
        carreux_entier_name = os.path.join(
            self.basin_structure.project_path,
            "results",
            "carreauxEntiers.csv",
        )
        _ = pd.read_csv(carreux_entier_name, index_col=0)
        self.sol = {
            "cin_s": values[0],
            "cvmar": values[1],
            "cvnb_s": values[2],
            "cvnh_s": values[3],
            "cvsb": values[4],
            "cvsi_s": values[5],
            "xinfma": values[6],
            "hinf_s": values[7],
            "hint_s": values[8],
            "hmar": values[9],
            "hnap_s": values[10],
            "hpot_s": values[11],
            "hsol_s": values[12],
            "hrimp_s": values[13],
            "tri_s": 0.0,
            "xla": float(self.compute_xla()),
        }

    def set_initial_soil_conditions(self, values: np.ndarray) -> None:
        """Populate the initial storage conditions in ``solInitial``."""
        self.solInitial = {
            "hsini": values[0],
            "hnini": values[1],
            "hmini": values[2],
            "q0": values[3],
            "tmur": values[4],
            "tstock": values[5],
        }
        self.snowmelt.apply_option_defaults(self.option, self.solInitial)

    def set_transfer(self, values: np.ndarray) -> None:
        """
        Populate the transfer block and derive the basin travel time.

        ``zn`` is exported as the average of the empirical time-of-concentration
        equations returned by :meth:`compute_tc`.
        """
        avg_tc, methods = self.compute_tc()
        self.time_of_concentrations = methods
        self.transfert = {
            "exxkt": values[0],
            "zn": avg_tc,
            "tc_struct": methods,
        }

    def set_snowmelt(self, values: np.ndarray) -> None:
        """Populate the hydrological snowmelt structure."""
        self.snowmelt.set_parameters(values)
        self.snowmelt.apply_option_defaults(self.option, self.solInitial)

    def set_evapotranspiration(self, values: np.ndarray) -> None:
        """Populate the evapotranspiration structure."""
        self.evapotranspiration.set_parameters(values)
        self.evapotranspiration.apply_option_defaults(self.option)

    @property
    def fonte(self) -> dict | None:
        return self.snowmelt.data

    @property
    def evapo(self) -> dict | None:
        return self.evapotranspiration.data

    def to_dict(self) -> dict:
        return {
            "option": self.option,
            "sol": self.sol,
            "solInitial": self.solInitial,
            "transfert": self.transfert,
            "fonte": self.snowmelt.to_dict(),
            "evapo": self.evapotranspiration.to_dict(),
        }
