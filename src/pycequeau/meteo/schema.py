from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class VariableSpec:
    canonical_name: str
    export_name: str
    canonical_unit: str
    accepted_names: tuple[str, ...]
    unsupported: bool = False
    unsupported_message: str | None = None


UNSUPPORTED_MESSAGE_WIND_COMPONENTS = (
    "Wind components 'u10' and 'v10' are not supported directly. "
    "Preprocess them into a scalar daily 'wind' variable before ingestion."
)


VARIABLE_SPECS: tuple[VariableSpec, ...] = (
    VariableSpec(
        canonical_name="temperature_max",
        export_name="tMax",
        canonical_unit="C",
        accepted_names=("tasmax", "tMax"),
    ),
    VariableSpec(
        canonical_name="temperature_min",
        export_name="tMin",
        canonical_unit="C",
        accepted_names=("tasmin", "tMin"),
    ),
    VariableSpec(
        canonical_name="precipitation",
        export_name="pTot",
        canonical_unit="mm d-1",
        accepted_names=("tp", "pr", "pTot"),
    ),
    VariableSpec(
        canonical_name="shortwave_radiation",
        export_name="rayonnement",
        canonical_unit="MJ m-2 d-1",
        accepted_names=("ssrd", "ssr", "rsds", "rayonnement"),
    ),
    VariableSpec(
        canonical_name="longwave_radiation",
        export_name="longwaveRad",
        canonical_unit="MJ m-2 d-1",
        accepted_names=("strd", "msdwlwrf", "longwaveRad"),
    ),
    VariableSpec(
        canonical_name="dewpoint_temperature",
        export_name="d2m",
        canonical_unit="C",
        accepted_names=("d2m",),
    ),
    VariableSpec(
        canonical_name="cloud_cover",
        export_name="nebulosite",
        canonical_unit="0-1",
        accepted_names=("tcc", "clt", "nebulosite"),
    ),
    VariableSpec(
        canonical_name="wind_speed",
        export_name="vitesseVent",
        canonical_unit="km h-1",
        accepted_names=("wind", "sfcWind", "vitesseVent"),
    ),
    VariableSpec(
        canonical_name="relative_humidity",
        export_name="humiditeRelative",
        canonical_unit="%",
        accepted_names=("hurs", "humiditeRelative"),
    ),
    VariableSpec(
        canonical_name="vapor_pressure",
        export_name="pression",
        canonical_unit="mmHg",
        accepted_names=("vp", "pression"),
    ),
    VariableSpec(
        canonical_name="surface_pressure",
        export_name="surfacePressure",
        canonical_unit="Pa",
        accepted_names=("sp", "surfacePressure"),
    ),
    VariableSpec(
        canonical_name="unsupported_wind_u10",
        export_name="u10",
        canonical_unit="m s-1",
        accepted_names=("u10",),
        unsupported=True,
        unsupported_message=UNSUPPORTED_MESSAGE_WIND_COMPONENTS,
    ),
    VariableSpec(
        canonical_name="unsupported_wind_v10",
        export_name="v10",
        canonical_unit="m s-1",
        accepted_names=("v10",),
        unsupported=True,
        unsupported_message=UNSUPPORTED_MESSAGE_WIND_COMPONENTS,
    ),
)


NAME_TO_SPEC: dict[str, VariableSpec] = {
    name: spec
    for spec in VARIABLE_SPECS
    for name in spec.accepted_names
}

CANONICAL_TO_SPEC: dict[str, VariableSpec] = {
    spec.canonical_name: spec
    for spec in VARIABLE_SPECS
    if not spec.unsupported
}


class MeteoSchema:
    """Registry of supported meteorological variables and aliases."""

    def __init__(self) -> None:
        self.variable_specs = VARIABLE_SPECS
        self.name_to_spec = NAME_TO_SPEC
        self.canonical_to_spec = CANONICAL_TO_SPEC

    def get_variable_spec(self, variable_name: str) -> VariableSpec:
        try:
            spec = self.name_to_spec[variable_name]
        except KeyError as exc:
            supported = sorted(
                name
                for name, candidate in self.name_to_spec.items()
                if not candidate.unsupported
            )
            supported_names = ", ".join(supported)
            raise ValueError(
                f"Unsupported meteorological variable '{variable_name}'. "
                f"Supported input names are: {supported_names}."
            ) from exc

        if spec.unsupported:
            raise ValueError(spec.unsupported_message or f"Unsupported variable '{variable_name}'.")
        return spec

    def get_export_name(self, variable_name: str) -> str:
        spec = self.canonical_to_spec.get(variable_name) or self.name_to_spec.get(variable_name)
        return spec.export_name if spec is not None else variable_name


DEFAULT_METEO_SCHEMA = MeteoSchema()


def get_variable_spec(variable_name: str) -> VariableSpec:
    return DEFAULT_METEO_SCHEMA.get_variable_spec(variable_name)
