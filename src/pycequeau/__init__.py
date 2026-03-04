"""Pycequeau - A library for the CEQUEAU hydrological-water temperature model."""
from __future__ import annotations

import re
from pathlib import Path

from pycequeau import core, meteo, physiographic, simulations

__author__ = "Eisinhower Rincon"
__email__ = "eisinhower.rincon@inrs.ca"


def _read_version_from_pyproject() -> str:
    """Read project version from pyproject.toml."""
    pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
    try:
        content = pyproject_path.read_text(encoding="utf-8")
    except OSError:
        return "0.0.0"

    match = re.search(r'^version\s*=\s*"([^"]+)"', content, flags=re.MULTILINE)
    return match.group(1) if match else "0.0.0"


__version__ = _read_version_from_pyproject()

__all__ = ["core", "meteo", "physiographic", "simulations"]
