# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('..'))


project = 'pycequeau'
copyright = '2023, Eisinhower Rincon'
author = 'Eisinhower Rincon'
release = 'beta'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "nbsphinx",
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
]
templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'build',
    'Thumbs.db',
    '.DS_Store',
]

# Read the Docs does not provide the full native GIS stack used by pycequeau.
# Mock these imports so autodoc can still import modules and render the API
# reference instead of generating empty pages.
autodoc_mock_imports = [
    "cftime",
    "geopandas",
    "matplotlib",
    "matplotlib.pyplot",
    "netCDF4",
    "ogr",
    "osgeo",
    "osr",
    "pandas",
    "pyproj",
    "rasterio",
    "rasterstats",
    "requests",
    "rioxarray",
    "scipy",
    "scipy.io",
    "scipy.io.matlab",
    "shapely",
    "shapely.geometry",
    "shapely.validation",
    "xarray",
    "xrspatial",
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'furo'

html_theme_options = {
    "navigation_with_keys": True,
}

# html_static_path = ['_static']

html_static_path = []

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

master_doc = 'index'

# -- Autodoc: avoid duplicate descriptions ---------------------------------
# ``pycequeau.meteo`` (and the thin ``meteo_calculator`` re-export module) expose
# calculator classes that are defined under ``pycequeau.meteo.calculators.*``.
# Documenting them again on those pages duplicates the canonical autoclass entries
# and triggers Sphinx "duplicate object description" warnings.

_METEO_CALC_REEXPORT_DOCNAMES = frozenset({
    'api/pycequeau.meteo',
})
_METEO_CALC_CLASS_NAMES = frozenset({
    'MeteoCalculator',
    'VaporPressureCalculator',
    'WindSpeedCalculator',
})
_SIMULATION_REEXPORT_DOCNAMES = frozenset({
    'api/pycequeau.simulations',
})
_SIMULATION_REEXPORT_CLASS_NAMES = frozenset({
    'SimulationParameterBase',
    'BasinParameterBase',
    'PycequeauParams',
    'HydrologicalParameters',
    'SnowmeltParameters',
    'EvapotranspirationParameters',
    'WaterQualityParameters',
    'WaterTemperatureParameters',
})


def _autodoc_skip_meteo_calculator_reexports(app, what, name, obj, skip, options):
    if app.env.docname not in _METEO_CALC_REEXPORT_DOCNAMES:
        return None
    if name not in _METEO_CALC_CLASS_NAMES:
        return None
    return True


def _autodoc_skip_simulations_reexports(app, what, name, obj, skip, options):
    if app.env.docname not in _SIMULATION_REEXPORT_DOCNAMES:
        return None
    if name not in _SIMULATION_REEXPORT_CLASS_NAMES:
        return None
    return True


def setup(app):
    app.connect('autodoc-skip-member', _autodoc_skip_meteo_calculator_reexports)
    app.connect('autodoc-skip-member', _autodoc_skip_simulations_reexports)
