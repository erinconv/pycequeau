# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
from datetime import datetime
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('..'))

import pycequeau

project = 'pycequeau'
copyright = f'{datetime.now().year}, Eisinhower Rincon'
author = 'Eisinhower Rincon'
release = pycequeau.__version__
version = release


html_title = "Pycequeau Official Documentation"
html_short_title = "pycequeau"

# The name of the Pygments (syntax highlighting) style to use for light and dark themes.
pygments_style = "sas"
pygments_dark_style = "lightbulb"

html_theme_options = {
    "navigation_with_keys": True,
    "dark_css_variables": {
        "color-admonition-background": "#202428",
        "color-admonition-title-background": "rgba(61, 148, 255, 0.30)",
        "color-admonition-title-background--important": "rgba(0, 191, 165, 0.35)",
        "color-admonition-title-background--warning": "rgba(255, 145, 0, 0.35)",
        "color-admonition-title-background--note": "rgba(0, 176, 255, 0.35)",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/erinconv/pycequeau",
            "html": """
                <svg viewBox="0 0 20 20" aria-hidden="true">
                    <path fill="currentColor" d="M10,0 C15.523,0 20,4.59 20,10.253 C20,14.782 17.138,18.624 13.167,19.981 C12.66,20.082 12.48,19.762 12.48,19.489 C12.48,19.151 12.492,18.047 12.492,16.675 C12.492,15.719 12.172,15.095 11.813,14.777 C14.04,14.523 16.38,13.656 16.38,9.718 C16.38,8.598 15.992,7.684 15.35,6.966 C15.454,6.707 15.797,5.664 15.252,4.252 C15.252,4.252 14.414,3.977 12.505,5.303 C11.706,5.076 10.85,4.962 10,4.958 C9.15,4.962 8.295,5.076 7.497,5.303 C5.586,3.977 4.746,4.252 4.746,4.252 C4.203,5.664 4.546,6.707 4.649,6.966 C4.01,7.684 3.619,8.598 3.619,9.718 C3.619,13.646 5.954,14.526 8.175,14.785 C7.889,15.041 7.63,15.493 7.54,16.156 C6.97,16.418 5.522,16.871 4.63,15.304 C4.63,15.304 4.101,14.319 3.097,14.247 C3.097,14.247 2.122,14.234 3.029,14.87 C3.029,14.87 3.684,15.185 4.139,16.37 C4.139,16.37 4.726,18.2 7.508,17.58 C7.513,18.437 7.522,19.245 7.522,19.489 C7.522,19.76 7.338,20.077 6.839,19.982 C2.865,18.627 0,14.783 0,10.253 C0,4.59 4.478,0 10,0" />
                </svg>
            """,
            "class": "",
        },
    ],
}

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

# Tutorial notebooks in this project are documentation examples with placeholder
# paths and environment-specific setup. Render the stored notebook content
# instead of executing them during docs builds.
nbsphinx_execute = 'never'

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

# html_static_path = ['_static']

html_static_path = ['_static']

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
    app.add_css_file('custom.css')
    app.connect('autodoc-skip-member', _autodoc_skip_meteo_calculator_reexports)
    app.connect('autodoc-skip-member', _autodoc_skip_simulations_reexports)
