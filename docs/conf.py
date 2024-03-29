# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
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
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# html_static_path = ['_static']

html_static_path = []

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

master_doc = 'index'
