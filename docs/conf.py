# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'pycequeau'
copyright = '2022, Eisinhower Rincon'
author = 'Eisinhower Rincon'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    "nbsphinx",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

templates_path = ['_templates']

source_suffix = ".rst"

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- html_theme_options = {}
html_theme_options = dict(
    repository_url="https://github.com/erinconv/pycequeau",
    repository_branch="main",
    path_to_docs="docs",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    home_page_in_toc=False,
    extra_navbar="",
    navbar_footer_text="",
)