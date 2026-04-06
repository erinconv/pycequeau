Contributing
============

``pycequeau`` is under active development, and contributions are welcome.
You can help by reporting bugs, improving the documentation, adding tutorials,
refining the API, or proposing new features.

If you find a problem, please open an issue on the project repository:
https://github.com/erinconv/pycequeau/issues

Ways to contribute
------------------

There are several useful ways to contribute to the project:

- Report bugs or unexpected behavior.
- Improve or correct the documentation.
- Add examples or tutorials for common workflows.
- Propose enhancements to the library.
- Submit code fixes or new features.

Before starting a larger change, it is a good idea to open an issue first so
the scope and implementation can be discussed with the maintainers.

Development setup
-----------------

Because the project depends on geospatial libraries, using Conda is
recommended.

1. Clone the repository:

.. code-block:: bash

   git clone https://github.com/erinconv/pycequeau.git
   cd pycequeau

2. Create and activate an environment:

.. code-block:: bash

   conda create -n pycequeau python=3.10 gdal=3.11.0
   conda activate pycequeau

   # Alternatively, use the provided environment file
   conda env create -f environment.yml
   conda activate pycequeau

3. Install the package in editable mode:

.. code-block:: bash

   pip install -e .

4. Install development extras for documentation work:

.. code-block:: bash

   pip install -e .[dev]

Branching and workflow
----------------------

Create a dedicated branch for each contribution. Use a short and descriptive
name, for example:

.. code-block:: bash

   git checkout -b fix-water-quality-docs

or

.. code-block:: bash

   git checkout -b add-netcdf-example

When preparing a contribution:

- Keep changes focused on a single topic whenever possible.
- Follow the existing project structure and naming patterns.
- Add or update docstrings and documentation when behavior changes.
- Include examples or notebooks when they make a feature easier to understand.

Code style
----------

There is no strict formatter enforced by the project yet, but contributors are
encouraged to follow standard Python best practices:

- Follow PEP 8 style conventions.
- Prefer clear function and variable names.
- Add docstrings to public classes and functions.
- Include type hints when they improve clarity.
- Keep new code consistent with the surrounding module.

The project currently includes ``autopep8`` as a development dependency if you
want to use it locally.

Testing your changes
--------------------

At a minimum, make sure the package still installs and imports correctly after
your changes.

Build and install from source:

.. code-block:: bash

   pip install -e .

Optional package build test:

.. code-block:: bash

   python -m build

If you add new functionality, include tests when practical and verify that the
affected modules still run as expected in your environment.

Building the documentation
--------------------------

If your contribution changes documentation, API docstrings, or notebooks, build
the Sphinx documentation locally before opening a pull request.

From the project root:

.. code-block:: bash

   cd docs
   make html

On Windows:

.. code-block:: bash

   cd docs
   make.bat html

The generated site will usually be available in ``docs/_build/html/`` when
using ``make`` or ``make.bat``. The repository also includes a
``docs/build_docs.bat`` helper that writes output to ``docs/build/html/``.

Pull requests
-------------

When your branch is ready:

1. Commit your changes with a clear message.
2. Push the branch to your fork or to the project branch you are using.
3. Open a pull request describing:

   - what changed,
   - why the change is needed,
   - how it was tested,
   - and whether documentation was updated.

A good pull request is easier to review when it is small, well described, and
limited to one main objective.

Documentation contributions
---------------------------

Documentation improvements are especially valuable for ``pycequeau``. Helpful
contributions include:

- clarifying installation steps,
- improving tutorial explanations,
- correcting outdated API examples,
- adding workflow notes for GIS preprocessing,
- and documenting common pitfalls.

If something was confusing while you were learning the library, improving that
part of the documentation is already a meaningful contribution.

Thank you
---------

Thank you for helping improve ``pycequeau``. Community feedback, bug reports,
and pull requests all help make the project more reliable and easier to use.
