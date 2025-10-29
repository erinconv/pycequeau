pycequeau
=======================================

Pycequeau is a library that obtains the necessary files to build and run the CEQUEAU hydrological-water temperature model.

To use this library it is highly recommended to clone the repository and create a fresh conda environment to have all the needed dependencies. To avoid compatibility issues, the new conda enviroment must contain a specific version of the GDAL library.

.. code-block:: bash
    
    conda create -n pycequeau python=3.10 gdal=3.11.0

Then, the package can be installed using the following command:

.. code-block:: bash

    pip install pycequeau

Alternatively, you can build and install the package from the source code by copying the git repository into your local machine and running the following command:

.. code-block:: bash

    git clone https://github.com/erinconv/pycequeau.git
    cd pycequeau
    pip install .

The tutorials and examples on how to use the library can be found in: https://pycequeau.readthedocs.io/en/latest/. The documentation is still under development, so any feedback/suggestions or contributions are welcome.