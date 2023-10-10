pycequeau
=======================================

Pycequeau is a library that obtains the necessary files to build and run the CEQUEAU hydrological-water temperature model.

To use this library, first clone the ithub repository

.. code-block:: bash
    
    git clone git@github.com:erinconv/pycequeau.git
    cd pycequeau

The pycequeau packge has some dependencies that may cause some  compatibilities troubles, specially in Windows systems. It is highly recommended to create a fresh conda environment to have all the needed dependencies. 

.. code-block:: bash

    conda install -c conda-forge mamba
    mamba create -n testenv gdal geopandas matplotlib xarray netcdf4 shapely rasterstats pyproj pytest rasterstats pandas=1.3.5 numpy=1.20.3

