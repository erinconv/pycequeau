pycequeau
=======================================

Pycequeau is a library that obtains the necessary files to build and run the CEQUEAU hydrological-water temperature model.

To use this library, first clone the ithub repository

.. code-block:: bash
    
    git clone git@github.com:erinconv/pycequeau.git
    cd pycequeau

This library was developed and tested in python 3.7, given incompatibilities with GDAL and newer python versions. To use it correctly, we recommend to use Anaconda as python interpreter. In the ``pycequeau`` folder, create the conda environment using the ``environment.yml`` file as follows:

.. code-block:: bash
    
    conda env create -f environment.yml